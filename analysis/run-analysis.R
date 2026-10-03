# All analyses of one or more results directories. Writes CSV tables to out_dir.
# Run from the repository root:
#   Rscript analysis/run-analysis.R <results_dir[,results_dir...]> <out_dir>
source("scripts/setup.R")
source("analysis/contrasts.R")
options(restart.cores = 3L)
args <- commandArgs(trailingOnly = TRUE)
out_dir <- args[2]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write <- \(x, name)
  readr::write_csv(x, file.path(out_dir, paste0(name, ".csv")))

cells <- design_table()
timed <- \(label, expr) {
  start <- Sys.time()
  value <- expr
  cat(sprintf(
    "%-28s %6.0f s\n",
    label,
    as.numeric(Sys.time() - start, units = "secs")
  ))
  value
}
# The per-cell summary is the slow step; reuse it if it was made from this directory.
study_file <- file.path(out_dir, "study.rds")
# The cache is valid only for the same set of result files (names and sizes).
dirs <- strsplit(args[1], ",")[[1]]
files <- unlist(lapply(
  dirs,
  list.files,
  pattern = "^task-.*\\.rds$",
  full.names = TRUE
))
# "terms": the cache format since per-term EDFs were added (2026-09-28)
fingerprint <- digest::digest(list(basename(files), file.size(files), "terms"))
cached <- if (file.exists(study_file)) readRDS(study_file)
study <- if (identical(cached$fingerprint, fingerprint)) {
  cached
} else {
  timed("summarise_study", summarise_study(dirs, cells))
}
study$fingerprint <- fingerprint
saveRDS(study, study_file)
tasks <- study$tasks
MAIN_FAMILIES <- c("gaussian", "poisson", "binomial")
PRIMARY <- c("mean", "beta")
ARMS_WITH_INTERVALS <- c(
  "reml_model",
  "reml_cl2",
  "ncv_model",
  "ncv_cl2",
  "ncv_cl2_bias",
  "ncv_cl2f_bias"
)

# Cell summaries, fit diagnostics, failures ----------------------------------------------

write(dplyr::left_join(cells, study$cells, by = "cell"), "cells")
write(study$failures, "failures")
fit_diagnostics <- study$fits |>
  dplyr::group_by(cell, fit) |>
  dplyr::summarise(
    n = dplyr::n(),
    converged = mean(converged),
    median_time = median(time),
    warnings = mean(n_warnings > 0),
    rho_mean = mean(rho),
    rho_boundary = mean(rho_boundary),
    .groups = "drop"
  )
write(fit_diagnostics, "fit-diagnostics")
# Per-term EDF (one value per fit; the terms table has one row per penalty) and log sp.
term_edf <- study$terms |>
  dplyr::group_by(cell, fit, term) |>
  dplyr::summarise(
    n_penalties = dplyr::n_distinct(penalty),
    mean_log_sp = mean(log_sp),
    .groups = "drop"
  ) |>
  dplyr::left_join(
    study$terms |>
      dplyr::distinct(cell, rep, fit, term, edf) |>
      dplyr::group_by(cell, fit, term) |>
      dplyr::summarise(
        n = dplyr::n(),
        edf_se = sd(edf) / sqrt(n),
        edf_mean = mean(edf),
        .groups = "drop"
      ),
    by = c("cell", "fit", "term")
  )
write(term_edf, "term-edf")

# Recommendation (§5) and question 1 --------------------------------------------------------

per_family_estimand <- \(f) {
  tidyr::expand_grid(family = MAIN_FAMILIES, estimand = PRIMARY) |>
    purrr::pmap(\(family, estimand) {
      dplyr::mutate(
        f(family, estimand),
        family = family,
        estimand = estimand,
        .before = 1
      )
    }) |>
    dplyr::bind_rows()
}
# The rule for the proposal and for the frequentist-CL2 variant.
PROPOSALS <- c(PROPOSAL, "ncv_cl2f_bias")
recs <- tidyr::expand_grid(
  proposal = PROPOSALS,
  family = MAIN_FAMILIES,
  estimand = PRIMARY
) |>
  purrr::pmap(\(proposal, family, estimand) {
    out <- recommend_with_uncertainty(
      tasks,
      study$cells,
      cells,
      family,
      estimand,
      proposal
    )
    list(
      decision = dplyr::mutate(out$decision, proposal = proposal, .before = 1),
      criteria = dplyr::mutate(
        out$criteria,
        proposal = proposal,
        family = family,
        estimand = estimand,
        .before = 1
      )
    )
  })
write(dplyr::bind_rows(purrr::map(recs, "decision")), "recommendation")
# The decision per signal level, with the same bootstrap flags as the pooled rule.
by_signal <- tidyr::expand_grid(
  proposal = PROPOSALS,
  family = MAIN_FAMILIES,
  estimand = PRIMARY
) |>
  purrr::pmap(\(proposal, family, estimand) {
    purrr::map(signal_levels(family), \(signal) {
      out <- recommend_with_uncertainty(
        tasks,
        study$cells,
        cells,
        family,
        estimand,
        proposal,
        signal
      )
      unresolved <- out$criteria$name[out$criteria$unresolved %in% TRUE]
      dplyr::mutate(
        out$decision,
        proposal = proposal,
        signal = signal,
        unresolved = paste(unresolved, collapse = ", "),
        .before = 1
      )
    })
  }) |>
  purrr::list_flatten() |>
  dplyr::bind_rows()
write(by_signal, "recommendation-by-signal")
write(dplyr::bind_rows(purrr::map(recs, "criteria")), "recommendation-criteria")
write(
  per_family_estimand(\(fam, est) question1(tasks, cells, fam, est)),
  "question1"
)

# Factor profiles of the core factorial (questions 1, 3, 4) --------------------------------

profile_values <- c(
  purrr::set_names(
    purrr::map(ARMS_WITH_INTERVALS, \(a) arm_value("coverage", a)),
    paste0("coverage_", ARMS_WITH_INTERVALS)
  ),
  log_mse_ratio_ncv_reml = log_ratio("sq_error", "ncv_cl2", "reml_cl2"),
  log_score_ratio_proposal_fallback = log_ratio("score", PROPOSAL, FALLBACK)
)
profiles <- per_family_estimand(\(fam, est) {
  tidyr::expand_grid(
    factor = c("error", "G", "signal"),
    value = names(profile_values)
  ) |>
    purrr::pmap(\(factor, value) {
      dplyr::mutate(
        factor_profile(tasks, cells, fam, est, factor, profile_values[[value]]),
        value = value
      )
    }) |>
    dplyr::bind_rows()
})
write(profiles, "profiles")

# Extensions vs their matching core cells (questions 2, 4, 5, 8, 9) ------------------------

# The term-type comparison sets beta (M1) against f (M2) as "the main term".
main_term <- dplyr::mutate(
  tasks,
  estimand = dplyr::if_else(estimand %in% c("beta", "f"), "main", estimand)
)
extension_values <- c(
  profile_values[c(
    "coverage_reml_model",
    "coverage_reml_cl2",
    "coverage_ncv_cl2",
    "coverage_ncv_cl2_bias",
    "coverage_ncv_cl2f_bias"
  )],
  profile_values["log_mse_ratio_ncv_reml"]
)
# Basis-size cells also get cross-cell log ratios (large vs default basis) of one arm's
# MSE and interval score, and the proposal/fallback interval-score ratio.
basis_values <- c(
  extension_values,
  profile_values["log_score_ratio_proposal_fallback"],
  log_mse_reml = log_value("sq_error", "reml_cl2"),
  log_mse_ncv = log_value("sq_error", "ncv_cl2"),
  log_score_reml_cl2 = log_value("score", "reml_cl2"),
  log_score_proposal = log_value("score", PROPOSAL)
)
extensions <- cells$cell[
  cells$block %in%
    c("dense_grid", "rough_truth", "term_type", "warp", "basis_size")
]
matched <- purrr::map(extensions, \(id) {
  core <- matching_core_cell(cells, id)
  if (length(core) != 1) return(NULL) # e.g. wiggly warp at G = 40: no iid comparator
  if (!id %in% tasks$cell) return(NULL) # extension not run (yet)
  cell <- as.list(cells[cells$cell == id, ])
  estimands <- if (cell$block == "term_type") "main" else
    scored_estimands(cell, PRIMARY)
  values <- if (cell$block == "basis_size") basis_values else extension_values
  tidyr::expand_grid(estimand = estimands, value = names(values)) |>
    purrr::pmap(\(estimand, value) {
      dplyr::mutate(
        matched_contrast(
          if (estimand == "main") main_term else tasks,
          id,
          core,
          estimand,
          values[[value]]
        ),
        cell = id,
        core_cell = core,
        estimand = estimand,
        value = value
      )
    }) |>
    dplyr::bind_rows()
}) |>
  dplyr::bind_rows()
write(matched, "matched-contrasts")

# Within-cell arm contrasts: AR(1) and pointwise NCV (questions 3, 7) ----------------------

within <- dplyr::bind_rows(
  dplyr::mutate(
    paired_contrast(tasks, "ar1_model", PROPOSAL, "score"),
    contrast = "ar1 vs proposal: interval score"
  ),
  dplyr::mutate(
    paired_contrast(tasks, "ar1_model", "ncv_cl2", "sq_error"),
    contrast = "ar1 vs ncv: squared error"
  ),
  dplyr::mutate(
    paired_contrast(tasks, "ncv_point", "reml_cl2", "sq_error"),
    contrast = "pointwise ncv vs reml: squared error"
  ),
  dplyr::mutate(
    paired_contrast(tasks, "ncv_cl2", "reml_cl2", "sq_error"),
    contrast = "ncv vs reml: squared error"
  )
) |>
  dplyr::mutate(
    ratio_lower = purrr::map_dbl(ratio_ci, 1),
    ratio_upper = purrr::map_dbl(ratio_ci, 2)
  ) |>
  dplyr::select(-ratio_ci)
write(within, "within-cell-contrasts")
cat("analysis written to", out_dir, "\n")
