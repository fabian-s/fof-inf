# Misspecified-truth stress test: cell summaries and SUMMARY.md from the task records of
# run-misspec.R. Run from the repository root (code copy):
#   Rscript misspec/analyse-misspec.R <results_dir> <out_dir> <ncv_xlarge_max_rep>
# Writes misspec-cells.csv, misspec-contrasts.csv, misspec-truths.csv, misspec-crn.csv,
# misspec-failures.csv, misspec-fits.csv and SUMMARY.md. No number in SUMMARY.md is typed
# by hand. Records of the reused T0 xlarge cells (study cells 113, 115) beyond the
# recomputed check replicates come from the study's own task files.

source("scripts/setup.R")
source("misspec/misspec-lib.R")

MS_ARMS <- tibble::tribble(
  ~arm,
  ~fit,
  ~se,
  "reml_model",
  "REML",
  "se_model_REML",
  "reml_cl2",
  "REML",
  "se_cl2_REML",
  "reml_cl2f",
  "REML",
  "se_cl2f_REML",
  "ncv_cl2",
  "NCV",
  "se_cl2_NCV"
)
Z <- qnorm(1 - (1 - LEVEL) / 2)

# Loading -------------------------------------------------------------------------------

# All records of a cell (list of task records), own files plus reused study files.
load_cell <- function(cell, results_dir, ncv_max) {
  own <- list.files(
    results_dir,
    sprintf("^task-c%03d-r[0-9]{3}\\.rds$", cell$cell),
    full.names = TRUE
  )
  recs <- lapply(own, readRDS)
  if (reused_cell(cell)) {
    have <- vapply(recs, `[[`, 0, "rep")
    extra <- setdiff(seq_len(N_REP), have)
    paths <- vapply(extra, \(r) study_record_path(cell$study_cell, r), "")
    study <- lapply(paths[file.exists(paths)], \(p) {
      x <- readRDS(p)
      x$source <- "study"
      x
    })
    recs <- c(recs, study)
  }
  # NCV in xlarge cells only for replicates <= ncv_max (also in reused study records)
  if (cell$basis == "xlarge") {
    recs <- lapply(recs, \(r) {
      if (r$rep <= ncv_max || !identical(r$status, "success")) return(r)
      r$records <- lapply(
        r$records,
        \(m) m[, !grepl("_NCV$", colnames(m)), drop = FALSE]
      )
      r
    })
  }
  recs
}

# Scoring --------------------------------------------------------------------------------

# Points x reps matrices of one arm (replicates in which the arm exists).
arm_matrices <- function(recs, estimand, spec) {
  est_col <- paste0("est_", spec$fit)
  ok <- vapply(
    recs,
    \(r) all(c(est_col, spec$se) %in% colnames(r$records[[estimand]])),
    TRUE
  )
  pick <- \(col)
    do.call(cbind, lapply(recs[ok], \(r) r$records[[estimand]][, col]))
  list(est = pick(est_col), se = pick(spec$se), truth = pick("truth"))
}

# Grid-average coverage (per replicate) with its MC SE, over a subset of points.
coverage_stats <- function(covered, rows = seq_len(nrow(covered))) {
  per_rep <- colMeans(covered[rows, , drop = FALSE])
  c(mean = mean(per_rep), se = sd(per_rep) / sqrt(length(per_rep)))
}

bump_region <- function() {
  p <- eval_points("beta")
  d2 <- ((p$s - MS_BUMP_CENTER[["s"]])^2 + (p$t - MS_BUMP_CENTER[["t"]])^2) /
    MS_BUMP_SD^2
  d2 <= MS_BUMP_REGION_SD^2
}

score_cell_arm <- function(m, estimand) {
  error <- m$est - m$truth
  covered <- abs(error) <= Z * m$se
  # bias-eliminated: interval re-centred by the Monte Carlo bias at each point
  bias <- rowMeans(error)
  covered_be <- abs(error - bias) <= Z * m$se
  sd_est <- apply(m$est, 1, sd)
  cov <- coverage_stats(covered)
  cov_be <- coverage_stats(covered_be)
  width <- colMeans(2 * Z * m$se)
  sq <- colMeans(error^2)
  out <- tibble::tibble(
    n_rep = ncol(m$est),
    coverage = cov[["mean"]],
    coverage_se = cov[["se"]],
    coverage_bias_elim = cov_be[["mean"]],
    coverage_bias_elim_se = cov_be[["se"]],
    width = mean(width),
    width_se = sd(width) / sqrt(length(width)),
    mse = mean(sq),
    mse_se = sd(sq) / sqrt(length(sq)),
    mean_absbias_sd = mean(abs(bias) / sd_est),
    se_sd_ratio = mean(sqrt(rowMeans(m$se^2)) / sd_est)
  )
  out <- dplyr::bind_cols(out, point_stats(error, covered))
  out$variance <- out$mse - out$bias2
  region <- if (estimand == "beta") bump_region() else NULL
  out$coverage_bump <- NA_real_
  out$coverage_bump_se <- NA_real_
  out$coverage_off <- NA_real_
  out$coverage_off_se <- NA_real_
  out$coverage_bump_bias_elim <- NA_real_
  if (!is.null(region)) {
    on <- coverage_stats(covered, which(region))
    off <- coverage_stats(covered, which(!region))
    out$coverage_bump <- on[["mean"]]
    out$coverage_bump_se <- on[["se"]]
    out$coverage_off <- off[["mean"]]
    out$coverage_off_se <- off[["se"]]
    out$coverage_bump_bias_elim <- coverage_stats(covered_be, which(region))[[
      "mean"
    ]]
  }
  list(row = out, per_rep = tibble::tibble(cov = colMeans(covered)))
}

# Per-replicate grid-average coverage of one arm, for paired contrasts.
rep_coverage <- function(recs, estimand, spec) {
  est_col <- paste0("est_", spec$fit)
  rows <- lapply(recs, \(r) {
    m <- r$records[[estimand]]
    if (!all(c(est_col, spec$se) %in% colnames(m))) return(NULL)
    tibble::tibble(
      rep = r$rep,
      coverage = mean(abs(m[, est_col] - m[, "truth"]) <= Z * m[, spec$se])
    )
  })
  dplyr::bind_rows(rows)
}

# Main ------------------------------------------------------------------------------------

analyse <- function(results_dir, out_dir, ncv_max) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  cells <- misspec_cells()
  key_cols <- c(
    "cell",
    "truth_label",
    "truth",
    "error",
    "G",
    "basis",
    "study_cell"
  )
  summaries <- list()
  per_rep <- list()
  failures <- list()
  crn <- list()
  fits <- list()
  for (i in seq_len(nrow(cells))) {
    cell <- cells[i, ]
    recs <- load_cell(cell, results_dir, ncv_max)
    ok <- vapply(recs, \(r) identical(r$status, "success"), TRUE)
    for (r in recs[!ok]) {
      failures[[length(failures) + 1]] <- cbind(
        cell = cell$cell,
        rep = r$rep,
        r$failure
      )
    }
    recs <- recs[ok]
    for (r in recs) {
      if (!is.null(r$crn) && "column" %in% names(r$crn)) {
        crn[[length(crn) + 1]] <- cbind(cell = cell$cell, rep = r$rep, r$crn)
      }
      if (is.null(r$source)) {
        fits[[length(fits) + 1]] <- cbind(cell = cell$cell, rep = r$rep, r$fits)
      }
    }
    if (!length(recs)) next
    for (estimand in model_estimands("M1")) {
      for (j in seq_len(nrow(MS_ARMS))) {
        spec <- MS_ARMS[j, ]
        m <- arm_matrices(recs, estimand, spec)
        if (is.null(m$est) || ncol(m$est) < 2) next
        s <- score_cell_arm(m, estimand)
        summaries[[length(summaries) + 1]] <- dplyr::bind_cols(
          cell[key_cols],
          tibble::tibble(estimand = estimand, arm = spec$arm),
          s$row
        )
        per_rep[[length(per_rep) + 1]] <- dplyr::bind_cols(
          cell[key_cols],
          tibble::tibble(estimand = estimand, arm = spec$arm),
          rep_coverage(recs, estimand, spec)
        )
      }
    }
  }
  out <- list(
    cells = dplyr::bind_rows(summaries),
    per_rep = dplyr::bind_rows(per_rep),
    failures = dplyr::bind_rows(failures),
    crn = dplyr::bind_rows(crn),
    fits = dplyr::bind_rows(fits),
    truths = describe_truths()
  )
  out$contrasts <- paired_contrasts(out$per_rep)
  write.csv(
    out$cells,
    file.path(out_dir, "misspec-cells.csv"),
    row.names = FALSE
  )
  write.csv(
    out$contrasts,
    file.path(out_dir, "misspec-contrasts.csv"),
    row.names = FALSE
  )
  write.csv(
    out$truths,
    file.path(out_dir, "misspec-truths.csv"),
    row.names = FALSE
  )
  write.csv(
    crn_summary(out$crn),
    file.path(out_dir, "misspec-crn.csv"),
    row.names = FALSE
  )
  write.csv(
    out$failures,
    file.path(out_dir, "misspec-failures.csv"),
    row.names = FALSE
  )
  write.csv(
    fit_summary(out$fits),
    file.path(out_dir, "misspec-fits.csv"),
    row.names = FALSE
  )
  writeLines(summary_md(out, ncv_max), file.path(out_dir, "SUMMARY.md"))
  invisible(out)
}

# Paired (common random numbers) coverage contrasts: each misspecified truth vs T0 in the
# same error/G/basis, and xlarge vs default for the same truth/error/G.
paired_contrasts <- function(per_rep) {
  drop <- c("cell", "truth", "study_cell")
  pair <- \(a, b, by, label) {
    b <- b[setdiff(
      names(b),
      setdiff(names(a), c(by, "estimand", "arm", "rep", "coverage"))
    )]
    dplyr::inner_join(a, b, by = c(by, "estimand", "arm", "rep")) |>
      dplyr::group_by(dplyr::across(dplyr::all_of(c(by, "estimand", "arm")))) |>
      dplyr::summarise(
        n_pairs = dplyr::n(),
        diff_se = sd(coverage.x - coverage.y) / sqrt(n_pairs),
        diff = mean(coverage.x - coverage.y),
        .groups = "drop"
      ) |>
      dplyr::mutate(contrast = label)
  }
  base <- dplyr::select(per_rep, -dplyr::all_of(drop))
  vs_t0 <- lapply(c("T1", "T2"), \(tl) {
    pair(
      dplyr::filter(base, truth_label == tl),
      dplyr::filter(base, truth_label == "T0") |> dplyr::select(-truth_label),
      c("error", "G", "basis"),
      paste(tl, "- T0")
    ) |>
      dplyr::mutate(truth_label = tl)
  })
  basis <- pair(
    dplyr::filter(base, basis == "xlarge") |> dplyr::select(-basis),
    dplyr::filter(base, basis == "default") |> dplyr::select(-basis),
    c("truth_label", "error", "G"),
    "xlarge - default"
  ) |>
    dplyr::mutate(basis = "xlarge vs default")
  dplyr::bind_rows(c(vs_t0, list(basis)))
}

crn_summary <- function(crn) {
  if (!nrow(crn)) return(data.frame())
  crn |>
    dplyr::group_by(cell, column) |>
    dplyr::summarise(
      n_rep = dplyr::n_distinct(rep),
      max_abs = max(max_abs),
      max_rel = max(max_rel),
      .groups = "drop"
    )
}

fit_summary <- function(fits) {
  if (!nrow(fits)) return(data.frame())
  fits |>
    dplyr::group_by(cell, fit) |>
    dplyr::summarise(
      n = dplyr::n(),
      converged = mean(converged),
      median_time = median(time),
      median_cl2_time = median(cl2_time),
      warnings = sum(n_warnings),
      .groups = "drop"
    )
}

# SUMMARY.md ------------------------------------------------------------------------------

fmt_cov <- \(m, se) sprintf("%.3f (%.3f)", m, se)

md_table <- function(df) {
  c(
    paste0("| ", paste(names(df), collapse = " | "), " |"),
    paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|"),
    apply(df, 1, \(r) paste0("| ", paste(r, collapse = " | "), " |"))
  )
}

summary_md <- function(out, ncv_max) {
  cells <- out$cells
  tr <- out$truths
  crn <- crn_summary(out$crn)
  n_fail <- nrow(out$failures)
  truth_tab <- data.frame(
    truth = paste(tr$truth, tr$name),
    `feature amplitude` = sprintf("%.2f", tr$feature_amplitude),
    `feature share of beta-term signal` = sprintf(
      "%.3f",
      tr$feature_signal_share
    ),
    `feature share of int beta^2` = sprintf("%.3f", tr$feature_l2_share),
    `rel. L2 proj. error default` = sprintf("%.3g", tr$rel_l2_default),
    `resid. signal share default` = sprintf(
      "%.3g",
      tr$resid_signal_share_default
    ),
    `rel. L2 proj. error xlarge` = sprintf("%.3g", tr$rel_l2_xlarge),
    `resid. signal share xlarge` = sprintf(
      "%.3g",
      tr$resid_signal_share_xlarge
    ),
    check.names = FALSE
  )
  cov_table <- \(est) {
    d <- dplyr::filter(cells, estimand == est)
    d$value <- fmt_cov(d$coverage, d$coverage_se)
    w <- tidyr::pivot_wider(
      d[c("truth_label", "error", "G", "basis", "arm", "value")],
      names_from = arm,
      values_from = value
    )
    as.data.frame(w)
  }
  be_table <- \(est) {
    d <- dplyr::filter(cells, estimand == est)
    d$value <- sprintf(
      "%.3f / %.3f / %.2f",
      d$coverage_bias_elim,
      d$pointwise_q05,
      d$mean_absbias_sd
    )
    as.data.frame(tidyr::pivot_wider(
      d[c("truth_label", "error", "G", "basis", "arm", "value")],
      names_from = arm,
      values_from = value
    ))
  }
  mse_table <- \(est) {
    d <- dplyr::filter(
      cells,
      estimand == est,
      arm %in% c("reml_cl2", "ncv_cl2")
    )
    d$value <- sprintf("%.4g (bias² %.2f)", d$mse, d$bias2 / d$mse)
    d$fit <- ifelse(d$arm == "reml_cl2", "REML", "NCV")
    as.data.frame(tidyr::pivot_wider(
      d[c("truth_label", "error", "G", "basis", "fit", "value")],
      names_from = fit,
      values_from = value
    ))
  }
  bump <- dplyr::filter(cells, estimand == "beta") |>
    dplyr::mutate(
      value = sprintf(
        "%s / %s",
        fmt_cov(coverage_bump, coverage_bump_se),
        fmt_cov(coverage_off, coverage_off_se)
      )
    )
  bump_tab <- as.data.frame(tidyr::pivot_wider(
    bump[c("truth_label", "error", "G", "basis", "arm", "value")],
    names_from = arm,
    values_from = value
  ))
  con <- out$contrasts |>
    dplyr::filter(arm %in% c("reml_cl2", "reml_model")) |>
    dplyr::mutate(value = sprintf("%+.3f (%.3f)", diff, diff_se))
  con_tab <- as.data.frame(tidyr::pivot_wider(
    con[c(
      "contrast",
      "truth_label",
      "error",
      "G",
      "basis",
      "arm",
      "estimand",
      "value"
    )],
    names_from = estimand,
    values_from = value
  ))
  n_rep <- range(cells$n_rep)
  crn_line <- if (nrow(crn)) {
    sprintf(
      "Common random numbers (T0 vs the study's records of the same cells, %d cell x column pairs, %d to %d replicates each): max relative difference %.2g, max absolute %.2g.",
      nrow(crn),
      min(crn$n_rep),
      max(crn$n_rep),
      max(crn$max_rel),
      max(crn$max_abs)
    )
  } else {
    "No common-random-numbers check was available."
  }
  c(
    "# Misspecified-truth stress test (generated by misspec/analyse-misspec.R)",
    "",
    sprintf(
      "Gaussian, R² 0.5, D = 61, M1. Replicates per cell x arm: %d to %d. NCV in xlarge cells: replicates <= %d. Nominal level %.2f, z intervals. Failed tasks: %d.",
      n_rep[1],
      n_rep[2],
      ncv_max,
      LEVEL,
      n_fail
    ),
    "",
    crn_line,
    "",
    "## Truths",
    "",
    md_table(truth_tab),
    "",
    "Signal share = share of int Var_X(int X beta ds) dt; resid. = residual of the L2 projection of the full beta onto the basis's ff tensor space.",
    "",
    "## Coverage of the exact truth, grid average (MC SE)",
    "",
    unlist(lapply(c("beta", "mean", "alpha", "gamma"), \(e) {
      c(paste0("### ", e), "", md_table(cov_table(e)), "")
    })),
    "## Bias-eliminated coverage / 5% pointwise quantile of coverage / mean |bias|/SD",
    "",
    unlist(lapply(
      c("beta", "mean"),
      \(e) c(paste0("### ", e), "", md_table(be_table(e)), "")
    )),
    "## MSE (share of squared bias)",
    "",
    unlist(lapply(
      c("beta", "mean"),
      \(e) c(paste0("### ", e), "", md_table(mse_table(e)), "")
    )),
    sprintf(
      "## beta coverage in the bump region (within %g sd of (%.1f, %.1f), %d grid points) / elsewhere",
      MS_BUMP_REGION_SD,
      MS_BUMP_CENTER[["s"]],
      MS_BUMP_CENTER[["t"]],
      sum(bump_region())
    ),
    "",
    md_table(bump_tab),
    "",
    "## Paired coverage contrasts (common random numbers), REML arms",
    "",
    md_table(con_tab),
    ""
  )
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3)
  stop("usage: analyse-misspec.R results_dir out_dir ncv_xlarge_max_rep")
analyse(args[1], args[2], as.numeric(args[3]))
