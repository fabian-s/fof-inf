# Covariance ablation: cell summaries, paired contrasts and SUMMARY.md from the task
# records of run-ablation.R. Run from the repository root (code copy):
#   Rscript ablation/analyse-ablation.R <results_dir> <out_dir>
# Writes ablation-cells.csv, ablation-contrasts.csv, ablation-failures.csv,
# ablation-crn.csv and SUMMARY.md to out_dir. No number in SUMMARY.md is typed by hand.

source("scripts/setup.R")
source("ablation/ablation-lib.R")

CRITS <- c("z", "t")

# Per-cell scoring ---------------------------------------------------------------------

# Critical value of a cell: z_{0.975} or t_{G-1, 0.975}.
critical_value <- function(crit, G) {
  q <- 1 - (1 - LEVEL) / 2
  if (crit == "z") qnorm(q) else qt(q, G - 1)
}

# Identifying columns of a cell x estimand (locals first: tibble() evaluates
# sequentially, so a new `cell` column would shadow the cell row).
cell_keys <- function(cell, estimand) {
  id <- cell$cell
  fam <- cell$family
  err <- cell$error
  g <- cell$G
  tibble::tibble(
    cell = id,
    family = fam,
    error = err,
    G = g,
    estimand = estimand
  )
}

# Records of the replicates in which an arm is available, as (points x reps) matrices.
arm_matrices <- function(results, estimand, spec) {
  est_col <- paste0("est_", spec$fit)
  se_col <- paste0("se_", spec$arm)
  ok <- vapply(
    results,
    \(r) all(c(est_col, se_col) %in% colnames(r$records[[estimand]])),
    TRUE
  )
  pick <- \(col)
    do.call(cbind, lapply(results[ok], \(r) r$records[[estimand]][, col]))
  list(
    reps = vapply(results[ok], `[[`, 0, "rep"),
    est = pick(est_col),
    se = pick(se_col),
    truth = pick("truth")
  )
}

# One cell x estimand x arm x critical value: summary row plus per-replicate scores.
score_arm <- function(m, cell, estimand, spec, crit, n_total) {
  linkinv <- FAMILIES[[cell$fit_family]]$linkinv
  s <- interval_scores(
    m$est,
    m$se,
    m$truth,
    critical_value(crit, cell$G),
    estimand,
    linkinv
  )
  n <- ncol(m$est)
  per_rep <- tibble::tibble(
    rep = m$reps,
    coverage = colMeans(s$covered),
    width = colMeans(s$width),
    score = colMeans(s$score)
  )
  # SE/SD ratio on the link scale, per point, then averaged over the grid
  se_sd <- sqrt(rowMeans(m$se^2)) / apply(m$est, 1, sd)
  row <- dplyr::bind_cols(
    cell_keys(cell, estimand),
    tibble::tibble(
      arm = spec$arm,
      fit = spec$fit,
      crit = crit,
      n_rep = n,
      n_missing = n_total - n,
      coverage_se = sd(per_rep$coverage) / sqrt(n),
      coverage = mean(per_rep$coverage),
      pointwise_q05 = unname(quantile(rowMeans(s$covered), 0.05)),
      width = mean(per_rep$width),
      score = mean(per_rep$score),
      se_sd_ratio = mean(se_sd)
    )
  )
  list(row = row, per_rep = dplyr::mutate(per_rep, arm = spec$arm, crit = crit))
}

# Paired coverage contrasts (per-replicate differences, so the MC SE is paired).
CONTRASTS <- tibble::tribble(
  ~contrast,
  ~arm_a,
  ~crit_a,
  ~arm_b,
  ~crit_b,
  "exact - shortcut",
  "cl2e_b",
  "z",
  "cl2s_b",
  "z",
  "CL2 - CR1",
  "cl2e_b",
  "z",
  "cr1_b",
  "z",
  "Bayes - freq (CL2)",
  "cl2e_b",
  "z",
  "cl2e_f",
  "z",
  "Bayes - freq (CR1)",
  "cr1_b",
  "z",
  "cr1_f",
  "z",
  "t - z (CL2)",
  "cl2e_b",
  "t",
  "cl2e_b",
  "z",
  "fixed - reselected (model)",
  "fix_model",
  "z",
  "model",
  "z",
  "fixed - reselected (CL2)",
  "fix_cl2e_b",
  "z",
  "cl2e_b",
  "z"
)

paired_contrasts <- function(per_rep, keys) {
  purrr::pmap(CONTRASTS, \(contrast, arm_a, crit_a, arm_b, crit_b) {
    a <- dplyr::filter(per_rep, arm == arm_a, crit == crit_a)
    b <- dplyr::filter(per_rep, arm == arm_b, crit == crit_b)
    d <- dplyr::inner_join(a, b, by = "rep", suffix = c("_a", "_b"))
    diff <- d$coverage_a - d$coverage_b
    dplyr::bind_cols(
      keys,
      tibble::tibble(
        contrast = contrast,
        n_rep = length(diff),
        diff_se = sd(diff) / sqrt(length(diff)),
        diff = mean(diff)
      )
    )
  }) |>
    dplyr::bind_rows()
}

summarise_cell <- function(files, cell) {
  all <- lapply(files, readRDS)
  ok <- vapply(all, \(r) identical(r$status, "success"), TRUE)
  results <- all[ok]
  failures <- dplyr::bind_rows(lapply(all, \(r) {
    if (is.null(r$failures)) return(NULL)
    dplyr::mutate(r$failures, cell = r$cell, rep = r$rep, status = r$status)
  }))
  crn <- dplyr::bind_rows(lapply(results, \(r) {
    if (is.null(r$crn)) return(NULL)
    dplyr::mutate(r$crn, cell = r$cell, rep = r$rep)
  }))
  rows <- list()
  contrasts <- list()
  for (estimand in model_estimands(cell$model)) {
    per_rep <- list()
    for (i in seq_len(nrow(COV_ARMS))) {
      spec <- COV_ARMS[i, ]
      m <- arm_matrices(results, estimand, spec)
      if (length(m$reps) == 0) next
      for (crit in CRITS) {
        out <- score_arm(m, cell, estimand, spec, crit, N_REP)
        rows[[length(rows) + 1]] <- out$row
        per_rep[[length(per_rep) + 1]] <- out$per_rep
      }
    }
    if (length(per_rep) == 0) next
    keys <- cell_keys(cell, estimand)
    contrasts[[estimand]] <- paired_contrasts(dplyr::bind_rows(per_rep), keys)
  }
  list(
    cells = dplyr::bind_rows(rows),
    contrasts = dplyr::bind_rows(contrasts),
    failures = failures,
    crn = crn,
    n_files = length(files),
    times = vapply(all, \(r) r$time %||% NA_real_, 0),
    software = dplyr::distinct(dplyr::bind_rows(lapply(all, \(r) {
      tibble::as_tibble(r$software)
    })))
  )
}

# Markdown --------------------------------------------------------------------------------

md_table <- function(df, digits = 1) {
  fmt <- \(x) {
    if (is.double(x)) formatC(x, format = "f", digits = digits) else
      as.character(x)
  }
  body <- vapply(
    seq_len(nrow(df)),
    \(i) {
      paste0(
        "| ",
        paste(vapply(df[i, ], \(x) fmt(x[[1]]), ""), collapse = " | "),
        " |"
      )
    },
    ""
  )
  c(
    paste0("| ", paste(names(df), collapse = " | "), " |"),
    paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|"),
    body,
    ""
  )
}

ARM_LABELS <- c(
  model = "model",
  cr1_b = "CR1 B",
  cr1_f = "CR1 F",
  cl2s_b = "CL2s B",
  cl2s_f = "CL2s F",
  cl2e_b = "CL2e B",
  cl2e_f = "CL2e F",
  fix_model = "fixλ model",
  fix_cl2e_b = "fixλ CL2e B"
)

# Coverage (%) per cell for the given arms and critical value, plus a mean row.
coverage_table <- function(cells_df, est, arms, crit_value) {
  wide <- cells_df |>
    dplyr::filter(estimand == est, crit == crit_value, arm %in% arms) |>
    dplyr::mutate(
      arm = factor(ARM_LABELS[arm], ARM_LABELS[arms]),
      cov = 100 * coverage
    ) |>
    dplyr::select(family, error, G, arm, cov) |>
    tidyr::pivot_wider(names_from = arm, values_from = cov) |>
    dplyr::arrange(
      factor(family, c("gaussian", "poisson", "binomial")),
      factor(error, c("iid", "ou", "fpc")),
      G
    )
  avg <- wide |>
    dplyr::summarise(dplyr::across(where(is.numeric) & !G, mean)) |>
    dplyr::mutate(family = "mean", error = "", G = NA)
  wide <- dplyr::bind_rows(wide, avg)
  wide$G <- ifelse(is.na(wide$G), "", as.character(wide$G))
  wide
}

contrast_table <- function(contrasts_df, est) {
  contrasts_df |>
    dplyr::filter(estimand == est) |>
    dplyr::mutate(
      value = sprintf("%+.1f (%.1f)", 100 * diff, 100 * diff_se),
      contrast = factor(contrast, CONTRASTS$contrast)
    ) |>
    dplyr::select(family, error, G, contrast, value) |>
    tidyr::pivot_wider(names_from = contrast, values_from = value) |>
    dplyr::arrange(
      factor(family, c("gaussian", "poisson", "binomial")),
      factor(error, c("iid", "ou", "fpc")),
      G
    ) |>
    dplyr::mutate(G = as.character(G))
}

# Mean contrast over cells, by family and overall (MC SE from the cell SEs).
contrast_means <- function(contrasts_df) {
  by_family <- contrasts_df |>
    dplyr::group_by(estimand, contrast, family) |>
    dplyr::summarise(
      diff = mean(diff),
      se = sqrt(sum(diff_se^2)) / dplyr::n(),
      .groups = "drop"
    )
  overall <- contrasts_df |>
    dplyr::group_by(estimand, contrast) |>
    dplyr::summarise(
      diff = mean(diff),
      se = sqrt(sum(diff_se^2)) / dplyr::n(),
      .groups = "drop"
    ) |>
    dplyr::mutate(family = "all")
  dplyr::bind_rows(by_family, overall) |>
    dplyr::mutate(
      value = sprintf("%+.1f (%.1f)", 100 * diff, 100 * se),
      family = factor(family, c("gaussian", "poisson", "binomial", "all")),
      contrast = factor(contrast, CONTRASTS$contrast)
    ) |>
    dplyr::select(estimand, contrast, family, value) |>
    tidyr::pivot_wider(names_from = family, values_from = value) |>
    dplyr::arrange(
      factor(estimand, c("beta", "mean", "alpha", "gamma")),
      contrast
    )
}

write_summary <- function(out, study_cells, path) {
  cells_df <- out$cells
  cov_arms <- c(
    "model",
    "cr1_b",
    "cr1_f",
    "cl2s_b",
    "cl2s_f",
    "cl2e_b",
    "cl2e_f"
  )
  # reproduction check of the paper's arm: CL2e B with z vs summaries/final/cells.csv
  check <- cells_df |>
    dplyr::filter(arm == "cl2e_b", crit == "z") |>
    dplyr::inner_join(
      dplyr::filter(study_cells, arm == "reml_cl2") |>
        dplyr::select(
          cell,
          estimand,
          study_coverage = coverage_mean,
          study_n = n_rep
        ),
      by = c("cell", "estimand")
    )
  crn <- out$crn
  crn_by <- crn |>
    dplyr::group_by(column) |>
    dplyr::summarise(
      max_abs = max(max_abs),
      max_rel = max(max_rel),
      .groups = "drop"
    )
  fixed <- cells_df |>
    dplyr::filter(
      crit == "z",
      arm %in% c("model", "fix_model", "cl2e_b", "fix_cl2e_b")
    )
  lines <- c(
    "# Covariance ablation (core cells, mid signal, D = 61, M1)",
    "",
    "Generated by `ablation/analyse-ablation.R`; all numbers computed from the task records.",
    "Coverage in % of nominal 95% pointwise intervals, grid-averaged, R = 200 replicates.",
    "Arms: model = Vc; CR1 = refund `sandwich = \"cluster\"` (dof_correction \"none\");",
    "CL2s / CL2e = CL2 with shortcut (diagonal) / exact (full-block) leverage adjustment;",
    "B = Bayesian (+ Vp - Ve), F = frequentist; fixλ = refit with smoothing parameters fixed at",
    "the cell's geometric-mean REML values. z = z_0.975, t = t_{G-1, 0.975}.",
    "",
    "## Run and checks",
    "",
    sprintf(
      "- Task files: %d of %d; failed tasks: %d; failure rows (incl. fixed-λ refits): %d.",
      out$n_files,
      18 * N_REP,
      sum(out$failures$status == "failure"),
      nrow(out$failures)
    ),
    sprintf(
      "- Replicates per arm: min %d, max %d.",
      min(cells_df$n_rep),
      max(cells_df$n_rep)
    ),
    sprintf(
      "- Software: refund %s, mgcv %s.",
      paste(unique(substr(out$software$refund_sha, 1, 8)), collapse = ", "),
      paste(unique(out$software$mgcv_version), collapse = ", ")
    ),
    sprintf(
      "- Worker time per task: median %.1f s, total %.1f core-hours.",
      median(out$times, na.rm = TRUE),
      sum(out$times, na.rm = TRUE) / 3600
    ),
    sprintf(
      "- CRN reproduction against the study's stored REML records: %d task records compared.",
      dplyr::n_distinct(paste(crn$cell, crn$rep))
    ),
    "",
    md_table(dplyr::mutate(
      crn_by,
      max_abs = formatC(max_abs, format = "e", digits = 1),
      max_rel = formatC(max_rel, format = "e", digits = 1)
    )),
    sprintf(
      "- REML + CL2e B (z) coverage vs `summaries/final/cells.csv` (reml_cl2): max |difference| %.2e over %d cell × estimand rows.",
      max(abs(check$coverage - check$study_coverage)),
      nrow(check)
    ),
    ""
  )
  for (est in c("beta", "mean", "alpha", "gamma")) {
    lines <- c(
      lines,
      sprintf("## %s", est),
      "",
      "Coverage (%), z critical value:",
      "",
      md_table(coverage_table(cells_df, est, cov_arms, "z")),
      "Coverage (%), t_{G-1} critical value:",
      "",
      md_table(coverage_table(cells_df, est, cov_arms, "t")),
      "Fixed vs reselected λ, coverage (%), z:",
      "",
      md_table(coverage_table(
        cells_df,
        est,
        c("model", "fix_model", "cl2e_b", "fix_cl2e_b"),
        "z"
      )),
      "SE/SD ratio (grid-averaged sqrt(mean SE²)/SD of the estimate, link scale):",
      "",
      md_table(
        coverage_table(
          dplyr::mutate(cells_df, coverage = se_sd_ratio / 100),
          est,
          c(cov_arms, "fix_model", "fix_cl2e_b"),
          "z"
        ),
        digits = 2
      ),
      "Paired coverage contrasts, percentage points (MC SE):",
      "",
      md_table(contrast_table(out$contrasts, est))
    )
  }
  lines <- c(
    lines,
    "## Mean contrasts over cells (percentage points, MC SE)",
    "",
    md_table(contrast_means(out$contrasts))
  )
  writeLines(lines, path)
}

# Main ---------------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: analyse-ablation.R results_dir out_dir")
results_dir <- args[1]
out_dir <- args[2]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cells <- ablation_cells()
files <- list.files(results_dir, "^task-c.*\\.rds$", full.names = TRUE)
file_cell <- as.integer(sub("^task-c(\\d+)-r.*", "\\1", basename(files)))
per_cell <- lapply(seq_len(nrow(cells)), \(i) {
  summarise_cell(files[file_cell == cells$cell[i]], cells[i, ])
})
parts <- c("cells", "contrasts", "failures", "crn", "software")
out <- lapply(
  setNames(nm = parts),
  \(p) dplyr::bind_rows(lapply(per_cell, `[[`, p))
)
out$software <- dplyr::distinct(out$software)
out$n_files <- sum(vapply(per_cell, `[[`, 0, "n_files"))
out$times <- unlist(lapply(per_cell, `[[`, "times"))
if (nrow(out$failures) == 0) {
  out$failures <- tibble::tibble(
    stage = character(),
    message = character(),
    cell = integer(),
    rep = integer(),
    status = character()
  )
}
write.csv(
  out$cells,
  file.path(out_dir, "ablation-cells.csv"),
  row.names = FALSE
)
write.csv(
  out$contrasts,
  file.path(out_dir, "ablation-contrasts.csv"),
  row.names = FALSE
)
write.csv(
  out$failures,
  file.path(out_dir, "ablation-failures.csv"),
  row.names = FALSE
)
write.csv(out$crn, file.path(out_dir, "ablation-crn.csv"), row.names = FALSE)
study_cells <- read.csv(Sys.getenv("STUDY_CELLS", "summaries/final/cells.csv"))
write_summary(out, study_cells, file.path(out_dir, "SUMMARY.md"))
cat("done:", nrow(out$cells), "rows\n")
