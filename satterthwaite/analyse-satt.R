# Satterthwaite df study: cell summaries, paired contrasts, df quantiles and SUMMARY.md
# from the task records of run-satt.R. Run from the repository root (code copy):
#   Rscript satterthwaite/analyse-satt.R <results_dir> <out_dir>
# Writes satt-cells.csv, satt-contrasts.csv, satt-df.csv, satt-failures.csv,
# satt-crn.csv and SUMMARY.md to out_dir. No number in SUMMARY.md is typed by hand.

source("scripts/setup.R")
source("ablation/ablation-lib.R")
source("satterthwaite/satt-lib.R")

CRIT_LABELS <- c(
  z = "z",
  t = "t_{G-1}",
  satt = "Satt",
  satt_diag = "Satt diag",
  satt_b2 = "Satt B2"
)

CONTRASTS <- tibble::tribble(
  ~contrast,
  ~crit_a,
  ~crit_b,
  "Satt - z",
  "satt",
  "z",
  "Satt - t",
  "satt",
  "t",
  "Satt diag - Satt",
  "satt_diag",
  "satt",
  "Satt B2 - t",
  "satt_b2",
  "t"
)

ESTIMANDS <- c("beta", "mean", "alpha", "gamma")

# Per-cell scoring ---------------------------------------------------------------------

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

# (points x reps) matrices of every record column.
record_matrices <- function(results, estimand) {
  cols <- colnames(results[[1]]$records[[estimand]])
  m <- lapply(setNames(nm = cols), \(col) {
    do.call(cbind, lapply(results, \(r) r$records[[estimand]][, col]))
  })
  m$reps <- vapply(results, `[[`, 0, "rep")
  m
}

score_crit <- function(m, cell, estimand, crit) {
  linkinv <- FAMILIES[[cell$fit_family]]$linkinv
  cv <- critical_values(crit, cell$G, m)
  s <- interval_scores(m$est_REML, m$se_cl2e_b, m$truth, cv, estimand, linkinv)
  n <- ncol(m$est_REML)
  per_rep <- tibble::tibble(
    rep = m$reps,
    coverage = colMeans(s$covered),
    width = colMeans(s$width)
  )
  n_df_na <- if (crit %in% c("satt", "satt_diag", "satt_b2")) {
    sum(is.na(if (crit == "satt_diag") m$df_diag else m$df_full))
  } else {
    0L
  }
  row <- dplyr::bind_cols(
    cell_keys(cell, estimand),
    tibble::tibble(
      crit = crit,
      n_rep = n,
      coverage_se = sd(per_rep$coverage) / sqrt(n),
      coverage = mean(per_rep$coverage),
      pointwise_q05 = unname(quantile(rowMeans(s$covered), 0.05)),
      width_se = sd(per_rep$width) / sqrt(n),
      width = mean(per_rep$width),
      score = mean(s$score),
      mean_crit = mean(cv),
      n_df_na = n_df_na
    )
  )
  list(row = row, per_rep = dplyr::mutate(per_rep, crit = crit))
}

paired_contrasts <- function(per_rep, keys) {
  purrr::pmap(CONTRASTS, \(contrast, crit_a, crit_b) {
    a <- dplyr::filter(per_rep, crit == crit_a)
    b <- dplyr::filter(per_rep, crit == crit_b)
    d <- dplyr::inner_join(a, b, by = "rep", suffix = c("_a", "_b"))
    diff <- d$coverage_a - d$coverage_b
    wratio <- d$width_a / d$width_b
    dplyr::bind_cols(
      keys,
      tibble::tibble(
        contrast = contrast,
        n_rep = length(diff),
        diff_se = sd(diff) / sqrt(length(diff)),
        diff = mean(diff),
        width_ratio = mean(wratio)
      )
    )
  }) |>
    dplyr::bind_rows()
}

# Distribution of the per-point df (pooled over points and replicates).
df_quantiles <- function(m, keys) {
  q <- \(x, p) unname(quantile(x, p, na.rm = TRUE))
  ratio_b2 <- (m$se_cl2e_b / m$se_cl2e_f)^4
  dplyr::bind_cols(
    keys,
    tibble::tibble(
      df_q05 = q(m$df_full, 0.05),
      df_median = q(m$df_full, 0.5),
      df_q95 = q(m$df_full, 0.95),
      df_median_rel = q(m$df_full, 0.5) / (keys$G - 1),
      df_diag_median = q(m$df_diag, 0.5),
      df_b2_median = q(m$df_full * ratio_b2, 0.5),
      # spread over replicates of each point's df, median over points
      df_rep_sd = median(apply(m$df_full, 1, sd, na.rm = TRUE)),
      n_na = sum(is.na(m$df_full))
    )
  )
}

summarise_cell <- function(files, cell) {
  all <- lapply(files, readRDS)
  ok <- vapply(all, \(r) identical(r$status, "success"), TRUE)
  results <- all[ok]
  failures <- dplyr::bind_rows(lapply(all[!ok], \(r) {
    dplyr::mutate(r$failures, cell = r$cell, rep = r$rep)
  }))
  crn <- dplyr::bind_rows(lapply(results, \(r) {
    if (is.null(r$crn)) return(NULL)
    dplyr::mutate(r$crn, cell = r$cell, rep = r$rep)
  }))
  rows <- contrasts <- dfq <- list()
  if (length(results) == 0) model_estimands <- \(model) character()
  for (estimand in model_estimands(cell$model)) {
    m <- record_matrices(results, estimand)
    keys <- cell_keys(cell, estimand)
    scored <- lapply(CRITS, \(crit) score_crit(m, cell, estimand, crit))
    rows[[estimand]] <- dplyr::bind_rows(lapply(scored, `[[`, "row"))
    per_rep <- dplyr::bind_rows(lapply(scored, `[[`, "per_rep"))
    contrasts[[estimand]] <- paired_contrasts(per_rep, keys)
    dfq[[estimand]] <- df_quantiles(m, keys)
  }
  list(
    cells = dplyr::bind_rows(rows),
    contrasts = dplyr::bind_rows(contrasts),
    df = dplyr::bind_rows(dfq),
    failures = failures,
    crn = crn,
    n_files = length(files),
    times = vapply(all, \(r) r$time %||% NA_real_, 0),
    df_times = vapply(all, \(r) r$df_time %||% NA_real_, 0),
    core_check = vapply(results, \(r) r$core_check, 0),
    G_check = all(vapply(results, \(r) r$G == cell$G, TRUE)),
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
    \(i)
      paste0(
        "| ",
        paste(vapply(df[i, ], \(x) fmt(x[[1]]), ""), collapse = " | "),
        " |"
      ),
    ""
  )
  c(
    paste0("| ", paste(names(df), collapse = " | "), " |"),
    paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|"),
    body,
    ""
  )
}

order_cells <- function(df) {
  dplyr::arrange(
    df,
    factor(family, c("gaussian", "poisson", "binomial")),
    factor(error, c("iid", "ou", "fpc")),
    G
  )
}

# One value per cell and critical value, wide over critical values, plus mean rows by G.
crit_table <- function(cells_df, est, value, scale = 100) {
  wide <- cells_df |>
    dplyr::filter(estimand == est) |>
    dplyr::mutate(
      crit = factor(CRIT_LABELS[crit], CRIT_LABELS),
      v = scale * .data[[value]]
    ) |>
    dplyr::select(family, error, G, crit, v) |>
    tidyr::pivot_wider(names_from = crit, values_from = v) |>
    order_cells()
  avg <- wide |>
    dplyr::group_by(G) |>
    dplyr::summarise(
      dplyr::across(where(is.numeric), mean),
      .groups = "drop"
    ) |>
    dplyr::mutate(family = "mean", error = "")
  dplyr::bind_rows(wide, avg) |> dplyr::mutate(G = as.character(G))
}

contrast_table <- function(contrasts_df, est) {
  contrasts_df |>
    dplyr::filter(estimand == est) |>
    dplyr::mutate(
      value = sprintf(
        "%+.2f (%.2f) ×%.3f",
        100 * diff,
        100 * diff_se,
        width_ratio
      ),
      contrast = factor(contrast, CONTRASTS$contrast)
    ) |>
    dplyr::select(family, error, G, contrast, value) |>
    tidyr::pivot_wider(names_from = contrast, values_from = value) |>
    order_cells() |>
    dplyr::mutate(G = as.character(G))
}

# Mean contrast over cells, by G and overall (MC SE from the cell SEs).
contrast_means <- function(contrasts_df) {
  summ <- \(d)
    dplyr::summarise(
      d,
      diff_m = mean(diff),
      se = sqrt(sum(diff_se^2)) / dplyr::n(),
      wr = mean(width_ratio),
      .groups = "drop"
    )
  by_G <- contrasts_df |>
    dplyr::group_by(estimand, contrast, G) |>
    summ() |>
    dplyr::mutate(G = paste0("G = ", G))
  overall <- contrasts_df |>
    dplyr::group_by(estimand, contrast) |>
    summ() |>
    dplyr::mutate(G = "all")
  dplyr::bind_rows(by_G, overall) |>
    dplyr::mutate(
      value = sprintf("%+.2f (%.2f) ×%.3f", 100 * diff_m, 100 * se, wr),
      contrast = factor(contrast, CONTRASTS$contrast)
    ) |>
    dplyr::select(estimand, contrast, G, value) |>
    tidyr::pivot_wider(names_from = G, values_from = value) |>
    dplyr::arrange(factor(estimand, ESTIMANDS), contrast)
}

df_table <- function(df_df) {
  df_df |>
    order_cells() |>
    dplyr::mutate(
      G = as.character(G),
      `median/(G-1)` = df_median_rel
    ) |>
    dplyr::select(
      estimand,
      family,
      error,
      G,
      q05 = df_q05,
      median = df_median,
      q95 = df_q95,
      `median/(G-1)`,
      `diag median` = df_diag_median,
      `B2 median` = df_b2_median,
      `rep SD` = df_rep_sd
    ) |>
    dplyr::arrange(factor(estimand, ESTIMANDS))
}

df_by_G <- function(df_df) {
  df_df |>
    dplyr::group_by(estimand, G) |>
    dplyr::summarise(
      cells = dplyr::n(),
      `min q05` = min(df_q05),
      `mean median` = mean(df_median),
      `max q95` = max(df_q95),
      `mean median/(G-1)` = mean(df_median_rel),
      .groups = "drop"
    ) |>
    dplyr::arrange(factor(estimand, ESTIMANDS), G) |>
    dplyr::mutate(G = as.character(G))
}

run_lines <- function(out, checks) {
  crn_by <- out$crn |>
    dplyr::group_by(column) |>
    dplyr::summarise(
      max_abs = max(max_abs),
      max_rel = max(max_rel),
      .groups = "drop"
    )
  c(
    "## Run and checks",
    "",
    sprintf(
      "- Task files: %d of %d; failed tasks: %d; replicates per cell × estimand: min %d, max %d.",
      out$n_files,
      18 * N_REP,
      nrow(out$failures),
      min(out$cells$n_rep),
      max(out$cells$n_rep)
    ),
    sprintf(
      "- Software: refund %s, mgcv %s.",
      paste(unique(substr(out$software$refund_sha, 1, 8)), collapse = ", "),
      paste(unique(out$software$mgcv_version), collapse = ", ")
    ),
    sprintf(
      "- Worker time per task: median %.1f s (df part median %.1f s), total %.1f core-hours.",
      median(out$times, na.rm = TRUE),
      median(out$df_times, na.rm = TRUE),
      sum(out$times, na.rm = TRUE) / 3600
    ),
    sprintf(
      "- df context = covariance: max relative difference between the frequentist CL2 covariance rebuilt from the df context's influence core and `pffr_vcov(sandwich = \"cl2\", cl2_adjustment = \"exact\", freq = TRUE)`: %.1e; df context G equals the cell's G in every replicate: %s.",
      max(out$core_check),
      all(out$G_check)
    ),
    sprintf(
      "- Undefined (NA) per-point df: %d (scored with z).",
      sum(out$df$n_na)
    ),
    sprintf(
      "- CRN reproduction against the study's stored REML records: %d task records compared.",
      dplyr::n_distinct(paste(out$crn$cell, out$crn$rep))
    ),
    "",
    md_table(dplyr::mutate(
      crn_by,
      max_abs = formatC(max_abs, format = "e", digits = 1),
      max_rel = formatC(max_rel, format = "e", digits = 1)
    )),
    checks,
    ""
  )
}

# Reproduction of the z (paper) and t arms against the study and the ablation.
reproduction_checks <- function(cells_df, study_cells, ablation_cells) {
  z_study <- cells_df |>
    dplyr::filter(crit == "z") |>
    dplyr::inner_join(
      dplyr::filter(study_cells, arm == "reml_cl2") |>
        dplyr::select(cell, estimand, ref = coverage_mean),
      by = c("cell", "estimand")
    )
  abl <- cells_df |>
    dplyr::filter(crit %in% c("z", "t")) |>
    dplyr::inner_join(
      dplyr::filter(ablation_cells, arm == "cl2e_b") |>
        dplyr::select(cell, estimand, crit, ref = coverage),
      by = c("cell", "estimand", "crit")
    )
  c(
    sprintf(
      "- z coverage vs `summaries/final/cells.csv` (reml_cl2): max |difference| %.2e over %d cell × estimand rows.",
      max(abs(z_study$coverage - z_study$ref)),
      nrow(z_study)
    ),
    sprintf(
      "- z and t_{G-1} coverage vs `ablation/summaries/ablation-cells.csv` (cl2e_b): max |difference| %.2e over %d rows.",
      max(abs(abl$coverage - abl$ref)),
      nrow(abl)
    )
  )
}

write_summary <- function(out, checks, path) {
  lines <- c(
    "# Satterthwaite / Bell-McCaffrey df for CL2 intervals (core cells, mid signal, D = 61, M1)",
    "",
    "Generated by `satterthwaite/analyse-satt.R`; all numbers computed from the task records.",
    "REML fit, CL2 exact Bayesian SE (as in the paper), nominal 95% pointwise intervals,",
    "grid-averaged, R = 200 replicates (common random numbers with the study).",
    "Critical values: z = z_0.975; t_{G-1}; Satt = per-point df of refund 79a346fb",
    "`coef.pffr(crit = \"satterthwaite\", df_gram = \"full\")` on the study's contrasts;",
    "Satt diag = `df_gram = \"diagonal\"`; Satt B2 = Satt df × (se_Bayes / se_freq)^4",
    "(the df of the Bayesian variance, which adds the fixed Vp - Ve to the sampled sandwich;",
    "refund's df describes the sandwich alone). q05 = 5% quantile over grid points of the",
    "pointwise coverage. Contrasts: paired coverage difference in percentage points (MC SE)",
    "and mean width ratio (×).",
    "",
    run_lines(out, checks)
  )
  for (est in ESTIMANDS) {
    lines <- c(
      lines,
      sprintf("## %s", est),
      "",
      "Coverage (%):",
      "",
      md_table(crit_table(out$cells, est, "coverage")),
      "5% pointwise coverage quantile (%):",
      "",
      md_table(crit_table(out$cells, est, "pointwise_q05")),
      "Mean width:",
      "",
      md_table(crit_table(out$cells, est, "width", scale = 1), digits = 3),
      "Paired contrasts, percentage points (MC SE) × width ratio:",
      "",
      md_table(contrast_table(out$contrasts, est))
    )
  }
  lines <- c(
    lines,
    "## Mean contrasts over cells (percentage points, MC SE, × width ratio)",
    "",
    md_table(contrast_means(out$contrasts)),
    "## Per-point df (full Gram), pooled over points and replicates",
    "",
    "By G (over the 9 cells of each G):",
    "",
    md_table(df_by_G(out$df), digits = 2),
    "Per cell (rep SD = median over points of the SD over replicates of the df):",
    "",
    md_table(df_table(out$df), digits = 1)
  )
  writeLines(lines, path)
}

# Main ---------------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: analyse-satt.R results_dir out_dir")
results_dir <- args[1]
out_dir <- args[2]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cells <- ablation_cells()
files <- list.files(results_dir, "^task-c.*\\.rds$", full.names = TRUE)
file_cell <- as.integer(sub("^task-c(\\d+)-r.*", "\\1", basename(files)))
per_cell <- lapply(seq_len(nrow(cells)), \(i) {
  summarise_cell(files[file_cell == cells$cell[i]], cells[i, ])
})
parts <- c("cells", "contrasts", "df", "failures", "crn", "software")
out <- lapply(
  setNames(nm = parts),
  \(p) dplyr::bind_rows(lapply(per_cell, `[[`, p))
)
out$software <- dplyr::distinct(out$software)
for (v in c("n_files")) out[[v]] <- sum(vapply(per_cell, `[[`, 0, v))
for (v in c("times", "df_times", "core_check", "G_check")) {
  out[[v]] <- unlist(lapply(per_cell, `[[`, v))
}
if (nrow(out$failures) == 0) {
  out$failures <- tibble::tibble(
    stage = character(),
    message = character(),
    cell = integer(),
    rep = integer()
  )
}
write.csv(out$cells, file.path(out_dir, "satt-cells.csv"), row.names = FALSE)
write.csv(
  out$contrasts,
  file.path(out_dir, "satt-contrasts.csv"),
  row.names = FALSE
)
write.csv(out$df, file.path(out_dir, "satt-df.csv"), row.names = FALSE)
write.csv(
  out$failures,
  file.path(out_dir, "satt-failures.csv"),
  row.names = FALSE
)
write.csv(out$crn, file.path(out_dir, "satt-crn.csv"), row.names = FALSE)
checks <- reproduction_checks(
  out$cells,
  read.csv(Sys.getenv("STUDY_CELLS", "summaries/final/cells.csv")),
  read.csv(Sys.getenv(
    "ABLATION_CELLS",
    "ablation/summaries/ablation-cells.csv"
  ))
)
write_summary(out, checks, file.path(out_dir, "SUMMARY.md"))
cat("done:", nrow(out$cells), "rows\n")
