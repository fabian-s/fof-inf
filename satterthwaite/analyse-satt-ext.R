# Satterthwaite df on the extension, small-G, null-effect and plasmode cells (2026-10-03):
# cell summaries, paired contrasts, df quantiles and a z-arm reproduction check against the
# study's summaries. Run from the repository root (on LRZ for the plasmode records):
#   Rscript satterthwaite/analyse-satt-ext.R <results_dir> <syn|plasmode> <out_dir>
# Writes satt-<mode>-cells.csv, -contrasts.csv, -df.csv, -failures.csv, -check.csv.
# Scoring, contrasts and df quantiles are those of analyse-satt.R (sourced up to its main).

source("scripts/setup.R")
source("ablation/ablation-lib.R")
source("satterthwaite/satt-lib.R")
src <- readLines("satterthwaite/analyse-satt.R")
eval(parse(text = src[seq_len(grep("^# Main -", src) - 1)]))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("usage: analyse-satt-ext.R results_dir syn|plasmode out_dir")
}
results_dir <- args[1]
mode <- args[2]
out_dir <- args[3]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (mode == "syn") {
  cells <- design_table()
  key_cols <- c("block", "signal", "truth", "D", "model", "covariate", "basis")
  study <- read.csv("summaries/final/cells.csv")
} else {
  source("plasmode/sources.R")
  source("plasmode/plasmode.R")
  cells <- plasmode_cells()
  key_cols <- c("app", "truth", "residual", "flip")
  study <- read.csv("summaries/plasmode/coverage.csv")
}
files <- list.files(results_dir, "^task-c.*\\.rds$", full.names = TRUE)
file_cell <- as.integer(sub("^task-c(\\d+)-r.*", "\\1", basename(files)))
cells <- cells[cells$cell %in% file_cell, ]
per_cell <- lapply(seq_len(nrow(cells)), \(i) {
  summarise_cell(files[file_cell == cells$cell[i]], cells[i, ])
})
add_keys <- \(d) {
  if (!nrow(d)) return(d)
  dplyr::left_join(d, cells[, c("cell", key_cols)], by = "cell")
}
parts <- c("cells", "contrasts", "df", "failures", "crn")
out <- lapply(setNames(nm = parts), \(p) {
  dplyr::bind_rows(lapply(per_cell, `[[`, p))
})
for (p in c("cells", "contrasts", "df")) out[[p]] <- add_keys(out[[p]])

# z arm vs the study's REML + CL2 (same replicates and seeds): coverage difference per
# cell x estimand where the study has the cell.
z <- dplyr::filter(out$cells, crit == "z")
if (mode == "syn") {
  ref <- dplyr::filter(study, arm == "reml_cl2") |>
    dplyr::transmute(
      cell,
      estimand,
      n_rep_study = n_rep,
      coverage_study = coverage_mean
    )
  check <- dplyr::inner_join(z, ref, by = c("cell", "estimand"))
} else {
  ref <- dplyr::filter(study, arm == "reml_cl2") |>
    dplyr::transmute(
      app,
      truth,
      residual,
      flip,
      error,
      G,
      estimand,
      n_rep_study = n_rep,
      coverage_study = coverage_mean
    )
  check <- dplyr::inner_join(
    z,
    ref,
    by = c("app", "truth", "residual", "flip", "error", "G", "estimand")
  )
}
check <- dplyr::transmute(
  check,
  cell,
  estimand,
  n_rep,
  n_rep_study,
  coverage,
  coverage_study,
  diff = coverage - coverage_study
)

info <- tibble::tibble(
  n_files = sum(vapply(per_cell, `[[`, 0, "n_files")),
  n_failed = nrow(out$failures),
  time_median = median(unlist(lapply(per_cell, `[[`, "times")), na.rm = TRUE),
  time_total_h = sum(unlist(lapply(per_cell, `[[`, "times")), na.rm = TRUE) /
    3600,
  df_time_median = median(
    unlist(lapply(per_cell, `[[`, "df_times")),
    na.rm = TRUE
  ),
  core_check_max = max(unlist(lapply(per_cell, `[[`, "core_check"))),
  G_check = all(unlist(lapply(per_cell, `[[`, "G_check"))),
  n_df_na = sum(out$df$n_na),
  crn_tasks = if (nrow(out$crn)) nrow(dplyr::distinct(out$crn, cell, rep)) else
    0L,
  crn_max_rel = if (nrow(out$crn)) max(out$crn$max_rel) else NA_real_,
  check_rows = nrow(check),
  check_max_abs_diff = if (nrow(check)) max(abs(check$diff)) else NA_real_,
  software = paste(
    unique(unlist(lapply(per_cell, \(p) {
      paste(p$software$refund_sha, p$software$mgcv_version)
    }))),
    collapse = "; "
  )
)
w <- \(d, name) {
  write.csv(
    d,
    file.path(out_dir, sprintf("satt-%s-%s.csv", mode, name)),
    row.names = FALSE
  )
}
w(out$cells, "cells")
w(out$contrasts, "contrasts")
w(out$df, "df")
w(out$failures, "failures")
w(check, "check")
if (nrow(out$crn)) {
  w(
    dplyr::summarise(
      dplyr::group_by(out$crn, column),
      tasks = dplyr::n_distinct(paste(cell, rep)),
      max_abs = max(max_abs),
      max_rel = max(max_rel)
    ),
    "crn"
  )
}
w(info, "info")
print(as.data.frame(info))
