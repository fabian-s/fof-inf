# Aggregate plasmode results into coverage, MSE and rule tables (aggregates only: no
# per-subject curves leave the cache). Run from the repository root:
#   Rscript plasmode/analyse-plasmode.R <app> [<app> ...]
# Results are read from PLASMODE_RESULTS/<app> (default PLASMODE_CACHE/results). Output
# goes to PLASMODE_SUMMARIES (default summaries/plasmode, which must be written from all
# applications at once; any other directory accepts a subset, e.g. for smoke runs).
source("scripts/setup.R")
source("analysis/contrasts.R")
source("plasmode/sources.R")
source("plasmode/plasmode.R")
source("plasmode/rule.R")
suppressPackageStartupMessages(library(dplyr))
options(restart.cores = as.integer(Sys.getenv("RESTART_CORES", "2")))

apps <- commandArgs(trailingOnly = TRUE)
cells <- plasmode_cells()
out_dir <- Sys.getenv("PLASMODE_SUMMARIES", "summaries/plasmode")
# The tables below are shared files: always write them from all applications at once.
if (out_dir == "summaries/plasmode") {
  stopifnot(
    "analyse all plasmode applications together" = setequal(apps, cells$app)
  )
}
results <- Sys.getenv("PLASMODE_RESULTS", file.path(PLASMODE_CACHE, "results"))
study <- summarise_study(file.path(results, apps), cells)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
design <- cells[c("cell", "app", "truth", "residual", "flip", "error", "G")]
design_keys <- c("app", "truth", "residual", "flip", "error", "G")

coverage <- study$cells |>
  left_join(design, by = "cell") |>
  filter(!is.na(coverage_mean)) |>
  select(
    all_of(design_keys),
    estimand,
    arm,
    n_rep,
    coverage_mean,
    coverage_se,
    width_mean,
    score_mean,
    mse,
    bias2,
    variance,
    pointwise_q05
  )
readr::write_csv(coverage, file.path(out_dir, "coverage.csv"))

# NCV/REML MSE ratio with a paired replicate-bootstrap interval (same data, both fits).
mse_ratio <- paired_contrast(study$tasks, "ncv_cl2", "reml_cl2", "sq_error") |>
  mutate(
    ratio_lower = vapply(ratio_ci, `[`, 1, 1),
    ratio_upper = vapply(ratio_ci, `[`, 1, 2)
  ) |>
  left_join(design, by = "cell") |>
  select(
    all_of(design_keys),
    estimand,
    n_pairs,
    mse_ratio_ncv_reml = ratio,
    ratio_lower,
    ratio_upper
  )
readr::write_csv(mse_ratio, file.path(out_dir, "mse-ratio.csv"))

# the original cells' records predate the covariance timing
if (!"cov_time" %in% names(study$fits)) study$fits$cov_time <- NA_real_
fits <- study$fits |>
  left_join(design, by = "cell") |>
  group_by(across(all_of(design_keys)), fit) |>
  summarise(
    n = n(),
    converged = mean(converged),
    median_time = median(time),
    median_cov_time = median(cov_time, na.rm = TRUE),
    rho = mean(rho),
    .groups = "drop"
  )
readr::write_csv(fits, file.path(out_dir, "fits.csv"))

# The per-dataset rule analysis on the new cells' apps.
rules <- rule_table(study$tasks, cells, apps)
readr::write_csv(rules, file.path(out_dir, "rule.csv"))
readr::write_csv(
  consistency_table(rules),
  file.path(out_dir, "consistency.csv")
)
readr::write_csv(
  attribution_table(study$tasks, cells, apps),
  file.path(out_dir, "attribution.csv")
)
readr::write_csv(flip_table(rules), file.path(out_dir, "flips.csv"))
readr::write_csv(
  block_arm(study$tasks, cells, apps),
  file.path(out_dir, "block-arm.csv")
)
cat("failures:", nrow(study$failures), "\n")
cat("plasmode summaries written to", out_dir, "\n")
