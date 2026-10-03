# Relative estimation error and interval informativeness per replicate, for the main
# interval arms. Aggregates only. Run from the repository root:
#   Rscript analysis/relative-error.R synthetic <out.csv>
#   Rscript analysis/relative-error.R plasmode <out.csv> <app> [<app> ...]
#
# Per replicate, estimand and arm (grid means over the evaluation points):
# - rel_rmse: sqrt(mean squared error / mean squared truth) of the arm's point estimate.
#   For the mean and the functional intercept the denominator is the truth's variation
#   around its grid mean (their level is not the effect of interest).
# - rel_halfwidth: sqrt(mean squared half-width / the same denominator).
# - detect (beta, gamma, f only): share of the points with |truth| > 25% of max |truth|
#   where the interval excludes zero on the side of the truth.
# - wrong_sign (beta, gamma, f only): share of all points where the interval excludes zero
#   on the side opposite to the truth (points with truth exactly 0 count as wrong).
source("scripts/setup.R")
suppressPackageStartupMessages(library(dplyr))
args <- commandArgs(trailingOnly = TRUE)
mode <- args[1]
out <- args[2]
cores <- as.integer(Sys.getenv("RESTART_CORES", "1"))

if (mode == "synthetic") {
  cells <- design_table()
  # RELERR_DIRS (colon-separated) restricts the synthetic run to some results directories
  dirs <- if (nzchar(Sys.getenv("RELERR_DIRS"))) {
    strsplit(Sys.getenv("RELERR_DIRS"), ":")[[1]]
  } else {
    list.dirs("results", recursive = FALSE)
  }
} else {
  source("plasmode/sources.R")
  source("plasmode/plasmode.R")
  cells <- plasmode_cells()
  results <- Sys.getenv(
    "PLASMODE_RESULTS",
    file.path(PLASMODE_CACHE, "results")
  )
  dirs <- file.path(results, args[-(1:2)])
}
files <- unlist(lapply(
  dirs,
  list.files,
  pattern = "^task-.*\\.rds$",
  full.names = TRUE
))
# smoke tests: RELERR_MAX_FILES limits the number of task files read
max_files <- as.integer(Sys.getenv("RELERR_MAX_FILES", "0"))
if (max_files > 0) files <- sample(files, max_files)
MAIN <- c("reml_model", "reml_cl2", "ncv_cl2", "ncv_cl2_bias", "ar1_model")
CENTRED <- c("mean", "alpha")
DETECT_THRESHOLD <- 0.25

# interval bounds as in point_scores() (R/summarise.R)
bounds <- function(record, arm, estimand, linkinv) {
  spec <- ARMS[ARMS$arm == arm, ]
  est <- record[, spec$est]
  se <- record[, spec$se]
  if (spec$bias_aware) se <- sqrt(se^2 + (est - record[, "est_REML"])^2)
  z <- qnorm(1 - (1 - LEVEL) / 2)
  b <- list(est = est, lower = est - z * se, upper = est + z * se)
  if (estimand == "mean") b <- lapply(b, linkinv)
  b
}

task_rows <- function(f) {
  r <- readRDS(f)
  if (!identical(r$status, "success")) return(NULL)
  cell <- cells[cells$cell == r$cell, ]
  stopifnot(nrow(cell) == 1)
  linkinv <- FAMILIES[[cell$fit_family]]$linkinv
  rows <- list()
  for (e in scored_estimands(cell, names(r$records))) {
    rec <- r$records[[e]]
    tr <- rec[, "truth"]
    denom <- if (e %in% CENTRED) mean((tr - mean(tr))^2) else mean(tr^2)
    big <- abs(tr) > DETECT_THRESHOLD * max(abs(tr))
    for (arm in MAIN) {
      spec <- ARMS[ARMS$arm == arm, ]
      if (!all(c(spec$est, spec$se) %in% colnames(rec))) next
      b <- bounds(rec, arm, e, linkinv)
      excl_pos <- b$lower > 0
      excl_neg <- b$upper < 0
      sided <- !(e %in% CENTRED)
      rows[[length(rows) + 1]] <- tibble(
        cell = r$cell,
        rep = r$rep,
        estimand = e,
        arm = arm,
        denom = denom,
        sq_error = mean((b$est - tr)^2),
        coverage = mean(b$lower <= tr & tr <= b$upper),
        sq_halfwidth = mean(((b$upper - b$lower) / 2)^2),
        detect = if (sided)
          mean((excl_pos & tr > 0 | excl_neg & tr < 0)[big]) else NA_real_,
        wrong_sign = if (sided)
          mean(excl_pos & tr <= 0 | excl_neg & tr >= 0) else NA_real_
      )
    }
  }
  bind_rows(rows)
}

per_rep <- bind_rows(parallel::mclapply(files, mc.cores = cores, task_rows))
stopifnot(nrow(per_rep) > 0)

rel <- per_rep |>
  mutate(
    rel_rmse = sqrt(sq_error / denom),
    rel_halfwidth = sqrt(sq_halfwidth / denom)
  ) |>
  group_by(cell, estimand, arm) |>
  summarise(
    n_rep = n(),
    denom = mean(denom),
    mse = mean(sq_error),
    coverage = mean(coverage),
    rel_rmse_cell = sqrt(mean(sq_error) / mean(denom)),
    rel_rmse_median = median(rel_rmse),
    rel_rmse_q90 = unname(quantile(rel_rmse, 0.9)),
    share_rel_gt_0.5 = mean(rel_rmse > 0.5),
    share_rel_gt_1 = mean(rel_rmse > 1),
    rel_halfwidth_median = median(rel_halfwidth),
    detect = mean(detect),
    wrong_sign = mean(wrong_sign),
    .groups = "drop"
  )
keys <- intersect(
  c("cell", "app", "truth", "residual", "flip", "error", "G"),
  names(cells)
)
if (mode == "plasmode") rel <- left_join(rel, cells[keys], by = "cell")
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
readr::write_csv(rel, out)
cat(sprintf(
  "%d task files, %d cells, %d rows -> %s\n",
  length(files),
  n_distinct(rel$cell),
  nrow(rel),
  out
))
