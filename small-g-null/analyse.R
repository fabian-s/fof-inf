# Core factorial at G = 20 and the null-effect cells (beta = 0), 2026-10-03. Aggregates
# only. Run from the repository root:
#   Rscript small-g-null/analyse.R <results_dir> <satt_results_dir> <out_dir>
# results_dir: task records of cells 124-152 (scripts/run-array.R); satt_results_dir: their
# Satterthwaite records (satterthwaite/run-satt-ext.R). Writes to out_dir:
# - cells.csv, failures.csv: summarise_study() of the new cells (all arms, estimands);
# - g-ladder.csv: core cells at G = 20 (new), 40 and 100 (summaries/final/cells.csv);
# - g-ladder-satt.csv: REML + CL2 with z, t_{G-1} and Satterthwaite critical values at
#   G = 20 (satterthwaite/summaries/satt-ext-*) and G = 40, 100 (core Satterthwaite study);
# - null.csv: beta = 0 cells, per arm: share of grid points whose interval excludes zero
#   (pointwise false-positive rate), share of replicates with any such point, and the size
#   of the estimate (RMS over the grid), also relative to the RMS of the non-null beta of the
#   matching core cell.
source("scripts/setup.R")
suppressPackageStartupMessages(library(dplyr))
args <- commandArgs(trailingOnly = TRUE)
results_dir <- args[1]
satt_dir <- args[2]
out_dir <- args[3]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
options(restart.cores = as.integer(Sys.getenv("RESTART_CORES", "1")))
cells <- design_table()
new_cells <- cells[cells$block %in% c("small_G", "null_beta"), ]
ARMS_MAIN <- c("reml_model", "reml_cl2", "ncv_model", "ncv_cl2", "ncv_cl2_bias")

s <- summarise_study(results_dir, new_cells)
s$cells <- left_join(
  new_cells,
  s$cells[, setdiff(names(s$cells), setdiff(names(new_cells), "cell"))],
  by = "cell"
)
readr::write_csv(s$cells, file.path(out_dir, "cells.csv"))
readr::write_csv(s$failures, file.path(out_dir, "failures.csv"))

# G ladder ---------------------------------------------------------------------------
ladder_cols <- c(
  "cell",
  "family",
  "signal",
  "error",
  "G",
  "estimand",
  "arm",
  "n_rep",
  "coverage_mean",
  "coverage_se",
  "pointwise_q05",
  "width_mean",
  "score_mean",
  "mse"
)
final <- readr::read_csv("summaries/final/cells.csv", show_col_types = FALSE)
ladder <- bind_rows(
  filter(final, block == "core"),
  filter(s$cells, block == "small_G")
) |>
  filter(arm %in% ARMS_MAIN) |>
  select(all_of(ladder_cols)) |>
  arrange(family, signal, error, estimand, arm, G)
readr::write_csv(ladder, file.path(out_dir, "g-ladder.csv"))

satt_new <- file.path("satterthwaite/summaries", "satt-syn-cells.csv")
if (file.exists(satt_new)) {
  sl <- bind_rows(
    readr::read_csv(
      "satterthwaite/summaries/satt-cells.csv",
      show_col_types = FALSE
    ),
    readr::read_csv(satt_new, show_col_types = FALSE) |>
      filter(
        cell %in% cells$cell[cells$block == "small_G" & cells$signal == "mid"]
      )
  ) |>
    select(
      cell,
      family,
      error,
      G,
      estimand,
      crit,
      n_rep,
      coverage,
      coverage_se,
      pointwise_q05,
      width,
      score,
      mean_crit
    ) |>
    arrange(family, error, estimand, crit, G)
  readr::write_csv(sl, file.path(out_dir, "g-ladder-satt.csv"))
}

# Null effect ------------------------------------------------------------------------
null_cells <- cells[cells$block == "null_beta", ]
# RMS of the non-null beta of the matching core cell (beta does not depend on the replicate)
ref_rms <- vapply(
  seq_len(nrow(null_cells)),
  \(i) {
    ref <- null_cells[i, ]
    ref$truth <- "smooth"
    sqrt(mean(simulate(ref, 1)$truth$beta^2))
  },
  0
)
null_cells$ref_rms <- ref_rms

z <- qnorm(1 - (1 - LEVEL) / 2)
null_rows <- function(est, se, crit, arm) {
  lower <- est - crit * se
  upper <- est + crit * se
  excl <- lower > 0 | upper < 0
  tibble(arm = arm, any = any(excl), excl = mean(excl), rms = sqrt(mean(est^2)))
}
files <- list.files(results_dir, "^task-c\\d+-r\\d+\\.rds$", full.names = TRUE)
fcell <- as.integer(sub("^task-c(\\d+)-r.*", "\\1", basename(files)))
per_rep <- bind_rows(parallel::mclapply(
  files[fcell %in% null_cells$cell],
  mc.cores = getOption("restart.cores"),
  \(f) {
    r <- readRDS(f)
    if (!identical(r$status, "success")) return(NULL)
    rec <- r$records$beta
    stopifnot(all(rec[, "truth"] == 0))
    rows <- lapply(ARMS_MAIN, \(arm) {
      spec <- ARMS[ARMS$arm == arm, ]
      est <- rec[, spec$est]
      se <- rec[, spec$se]
      if (spec$bias_aware) se <- sqrt(se^2 + (est - rec[, "est_REML"])^2)
      null_rows(est, se, z, arm)
    })
    sf <- file.path(satt_dir, basename(f))
    if (file.exists(sf)) {
      sr <- readRDS(sf)
      if (identical(sr$status, "success")) {
        m <- sr$records$beta
        G <- r$cell_spec$G
        q <- 1 - (1 - LEVEL) / 2
        rows <- c(
          rows,
          list(
            null_rows(
              m[, "est_REML"],
              m[, "se_cl2e_b"],
              qt(q, G - 1),
              "reml_cl2_t"
            ),
            null_rows(
              m[, "est_REML"],
              m[, "se_cl2e_b"],
              qt(q, pmax(m[, "df_full"], 1)),
              "reml_cl2_satt"
            )
          )
        )
      }
    }
    mutate(bind_rows(rows), cell = r$cell, rep = r$rep)
  }
))
null <- per_rep |>
  group_by(cell, arm) |>
  summarise(
    n_rep = n(),
    fpr_se = sd(excl) / sqrt(n()),
    fpr = mean(excl),
    fpr_q90 = unname(quantile(excl, 0.9)),
    any_excl = mean(any),
    share_gt_5pct = mean(excl > 0.05),
    rms = mean(rms),
    .groups = "drop"
  ) |>
  left_join(
    null_cells[, c("cell", "family", "error", "G", "ref_rms")],
    by = "cell"
  ) |>
  mutate(rms_rel = rms / ref_rms) |>
  relocate(family, error, G, .after = cell)
readr::write_csv(null, file.path(out_dir, "null.csv"))
cat(sprintf(
  "%d new cells summarised, %d failures; null: %d rows\n",
  n_distinct(s$cells$cell),
  nrow(s$failures),
  nrow(null)
))
