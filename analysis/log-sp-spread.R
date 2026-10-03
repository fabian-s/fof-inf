# Spread of the selected log smoothing parameters and term EDFs over replicates, per cell,
# fit, term and penalty, from the study's per-replicate term table.
#   Rscript analysis/log-sp-spread.R [study.rds] [out.csv]
suppressPackageStartupMessages(library(dplyr))
args <- commandArgs(trailingOnly = TRUE)
study <- readRDS(
  if (length(args) >= 1) args[1] else "summaries/final/study.rds"
)
out <- if (length(args) >= 2) args[2] else "summaries/final/log-sp-spread.csv"
spread <- study$terms |>
  group_by(cell, fit, term, penalty) |>
  summarise(
    n_rep = n(),
    log_sp_mean = mean(log_sp),
    log_sp_sd = sd(log_sp),
    log_sp_q10 = unname(quantile(log_sp, 0.1)),
    log_sp_q90 = unname(quantile(log_sp, 0.9)),
    edf_sd = sd(edf),
    .groups = "drop"
  )
readr::write_csv(spread, out)
cat(nrow(spread), "rows ->", out, "\n")
