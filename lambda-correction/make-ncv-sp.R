# Fixed smoothing parameters for the NCV fixed-lambda oracle: per core cell (mid signal)
# and penalty, the geometric mean of the study's NCV smoothing parameters over its 200
# replicates. Run from the repository root: Rscript lambda-correction/make-ncv-sp.R
source("scripts/setup.R")
source("ablation/ablation-lib.R")
study <- readRDS("summaries/final/study.rds")
core <- ablation_cells()
sp <- study$terms |>
  dplyr::filter(fit == "NCV", cell %in% core$cell) |>
  dplyr::group_by(cell, term, penalty) |>
  dplyr::summarise(
    n = dplyr::n(),
    mean_log_sp = mean(log_sp),
    sd_log_sp = sd(log_sp),
    .groups = "drop"
  )
stopifnot(all(sp$n == N_REP))
write.csv(sp, "lambda-correction/ncv-sp-fixed.csv", row.names = FALSE)
print(as.data.frame(sp), digits = 3)
