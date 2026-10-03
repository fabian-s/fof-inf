# Fixed smoothing parameters for the fixed-lambda arm: per core cell (mid signal) and
# per penalty, the geometric mean of the study's REML smoothing parameters over its 200
# replicates (exp of the mean log sp). Per penalty, not per term: the ff term has two
# penalties, which summaries/final/term-edf.csv pools into one mean_log_sp.
# Run from the repository root: Rscript ablation/make-sp.R  ->  ablation/sp-fixed.csv

source("scripts/setup.R")

core_cells <- function(cells = design_table()) {
  cells[cells$block == "core" & cells$signal == "mid", ]
}

study <- readRDS("summaries/final/study.rds")
core <- core_cells()
sp <- study$terms |>
  dplyr::filter(fit == "REML", cell %in% core$cell) |>
  dplyr::group_by(cell, term, penalty) |>
  dplyr::summarise(n = dplyr::n(), mean_log_sp = mean(log_sp), .groups = "drop")
stopifnot(all(sp$n == N_REP))

# Cross-check: for single-penalty terms the per-penalty mean is term-edf.csv's value.
edf <- read.csv("summaries/final/term-edf.csv", check.names = FALSE) |>
  dplyr::filter(fit == "REML", n_penalties == 1)
check <- dplyr::inner_join(sp, edf, by = c("cell", "term"))
stopifnot(
  nrow(check) == 2 * nrow(core),
  max(abs(check$mean_log_sp.x - check$mean_log_sp.y)) < 1e-8
)
write.csv(sp, "ablation/sp-fixed.csv", row.names = FALSE)
cat("wrote", nrow(sp), "rows\n")
