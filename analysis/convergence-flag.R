# Scaled-t fits: accuracy and coverage in the replicates whose NCV fit mgcv flags as not converged
# against the others, per interval, error process and estimand. Aggregates only. Run from the repository root:
#   Rscript analysis/convergence-flag.R <results_dir> summaries/final/convergence-flag.csv
# results_dir holds the task records of the study's first run (task-<id>.rds, 95-cell table).
source("scripts/setup.R")
suppressPackageStartupMessages(library(dplyr))
args <- commandArgs(trailingOnly = TRUE)
results <- args[1]
cells <- design_table()
scat <- cells[cells$family == "scat", ]
stopifnot(nrow(scat) == 2, all(scat$cell <= 95))
ARMS_FLAG <- c(reml_cl2 = "REML", ncv_cl2 = "NCV", ncv_model = "NCV")
CENTRED <- c("mean", "alpha")

# `flagged` is the NCV fit's flag in this replicate, for every interval: REML + CL2 split by it
# shows whether the flagged replicates are harder for any fit.
task_rows <- function(cell, rep) {
  path <- file.path(results, sprintf("task-%d.rds", (rep - 1) * 95 + cell$cell))
  r <- readRDS(path)
  stopifnot(identical(r$status, "success"), r$cell == cell$cell, r$rep == rep)
  linkinv <- FAMILIES[[cell$fit_family]]$linkinv
  flagged <- !r$fits$converged[r$fits$fit == "NCV"]
  stopifnot(r$fits$converged[r$fits$fit == "REML"])
  rows <- list()
  for (arm in names(ARMS_FLAG)) {
    fit <- ARMS_FLAG[[arm]]
    for (e in names(r$records)) {
      rec <- r$records[[e]]
      s <- point_scores(rec, arm, e, linkinv)
      tr <- rec[, "truth"]
      est <- rec[, paste0("est_", fit)]
      if (e == "mean") est <- linkinv(est)
      denom <- if (e %in% CENTRED) mean((tr - mean(tr))^2) else mean(tr^2)
      rows[[length(rows) + 1]] <- tibble(
        error = cell$error,
        arm = arm,
        estimand = e,
        flagged = flagged,
        coverage = mean(s$covered),
        rel = sqrt(mean((est - tr)^2) / denom)
      )
    }
  }
  bind_rows(rows)
}

per_rep <- bind_rows(lapply(seq_len(nrow(scat)), \(i) {
  bind_rows(lapply(seq_len(N_REP), \(rep) task_rows(scat[i, ], rep)))
}))
out <- per_rep |>
  group_by(error, arm, estimand, flagged) |>
  summarise(
    n_rep = n(),
    coverage_se = sd(coverage) / sqrt(n()),
    coverage = mean(coverage),
    rel = median(rel),
    .groups = "drop"
  )
readr::write_csv(out, args[2])
print(as.data.frame(out))
