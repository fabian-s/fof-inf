# Gaussian NCV lambda-correction worker (core Gaussian cells, mid signal, reps 1-200,
# rep-major, lock-dir claiming). Run from the repository root (code copy):
#   Rscript lambda-correction/run-lc-ncv.R <out_dir> <ref_dir> [max_tasks] [reps]
source("scripts/setup.R")
source("ablation/ablation-lib.R")
source("lambda-correction/lc-lib.R")
source("lambda-correction/ncv-lib.R")
source("lambda-correction/lc-run-lib.R")

STUDY_N_CELLS <- 95
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2)
  stop("usage: run-lc-ncv.R out_dir ref_dir [max_tasks] [reps]")
reps <- if (length(args) >= 4) eval(parse(text = args[4])) else seq_len(N_REP)
cells <- ablation_cells()
cells <- cells[cells$family == "gaussian", ]
sp_table <- read.csv("lambda-correction/ncv-sp-fixed.csv")
run_worker(
  args[1],
  cells,
  ablation_tasks(cells, reps),
  ref_path = \(task)
    file.path(
      args[2],
      sprintf("task-%d.rds", (task$rep - 1) * STUDY_N_CELLS + task$cell)
    ),
  max_tasks = if (length(args) >= 3) as.numeric(args[3]) else Inf,
  task_fun = \(task, cell, ref) ncv_task(task, cell, ref, sp_table)
)
