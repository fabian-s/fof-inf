# Lambda-correction study worker. Run from the repository root (code copy):
#   Rscript lambda-correction/run-lc.R <out_dir> <ref_dir> [max_tasks] [reps] [nojk|jk] [cells]
# cells: comma-separated cell ids (default: all 18 core cells).
# Same task table, data and REML fits as the covariance ablation (ablation/): the 18 core
# cells at the default signal, replicates 1-200, rep-major. Workers claim tasks through
# atomic lock directories. ref_dir holds the study's task records (task-<id>.rds) for the
# common-random-numbers check; `reps` (e.g. "1:2") restricts the replicates (pilot).

source("scripts/setup.R")
source("ablation/ablation-lib.R") # ablation_cells(), ablation_tasks()
source("lambda-correction/lc-lib.R")
source("lambda-correction/lc-run-lib.R")

STUDY_N_CELLS <- 95

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: run-lc.R out_dir ref_dir [max_tasks] [reps]")
reps <- if (length(args) >= 4) eval(parse(text = args[4])) else seq_len(N_REP)
cells <- ablation_cells()
if (length(args) >= 6)
  cells <- cells[cells$cell %in% as.integer(strsplit(args[6], ",")[[1]]), ]
run_worker(
  args[1],
  cells,
  ablation_tasks(cells, reps),
  ref_path = \(task) {
    file.path(
      args[2],
      sprintf("task-%d.rds", (task$rep - 1) * STUDY_N_CELLS + task$cell)
    )
  },
  jackknife = !(length(args) >= 5 && args[5] == "nojk"),
  max_tasks = if (length(args) >= 3) as.numeric(args[3]) else Inf
)
