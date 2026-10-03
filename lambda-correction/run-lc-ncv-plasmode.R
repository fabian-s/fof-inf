# Gaussian NCV lambda correction on the plasmode cells (IJ arms; no one-step, no fixed
# lambda). Run from the repository root (code copy):
#   Rscript lambda-correction/run-lc-ncv-plasmode.R <out_dir> <ref_dir> [max_tasks] [reps] [apps]
# ref_dir/<app>/task-cXXX-rYYY.rds: the plasmode study's records (CRN on est_NCV, se_cl2_NCV).
# Cardiology apps need PLASMODE_CARDIO=deid on a cluster; truth caches are read, never built.
source("scripts/setup.R")
source("plasmode/sources.R")
source("plasmode/plasmode.R")
source("lambda-correction/lc-lib.R")
source("lambda-correction/ncv-lib.R")
source("lambda-correction/lc-run-lib.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2)
  stop("usage: run-lc-ncv-plasmode.R out_dir ref_dir [max_tasks] [reps] [apps]")
reps <- if (length(args) >= 4 && args[4] != "all")
  eval(parse(text = args[4])) else seq_len(N_REP)
cells <- plasmode_cells()
cells <- cells[cells$app != "ecg_gls4", ]
if (length(args) >= 5)
  cells <- cells[cells$app %in% strsplit(args[5], ",")[[1]], ]
needed <- unique(truth_file(cells$app, cells$truth))
if (!all(file.exists(needed)))
  stop(
    "missing truth caches: ",
    paste(needed[!file.exists(needed)], collapse = ", ")
  )
tasks <- tidyr::expand_grid(rep = reps, cell = cells$cell)
tasks$task <- seq_len(nrow(tasks))
run_worker(
  args[1],
  cells,
  tasks,
  ref_path = \(task)
    file.path(
      args[2],
      cells$app[cells$cell == task$cell],
      sprintf("task-c%03d-r%03d.rds", task$cell, task$rep)
    ),
  max_tasks = if (length(args) >= 3) as.numeric(args[3]) else Inf,
  task_fun = \(task, cell, ref)
    ncv_task(task, cell, ref, NULL, plasmode_simulate, onestep = FALSE)
)
