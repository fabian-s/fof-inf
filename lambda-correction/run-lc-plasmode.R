# Lambda-correction on the plasmode cells (IJ and one-step arms only, no jackknife refits).
# Run from the repository root (code copy):
#   Rscript lambda-correction/run-lc-plasmode.R <out_dir> <ref_dir> [max_tasks] [reps] [apps]
# Same cells, replicates, data and REML fits as the plasmode study (plasmode/): every cell
# except the original patient-data app ecg_gls4 (local only), replicates 1-200, rep-major.
# ref_dir/<app>/task-cXXX-rYYY.rds are the plasmode study's records (CRN check). The
# cardiology apps ecg_gls4b and af need PLASMODE_CARDIO=deid on a cluster (plasmode/sources.R);
# PLASMODE_CACHE must hold the study's truth caches (they are read, never rebuilt here).

source("scripts/setup.R")
source("plasmode/sources.R")
source("plasmode/plasmode.R")
source("lambda-correction/lc-lib.R")
source("lambda-correction/lc-run-lib.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2)
  stop("usage: run-lc-plasmode.R out_dir ref_dir [max_tasks] [reps] [apps]")
reps <- if (length(args) >= 4 && args[4] != "all")
  eval(parse(text = args[4])) else seq_len(N_REP)
cells <- plasmode_cells()
cells <- cells[!cells$app %in% "ecg_gls4", ]
if (length(args) >= 5)
  cells <- cells[cells$app %in% strsplit(args[5], ",")[[1]], ]
# never build a truth cache from a worker: every needed cache must exist
needed <- unique(c(truth_file(cells$app, cells$truth)))
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
  ref_path = \(task) {
    app <- cells$app[cells$cell == task$cell]
    file.path(
      args[2],
      app,
      sprintf("task-c%03d-r%03d.rds", task$cell, task$rep)
    )
  },
  simulate_fun = plasmode_simulate,
  jackknife = FALSE,
  max_tasks = if (length(args) >= 3) as.numeric(args[3]) else Inf
)
