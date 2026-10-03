# Run the plasmode tasks of one app on one of several worker processes.
# Run from the repository root:
#   Rscript plasmode/run-plasmode.R <app> <worker> <n_workers> [max_rep]
# With worker "prepare", the app's truth (and beat) caches are built and nothing else:
# run it once before starting parallel workers, which otherwise all build them.
# Results go to PLASMODE_RESULTS/<app> (default PLASMODE_CACHE/results; outside the
# repository: task files hold per-subject curves).
source("scripts/setup.R")
source("plasmode/sources.R")
source("plasmode/plasmode.R")
suppressPackageStartupMessages(library(refund))

# Peak resident memory of this process (Linux), for sizing --mem.
peak_memory_mb <- function() {
  status <- readLines("/proc/self/status")
  kb <- as.numeric(gsub("\\D", "", grep("^VmHWM", status, value = TRUE)))
  round(kb / 1024)
}

args <- commandArgs(trailingOnly = TRUE)
app <- args[1]
prepare <- identical(args[2], "prepare")
worker <- if (prepare) NA_integer_ else as.integer(args[2])
n_workers <- as.integer(args[3])
max_rep <- if (length(args) >= 4) as.integer(args[4]) else N_REP
# The cardiology data are patient data: this machine only, never a cluster job.
# On a cluster, the extension's cardiology apps run from the de-identified input only.
on_cluster <- nzchar(Sys.getenv("SLURM_JOB_ID"))
deid_ok <- app %in% PLASMODE_DEID_APPS && identical(PLASMODE_CARDIO, "deid")
if (app %in% PLASMODE_LOCAL_ONLY && on_cluster && !deid_ok) {
  stop(
    app,
    " uses patient data: on a cluster only de-identified (PLASMODE_CARDIO=deid)"
  )
}

cells <- plasmode_cells()
if (!app %in% cells$app) stop("unknown plasmode app: ", app)
if (prepare) {
  mine <- cells[cells$app == app, ]
  for (source in unique(c(mine$truth, mine$residual))) {
    start <- proc.time()[["elapsed"]]
    if (startsWith(source, "beat")) {
      plasmode_beats(app)
    } else {
      plasmode_truth(app, source)
    }
    cat(sprintf(
      "%s %s: %.1f s\n",
      app,
      source,
      proc.time()[["elapsed"]] - start
    ))
  }
  mid <- file.path(PLASMODE_CACHE, paste0("truth-", app, "-MID.rds"))
  if (file.exists(mid)) cat("MID interior:", readRDS(mid)$mid$interior, "\n")
  cat("peak memory (MB):", peak_memory_mb(), "\n")
  quit(save = "no")
}
tasks <- task_table(cells, reps = seq_len(max_rep))
mine <- tasks[tasks$cell %in% cells$cell[cells$app == app], ]
mine <- mine[seq_len(nrow(mine)) %% n_workers == worker, ]
results <- Sys.getenv("PLASMODE_RESULTS", file.path(PLASMODE_CACHE, "results"))
out_dir <- file.path(results, app)
for (id in mine$task) {
  start <- proc.time()[["elapsed"]]
  run_task(id, out_dir, cells, tasks, simulate_fun = plasmode_simulate)
  task <- mine[mine$task == id, ]
  cat(sprintf(
    "task %d cell %d rep %d: %.1f s\n",
    id,
    task$cell,
    task$rep,
    proc.time()[["elapsed"]] - start
  ))
}
cat("worker", worker, "done:", nrow(mine), "tasks\n")
cat("peak memory (MB):", peak_memory_mb(), "\n")
