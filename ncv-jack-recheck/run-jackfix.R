# LRZ worker for the corrected NCV covariance (jackfix_vcov) on the study's replicates or
# the plasmode cells: walks the rep-major task table and claims tasks via atomic lock dirs.
#   Rscript ncv-jack-recheck/run-jackfix.R <out_dir> syn <cell,cell,...> <rep_from> <rep_to>
#   Rscript ncv-jack-recheck/run-jackfix.R <out_dir> plasmode <app,app,...> <rep_from> <rep_to>
# Plasmode task files hold values derived from real covariates: they stay on LRZ, only
# aggregates leave it.
source("scripts/setup.R")
source("ncv-jack-recheck/jackfix-lib.R")
suppressMessages(library(refund))
args <- commandArgs(trailingOnly = TRUE)
out <- args[1]
mode <- args[2]
which <- strsplit(args[3], ",")[[1]]
reps <- seq(as.integer(args[4]), as.integer(args[5]))
if (mode == "syn") {
  cells <- design_table()
  cells <- cells[cells$cell %in% as.integer(which), ]
  sim_fun <- simulate
} else {
  source("plasmode/sources.R")
  source("plasmode/plasmode.R")
  if (
    any(which %in% PLASMODE_LOCAL_ONLY & !which %in% PLASMODE_DEID_APPS) ||
      (any(which %in% PLASMODE_DEID_APPS) &&
        !identical(PLASMODE_CARDIO, "deid"))
  )
    stop(
      "cardiology apps on a cluster only from the de-identified input (PLASMODE_CARDIO=deid)"
    )
  cells <- plasmode_cells()
  cells <- cells[
    cells$app %in%
      which &
      cells$flip == "curve" &
      cells$G > 40 &
      cells$error == "attached",
  ]
  sim_fun <- plasmode_simulate
}
dir.create(file.path(out, "locks"), recursive = TRUE, showWarnings = FALSE)
tasks <- expand.grid(cell = cells$cell, rep = reps)
tasks <- tasks[order(tasks$rep, tasks$cell), ]
for (i in seq_len(nrow(tasks))) {
  id <- tasks$cell[i]
  rep <- tasks$rep[i]
  f <- file.path(out, sprintf("c%03d-r%03d.rds", id, rep))
  if (file.exists(f)) next
  if (!dir.create(file.path(out, "locks", basename(f)), showWarnings = FALSE))
    next
  res <- tryCatch(
    jackfix_task(cells[cells$cell == id, ], rep, sim_fun),
    error = \(e) list(error = conditionMessage(e))
  )
  res$mgcv <- paste(
    format(packageVersion("mgcv")),
    packageDescription("mgcv")$Patched
  )
  tmp <- tempfile(tmpdir = out)
  saveRDS(res, tmp)
  file.rename(tmp, f)
  cat(sprintf(
    "cell %d rep %d %s\n",
    id,
    rep,
    if (is.null(res$error)) sprintf("%.1fs", res$time) else res$error
  ))
}
