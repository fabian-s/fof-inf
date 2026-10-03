# Satterthwaite df worker for the extension cells, the new small-G / null-effect cells and
# the plasmode cells (2026-10-03). Run from the repository root (code copy):
#   Rscript satterthwaite/run-satt-ext.R <out_dir> syn <cells, R expression> <rep_from> <rep_to>
#   Rscript satterthwaite/run-satt-ext.R <out_dir> plasmode <app,app,...> <rep_from> <rep_to>
# Walks the rep-major task table and claims tasks through atomic lock directories, so any
# number of workers can run or be restarted. Records as run-satt.R (REML fit, exact CL2,
# per-point df). The z arm is checked against the study's summaries in the analysis.
# Plasmode task files hold values derived from real covariates: they stay on LRZ, only
# aggregates leave it.

source("scripts/setup.R")
source("ablation/ablation-lib.R")
source("satterthwaite/satt-lib.R")
suppressMessages(library(refund))

# Max absolute and relative differences between the REML columns and the study's stored
# record of the same cell and replicate (common random numbers), searched in the
# colon-separated directories of SATT_REF_DIRS (plasmode: <dir>/<app>/); NULL if none.
crn_check <- function(records, cell, rep) {
  dirs <- strsplit(Sys.getenv("SATT_REF_DIRS"), ":")[[1]]
  if (!is.null(cell$app)) dirs <- file.path(dirs, cell$app)
  name <- sprintf("task-c%03d-r%03d.rds", cell$cell, rep)
  # the first study run (cells 1-95) names files by the id of its 95-cell task table
  if (is.null(cell$app) && cell$cell <= 95) {
    name <- c(name, sprintf("task-%d.rds", (rep - 1) * 95 + cell$cell))
  }
  path <- Filter(file.exists, as.vector(outer(dirs, name, file.path)))
  if (!length(path)) return(NULL)
  ref <- readRDS(path[[1]])
  if (
    !identical(ref$status, "success") || ref$cell != cell$cell || ref$rep != rep
  ) {
    return(NULL)
  }
  rows <- lapply(names(records), \(estimand) {
    a <- records[[estimand]]
    b <- ref$records[[estimand]]
    do.call(
      rbind,
      lapply(names(CRN_COLUMNS), \(col) {
        x <- a[, col]
        y <- b[, CRN_COLUMNS[[col]]]
        data.frame(
          estimand = estimand,
          column = col,
          max_abs = max(abs(x - y)),
          max_rel = max(abs(x - y) / pmax(abs(y), 1e-8))
        )
      })
    )
  })
  do.call(rbind, rows)
}

run_task_ext <- function(cell, rep, sim_fun) {
  stage <- "simulate"
  out <- tryCatch(
    {
      sim <- sim_fun(cell, rep)
      stage <- "fit:REML"
      fit <- fit_model(sim, cell, "REML")
      stage <- "df"
      res <- satt_records(fit, sim, cell)
      c(
        list(status = "success"),
        res,
        list(
          fit_time = fit$restart$time,
          crn = crn_check(res$records, cell, rep)
        )
      )
    },
    error = \(e) {
      list(
        status = "failure",
        failures = data.frame(stage = stage, message = conditionMessage(e))
      )
    }
  )
  c(out, list(cell = cell$cell, rep = rep))
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("usage: run-satt-ext.R out_dir syn|plasmode cells|apps rep_from rep_to")
}
out_dir <- args[1]
mode <- args[2]
reps <- seq(as.integer(args[4]), as.integer(args[5]))
software <- software_versions()
check_software(software)
stopifnot(exists("pffr_influence_df", asNamespace("refund")))
if (mode == "syn") {
  cells <- design_table()
  cells <- cells[
    cells$cell %in% eval(parse(text = paste0("c(", args[3], ")"))),
  ]
  sim_fun <- simulate
} else if (mode == "plasmode") {
  source("plasmode/sources.R")
  source("plasmode/plasmode.R")
  apps <- strsplit(args[3], ",")[[1]]
  on_cluster <- nzchar(Sys.getenv("SLURM_JOB_ID"))
  cardio <- apps[apps %in% PLASMODE_LOCAL_ONLY]
  if (
    on_cluster &&
      length(cardio) &&
      (!all(cardio %in% PLASMODE_DEID_APPS) ||
        !identical(PLASMODE_CARDIO, "deid"))
  ) {
    stop(
      "cardiology apps on a cluster only from the de-identified input (PLASMODE_CARDIO=deid)"
    )
  }
  cells <- plasmode_cells()
  cells <- cells[
    cells$app %in% apps & cells$flip == "curve" & cells$error == "attached",
  ]
  sim_fun <- plasmode_simulate
} else {
  stop("mode must be syn or plasmode")
}
if (nrow(cells) == 0) stop("no cells selected")
dir.create(file.path(out_dir, "locks"), recursive = TRUE, showWarnings = FALSE)
tasks <- expand.grid(cell = cells$cell, rep = reps)
tasks <- tasks[order(tasks$rep, tasks$cell), ]
for (i in seq_len(nrow(tasks))) {
  id <- tasks$cell[i]
  rep <- tasks$rep[i]
  path <- file.path(out_dir, sprintf("task-c%03d-r%03d.rds", id, rep))
  if (file.exists(path)) next
  lock <- file.path(out_dir, "locks", sub("\\.rds$", ".lock", basename(path)))
  if (!dir.create(lock, showWarnings = FALSE)) next
  start <- proc.time()[["elapsed"]]
  record <- run_task_ext(cells[cells$cell == id, ], rep, sim_fun)
  record$mode <- mode
  record$software <- software
  record$time <- proc.time()[["elapsed"]] - start
  atomic_save(record, path)
  cat(sprintf(
    "cell %d rep %d: %s, %.1f s (df %.1f s)\n",
    id,
    rep,
    record$status,
    record$time,
    record$df_time %||% NA_real_
  ))
}
