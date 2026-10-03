# Misspecified-truth stress test: worker. Run from the repository root (code copy):
#   Rscript misspec/run-misspec.R <out_dir> <ncv_xlarge_max_rep> [rep_max] [hours]
# Walks the rep-major task table (reps <= rep_max) and claims tasks through atomic lock
# directories (<out_dir>/locks), so any number of workers can run or be restarted. Stops
# claiming new tasks after `hours` (default 22) so no task is killed by the time limit.
# xlarge cells get the NCV fit only in replicates <= ncv_xlarge_max_rep (cost).

source("scripts/setup.R")
source("misspec/misspec-lib.R")

XLARGE_TASK_HOURS <- 0.6 # head room before the deadline for one xlarge task

# One replicate -----------------------------------------------------------------------

# Per estimand: truth plus est/se columns of every fit (as R/score.R's score_task, which
# requires the NCV fit).
assemble_records <- function(inputs, sim, cell) {
  lapply(setNames(nm = model_estimands(cell$model)), \(estimand) {
    columns <- list(truth = as.vector(sim$truth[[estimand]]))
    for (method in names(inputs)) {
      value <- inputs[[method]]$estimands[[estimand]]
      columns[paste0(names(value), "_", method)] <- value
    }
    out <- do.call(cbind, columns)
    if (!all(is.finite(out))) stop("non-finite estimate or SE for ", estimand)
    out
  })
}

# Max abs/rel differences to the study's record of the same data (T0 cells only).
crn_check <- function(records, cell, rep) {
  if (is.na(cell$study_cell)) return(NULL)
  path <- study_record_path(cell$study_cell, rep)
  if (!file.exists(path)) return(data.frame(note = paste("missing", path)))
  ref <- readRDS(path)
  if (ref$cell != cell$study_cell || ref$rep != rep)
    stop("wrong study record ", path)
  rows <- lapply(names(records), \(estimand) {
    a <- records[[estimand]]
    b <- ref$records[[estimand]]
    cols <- intersect(colnames(a), colnames(b))
    do.call(
      rbind,
      lapply(cols, \(col) {
        data.frame(
          estimand = estimand,
          column = col,
          max_abs = max(abs(a[, col] - b[, col])),
          max_rel = max(abs(a[, col] - b[, col]) / pmax(abs(b[, col]), 1e-8))
        )
      })
    )
  })
  do.call(rbind, rows)
}

run_misspec_task <- function(task, cell, ncv_xlarge_max_rep) {
  stage <- "simulate"
  tryCatch(
    {
      sim <- simulate(cell, task$rep)
      methods <- "REML"
      if (cell$basis == "default" || task$rep <= ncv_xlarge_max_rep) {
        methods <- c(methods, "NCV")
      }
      inputs <- list()
      for (method in methods) {
        stage <- paste0("fit:", method)
        inputs[[method]] <- run_method(sim, cell, method)
      }
      stage <- "assemble"
      records <- assemble_records(inputs, sim, cell)
      list(
        status = "success",
        records = records,
        fits = do.call(rbind, lapply(inputs, `[[`, "fit_row")),
        warnings = lapply(inputs, `[[`, "warnings"),
        diag = sim$diag,
        crn = crn_check(records, cell, task$rep)
      )
    },
    error = \(e) {
      list(
        status = "failure",
        failure = data.frame(stage = stage, message = conditionMessage(e))
      )
    }
  )
}

# Worker loop -------------------------------------------------------------------------

task_path <- function(out_dir, task) {
  file.path(out_dir, sprintf("task-c%03d-r%03d.rds", task$cell, task$rep))
}

run_worker <- function(out_dir, ncv_xlarge_max_rep, rep_max, hours) {
  t0 <- proc.time()[["elapsed"]]
  software <- software_versions()
  check_software(software)
  cells <- misspec_cells()
  tasks <- misspec_tasks(cells)
  tasks <- tasks[tasks$rep <= rep_max, ]
  dir.create(
    file.path(out_dir, "locks"),
    recursive = TRUE,
    showWarnings = FALSE
  )
  for (id in tasks$task) {
    task <- tasks[tasks$task == id, ]
    cell <- cells[cells$cell == task$cell, ]
    left <- hours - (proc.time()[["elapsed"]] - t0) / 3600
    if (left < (if (cell$basis == "xlarge") XLARGE_TASK_HOURS else 0.05)) next
    path <- task_path(out_dir, task)
    if (file.exists(path)) next
    lock <- file.path(out_dir, "locks", sub("\\.rds$", ".lock", basename(path)))
    if (!dir.create(lock, showWarnings = FALSE)) next
    writeLines(Sys.getenv("SLURM_JOB_ID", "local"), file.path(lock, "owner"))
    start <- proc.time()[["elapsed"]]
    record <- run_misspec_task(task, cell, ncv_xlarge_max_rep)
    record <- c(
      record,
      list(
        task = task$task,
        cell = task$cell,
        rep = task$rep,
        software = software,
        time = proc.time()[["elapsed"]] - start
      )
    )
    atomic_save(record, path)
    cat(sprintf(
      "task %d (cell %d rep %d): %s, %.1f s\n",
      id,
      task$cell,
      task$rep,
      record$status,
      record$time
    ))
  }
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("usage: run-misspec.R out_dir ncv_xlarge_max_rep [rep_max] [hours]")
}
run_worker(
  args[1],
  as.numeric(args[2]),
  if (length(args) >= 3) as.numeric(args[3]) else N_REP,
  if (length(args) >= 4) as.numeric(args[4]) else 22
)
