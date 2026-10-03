# Satterthwaite df worker. Run from the repository root (code copy):
#   Rscript satterthwaite/run-satt.R <out_dir> <ref_dir> [max_tasks]
# Walks the rep-major task table of the 18 core cells and claims tasks through atomic
# lock directories (<out_dir>/locks), so any number of workers can run or be restarted.
# ref_dir holds the study's task records (task-<id>.rds, 95-cell task table) against
# which the REML fit is checked (common random numbers); a missing record gives NA.

source("scripts/setup.R")
source("ablation/ablation-lib.R")
source("satterthwaite/satt-lib.R")

STUDY_N_CELLS <- 95 # task table of the study run whose records are in ref_dir

# One replicate: satt_records() (satt-lib.R) ------------------------------------------

run_satt_task <- function(task, cell, ref_dir) {
  stage <- "simulate"
  out <- tryCatch(
    {
      sim <- simulate(cell, task$rep)
      stage <- "fit:REML"
      fit <- fit_model(sim, cell, "REML")
      stage <- "df"
      res <- satt_records(fit, sim, cell)
      ref_path <- file.path(
        ref_dir,
        sprintf("task-%d.rds", (task$rep - 1) * STUDY_N_CELLS + task$cell)
      )
      c(
        list(status = "success"),
        res,
        list(
          fit_time = fit$restart$time,
          crn = crn_check(res$records, ref_path, task)
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
  c(out, list(task = task$task, cell = task$cell, rep = task$rep))
}

# Max absolute and relative differences between the REML columns and the study record
# (as in ablation/run-ablation.R).
crn_check <- function(records, ref_path, task) {
  if (!file.exists(ref_path)) return(NULL)
  ref <- readRDS(ref_path)
  if (ref$cell != task$cell || ref$rep != task$rep) {
    stop(
      "reference record ",
      ref_path,
      " is not cell/rep ",
      task$cell,
      "/",
      task$rep
    )
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

# Worker loop -------------------------------------------------------------------------

task_path <- function(out_dir, task) {
  file.path(out_dir, sprintf("task-c%03d-r%03d.rds", task$cell, task$rep))
}

run_worker <- function(out_dir, ref_dir, max_tasks = Inf) {
  software <- software_versions()
  check_software(software)
  stopifnot(exists("pffr_influence_df", asNamespace("refund")))
  cells <- ablation_cells()
  tasks <- ablation_tasks(cells)
  dir.create(
    file.path(out_dir, "locks"),
    recursive = TRUE,
    showWarnings = FALSE
  )
  done <- 0
  for (id in tasks$task) {
    if (done >= max_tasks) break
    task <- tasks[tasks$task == id, ]
    path <- task_path(out_dir, task)
    if (file.exists(path)) next
    lock <- file.path(out_dir, "locks", sub("\\.rds$", ".lock", basename(path)))
    if (!dir.create(lock, showWarnings = FALSE)) next
    start <- proc.time()[["elapsed"]]
    record <- run_satt_task(task, cells[cells$cell == task$cell, ], ref_dir)
    record$software <- software
    record$time <- proc.time()[["elapsed"]] - start
    atomic_save(record, path)
    cat(sprintf(
      "task %d (cell %d rep %d): %s, %.1f s (df %.1f s)\n",
      id,
      task$cell,
      task$rep,
      record$status,
      record$time,
      record$df_time %||% NA_real_
    ))
    done <- done + 1
  }
  invisible(done)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: run-satt.R out_dir ref_dir [max_tasks]")
run_worker(
  args[1],
  args[2],
  if (length(args) >= 3) as.numeric(args[3]) else Inf
)
