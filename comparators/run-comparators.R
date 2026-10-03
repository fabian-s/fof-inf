# Comparator study worker. Run from the repository root (code copy), BLAS pinned:
#   Rscript comparators/run-comparators.R <out_dir> <ref_syn> <ref_plas> [max_tasks]
# Workers walk the rep-major task table and claim tasks through atomic lock directories
# (<out_dir>/locks), so any number of workers can run or be restarted.
# ref_syn: the study's synthetic task records (task-<id>.rds, 95-cell task table);
# ref_plas: the plasmode results (<app>/task-cXXX-rYYY.rds). Missing records give no
# CRN rows. Environment:
#   COMP_KEYS      comma-separated cell keys to run (default all 18; e.g. s012,p081)
#   COMP_MAX_REP   replicates 1..COMP_MAX_REP (default N_REP = 200)
#   COMP_BOOT_B    bootstrap size (default 199)
#   COMP_BOOT_REPS bootstrap arm only for replicates <= this (default N_REP); per source
#                  via COMP_BOOT_REPS_SYN / COMP_BOOT_REPS_PLAS
# Cardiology apps (ecg_gls4b, af) run on a cluster from the de-identified input only
# (PLASMODE_CARDIO=deid); their task files stay on the cluster.

source("scripts/setup.R")
source("plasmode/sources.R")
source("plasmode/plasmode.R")
source("comparators/comparators-lib.R")
suppressPackageStartupMessages(library(refund))

env_int <- function(name, default) {
  value <- Sys.getenv(name)
  if (nzchar(value)) as.integer(value) else as.integer(default)
}

task_path <- function(out_dir, task) {
  file.path(out_dir, sprintf("task-%s-r%03d.rds", task$key, task$rep))
}

peak_memory_mb <- function() {
  status <- readLines("/proc/self/status")
  round(
    as.numeric(gsub("\\D", "", grep("^VmHWM", status, value = TRUE))) / 1024
  )
}

run_one <- function(task, cell, ref_syn, ref_plas, boot_reps, B) {
  stage <- "simulate"
  boot <- task$rep <= boot_reps[[cell$source]]
  tryCatch(
    {
      sim <- if (cell$source == "synthetic") {
        simulate(cell, task$rep)
      } else {
        plasmode_simulate(cell, task$rep)
      }
      stage <- "fit"
      out <- run_comparators(sim, cell, task$rep, boot = boot, B = B)
      stage <- "crn"
      out$crn <- crn_check(
        out$records,
        reference_path(task, cell, ref_syn, ref_plas),
        cell,
        task$rep
      )
      c(list(status = "success", boot = boot, B = B), out)
    },
    error = \(e) {
      list(
        status = "failure",
        boot = boot,
        B = B,
        failures = data.frame(stage = stage, message = conditionMessage(e))
      )
    }
  )
}

run_worker <- function(out_dir, ref_syn, ref_plas, max_tasks = Inf) {
  software <- software_versions()
  check_software(software)
  cells <- comparator_cells()
  keys <- strsplit(
    Sys.getenv("COMP_KEYS", paste(cells$key, collapse = ",")),
    ","
  )[[
    1
  ]]
  if (!all(keys %in% cells$key)) stop("unknown cell key in COMP_KEYS")
  on_cluster <- nzchar(Sys.getenv("SLURM_JOB_ID"))
  cardio <- cells$app[cells$key %in% keys] %in% PLASMODE_LOCAL_ONLY
  if (any(cardio) && on_cluster && !identical(PLASMODE_CARDIO, "deid")) {
    stop(
      "cardiology apps on a cluster only from de-identified input (PLASMODE_CARDIO=deid)"
    )
  }
  tasks <- comparator_tasks(keys, seq_len(env_int("COMP_MAX_REP", N_REP)))
  B <- env_int("COMP_BOOT_B", BOOT_B)
  all_reps <- env_int("COMP_BOOT_REPS", N_REP)
  boot_reps <- list(
    synthetic = env_int("COMP_BOOT_REPS_SYN", all_reps),
    plasmode = env_int("COMP_BOOT_REPS_PLAS", all_reps)
  )
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
    cell <- cell_spec(task$key)
    start <- elapsed()
    record <- run_one(task, cell, ref_syn, ref_plas, boot_reps, B)
    record[c("task", "key", "cell", "rep", "source", "app")] <- list(
      id,
      task$key,
      cell$cell,
      task$rep,
      cell$source,
      if ("app" %in% names(cell)) cell$app else NA_character_
    )
    record$software <- software
    record$time <- elapsed() - start
    record$peak_mb <- peak_memory_mb()
    atomic_save(record, path)
    cat(sprintf(
      "[%s] task %d (%s rep %d): %s, %.1f s, %d MB\n",
      format(Sys.time(), "%H:%M:%S"),
      id,
      task$key,
      task$rep,
      record$status,
      record$time,
      record$peak_mb
    ))
    done <- done + 1
  }
  invisible(done)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("usage: run-comparators.R out_dir ref_syn ref_plas [max_tasks]")
}
run_worker(
  args[1],
  args[2],
  args[3],
  if (length(args) >= 4) as.numeric(args[4]) else Inf
)
