# CL2-arms worker (GLS raw / fpc + CL2, pcre + CL2) on the comparator cells, common
# random numbers with the study and the comparator run. Run from the repository root
# (code copy), BLAS pinned:
#   Rscript comparators/run-cl2-arms.R <out_dir> <ref_syn> <ref_plas> <comp_dir> [max_tasks]
# Workers walk the rep-major task table (comparator_tasks) and claim tasks through atomic
# lock directories (<out_dir>/locks). Environment:
#   CL2_KEYS       comma-separated cell keys (default all 18)
#   CL2_MAX_REP    replicates 1..CL2_MAX_REP (default N_REP = 200)
#   CL2_PCRE_REPS  pcre (refit) only for replicates <= this (default 100)
# Cardiology apps only from the de-identified input (PLASMODE_CARDIO=deid); their task
# files stay on the cluster.

source("scripts/setup.R")
source("plasmode/sources.R")
source("plasmode/plasmode.R")
source("comparators/comparators-lib.R")
source("comparators/cl2-arms-lib.R")
suppressPackageStartupMessages(library(refund))

env_int <- function(name, default) {
  value <- Sys.getenv(name)
  if (nzchar(value)) as.integer(value) else as.integer(default)
}

task_file <- function(dir, task) {
  file.path(dir, sprintf("task-%s-r%03d.rds", task$key, task$rep))
}

peak_memory_mb <- function() {
  status <- readLines("/proc/self/status")
  round(
    as.numeric(gsub("\\D", "", grep("^VmHWM", status, value = TRUE))) / 1024
  )
}

CRN_STUDY <- c("truth", "est_REML", "se_model_REML", "se_cl2_REML")
CRN_COMP <- c(
  CRN_STUDY,
  "est_GLSraw",
  "se_GLSraw",
  "est_GLSfpc",
  "se_GLSfpc",
  "est_PCRE",
  "se_PCRE"
)

run_one <- function(task, cell, ref_syn, ref_plas, comp_dir, pcre) {
  stage <- "simulate"
  tryCatch(
    {
      sim <- if (cell$source == "synthetic") {
        simulate(cell, task$rep)
      } else {
        plasmode_simulate(cell, task$rep)
      }
      stage <- "fit"
      out <- run_cl2_arms(sim, cell, pcre = pcre)
      stage <- "crn"
      out$crn <- rbind(
        crn_against(
          out$records,
          reference_path(task, cell, ref_syn, ref_plas),
          cell,
          task$rep,
          CRN_STUDY,
          "study"
        ),
        crn_against(
          out$records,
          task_file(comp_dir, task),
          cell,
          task$rep,
          CRN_COMP,
          "comparators"
        )
      )
      c(list(status = "success", pcre = pcre), out)
    },
    error = \(e) {
      list(
        status = "failure",
        pcre = pcre,
        failures = data.frame(stage = stage, message = conditionMessage(e))
      )
    }
  )
}

run_worker <- function(out_dir, ref_syn, ref_plas, comp_dir, max_tasks = Inf) {
  software <- software_versions()
  check_software(software)
  cells <- comparator_cells()
  keys <- strsplit(
    Sys.getenv("CL2_KEYS", paste(cells$key, collapse = ",")),
    ","
  )[[1]]
  if (!all(keys %in% cells$key)) stop("unknown cell key in CL2_KEYS")
  on_cluster <- nzchar(Sys.getenv("SLURM_JOB_ID"))
  cardio <- cells$app[cells$key %in% keys] %in% PLASMODE_LOCAL_ONLY
  if (any(cardio) && on_cluster && !identical(PLASMODE_CARDIO, "deid")) {
    stop("cardiology apps on a cluster only from de-identified input")
  }
  tasks <- comparator_tasks(keys, seq_len(env_int("CL2_MAX_REP", N_REP)))
  pcre_reps <- env_int("CL2_PCRE_REPS", 100)
  dir.create(
    file.path(out_dir, "locks"),
    recursive = TRUE,
    showWarnings = FALSE
  )
  done <- 0
  for (id in tasks$task) {
    if (done >= max_tasks) break
    task <- tasks[tasks$task == id, ]
    path <- task_file(out_dir, task)
    if (file.exists(path)) next
    lock <- file.path(out_dir, "locks", sub("\\.rds$", ".lock", basename(path)))
    if (!dir.create(lock, showWarnings = FALSE)) next
    cell <- cell_spec(task$key)
    start <- elapsed()
    record <- run_one(
      task,
      cell,
      ref_syn,
      ref_plas,
      comp_dir,
      task$rep <= pcre_reps
    )
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
if (length(args) < 4) {
  stop("usage: run-cl2-arms.R out_dir ref_syn ref_plas comp_dir [max_tasks]")
}
run_worker(
  args[1],
  args[2],
  args[3],
  args[4],
  if (length(args) >= 5) as.numeric(args[5]) else Inf
)
