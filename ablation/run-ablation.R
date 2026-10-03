# Covariance ablation worker. Run from the repository root (code copy):
#   Rscript ablation/run-ablation.R <out_dir> <ref_dir> [max_tasks]
# Each worker walks the rep-major task table and claims tasks through atomic lock
# directories (<out_dir>/locks), so any number of workers can run or be restarted.
# ref_dir holds the study's task records (task-<id>.rds, 95-cell task table) against which
# the REML fit is checked (common random numbers); a missing record gives NA checks.

source("scripts/setup.R")
source("ablation/ablation-lib.R")

STUDY_N_CELLS <- 95 # task table of the study run whose records are in ref_dir

# One replicate -----------------------------------------------------------------------

# Estimates and SEs of every arm of one fit, per estimand.
fit_columns <- function(fit, label, arms, sim, cell) {
  covs <- lapply(seq_len(nrow(arms)), \(i) arm_covariance(fit, arms[i, ]))
  names(covs) <- arms$arm
  lapply(setNames(nm = model_estimands(cell$model)), \(estimand) {
    L <- evaluation_matrix(fit, sim, estimand)
    se <- lapply(covs, \(V) sqrt(rowSums((L %*% V) * L)))
    c(
      setNames(list(drop(L %*% fit$coefficients)), paste0("est_", label)),
      setNames(se, paste0("se_", names(se)))
    )
  })
}

fit_row <- function(fit, label) {
  data.frame(
    fit = label,
    time = fit$restart$time,
    converged = isTRUE(fit$converged) &&
      (is.null(fit$outer.info$conv) ||
        fit$outer.info$conv == "full convergence"),
    n_warnings = length(fit$restart$warnings),
    model_covariance = if (!is.null(fit$Vc)) "Vc" else "Vp"
  )
}

# Fixed-lambda refit with the cell's geometric-mean REML smoothing parameters.
fit_fixed <- function(sim, cell, sp) {
  start <- proc.time()[["elapsed"]]
  result <- collect_warnings(fit_pffr(sim, cell, "REML", sp = sp))
  fit <- result$value
  # with every sp fixed, mgcv leaves fit$sp empty and reports them in full.sp
  if (!isTRUE(all.equal(unname(fit$full.sp), unname(sp), tolerance = 1e-8))) {
    stop("fixed-lambda fit did not keep the supplied smoothing parameters")
  }
  fit$restart <- list(
    time = proc.time()[["elapsed"]] - start,
    warnings = result$warnings
  )
  fit
}

# Max absolute and relative differences between the REML columns and the study record.
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

run_ablation_task <- function(task, cell, sp_table, ref_dir) {
  stage <- "simulate"
  failures <- list()
  out <- tryCatch(
    {
      sim <- simulate(cell, task$rep)
      stage <- "fit:REML"
      fit <- fit_model(sim, cell, "REML")
      stage <- "cov:REML"
      cols <- fit_columns(
        fit,
        "REML",
        COV_ARMS[COV_ARMS$fit == "REML", ],
        sim,
        cell
      )
      fits <- fit_row(fit, "REML")
      # A failing fixed-lambda refit must not drop the replicate's REML arms.
      fix <- tryCatch(
        {
          fit_fix <- fit_fixed(
            sim,
            cell,
            fixed_sp(sp_table, cell$cell, names(fit$sp))
          )
          list(
            cols = fit_columns(
              fit_fix,
              "FIX",
              COV_ARMS[COV_ARMS$fit == "FIX", ],
              sim,
              cell
            ),
            row = fit_row(fit_fix, "FIX")
          )
        },
        error = \(e) {
          failures[[length(failures) + 1]] <<- data.frame(
            stage = "FIX",
            message = conditionMessage(e)
          )
          NULL
        }
      )
      stage <- "assemble"
      records <- lapply(setNames(nm = names(cols)), \(estimand) {
        do.call(
          cbind,
          c(
            list(truth = sim$truth[[estimand]]),
            cols[[estimand]],
            fix$cols[[estimand]]
          )
        )
      })
      ref_path <- file.path(
        ref_dir,
        sprintf("task-%d.rds", (task$rep - 1) * STUDY_N_CELLS + task$cell)
      )
      list(
        status = "success",
        records = records,
        fits = rbind(fits, fix$row),
        crn = crn_check(records, ref_path, task),
        failures = do.call(rbind, failures)
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

# Worker loop -------------------------------------------------------------------------

task_path <- function(out_dir, task) {
  file.path(out_dir, sprintf("task-c%03d-r%03d.rds", task$cell, task$rep))
}

run_worker <- function(out_dir, ref_dir, max_tasks = Inf) {
  software <- software_versions()
  check_software(software)
  cells <- ablation_cells()
  tasks <- ablation_tasks(cells)
  sp_table <- read.csv("ablation/sp-fixed.csv")
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
    record <- run_ablation_task(
      task,
      cells[cells$cell == task$cell, ],
      sp_table,
      ref_dir
    )
    record$software <- software
    record$time <- proc.time()[["elapsed"]] - start
    atomic_save(record, path)
    cat(sprintf(
      "task %d (cell %d rep %d): %s, %.1f s\n",
      id,
      task$cell,
      task$rep,
      record$status,
      record$time
    ))
    done <- done + 1
  }
  invisible(done)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: run-ablation.R out_dir ref_dir [max_tasks]")
run_worker(
  args[1],
  args[2],
  if (length(args) >= 3) as.numeric(args[3]) else Inf
)
