# Lambda-correction workers: one task and the lock-dir worker loop (shared by run-lc.R
# and run-lc-plasmode.R). Sourced after scripts/setup.R and lc-lib.R.

lc_task <- function(
  task,
  cell,
  ref_path,
  simulate_fun = simulate,
  jackknife = TRUE
) {
  stage <- "simulate"
  tryCatch(
    {
      sim <- simulate_fun(cell, task$rep)
      stage <- "fit"
      fit <- fit_model(sim, cell, "REML")
      stage <- "covariances"
      cov <- collect_warnings(lc_covariances(fit, jackknife = jackknife))
      stage <- "assemble"
      records <- lapply(
        setNames(nm = model_estimands(cell$model)),
        \(estimand) {
          L <- evaluation_matrix(fit, sim, estimand)
          se <- lapply(cov$value$V, \(V) sqrt(pmax(rowSums((L %*% V) * L), 0)))
          do.call(
            cbind,
            c(
              list(
                truth = sim$truth[[estimand]],
                est_REML = drop(L %*% fit$coefficients)
              ),
              setNames(se, paste0("se_", names(se)))
            )
          )
        }
      )
      list(
        status = "success",
        records = records,
        crn = crn_check_lc(records, ref_path, task),
        omega = cov$value$omega,
        time_parts = cov$value$time,
        fit_time = fit$restart$time,
        n_jk_nonconv = cov$value$n_jk_nonconv,
        converged = isTRUE(fit$converged) &&
          fit$outer.info$conv == "full convergence",
        warnings = c(fit$restart$warnings, cov$warnings)
      )
    },
    error = \(e)
      list(status = "failure", stage = stage, message = conditionMessage(e))
  )
}

# CRN check on truth, estimate and CL2 SE only.
crn_check_lc <- function(records, ref_path, task) {
  if (!file.exists(ref_path)) return(NULL)
  ref <- readRDS(ref_path)
  if (ref$cell != task$cell || ref$rep != task$rep)
    stop("reference record mismatch: ", ref_path)
  cols <- c(truth = "truth", est_REML = "est_REML", se_cl2 = "se_cl2_REML")
  do.call(
    rbind,
    lapply(names(records), \(estimand) {
      do.call(
        rbind,
        lapply(names(cols), \(col) {
          x <- records[[estimand]][, col]
          y <- ref$records[[estimand]][, cols[[col]]]
          data.frame(
            estimand = estimand,
            column = col,
            max_abs = max(abs(x - y)),
            max_rel = max(abs(x - y) / pmax(abs(y), 1e-8))
          )
        })
      )
    })
  )
}

task_path <- function(out_dir, task) {
  file.path(out_dir, sprintf("task-c%03d-r%03d.rds", task$cell, task$rep))
}

# Walk the (rep-major) task table, claiming tasks through lock directories. ref_path(task)
# gives the study record for the CRN check.
run_worker <- function(
  out_dir,
  cells,
  tasks,
  ref_path,
  simulate_fun = simulate,
  jackknife = TRUE,
  max_tasks = Inf,
  task_fun = NULL
) {
  task_fun <- task_fun %||%
    \(task, cell, ref) lc_task(task, cell, ref, simulate_fun, jackknife)
  software <- software_versions()
  check_software(software)
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
    record <- task_fun(task, cells[cells$cell == task$cell, ], ref_path(task))
    record <- c(
      record,
      list(task = task$task, cell = task$cell, rep = task$rep)
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
