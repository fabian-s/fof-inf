# Provenance and storage -------------------------------------------------------

# Software versions recorded in every result file (checked at start, see below).
software_versions <- function() {
  list(
    refund_sha = loaded_refund_sha(),
    mgcv_version = as.character(getNamespaceVersion("mgcv"))
  )
}

# The refund build actually loaded: SHA file next to the installed package (LRZ
# library), else env REFUND_SHA (local installs). Must match REFUND_SHA_EXPECTED.
loaded_refund_sha <- function() {
  sha_file <- file.path(dirname(find.package("refund")), "REFUND_SHA")
  sha <- if (file.exists(sha_file)) trimws(readLines(sha_file, 1)) else
    Sys.getenv("REFUND_SHA")
  sha
}

check_software <- function(software) {
  sha <- software$refund_sha
  if (nchar(sha) < 7 || !startsWith(REFUND_SHA_EXPECTED, sha)) {
    stop("Loaded refund build '", sha, "' is not ", REFUND_SHA_EXPECTED)
  }
  if (package_version(software$mgcv_version) != package_version(MGCV_VERSION)) {
    stop("Loaded mgcv ", software$mgcv_version, "; expected ", MGCV_VERSION)
  }
  if (!exists("pffr_vcov", asNamespace("refund")))
    stop("refund lacks pffr_vcov()")
}

atomic_save <- function(record, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- tempfile(pattern = ".task-", tmpdir = dirname(path))
  on.exit(unlink(tmp))
  saveRDS(record, tmp)
  if (!file.rename(tmp, path)) stop("Could not atomically save ", path)
  invisible(path)
}

# One dataset, all its fits -----------------------------------------------------

# Fit one method and extract everything scoring needs from it (plus, on request, the
# fit's concurvity).
run_method <- function(sim, cell, method, concurvity = FALSE) {
  fit <- fit_model(sim, cell, method)
  inputs <- method_inputs(fit, method, sim, cell)
  if (concurvity) inputs$concurvity <- check_concurvity(fit, cell)
  inputs
}

# A saved task counts as done only if it succeeded completely and is the same task.
saved_complete <- function(path, task) {
  if (!file.exists(path)) return(FALSE)
  saved <- readRDS(path)
  if (!identical(saved$status, "success")) return(FALSE)
  if (saved$cell != task$cell || saved$rep != task$rep) {
    stop(
      "Saved task ",
      task$task,
      " has cell/rep ",
      saved$cell,
      "/",
      saved$rep,
      ", but the task table says ",
      task$cell,
      "/",
      task$rep
    )
  }
  length(saved$extra_failures) == 0
}

requested_methods <- function(cell) {
  extras <- strsplit(trimws(cell$extra), " +")[[1]]
  c(
    "REML",
    "NCV",
    c(
      ar1 = "AR1",
      pointwise = "NCV_point",
      block = "NCV_block",
      block_alt = "NCV_alt"
    )[extras[extras != ""]]
  )
}

run_task <- function(
  task_id,
  out_dir,
  cells = design_table(),
  tasks = task_table(cells),
  simulate_fun = simulate # the plasmode applications pass their own generator
) {
  task <- tasks[tasks$task == task_id, , drop = FALSE]
  if (nrow(task) != 1) stop("Task id must identify exactly one row: ", task_id)
  cell <- cells[cells$cell == task$cell, , drop = FALSE]
  if (nrow(cell) != 1) stop("Unknown or duplicated cell: ", task$cell)
  # Files are named by cell and replicate, so they stay valid when cells are added.
  path <- file.path(
    out_dir,
    sprintf("task-c%03d-r%03d.rds", task$cell, task$rep)
  )
  software <- software_versions()
  stage <- "version"
  record <- tryCatch(
    {
      check_software(software)
      if (saved_complete(path, task)) return(invisible(path))
      stage <- "simulate"
      sim <- simulate_fun(cell, task$rep)
      inputs <- list()
      extra_failures <- list()
      for (method in requested_methods(cell)) {
        stage <- paste0("fit:", method)
        if (method %in% c("REML", "NCV")) {
          # cheap diagnostic: only the REML fit of a few replicates per cell
          check <- method == "REML" && task$rep %in% CONCURVITY_REPS
          inputs[[method]] <- run_method(sim, cell, method, check)
          next
        }
        # A failing extra fit must not drop the replicate's REML and NCV results.
        inputs[[method]] <- tryCatch(
          run_method(sim, cell, method),
          error = \(e) {
            extra_failures[[method]] <<- conditionMessage(e)
            NULL
          }
        )
      }
      stage <- "score"
      c(
        list(
          status = "success",
          extra_failures = extra_failures,
          cell_spec = cell,
          concurvity = inputs$REML$concurvity # NULL except for CONCURVITY_REPS
        ),
        score_task(inputs, sim, cell, task)
      )
    },
    error = function(e) {
      # A failure must never replace an existing successful task.
      if (file.exists(path) && identical(readRDS(path)$status, "success"))
        stop(e)
      list(
        status = "failure",
        task = task_id,
        cell = task$cell,
        rep = task$rep,
        failure = data.frame(
          task = task_id,
          cell = task$cell,
          rep = task$rep,
          stage = stage,
          message = conditionMessage(e)
        )
      )
    }
  )
  atomic_save(c(record, list(software = software)), path)
}


# Array chunks use the full task table, so run 2 retains the original task ids.
array_task_ids <- function(tasks, array_id, chunk_size, run) {
  args <- c(array_id, chunk_size, run)
  if (
    anyNA(args) ||
      any(!is.finite(args)) ||
      any(args < 1 | args != floor(args)) ||
      !run %in% c(1, 2)
  )
    stop("Array id/chunk size must be positive integers; run must be 1 or 2.")
  selected <- tasks[
    tasks$rep > (run - 1) * N_REP / 2 & tasks$rep <= run * N_REP / 2,
  ]
  chunk <- (seq_len(nrow(selected)) - 1) %/% chunk_size + 1
  selected$task[chunk == array_id]
}
