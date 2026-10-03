testthat::test_that("tasks are looked up by id, saved atomically and skipped", {
  withr::local_dir("..")
  withr::local_envvar(REFUND_SHA = substr(REFUND_SHA_EXPECTED, 1, 8))
  env <- new.env(parent = environment())
  sys.source("R/run.R", env)
  cell <- fixture_cell()
  cell$extra <- "pointwise ar1"
  tasks <- data.frame(task = c(83, 7), cell = 17, rep = c(2, 1))
  seen <- character()
  env$simulate <- function(cell, rep) {
    testthat::expect_equal(rep, 1)
    fixture_sim()
  }
  env$fit_model <- function(sim, cell, method) {
    seen <<- c(seen, method)
    method
  }
  env$method_inputs <- function(fit, method, sim, cell) {
    force(fit) # evaluates the (mocked) fit, as the real method_inputs does
    fixture_inputs(sim, cell, method)[[method]]
  }
  out <- withr::local_tempdir()
  path <- env$run_task(7, out, cell, tasks)
  record <- readRDS(path)
  testthat::expect_identical(record$status, "success")
  testthat::expect_equal(record$task, 7)
  testthat::expect_equal(record$rep, 1)
  testthat::expect_setequal(seen, c("REML", "NCV", "AR1", "NCV_point"))
  # rep 1 runs the concurvity check; on the mocked fit it records NA, not a failure
  testthat::expect_true(1 %in% CONCURVITY_REPS)
  testthat::expect_true(is.na(record$concurvity$concurvity))
  testthat::expect_length(
    list.files(out, recursive = TRUE, all.files = TRUE),
    1
  )
  before <- tools::md5sum(path)
  env$fit_model <- function(...) stop("must skip")
  testthat::expect_identical(env$run_task(7, out, cell, tasks), path)
  testthat::expect_identical(tools::md5sum(path), before)
  env$MGCV_VERSION <- "0.0.0"
  testthat::expect_error(env$run_task(7, out, cell, tasks), "expected 0.0.0")
  env$MGCV_VERSION <- MGCV_VERSION
  testthat::expect_identical(tools::md5sum(path), before)
  testthat::expect_named(record$software, c("refund_sha", "mgcv_version"))
  testthat::expect_match(basename(path), "^task-c017-r001\\.rds$")
})

testthat::test_that("a failed fit is saved with task, stage and message", {
  withr::local_dir("..")
  withr::local_envvar(REFUND_SHA = substr(REFUND_SHA_EXPECTED, 1, 8))
  env <- new.env(parent = environment())
  sys.source("R/run.R", env)
  env$simulate <- function(...) fixture_sim()
  env$fit_model <- function(...) stop("deliberate fit error")
  cell <- fixture_cell()
  tasks <- data.frame(task = 91:93, cell = 17, rep = 1:3)
  out <- withr::local_tempdir()
  failure <- readRDS(env$run_task(91, out, cell, tasks))
  testthat::expect_identical(failure$status, "failure")
  testthat::expect_identical(failure$failure$stage, "fit:REML")
  testthat::expect_identical(failure$failure$message, "deliberate fit error")
  testthat::expect_equal(failure$failure$task, 91)
  testthat::expect_length(
    list.files(out, recursive = TRUE, all.files = TRUE),
    1
  )
  env$MGCV_VERSION <- "0.0.0"
  mismatch <- readRDS(env$run_task(92, out, cell, tasks))
  testthat::expect_identical(mismatch$failure$stage, "version")
  testthat::expect_match(mismatch$failure$message, "expected 0.0.0")
  env$MGCV_VERSION <- MGCV_VERSION
  withr::local_envvar(REFUND_SHA = "0000000")
  wrong_sha <- readRDS(env$run_task(93, out, cell, tasks))
  testthat::expect_identical(wrong_sha$failure$stage, "version")
  testthat::expect_match(wrong_sha$failure$message, "is not")
})

testthat::test_that("array chunks preserve ids and split the two replicate runs", {
  tasks <- task_table(data.frame(cell = c(19, 3)))
  first <- array_task_ids(tasks, 1, 3, 1)
  second <- array_task_ids(tasks, 2, 3, 1)
  testthat::expect_equal(first, 1:3)
  testthat::expect_equal(second, 4:6)
  testthat::expect_equal(array_task_ids(tasks, 1, 3, 2), 201:203)
  testthat::expect_length(array_task_ids(tasks, 1000, 3, 1), 0)
  testthat::expect_error(array_task_ids(tasks, NA, 3, 1), "positive integers")
  testthat::expect_error(array_task_ids(tasks, 1, 0, 1), "positive integers")
  testthat::expect_error(array_task_ids(tasks, 1, 3, 3), "run must be")
})

testthat::test_that("a failing extra fit keeps the task; saved failures are retried", {
  withr::local_dir("..")
  withr::local_envvar(REFUND_SHA = substr(REFUND_SHA_EXPECTED, 1, 8))
  env <- new.env(parent = environment())
  sys.source("R/run.R", env)
  cell <- fixture_cell()
  cell$extra <- "ar1"
  tasks <- data.frame(task = 5, cell = 17, rep = 1)
  env$simulate <- function(...) fixture_sim()
  env$method_inputs <- function(fit, method, sim, cell) {
    force(fit) # evaluates the (mocked) fit, as the real method_inputs does
    fixture_inputs(sim, cell, method)[[method]]
  }
  env$fit_model <- function(sim, cell, method) {
    if (method == "AR1") stop("ar1 broke")
    method
  }
  out <- withr::local_tempdir()
  record <- readRDS(env$run_task(5, out, cell, tasks))
  testthat::expect_identical(record$status, "success")
  testthat::expect_identical(record$extra_failures$AR1, "ar1 broke")

  env$fit_model <- function(...) stop("transient")
  tasks$task <- 6
  tasks$rep <- 2 # a different replicate, hence a different file
  tasks$rep <- 2 # a different replicate, hence a different file
  testthat::expect_identical(
    readRDS(env$run_task(6, out, cell, tasks))$status,
    "failure"
  )
  env$fit_model <- function(sim, cell, method) method
  testthat::expect_identical(
    readRDS(env$run_task(6, out, cell, tasks))$status,
    "success"
  )
})
