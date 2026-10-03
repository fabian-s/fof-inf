testthat::test_that("scoring packs aligned numeric records and keeps response truth", {
  sim <- fixture_sim(family = "poisson")
  cell <- fixture_cell(family = "poisson")
  methods <- c("REML", "NCV", "AR1", "NCV_point")
  inputs <- fixture_inputs(sim, cell, methods)
  task <- data.frame(task = 23, cell = 17, rep = 2)
  result <- score_task(inputs, sim, cell, task)
  expected <- c(
    "truth",
    "est_REML",
    "se_model_REML",
    "se_cl2_REML",
    "est_NCV",
    "se_model_NCV",
    "se_cl2_NCV",
    "est_AR1",
    "se_model_AR1",
    "est_NCV_point",
    "se_model_NCV_point"
  )
  for (e in names(sim$truth)) {
    testthat::expect_true(is.numeric(result$records[[e]]))
    testthat::expect_identical(colnames(result$records[[e]]), expected)
    testthat::expect_equal(result$records[[e]][, "truth"], sim$truth[[e]])
    testthat::expect_equal(
      result$records[[e]][, "est_REML"],
      inputs$REML$estimands[[e]]$est
    )
  }
  testthat::expect_identical(result$task, task$task)
  testthat::expect_equal(result$fits$n_warnings, rep(1, length(methods)))
  testthat::expect_identical(result$warnings$REML, "fixture warning")
  broken <- inputs
  broken$NCV$estimands$gamma$se_cl2 <- NULL
  testthat::expect_error(
    score_task(broken, sim, cell, task),
    "CL2 standard errors missing"
  )
  sim$truth$beta <- 0
  testthat::expect_error(
    score_task(inputs, sim, cell, task),
    "Truth does not match"
  )
})

testthat::test_that("M2 scores only its defined estimands", {
  sim <- fixture_sim("M2")
  cell <- fixture_cell("M2")
  result <- score_task(
    fixture_inputs(sim, cell),
    sim,
    cell,
    data.frame(task = 1)
  )
  testthat::expect_named(result$records, c("mean", "f", "gamma"))
  testthat::expect_error(score_task(list(), sim, cell, NULL), "REML and NCV")
})
