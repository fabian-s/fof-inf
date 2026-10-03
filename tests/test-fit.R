testthat::test_that("evaluation matrices reproduce refund coefficients and predictions", {
  for (family in c("gaussian", "poisson")) {
    sim <- fixture_sim(family = family)
    cell <- fixture_cell(family = family)
    fits <- lapply(c("REML", "NCV"), \(m) fit_model(sim, cell, m)) |>
      setNames(c("REML", "NCV"))
    inputs <- interval_inputs(fits, sim, cell)
    for (method in names(fits)) {
      fit <- fits[[method]]
      input <- inputs[[method]]
      grids <- lapply(c("alpha", "beta", "gamma"), \(e) {
        sm <- estimand_smooth(fit, e)
        estimand_grid(sm, e)
      })
      co <- coef(fit, eval_grid = grids, sandwich = "none", seWithMean = FALSE)
      for (i in seq_along(grids)) {
        e <- c("alpha", "beta", "gamma")[i]
        expected <- as.vector(co$smterms[[i]]$value)
        if (e == "alpha") expected <- expected + co$pterms[1, "value"]
        testthat::expect_equal(
          unname(input$estimands[[e]]$est),
          expected,
          tolerance = 1e-10
        )
      }
      testthat::expect_equal(
        unname(input$estimands$mean$est),
        as.vector(predict(fit, sim$test)),
        tolerance = 1e-10
      )
      # refund's model-based SEs use Vc, which is invalid for NCV fits (see fit.R)
      if (method == "REML") {
        testthat::expect_equal(
          input$estimands$beta$se_model,
          as.vector(co$smterms[[2]]$se),
          tolerance = 1e-10
        )
      }
      testthat::expect_identical(
        input$fit_row$model_covariance,
        if (method == "NCV" || is.null(fit$Vc)) "Vp" else "Vc"
      )
      testthat::expect_true(all(is.finite(input$estimands$beta$se_cl2)))
      testthat::expect_identical(fit$pffr$sandwich, "none")
    }
    scored <- score_task(
      inputs,
      sim,
      cell,
      data.frame(task = 1, cell = 17, rep = 1)
    )
    testthat::expect_equal(nrow(scored$records$mean), N_TEST * length(EVAL_T))
  }
})

testthat::test_that("M2 keeps the fixed x knots and constrained basis", {
  for (family in c("gaussian", "poisson")) {
    sim <- fixture_sim("M2", family)
    fit <- fit_model(sim, fixture_cell("M2", family), "REML")
    sm <- estimand_smooth(fit, "f")
    margin <- sm$margin[[1]]
    testthat::expect_equal(margin$knots, x_knots())
    # Remove the sample centring constraint to compare the underlying [0,1] basis.
    testthat::expect_equal(margin$m, BASIS_X$m)
    testthat::expect_equal(sm$margin[[2]]$m, BASIS_T$m)
    attr(margin, "qrc") <- NULL
    attr(margin, "nCons") <- 0
    testthat::expect_equal(
      unname(mgcv::PredictMat(margin, data.frame(x = EVAL_X))),
      unname(ps_eval(ps_basis(BASIS_X), EVAL_X)),
      tolerance = 1e-10
    )
    grids <- list(NULL, estimand_grid(sm, "f"), NULL)
    co <- coef(fit, eval_grid = grids, se = FALSE)
    L <- evaluation_matrix(fit, sim, "f")
    testthat::expect_equal(
      drop(L %*% fit$coefficients),
      as.vector(co$smterms[[2]]$value),
      tolerance = 1e-10
    )
    p <- expand.grid(x = sim$data$x, t = EVAL_T)
    names(p) <- sm$term
    values <- mgcv::PredictMat(sm, p) %*%
      fit$coefficients[sm$first.para:sm$last.para]
    testthat::expect_equal(
      colSums(matrix(values, nrow(sim$data))),
      rep(0, length(EVAL_T)),
      tolerance = 1e-8
    )
  }
})

testthat::test_that("AR1 profiling recovers strong correlation and omits CL2", {
  sim <- fixture_sim(rho = 0.9)
  cell <- fixture_cell()
  fit <- fit_model(sim, cell, "AR1")
  testthat::expect_true(fit$AR1.rho >= 0 && fit$AR1.rho <= 0.99)
  testthat::expect_gt(fit$AR1.rho, 0.75)
  testthat::expect_identical(
    fit$rho_boundary,
    fit$AR1.rho %in% range(AR_RHO_GRID)
  )
  testthat::expect_equal(sum(fit$model[["(AR.start)"]]), nrow(sim$data))
  input <- interval_inputs(list(AR1 = fit), sim, cell)$AR1
  testthat::expect_null(input$cl2_adjustment)
  testthat::expect_null(input$estimands$beta$se_cl2)
  testthat::expect_true(all(is.finite(input$estimands$mean$se_model)))
  testthat::expect_error(
    fit_model(sim, fixture_cell(family = "poisson"), "AR1"),
    "Gaussian"
  )
})

testthat::test_that("pointwise NCV runs without a CL2 covariance", {
  sim <- fixture_sim()
  fit <- fit_model(sim, fixture_cell(), "NCV_point")
  input <- interval_inputs(list(NCV_point = fit), sim, fixture_cell())$NCV_point
  testthat::expect_null(input$cl2_adjustment)
  testthat::expect_true(all(is.finite(input$estimands$beta$est)))
})

testthat::test_that("covariance warnings and optimizer failures enter diagnostics", {
  sim <- fixture_sim("M2")
  cell <- fixture_cell("M2")
  fit <- fit_model(sim, cell, "REML")
  original <- refund:::pffr_vcov
  testthat::local_mocked_bindings(
    pffr_vcov = function(object, sandwich, ...) {
      if (sandwich == "cl2") warning("covariance test warning")
      original(object, sandwich = sandwich, ...)
    },
    .package = "refund"
  )
  input <- interval_inputs(list(REML = fit), sim, cell)$REML
  testthat::expect_true(
    "covariance test warning" %in% input$warnings
  )
  fit$outer.info$conv <- "iteration limit reached"
  testthat::expect_false(
    method_inputs(fit, "REML", sim, cell)$fit_row$converged
  )
})

testthat::test_that("concurvity is mgcv's worst case per smooth, warns and never fails", {
  sim <- fixture_sim()
  cell <- fixture_cell()
  fit <- fit_model(sim, cell, "REML")
  conc <- testthat::expect_silent(check_concurvity(fit, cell))
  labels <- vapply(fit$smooth, `[[`, "", "label", USE.NAMES = FALSE)
  testthat::expect_identical(conc$term, labels)
  worst <- suppressWarnings(mgcv::concurvity(fit, full = TRUE))["worst", ]
  testthat::expect_equal(conc$concurvity, unname(worst[labels]))
  testthat::expect_true(all(conc$concurvity >= 0 & conc$concurvity <= 1))
  flagged <- collect_warnings(check_concurvity(fit, cell, threshold = 0))
  testthat::expect_length(flagged$warnings, length(labels))
  testthat::expect_match(
    flagged$warnings[1],
    "cell 17: worst-case concurvity .* of s\\(yindex"
  )
  broken <- testthat::expect_silent(check_concurvity("not a fit", cell))
  testthat::expect_true(is.na(broken$concurvity))
})
