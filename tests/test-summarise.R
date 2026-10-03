test_that("bias-aware interval adds the NCV-REML difference in quadrature", {
  r <- cbind(
    truth = 0,
    est_REML = 1,
    se_model_REML = 1,
    se_cl2_REML = 1,
    est_NCV = 4,
    se_model_NCV = 1,
    se_cl2_NCV = 3
  )
  s <- point_scores(r, "ncv_cl2_bias", "beta", identity)
  z <- qnorm(1 - (1 - LEVEL) / 2)
  # delta = est_NCV - est_REML = 3, se = sqrt(3^2 + 3^2)
  expect_equal(unname(s$width), 2 * z * sqrt(18))
})

test_that("mean intervals are transformed to the response scale", {
  r <- cbind(
    truth = exp(0.1),
    est_REML = 0,
    se_model_REML = 0.1,
    se_cl2_REML = 0.1,
    est_NCV = 0,
    se_model_NCV = 0.1,
    se_cl2_NCV = 0.1
  )
  s <- point_scores(r, "reml_cl2", "mean", exp)
  z <- qnorm(1 - (1 - LEVEL) / 2)
  expect_equal(unname(s$width), exp(0.1 * z) - exp(-0.1 * z))
  expect_true(s$covered)
})

test_that("squared bias is unbiased and calibration error discounts MC noise", {
  set.seed(1)
  errors <- matrix(rnorm(200 * 400, mean = 0), 200) # zero bias
  expect_lt(abs(unbiased_bias2(errors)), 0.002)
  expect_gt(mean(rowMeans(errors)^2), 0.002) # the naive version is biased up
  expect_equal(calibration_error(c(0.95, 0.95), c(0.01, 0.01)), 0)
  expect_equal(calibration_error(c(0.90, 0.90), c(0, 0)), 0.05)
  expect_equal(calibration_error(c(0.99, 0.99), c(0, 0)), 0) # overcoverage does not count
  expect_equal(calibration_error(c(0.90, 0.99), c(0, 0)), sqrt(0.05^2 / 2))
})

test_that("the recommendation needs its complete pool and applies the rule", {
  cells <- design_table()
  pool <- recommendation_pool(cells, "gaussian")
  rows <- \(arm, coverage, score) {
    tibble::tibble(
      cell = pool,
      estimand = "beta",
      arm = arm,
      coverage_mean = coverage,
      coverage_se = 0.01,
      score_mean = score
    )
  }
  full <- dplyr::bind_rows(rows(PROPOSAL, 0.95, 0.9), rows(FALLBACK, 0.95, 1))
  expect_identical(
    recommend(full, cells, "gaussian", "beta")$decision,
    "proposal"
  )
  worse <- dplyr::bind_rows(rows(PROPOSAL, 0.90, 0.9), rows(FALLBACK, 0.95, 1))
  expect_identical(
    recommend(worse, cells, "gaussian", "beta")$decision,
    "fallback"
  )
  expect_identical(
    recommend(full[-1, ], cells, "gaussian", "beta")$decision,
    "incomplete pool"
  )
})
