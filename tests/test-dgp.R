design_effect <- function(R) sum(R) / sum(diag(R))
t61 <- t_grid(61)

test_that("kernels have unit diagonal, are positive definite and hit their DE targets", {
  targets <- c(ou = 8, ou0 = 8, fpc = 11.8) # fpc: frozen DTI kernel, measured
  for (error in names(targets)) {
    R <- correlation_kernel(error, t61)
    expect_equal(diag(R), rep(1, 61), tolerance = 1e-12)
    expect_gt(min(eigen(R, only.values = TRUE)$values), 0)
    expect_equal(design_effect(R), targets[[error]], tolerance = 0.01)
  }
  R <- correlation_kernel("osc", t61) # sign-changing correlations
  expect_equal(diag(R), rep(1, 61), tolerance = 1e-12)
  expect_gt(min(eigen(R, only.values = TRUE)$values), 0)
  expect_lt(min(R), -0.4)
  # the dense-grid kernel is positive definite too (used for every draw)
  expect_gt(min(eigen(correlation_kernel("ou0", T_DENSE))$values), 0)
})

test_that("copula draws have the family's marginal mean (and variance for scat)", {
  set.seed(1)
  u <- runif(2e5)
  eta <- c(
    gaussian = 0.3,
    scat = 0.3,
    poisson = log(3),
    binomial = qlogis(0.2),
    betar = qlogis(0.3),
    nb = log(3)
  )
  for (fam in names(eta)) {
    y <- FAMILIES[[fam]]$draw(qnorm(u), eta[[fam]], 1.5)
    mu <- FAMILIES[[fam]]$linkinv(eta[[fam]])
    expect_equal(mean(y), mu, tolerance = 0.01 * max(1, mu), label = fam)
  }
  expect_equal(
    var(FAMILIES$scat$draw(qnorm(u), 0, 1.5)),
    1.5^2,
    tolerance = 0.05
  )
})

test_that("truths have the stated signal split and GLM levels", {
  pop <- population_sample()
  for (model in c("M1", "M2")) {
    truth <- make_truth(model, "smooth")
    main <- main_contribution(truth, pop, T_DENSE)
    zg <- outer(pop$z, truth$gamma(T_DENSE))
    expect_equal(mean(apply(main, 2, var)), SIGNAL_SHARE_MAIN, tolerance = 0.1)
    expect_equal(mean(apply(main + zg, 2, var)), 1, tolerance = 0.1)
  }
  cell <- as.list(design_table()[1, ])
  cell$family <- "binomial"
  truth <- make_truth("M1", "smooth")
  c0 <- family_intercept(cell, truth, 1)
  eta <- linear_predictor(truth, pop, t61, 1, c0)
  expect_equal(mean(plogis(eta)), MEAN_LEVEL[["binomial"]], tolerance = 1e-6)
  expect_lt(max(abs(eta - c0)), 8) # bounded linear predictor at the high GLM level
})

test_that("warps are monotone and targets match brute-force Monte Carlo", {
  a <- c(-WARP_MAX, 0, WARP_MAX)
  h <- warp_times(a, T_DENSE)
  expect_true(all(diff(t(h)) > 0))
  expect_equal(h[, c(1, 241)], cbind(c(0, 0, 0), c(1, 1, 1)))

  cell <- as.list(design_table()[design_table()$block == "warp", ][1, ])
  truth <- make_truth("M1", "smooth")
  sc <- signal_scale(cell)
  target <- estimand_targets(cell, truth, sc$scale, 0, test_cohort())$mean
  test <- test_cohort()
  one <- lapply(test, \(v) if (is.matrix(v)) v[1, , drop = FALSE] else v[1])
  set.seed(2)
  times <- warp_times(warp_amplitude(runif(2e4)), EVAL_T)
  mc <- colMeans(linear_predictor(
    truth,
    lapply(one, \(v) {
      if (is.matrix(v)) v[rep(1, nrow(times)), ] else rep(v, nrow(times))
    }),
    times,
    sc$scale,
    0
  ))
  first_curve <- matrix(target, N_TEST)[1, ]
  expect_lt(max(abs(mc - first_curve)), 0.02)
})

test_that("common random numbers: G and D subsets reuse the same draws", {
  d <- design_table()
  big <- as.list(d[
    d$block == "core" &
      d$family == "gaussian" &
      d$error == "fpc" &
      d$G == 100 &
      d$signal == "mid",
  ])
  small <- modifyList(big, list(G = 40L))
  dense <- modifyList(big, list(D = 241L))
  sb <- simulate(big, 3)
  expect_equal(simulate(small, 3)$data$Y, sb$data$Y[1:40, ])
  expect_equal(simulate(dense, 3)$data$Y[, T_SUBSAMPLE[["61"]]], sb$data$Y)
})

test_that("truths are representable: a near-noise-free fit recovers them", {
  skip_if(Sys.getenv("RESTART_FAST") == "true")
  d <- design_table()
  cell <- as.list(d[
    d$block == "core" &
      d$family == "gaussian" &
      d$error == "iid" &
      d$G == 100 &
      d$signal == "high",
  ])
  sim <- simulate(cell, 1)
  truth <- make_truth("M1", "smooth")
  cov <- list(X = unclass(sim$data$X), z = sim$data$z)
  sim$data$Y <- I(
    linear_predictor(truth, cov, sim$yind, 1, 0) +
      matrix(rnorm(length(sim$data$Y), sd = 1e-3), nrow(sim$data$Y))
  )
  fit <- refund::pffr(
    Y ~ ff(X, xind = S_GRID, splinepars = BASIS_FF) + z,
    data = sim$data,
    yind = sim$yind,
    bs.yindex = BASIS_T,
    bs.int = BASIS_T,
    sandwich = "none"
  )
  est <- coef(fit, n1 = length(EVAL_T), n2 = length(EVAL_S))$smterms
  gamma_hat <- est[["z(yindex)"]]$coef
  expect_lt(max(abs(gamma_hat$value - truth$gamma(gamma_hat$yindex.vec))), 0.01)
  ff_hat <- est[["ff(X)"]]$coef
  beta_true <- diag(truth$beta(ff_hat$X.smat, ff_hat$X.tmat))
  expect_lt(max(abs(ff_hat$value - beta_true)) / max(abs(beta_true)), 0.02)
})

test_that("heteroskedastic processes keep the average error variance at 1", {
  R <- correlation_kernel("fpc_raw", t61)
  expect_equal(mean(diag(R)), 1, tolerance = 1e-12)
  expect_gt(max(diag(R)) / min(diag(R)), 3) # the DTI variance profile along t
  set.seed(4)
  expect_equal(
    mean(subject_scale("fpc_het", rnorm(1e6))^2),
    1,
    tolerance = 0.01
  )
  expect_identical(subject_scale("fpc", c(-1, 0, 2)), c(1, 1, 1))
  cell <- modifyList(
    as.list(design_table()[1, ]),
    list(family = "poisson", error = "fpc_het", signal = "mid")
  )
  expect_error(simulate(cell, 1), "Gaussian responses only")
})
