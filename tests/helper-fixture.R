# Small data independent of dgp.R ----------------------------------------------

fixture_sim <- function(model = "M1", family = "gaussian", rho = 0, seed = 42) {
  withr::local_seed(seed)
  G <- 20L
  t <- seq(0, 1, length.out = 31) # smallest grid containing the production EVAL_T
  covariates <- function(n) {
    z <- rnorm(n)
    if (model == "M2") return(data.frame(x = runif(n, 0.15, 0.85), z = z))
    scores <- matrix(rnorm(n * 12), n) |> sweep(2, seq_len(12), "/")
    basis <- cos(outer(0:11, S_GRID * pi))
    data.frame(X = I(scores %*% basis), z = z)
  }
  data <- covariates(G)
  test <- covariates(N_TEST)
  eta <- function(d) {
    main <- if (model == "M1") {
      (rowSums(d$X) - (d$X[, 1] + d$X[, ncol(d$X)]) / 2) / (ncol(d$X) - 1)
    } else d$x - mean(data$x)
    outer(main, 0.3 * sin(pi * t)) + outer(d$z, 0.2 * cos(pi * t)) + 1
  }
  errors <- matrix(rnorm(G * length(t)), G) %*%
    chol(rho^abs(outer(seq_along(t), seq_along(t), "-")))
  data$Y <- I(
    if (family == "poisson") {
      matrix(rpois(length(errors), exp(eta(data))), G)
    } else eta(data) + errors
  )
  truth <- lapply(model_estimands(model), \(e) rep(0, nrow(eval_points(e)))) |>
    setNames(model_estimands(model))
  truth$mean <- as.vector(FAMILIES[[family]]$linkinv(eta(test)))
  truth$gamma <- 0.2 * cos(pi * EVAL_T)
  if (model == "M1") {
    truth$alpha <- rep(1, length(EVAL_T))
    truth$beta <- 0.3 * sin(pi * eval_points("beta")$t)
  } else {
    p <- eval_points("f")
    truth$f <- 0.3 * (p$x - mean(data$x)) * sin(pi * p$t)
  }
  list(
    data = data,
    yind = t,
    xind = S_GRID,
    test = test,
    truth = truth,
    diag = list()
  )
}

# Synthetic inputs isolate packing and orchestration tests from fitting.
fixture_inputs <- function(sim, cell, methods = c("REML", "NCV")) {
  lapply(methods, \(method) {
    estimands <- lapply(
      sim$truth,
      \(truth)
        list(
          est = seq_along(truth) / length(truth),
          se_model = rep(0.1, length(truth)),
          # as interval_inputs: CL2 only for REML and NCV
          se_cl2 = if (method %in% c("REML", "NCV")) rep(0.2, length(truth))
        )
    )
    list(
      estimands = estimands,
      fit_row = data.frame(
        fit = method,
        time = 0.1,
        cl2_time = 0.05,
        converged = TRUE,
        n_warnings = 1,
        criterion = 1,
        rho = NA_real_,
        rho_boundary = FALSE,
        model_covariance = "Vp"
      ),
      terms = data.frame(
        fit = method,
        term = "test",
        edf = 2,
        penalty = "test",
        log_sp = 0
      ),
      warnings = "fixture warning"
    )
  }) |>
    setNames(methods)
}

fixture_cell <- function(model = "M1", family = "gaussian") {
  cell <- as.data.frame(cell_defaults())
  cell$model <- model
  cell$family <- cell$fit_family <- family
  cell$G <- 20L
  cell$D <- 31L
  cell$cell <- 17L
  cell
}
