# Data generation: covariates, error processes, warps, simulate(), estimand targets.

# Covariates --------------------------------------------------------------------------

# SDs of the covariate's basis coefficients for a cell's `covariate` level.
covariate_sd <- function(covariate = "rich") {
  switch(
    covariate,
    rich = X_SD,
    lowrank = X_SD_LOWRANK,
    stop("unknown covariate: ", covariate)
  )
}

# Covariates from standard-normal/uniform draws: X on S_GRID (M1), x (M2), z.
make_covariates <- function(xi, z, x, x_sd = X_SD) {
  list(
    X = (xi * rep(x_sd, each = nrow(xi))) %*% t(covariate_basis(S_GRID)),
    z = z,
    x = x
  )
}

draw_covariates <- function(n, x_sd = X_SD) {
  make_covariates(matrix(rnorm(n * X_N_TERMS), n), rnorm(n), runif(n), x_sd)
}

# The frozen test cohort (same in every cell and replicate); RNG state untouched.
test_cohort <- function(x_sd = X_SD) {
  withr::with_seed(SEED_TEST, draw_covariates(N_TEST, x_sd))
}

# Large frozen population sample, used only to set GLM intercepts.
population_sample <- function(x_sd = X_SD) {
  withr::with_seed(SEED_POPULATION, draw_covariates(N_POPULATION, x_sd))
}

# Signal level and intercept ------------------------------------------------------------

# sqrt(S) and error SD of a cell. Gaussian/scat: S = 1 and sigma from R^2;
# other families: sqrt(S) from the signal level, no additive error.
signal_scale <- function(cell) {
  level <- SIGNAL_LEVELS[[cell$family]][[cell$signal]]
  if (cell$family %in% c("gaussian", "scat")) {
    return(list(scale = 1, sigma = sqrt((1 - level) / level)))
  }
  list(scale = level, sigma = NA_real_)
}

# Intercept shift that gives the family's average response level MEAN_LEVEL
# (averaged over the population and t); 0 for Gaussian and scat.
family_intercept <- function(cell, truth, scale) {
  if (cell$family %in% c("gaussian", "scat")) return(0)
  pop <- population_sample(covariate_sd(cell$covariate %||% "rich"))
  eta0 <- linear_predictor(truth, pop, t_grid(61), scale, 0)
  linkinv <- FAMILIES[[cell$family]]$linkinv
  target <- MEAN_LEVEL[[cell$family]]
  uniroot(
    \(c) mean(linkinv(c + eta0)) - target,
    INTERCEPT_RANGE,
    tol = 1e-10
  )$root
}


# Error processes -------------------------------------------------------------------------

# Kernel of an error process on grid t: a correlation matrix (unit diagonal), except for
# fpc_raw*, whose variance varies along t (average 1; Gaussian responses only).
correlation_kernel <- function(error, t) {
  lag <- abs(outer(t, t, "-"))
  switch(
    error,
    iid = diag(length(t)),
    warp = diag(length(t)), # all dependence comes from the warps
    ou = (1 - OU_NUGGET) * exp(-lag / OU_RANGE) + OU_NUGGET * diag(length(t)),
    ou0 = exp(-lag / OU0_RANGE),
    osc = (1 - OSC_NUGGET) *
      exp(-lag / OSC_RANGE) *
      cos(2 * pi * lag / OSC_PERIOD) +
      OSC_NUGGET * diag(length(t)),
    fpc = fpc_kernel(t),
    fpc_het = fpc_kernel(t),
    fpc_raw = fpc_kernel(t, standardize = FALSE),
    fpc_raw_het = fpc_kernel(t, standardize = FALSE),
    stop("unknown error process: ", error)
  )
}

# Smooth random curves with a realistic, irregular eigenstructure: the residual FPCs of
# the DTI application (eigenfunctions interpolated to t) plus its white noise,
# standardized to unit variance at every t (the residual variance is not constant).
fpc_kernel <- function(t, standardize = TRUE) {
  k <- readRDS(FPC_KERNEL_FILE)
  phi <- apply(
    k$efunctions,
    2,
    \(e) stats::spline(k$grid, e, xout = t, method = "natural")$y
  )
  raw <- phi %*% (k$evalues * t(phi)) + k$sigma2 * diag(length(t))
  if (standardize) cov2cor(raw) else raw / mean(diag(raw))
}

# Per-subject error scale: exp(HET_SD z - HET_SD^2) for the *_het processes (so that
# E[scale^2] = 1), else 1.
subject_scale <- function(error, z) {
  if (!grepl("_het$", error)) return(rep(1, length(z)))
  exp(HET_SD * z - HET_SD^2)
}

# Latent Gaussian curves (rows) with the error's correlation on T_DENSE, from
# standard normals `n01` (G x 241).
latent_curves <- function(error, n01) {
  if (error %in% c("iid", "warp")) return(n01)
  n01 %*% chol(correlation_kernel(error, T_DENSE))
}

# Warps ----------------------------------------------------------------------------------

# a ~ N(0, WARP_SD^2) truncated to |a| <= WARP_MAX, by inverse CDF from uniforms.
warp_amplitude <- function(u) {
  p <- pnorm(WARP_MAX / WARP_SD)
  WARP_SD * qnorm(1 - p + u * (2 * p - 1))
}

# h(t) = t + a sin(pi t) / pi: one row of warped times per amplitude.
warp_times <- function(a, t)
  outer(rep(1, length(a)), t) + outer(a, sin(pi * t) / pi)

# Quadrature nodes and weights for expectations over the truncated warp law.
warp_quadrature <- function() {
  a <- seq(-WARP_MAX, WARP_MAX, length.out = WARP_QUAD_N)
  w <- dnorm(a, 0, WARP_SD)
  list(a = a, w = w / sum(w))
}

# Simulation ---------------------------------------------------------------------------

# All random draws for replicate `rep`, at full size (G_MAX curves, dense grid) and in
# a fixed order, so every cell sees the same draws (common random numbers).
draw_replicate <- function(rep) {
  withr::with_seed(
    SEED_BASE + rep,
    list(
      xi = matrix(rnorm(G_MAX * X_N_TERMS), G_MAX),
      z = rnorm(G_MAX),
      x = runif(G_MAX),
      warp_u = runif(G_MAX),
      n01 = matrix(rnorm(G_MAX * length(T_DENSE)), G_MAX)
    )
  )
}

simulate <- function(cell, rep) {
  draws <- draw_replicate(rep)
  keep <- seq_len(cell$G)
  t_index <- T_SUBSAMPLE[[as.character(cell$D)]]
  t <- T_DENSE[t_index]
  x_sd <- covariate_sd(cell$covariate %||% "rich")
  cov <- make_covariates(
    draws$xi[keep, , drop = FALSE],
    draws$z[keep],
    draws$x[keep],
    x_sd
  )
  truth <- make_truth(cell$model, cell$truth, x_sd)
  sc <- signal_scale(cell)
  intercept <- family_intercept(cell, truth, sc$scale)

  times <- t
  if (cell$error == "warp")
    times <- warp_times(warp_amplitude(draws$warp_u[keep]), t)
  eta <- linear_predictor(truth, cov, times, sc$scale, intercept)
  latent <- latent_curves(cell$error, draws$n01[keep, , drop = FALSE]) *
    subject_scale(cell$error, cov$z)
  if (cell$family != "gaussian" && grepl("^fpc_raw|_het$", cell$error)) {
    stop("error process ", cell$error, " is for Gaussian responses only")
  }
  y <- matrix(
    FAMILIES[[cell$family]]$draw(latent[, t_index], eta, sc$sigma),
    nrow = cell$G
  )

  list(
    data = model_frame(cell$model, cov, y),
    yind = t,
    xind = S_GRID,
    test = model_frame(cell$model, test_cohort(x_sd)),
    truth = estimand_targets(cell, truth, sc$scale, intercept, cov),
    diag = list(x_rank = centred_rank(cov$X), intercept = intercept)
  )
}

model_frame <- function(model, cov, y = NULL) {
  frame <- if (model == "M1") {
    data.frame(X = I(cov$X), z = cov$z)
  } else {
    data.frame(x = cov$x, z = cov$z)
  }
  if (!is.null(y)) frame$Y <- I(y)
  frame
}

# Number of principal components explaining 99.5% of the centred covariate variance.
centred_rank <- function(X) {
  ev <- svd(scale(X, scale = FALSE), nu = 0, nv = 0)$d^2
  which(cumsum(ev) / sum(ev) >= RANK_PVE)[1]
}

# Estimand targets -----------------------------------------------------------------------

# Truth of every estimand of the cell's model, in eval_points() order. `mean` is on
# the response scale, all others on the link scale. For warp cells, the targets are
# expectations over the warp law (observed clock): exact for Gaussian responses; for
# other families the coefficient targets are the link-scale approximation.
estimand_targets <- function(cell, truth, scale, intercept, cov) {
  at <- \(t) targets_at(t, cell, truth, scale, intercept, cov)
  if (cell$error != "warp") return(at(EVAL_T))
  quad <- warp_quadrature()
  per_node <- lapply(quad$a, \(a) at(drop(warp_times(a, EVAL_T))))
  Map(
    \(name) Reduce(`+`, Map(\(tg, w) w * tg[[name]], per_node, quad$w)),
    names(per_node[[1]])
  )
}

# Targets when every curve's clock reads t (a vector of length EVAL_T).
targets_at <- function(t, cell, truth, scale, intercept, cov) {
  linkinv <- FAMILIES[[cell$family]]$linkinv
  out <- list(
    mean = as.vector(linkinv(linear_predictor(
      truth,
      test_cohort(covariate_sd(cell$covariate %||% "rich")),
      t,
      scale,
      intercept
    ))),
    gamma = scale * truth$gamma(t)
  )
  if (cell$model == "M1") {
    out$beta <- as.vector(scale * truth$beta(EVAL_S, t))
    out$alpha <- intercept + scale * truth$bump(t)
  } else {
    # pffr centres f on the sample: sum_i f(x_i, t) = 0 at every t
    centre <- colMeans(truth$f(cov$x, t))
    out$f <- as.vector(scale * sweep(truth$f(EVAL_X, t), 2, centre))
  }
  out[model_estimands(cell$model)]
}
