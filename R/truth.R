# Frozen truths, evaluation points and observed-clock (warp) targets.

# Evaluation points ----------------------------------------------------------------

# Points at which each estimand is scored, in vector order (first column fastest). The
# mean's times are `eval_t` (plasmode apps on native grids pass their own).
eval_points <- function(estimand, eval_t = EVAL_T) {
  switch(
    estimand,
    mean = expand.grid(curve = seq_len(N_TEST), t = eval_t),
    beta = expand.grid(s = EVAL_S, t = EVAL_T),
    f = expand.grid(x = EVAL_X, t = EVAL_T),
    alpha = data.frame(t = EVAL_T),
    gamma = data.frame(t = EVAL_T),
    stop("unknown estimand: ", estimand)
  )
}

model_estimands <- function(model) {
  switch(
    model,
    M1 = c("mean", "beta", "alpha", "gamma"),
    M2 = c("mean", "f", "gamma"),
    stop("unknown model: ", model)
  )
}

# Spline bases on [0, 1] -------------------------------------------------------------

# mgcv P-spline basis on [0, 1] with the same knots pffr uses for a grid spanning
# [0, 1]. Returns the smooth object; evaluate with ps_eval().
ps_basis <- function(basis) {
  # knots depend only on the data range [0, 1]; the points just have to cover it
  mgcv::smoothCon(
    mgcv::s(u, bs = basis$bs, k = basis$k, m = basis$m),
    data = data.frame(u = seq(0, 1, length.out = PROJECTION_N)),
    absorb.cons = FALSE
  )[[1]]
}

ps_eval <- function(sm, u) mgcv::PredictMat(sm, data.frame(u = u))

# Knots for M2's s(x): fixed to [0, 1] so the fitted x-basis does not depend on the
# sample range of x (pffr/gam `knots` argument).
x_knots <- function() ps_basis(BASIS_X)$knots

# Composite Simpson weights on an equidistant grid with an odd number of points
# (pffr's default ff integration).
simpson_weights <- function(grid) {
  n <- length(grid)
  stopifnot("Simpson needs an odd number of points" = n %% 2 == 1)
  w <- rep(c(2, 4), length.out = n)
  w[c(1, n)] <- 1
  w * diff(grid[1:2]) / 3
}

# Least-squares projection of f(u, v) onto the tensor basis bu x bv; returns the
# coefficient matrix Theta with f(u, v) ~ Bu(u) Theta Bv(v)'.
project_tensor <- function(f, bu, bv, n = PROJECTION_N) {
  u <- seq(0, 1, length.out = n)
  Bu <- ps_eval(bu, u)
  Bv <- ps_eval(bv, u)
  Fuv <- outer(u, u, f)
  solve(crossprod(Bu), t(Bu)) %*% Fuv %*% t(solve(crossprod(Bv), t(Bv)))
}

project_curve <- function(f, b, n = PROJECTION_N) {
  u <- seq(0, 1, length.out = n)
  B <- ps_eval(b, u)
  drop(solve(crossprod(B), crossprod(B, f(u))))
}

# Raw truth formulas --------------------------------------------------------------------

# `wiggly` adds a component with 1.5 periods in each direction; in s it loads on the
# covariate's weakly identified high-frequency directions.
TRUTH_FORMULAS <- list(
  beta = list(
    smooth = \(s, t) cos(pi * s) * sin(pi * t) + (s - 0.5) * cos(2 * pi * t),
    wiggly = \(s, t) {
      cos(pi * s) *
        sin(pi * t) +
        (s - 0.5) * cos(2 * pi * t) +
        0.6 * sin(3 * pi * s) * sin(3 * pi * t)
    }
  ),
  f = list(
    smooth = \(x, t) 2 * (x - 0.5) * cos(pi * t) + sin(pi * x) * t,
    wiggly = \(x, t) {
      2 *
        (x - 0.5) *
        cos(pi * t) +
        sin(pi * x) * t +
        0.6 * sin(3 * pi * x) * sin(3 * pi * t)
    }
  ),
  gamma = list(
    smooth = \(t) sin(2 * pi * t),
    wiggly = \(t) sin(2 * pi * t) + 0.5 * sin(5 * pi * t)
  )
)

raw_bump <- function(t)
  BUMP_HEIGHT * exp(-(t - BUMP_CENTER)^2 / (2 * BUMP_SD^2))

# Covariate basis ------------------------------------------------------------------------

# psi_0(s) = 1, psi_j(s) = sqrt(2) cos(j pi s); columns scaled by X_SD (or X_SD_LOWRANK) give
# X's law.
covariate_basis <- function(s) {
  j <- seq_len(X_N_TERMS - 1)
  cbind(1, sqrt(2) * cos(outer(s, j) * pi))
}

# Truth object -------------------------------------------------------------------------

truth_bases <- function(model) {
  bt <- ps_basis(BASIS_T)
  if (model == "M1") {
    ff_margin <- \(i)
      list(bs = BASIS_FF$bs, k = BASIS_FF$k[i], m = BASIS_FF$m[[i]])
    return(list(u = ps_basis(ff_margin(1)), v = ps_basis(ff_margin(2)), t = bt))
  }
  list(u = ps_basis(BASIS_X), v = bt, t = bt)
}

# Population variance of the main term's contribution, averaged over t in [0, 1].
main_term_variance <- function(main, model, x_sd = X_SD) {
  if (model == "M1") {
    # Var_X(sum_s w_s X(s) beta(s, t)) = sum_j sd_j^2 (sum_s w_s psi_j(s) beta(s, t))^2
    loads <- t(covariate_basis(S_GRID) * simpson_weights(S_GRID)) %*%
      main(S_GRID, T_DENSE)
    return(mean(colSums((x_sd * loads)^2)))
  }
  f_values <- main(seq(0, 1, length.out = PROJECTION_N), T_DENSE) # x ~ U(0, 1)
  mean(apply(f_values, 2, \(col) mean(col^2) - mean(col)^2))
}

# Truth functions for one model and roughness. They are splines in the fitted bases
# (exactly representable), scaled so the non-intercept signal S = 1, of which the
# main term (beta or f) contributes SIGNAL_SHARE_MAIN and gamma the rest.
# `null` (extension 2026-10-03): the smooth truth with the main term set to zero; gamma,
# the bump and (Gaussian) the error SD stay those of the smooth truth.
make_truth <- function(model, roughness, x_sd = X_SD) {
  b <- truth_bases(model)
  main_name <- if (model == "M1") "beta" else "f"
  null_main <- roughness == "null"
  if (null_main) roughness <- "smooth"
  theta <- project_tensor(TRUTH_FORMULAS[[main_name]][[roughness]], b$u, b$v)
  gamma_coef <- project_curve(TRUTH_FORMULAS$gamma[[roughness]], b$t)
  main <- \(u, t) ps_eval(b$u, u) %*% theta %*% t(ps_eval(b$v, t))
  gamma <- \(t) drop(ps_eval(b$t, t) %*% gamma_coef)

  main_scale <- sqrt(SIGNAL_SHARE_MAIN / main_term_variance(main, model, x_sd))
  if (null_main) main_scale <- 0
  gamma_scale <- sqrt((1 - SIGNAL_SHARE_MAIN) / mean(gamma(T_DENSE)^2))
  # The bump is projected onto the intercept basis too, so alpha is representable.
  bump_coef <- project_curve(raw_bump, b$t)
  truth <- list(
    model = model,
    main_name = main_name,
    gamma = \(t) gamma_scale * gamma(t),
    bump = \(t) drop(ps_eval(b$t, t) %*% bump_coef)
  )
  truth[[main_name]] <- \(u, t) main_scale * main(u, t)
  truth
}

# Linear predictor -----------------------------------------------------------------------

# Main-term contribution for covariates `cov` (list with X or x) at times t:
# a (curves x t) matrix. For M1 the integral uses pffr's Simpson weights.
main_contribution <- function(truth, cov, t) {
  if (truth$model == "M1") {
    return(
      (cov$X * rep(simpson_weights(S_GRID), each = nrow(cov$X))) %*%
        truth$beta(S_GRID, t)
    )
  }
  truth$f(cov$x, t)
}

# eta_i(t) = intercept + scale * (bump(t) + main_i(t) + z_i gamma(t)); `scale` is sqrt(S).
# t may be a matrix of per-curve times (warped clocks), one row per curve.
linear_predictor <- function(truth, cov, t, scale, intercept) {
  if (is.matrix(t)) {
    rows <- lapply(seq_len(nrow(t)), \(i) {
      one <- lapply(cov, \(v) if (is.matrix(v)) v[i, , drop = FALSE] else v[i])
      linear_predictor(truth, one, t[i, ], scale, intercept)
    })
    return(do.call(rbind, rows))
  }
  signal <- outer(rep(1, length(cov$z)), truth$bump(t)) +
    main_contribution(truth, cov, t) +
    outer(cov$z, truth$gamma(t))
  intercept + scale * signal
}
