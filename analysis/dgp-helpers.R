# Shared helpers for describing the data-generating process (used by describe-cells.R and
# dgp-visualizations.qmd). Source after scripts/setup.R.

# Residuals Y - E[Y | X, z] of one replicate. Recomputes the conditional mean with the
# functions simulate() uses (R/dgp.R is frozen during the runs, so it is not changed to
# return it); for warp cells the conditional mean is the expectation over the warp law.
residual_curves <- function(cell, rep) {
  draws <- draw_replicate(rep)
  keep <- seq_len(cell$G)
  t_index <- T_SUBSAMPLE[[as.character(cell$D)]]
  t <- T_DENSE[t_index]
  cov <- make_covariates(
    draws$xi[keep, , drop = FALSE],
    draws$z[keep],
    draws$x[keep]
  )
  truth <- make_truth(cell$model, cell$truth)
  sc <- signal_scale(cell)
  intercept <- family_intercept(cell, truth, sc$scale)
  linkinv <- FAMILIES[[cell$family]]$linkinv
  mean_at <- \(times)
    linkinv(linear_predictor(truth, cov, times, sc$scale, intercept))
  conditional_mean <- mean_at(t)
  if (cell$error == "warp") {
    quad <- warp_quadrature()
    conditional_mean <- Reduce(
      `+`,
      Map(
        \(a, w) {
          w * mean_at(drop(warp_times(a, t))) # same clock for every curve
        },
        quad$a,
        quad$w
      )
    )
  }
  sim <- simulate(cell, rep)
  unclass(sim$data$Y) - conditional_mean
}
