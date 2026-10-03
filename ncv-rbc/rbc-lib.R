# Bias-corrected intervals around the NCV fit (exploratory, 2026-10-01 night).
#
# For a penalized fit at fixed lambda, E(theta-hat) = A theta with A = F^-1 X'WX,
# F = X'WX + S_lambda; mgcv's Ve = F^-1 X'WX F^-1 phi and Vp = F^-1 phi give
# A = Ve Vp^-1. One-step bias correction ("twicing"): theta-tilde = (2I - A) theta-hat,
# whose bias is -(I - A)^2 theta instead of -(I - A) theta. Its covariance is propagated
# through the same linear map (robust bias correction in the sense of Calonico, Cattaneo
# and Farrell 2018): V(theta-tilde) = M V M', M = 2I - A, V = CL2 of theta-hat.
#
# Arms (z intervals, all from the one NCV fit):
#   ncv_cl2      NCV estimate, Bayesian CL2 (the study's arm; reproduced as a check)
#   ncv_cl2f     NCV estimate, frequentist CL2
#   rbc_f        twiced estimate, M V_cl2f M'
#   rbc_b        twiced estimate, M V_cl2 M' (Bayesian CL2 inside)
#   rbc2_f       two correction steps: M2 = 3I - 3A + A^2 (bias -(I - A)^3 theta)
#   rbc_p        twiced estimate, M V_cl2f M' + (I - A)(Vp - Ve)(I - A)': the Bayesian bias
#                allowance propagated correctly. Vp - Ve = F^-1 S F^-1 phi is the prior
#                covariance of the bias -(I - A) theta under theta ~ N(0, phi S^-); the twiced
#                bias -(I - A)^2 theta has prior covariance (I - A)(Vp - Ve)(I - A)'.
#   reml_cl2     REML fit, Bayesian CL2 (the study's recipe) on the same data

rbc_task <- function(cell, rep, simulate_fun = simulate, with_reml = FALSE) {
  sim <- simulate_fun(cell, rep)
  fit <- fit_model(sim, cell, "NCV")
  V <- list(
    b = refund:::pffr_vcov(
      fit,
      sandwich = "cl2",
      cl2_adjustment = "exact",
      freq = FALSE
    ),
    f = refund:::pffr_vcov(
      fit,
      sandwich = "cl2",
      cl2_adjustment = "exact",
      freq = TRUE
    )
  )
  p <- length(fit$coefficients)
  A <- t(solve(fit$Vp, fit$Ve)) # Ve Vp^-1 (both symmetric)
  I <- diag(p)
  M1 <- 2 * I - A
  M2 <- 3 * I - 3 * A + A %*% A
  theta <- fit$coefficients
  arms <- list(
    ncv_cl2 = list(theta = theta, V = V$b),
    ncv_cl2f = list(theta = theta, V = V$f),
    rbc_f = list(theta = drop(M1 %*% theta), V = M1 %*% V$f %*% t(M1)),
    rbc_b = list(theta = drop(M1 %*% theta), V = M1 %*% V$b %*% t(M1)),
    rbc2_f = list(theta = drop(M2 %*% theta), V = M2 %*% V$f %*% t(M2)),
    rbc_p = list(
      theta = drop(M1 %*% theta),
      V = M1 %*% V$f %*% t(M1) + (I - A) %*% (fit$Vp - fit$Ve) %*% t(I - A)
    )
  )
  if (with_reml) {
    fitR <- fit_model(sim, cell, "REML")
    VR <- refund:::pffr_vcov(
      fitR,
      sandwich = "cl2",
      cl2_adjustment = "exact",
      freq = FALSE
    )
  }
  linkinv <- FAMILIES[[cell$fit_family]]$linkinv
  z <- qnorm(1 - (1 - LEVEL) / 2)
  rows <- list()
  for (e in model_estimands(cell$model)) {
    L <- evaluation_matrix(fit, sim, e)
    truth <- as.vector(sim$truth[[e]])
    tf <- if (e == "mean") linkinv else identity
    LR <- if (with_reml) evaluation_matrix(fitR, sim, e)
    for (nm in c(names(arms), if (with_reml) "reml_cl2")) {
      a <- if (nm == "reml_cl2") list(theta = fitR$coefficients, V = VR) else
        arms[[nm]]
      if (nm == "reml_cl2") L <- LR
      est <- drop(L %*% a$theta)
      se <- sqrt(rowSums((L %*% a$V) * L))
      lo <- tf(est - z * se)
      hi <- tf(est + z * se)
      rows[[length(rows) + 1]] <- data.frame(
        cell = cell$cell,
        rep = rep,
        estimand = e,
        arm = nm,
        point = seq_along(truth),
        truth = truth,
        est = tf(est),
        lower = lo,
        upper = hi
      )
    }
  }
  list(points = do.call(rbind, rows), time = fit$restart$time)
}
