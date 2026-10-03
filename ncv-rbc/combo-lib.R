# Bias correction (twicing) combined with the cluster-robust
# λ̂-uncertainty term for Gaussian NCV fits (lambda-correction/ncv-lib.R, ncv_covariances()).
# Needs lambda-correction/lc-lib.R and ncv-lib.R sourced.
#
# From ncv_covariances(): cl2 = V_f + B2 and ij_linx = V_f,ρ + B2, where V_f is the frequentist
# CL2, B2 = Vp − Ve the Bayesian allowance and V_f,ρ the CL2 with the IJ λ̂ term (cross form).
# Arms (z intervals; A = Ve Vp^-1, M1 = 2I − A, M2 = 3I − 3A + A²):
#   ncv_cl2      NCV estimate, Bayesian CL2 (study arm)
#   ncv_ijx      NCV estimate, λ-corrected CL2 (their arm)
#   rbc_p        M1 θ̂, M1 V_f M1' + (I − A) B2 (I − A)'
#   rbc_p_ijx    M1 θ̂, M1 V_f,ρ M1' + (I − A) B2 (I − A)'
#   rbc2_p       M2 θ̂, M2 V_f M2' + (I − A)² B2 (I − A)²'
#   rbc2_p_ijx   M2 θ̂, M2 V_f,ρ M2' + (I − A)² B2 (I − A)²'
#   ncv_ij, rbc_p_ij, rbc2_p_ij   as above with their IJ term without the cross term (ij_lin)
# The λ̂ term enters through M at first order (∂A/∂ρ θ̂ is ignored).

combo_task <- function(cell, rep, simulate_fun = simulate) {
  if (cell$fit_family != "gaussian") stop("Gaussian cells only")
  sim <- simulate_fun(cell, rep)
  fit <- fit_model(sim, cell, "NCV")
  cov <- ncv_covariances(fit, onestep = FALSE)$V
  B2 <- fit$Vp - fit$Ve
  Vf <- cov$cl2 - B2
  Vfr <- cov$ij_linx - B2
  Vfr0 <- cov$ij_lin - B2
  p <- length(fit$coefficients)
  A <- t(solve(fit$Vp, fit$Ve))
  I <- diag(p)
  IA <- I - A
  M1 <- 2 * I - A
  M2 <- 3 * I - 3 * A + A %*% A
  th <- fit$coefficients
  B2_1 <- IA %*% B2 %*% t(IA)
  B2_2 <- IA %*% B2_1 %*% t(IA)
  arms <- list(
    ncv_cl2 = list(theta = th, V = cov$cl2),
    ncv_ijx = list(theta = th, V = cov$ij_linx),
    rbc_p = list(theta = drop(M1 %*% th), V = M1 %*% Vf %*% t(M1) + B2_1),
    rbc_p_ijx = list(theta = drop(M1 %*% th), V = M1 %*% Vfr %*% t(M1) + B2_1),
    rbc2_p = list(theta = drop(M2 %*% th), V = M2 %*% Vf %*% t(M2) + B2_2),
    rbc2_p_ijx = list(theta = drop(M2 %*% th), V = M2 %*% Vfr %*% t(M2) + B2_2),
    ncv_ij = list(theta = th, V = cov$ij_lin),
    rbc_p_ij = list(theta = drop(M1 %*% th), V = M1 %*% Vfr0 %*% t(M1) + B2_1),
    rbc2_p_ij = list(theta = drop(M2 %*% th), V = M2 %*% Vfr0 %*% t(M2) + B2_2)
  )
  z <- qnorm(1 - (1 - LEVEL) / 2)
  rows <- list()
  for (e in model_estimands(cell$model)) {
    L <- evaluation_matrix(fit, sim, e)
    truth <- as.vector(sim$truth[[e]])
    for (nm in names(arms)) {
      a <- arms[[nm]]
      est <- drop(L %*% a$theta)
      se <- sqrt(pmax(rowSums((L %*% a$V) * L), 0))
      rows[[length(rows) + 1]] <- data.frame(
        cell = cell$cell,
        rep = rep,
        estimand = e,
        arm = nm,
        point = seq_along(truth),
        truth = truth,
        est = est,
        lower = est - z * se,
        upper = est + z * se
      )
    }
  }
  list(points = do.call(rbind, rows), time = fit$restart$time)
}
