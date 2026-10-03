# Checks of the whitened-fit CL2 on one synthetic replicate (default: ou, G = 40, rep 1):
# (1) cl2_from_design() on the REML fit reproduces refund:::pffr_vcov(cl2, exact, Bayesian);
# (2) on the GLS-fpc fit, a direct dense computation (H = Xw Vp Xw', residual-maker blocks
#     ((I - H)(I - H)')_gg, their inverse square roots, the G/(G-1) sum of the leverage-adjusted
#     score outer products plus Vp - Ve) reproduces cl2_from_design();
# (3) eigenvalues of the whitened curve blocks H_gg lie in [0, 1] and sum to the EDF.
# Run from the repository root: Rscript comparators/check-cl2-arms.R [key] [rep]

source("scripts/setup.R")
source("plasmode/sources.R")
source("plasmode/plasmode.R")
source("comparators/comparators-lib.R")
source("comparators/cl2-arms-lib.R")

args <- commandArgs(trailingOnly = TRUE)
key <- if (length(args) >= 1) args[1] else "s009"
rep <- if (length(args) >= 2) as.integer(args[2]) else 1L
cell <- cell_spec(key)
sim <- if (cell$source == "synthetic") simulate(cell, rep) else
  plasmode_simulate(cell, rep)
fit <- fit_model(sim, cell, "REML")
G <- nrow(sim$data)
E <- residual_curves(fit, G)
D <- ncol(E)
cid <- rep(seq_len(G), each = D)

# (1)
V_ref <- refund:::pffr_vcov(
  fit,
  sandwich = "cl2",
  cl2_adjustment = "exact",
  freq = FALSE
)
V_own <- cl2_from_design(
  gam_lpmatrix(fit),
  as.vector(fit$y),
  fit$coefficients,
  fit$Vp,
  fit$Ve,
  fit$sig2,
  cid
)
cat(sprintf(
  "(1) REML: max|V_own - V_refund| / max|V_refund| = %.2e\n",
  max(abs(V_own - V_ref)) / max(abs(V_ref))
))

# (2)
fpc <- residual_fpca(E, sim$yind)
Sigma <- fpc$efunctions %*%
  (fpc$evalues * t(fpc$efunctions)) +
  fpc$sigma2 * diag(D)
gls <- fit_gls_setup(fit, Sigma, G, D)
m <- gls$fit
V_core <- cl2_from_design(gls$X, gls$y, m$coefficients, m$Vp, m$Ve, m$sig2, cid)
s <- 1 / sqrt(m$sig2)
Xw <- gls$X * s
z <- drop(gls$y - gls$X %*% m$coefficients) * s
H <- Xw %*% m$Vp %*% t(Xw)
R <- diag(nrow(H)) - H
meat <- matrix(0, ncol(Xw), ncol(Xw))
hat_eig <- numeric()
for (g in seq_len(G)) {
  idx <- which(cid == g)
  M <- tcrossprod(R[idx, , drop = FALSE]) # ((I - H)(I - H)')_gg
  e <- eigen(M, symmetric = TRUE)
  A <- e$vectors %*% (t(e$vectors) / sqrt(e$values))
  u <- crossprod(Xw[idx, ], A %*% z[idx])
  meat <- meat + tcrossprod(u)
  hat_eig <- c(
    hat_eig,
    eigen(H[idx, idx], symmetric = TRUE, only.values = TRUE)$values
  )
}
V_direct <- G / (G - 1) * m$Vp %*% meat %*% m$Vp + (m$Vp - m$Ve)
cat(sprintf(
  "(2) GLS-fpc: max|V_direct - V_core| / max|V_core| = %.2e\n",
  max(abs(V_direct - V_core)) / max(abs(V_core))
))
# (3)
cat(sprintf(
  "(3) whitened H_gg eigenvalues in [%.3g, %.3g]; sum = %.3f, EDF = %.3f\n",
  min(hat_eig),
  max(hat_eig),
  sum(hat_eig),
  sum(m$edf)
))
lev_reml <- attr(V_own, "diag")
cat(sprintf(
  "    REML H_gg eigenvalues in [%.3g, %.3g]; whitened max h_ii %.3g (REML %.3g)\n",
  min(lev_reml$min_hat_eig),
  max(lev_reml$max_leverage),
  max(diag(H)),
  max(lev_reml$max_obs_leverage)
))
L <- reml_matrices(fit, sim, cell)$beta
ratio <- se_rows(L, V_core) / se_rows(L, model_cov(m))
cat(sprintf(
  "    beta: SE(GLS CL2) / SE(GLS Vc) median %.3f; SE(GLS CL2) / SE(REML CL2) median %.3f\n",
  median(ratio),
  median(se_rows(L, V_core) / se_rows(L, V_ref))
))
