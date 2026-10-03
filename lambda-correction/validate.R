# Numerical checks of the lambda correction on single replicates (local, ~minutes each):
#   Rscript lambda-correction/validate.R <cell> <rep> [n_fd_curves]
# (1) analytic curve scores s_g vs central finite differences of omega-hat in the curve's
#     prior weight; (2) IJ and one-step leave-curve-out vs exact leave-curve-out refits
#     (zero prior weights, same parameterisation), for omega-hat and theta-hat.
source("scripts/setup.R")
source("lambda-correction/lc-lib.R")

args <- commandArgs(trailingOnly = TRUE)
cell_id <- as.integer(args[1])
rep_id <- as.integer(args[2])
n_fd <- if (length(args) >= 3) as.integer(args[3]) else 3
out_dir <- "lambda-correction/validation"
dir.create(out_dir, showWarnings = FALSE)

cell <- design_table()
cell <- cell[cell$cell == cell_id, ]
sim <- simulate(cell, rep_id)
fit <- fit_model(sim, cell, "REML")
cat(sprintf(
  "cell %d (%s/%s/G=%d) rep %d, fit %.1f s\n",
  cell_id,
  cell$fit_family,
  cell$error,
  cell$G,
  rep_id,
  fit$restart$time
))

t0 <- proc.time()[[3]]
P <- lc_pieces(fit)
s_ij <- ij_scores(P)
t_ij <- proc.time()[[3]] - t0
t0 <- proc.time()[[3]]
os <- onestep_loo(P)
t_os <- proc.time()[[3]] - t0
cat(sprintf("IJ scores %.2f s, one-step %.2f s\n", t_ij, t_os))
Hinv <- pinv_sym(P$H)

# (1) finite differences in the curve weight (Gaussian: mgcv weights keep n fixed) -------
G0 <- gam_setup(fit)
eps <- 1e-3
s_now <- ij_scores(P, n_term = FALSE)
fd <- lapply(P$groups[seq_len(n_fd)], \(g) {
  w <- rep(1, length(P$y))
  ig <- P$cluster == g
  w[ig] <- 1 + eps
  up <- refit_weighted(fit, G0, w)
  w[ig] <- 1 - eps
  dn <- refit_weighted(fit, G0, w)
  data.frame(
    curve = g,
    par = seq_along(up$omega),
    fd = (up$omega - dn$omega) / (2 * eps),
    ij = -drop(Hinv %*% s_now[g, ])
  )
})
fd <- do.call(rbind, fd)
cat("\n(1) d omega-hat / d w_g: finite difference vs analytic\n")
print(fd, digits = 4)

# (2) exact jackknife ----------------------------------------------------------------
t0 <- proc.time()[[3]]
jk <- jackknife_refits(fit, P, G0)
t_jk <- proc.time()[[3]] - t0
cat(sprintf(
  "\njackknife: %d refits, %.1f s, converged %d\n",
  length(jk),
  t_jk,
  sum(sapply(jk, `[[`, "converged"))
))
d_om_jk <- t(sapply(jk, \(m) m$omega)) -
  matrix(
    c(log(fit$sp), if (P$gaussian) log(fit$sig2)),
    length(jk),
    ncol(s_ij),
    byrow = TRUE
  )
d_om_ij <- s_ij %*% Hinv # leave-out: delta w = -1 => +H^-1 s_g (uncentred)
d_om_os <- os$d_omega
d_th_jk <- sapply(jk, \(m) m$theta) - P$theta
K <- length(P$sp)
d_th_os <- os$d_theta + P$J %*% t(d_om_os[, seq_len(K), drop = FALSE])
d_th_ij <- -P$Ainv %*%
  sapply(P$groups, \(g) {
    ig <- P$cluster == g
    crossprod(P$X[ig, , drop = FALSE], P$y[ig] - P$mu[ig])
  }) +
  P$J %*% t(d_om_ij[, seq_len(K), drop = FALSE])

cmp <- \(a, b, label) {
  data.frame(
    what = label,
    par = seq_len(ncol(a)),
    cor = diag(cor(a, b)),
    rel_rmse = sqrt(colMeans((a - b)^2)) / apply(b, 2, sd),
    sd_ratio = apply(a, 2, sd) / apply(b, 2, sd)
  )
}
cat(
  "\n(2) leave-curve-out change of omega-hat, vs exact jackknife (per parameter)\n"
)
print(
  rbind(cmp(d_om_ij, d_om_jk, "IJ"), cmp(d_om_os, d_om_jk, "one-step")),
  digits = 3
)

# variance of omega-hat: naive (REML Hessian), IJ, one-step JK, exact JK
Vom <- list(
  naive = Hinv,
  ij = omega_cov(-d_om_ij, "ij"),
  onestep = omega_cov(d_om_os, "jk"),
  jackknife = omega_cov(d_om_jk, "jk")
)
cat("\nSD of omega-hat components\n")
print(sapply(Vom, \(V) sqrt(diag(V))), digits = 3)

# effect on the estimands' SEs (grid-averaged ratio to CL2)
core <- refund:::pffr_influence(fit, sandwich = "cl2", cl2_adjustment = "exact")
V_cl2 <- core$correction * tcrossprod(core$K) + core$B2
th_jk_var <- (length(jk) - 1) /
  length(jk) *
  tcrossprod(d_th_jk - rowMeans(d_th_jk))
covs <- list(
  cl2 = V_cl2,
  naive = V_cl2 + naive_correction(P),
  vc = V_cl2 + fit$Vc - fit$Vp,
  ij_lin = add_rho_part(core, rho_part_linear(P, -d_om_ij), "ij"),
  ij_linx = add_rho_part(
    core,
    rho_part_linear(P, -d_om_ij),
    "ij",
    cross = TRUE
  ),
  os_lin = add_rho_part(core, rho_part_linear(P, -d_om_os), "jk"),
  jk_lin = add_rho_part(core, rho_part_linear(P, -d_om_jk), "jk"),
  jk_full = th_jk_var + core$B2
)
se_ratio <- sapply(model_estimands(cell$model), \(e) {
  L <- evaluation_matrix(fit, sim, e)
  se <- sapply(covs, \(V) sqrt(rowSums((L %*% V) * L)))
  colMeans(se / se[, "cl2"])
})
cat("\nSE / SE(CL2), grid-averaged\n")
print(se_ratio, digits = 4)
cat("\ntheta influence vs jackknife (all coefficients pooled)\n")
print(
  rbind(
    cmp(matrix(d_th_ij), matrix(d_th_jk), "IJ total"),
    cmp(matrix(d_th_os), matrix(d_th_jk), "one-step total")
  ),
  digits = 3
)

saveRDS(
  list(
    cell = cell,
    rep = rep_id,
    fd = fd,
    d_om = list(ij = d_om_ij, os = d_om_os, jk = d_om_jk),
    d_th = list(ij = d_th_ij, os = d_th_os, jk = d_th_jk),
    Vom = Vom,
    se_ratio = se_ratio,
    time = c(
      fit = fit$restart$time,
      ij = t_ij,
      onestep = t_os,
      jackknife = t_jk
    )
  ),
  file.path(out_dir, sprintf("validate-c%03d-r%03d.rds", cell_id, rep_id))
)
