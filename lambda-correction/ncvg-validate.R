# Checks of the GLM NCV lambda correction on one replicate (local):
#   Rscript lambda-correction/ncvg-validate.R <cell> <rep> [n_fd_curves] [exact]
# exact = "exact": leave-curve-out NCV minimisation (zero curve weight) for every curve.
source("scripts/setup.R")
source("lambda-correction/lc-lib.R")
source("lambda-correction/ncv-lib.R")
source("lambda-correction/ncvg-lib.R")

args <- commandArgs(trailingOnly = TRUE)
cell_id <- as.integer(args[1])
rep_id <- as.integer(args[2])
n_fd <- if (length(args) >= 3) as.integer(args[3]) else 3
exact <- length(args) >= 4 && args[4] == "exact"
cells <- design_table()
cell <- cells[cells$cell == cell_id, ]
sim <- simulate(cell, rep_id)
fit <- fit_model(sim, cell, "NCV")
cat(sprintf(
  "cell %d (%s/%s/G=%d) rep %d, NCV fit %.1f s, conv: %s\n",
  cell_id,
  cell$fit_family,
  cell$error,
  cell$G,
  rep_id,
  fit$restart$time,
  fit$outer.info$conv
))
P <- ncvg_pieces(fit)
K <- length(P$sp)
G <- length(P$groups)
cat("n =", nrow(P$X), " p =", ncol(P$X), " K =", K, " G =", G, "\n")
print(P$check)
t0 <- proc.time()[[3]]
st <- ncvg_state(P)
cat(
  "state time",
  proc.time()[[3]] - t0,
  " max |beta - fit coef|",
  max(abs(st$beta - P$theta)),
  "\n"
)

cat("\n(1) value and gradient vs mgcv\n")
print(c(mine = ncvg_value(st), mgcv = unname(fit$gcv.ubre)), digits = 12)
t0 <- proc.time()[[3]]
g <- ncvg_grad(P, st)
cat("grad time", proc.time()[[3]] - t0, "\n")
print(rbind(mine = g, mgcv = drop(fit$outer.info$grad)))

cat("\n(2) gradient vs central differences of the value (refits of beta)\n")
fd_g <- vapply(
  seq_len(K),
  \(k) {
    e <- replace(rep(0, K), k, 1e-5)
    (ncvg_value(ncvg_state(P, P$sp * exp(e))) -
      ncvg_value(ncvg_state(P, P$sp * exp(-e)))) /
      2e-5
  },
  0
)
print(rbind(analytic = g, fd = fd_g))

cat("\n(3) curve-weight derivatives vs central differences of the value\n")
dw <- ncvg_dderiv(
  P,
  st,
  matrix(0, K, n_fd),
  diag(G)[, seq_len(n_fd), drop = FALSE]
)
fd_w <- vapply(
  seq_len(n_fd),
  \(h) {
    wp <- replace(rep(1, G), h, 1 + 1e-5)
    wm <- replace(rep(1, G), h, 1 - 1e-5)
    (ncvg_value(ncvg_state(P, w = wp)) - ncvg_value(ncvg_state(P, w = wm))) /
      2e-5
  },
  0
)
print(rbind(analytic = dw, fd = fd_w))

cat(
  "\n(4) Hessian and scores (FD in rho of exact derivatives) vs the other route\n"
)
t0 <- proc.time()[[3]]
hs <- ncvg_hessian_scores(P)
cat("hessian+scores time", proc.time()[[3]] - t0, "\n")
print(round(hs$H, 4))
print(round(P$H_mgcv, 4))
# scores by central differences in w_h of the exact rho-gradient
fd_s <- t(vapply(
  seq_len(n_fd),
  \(h) {
    wp <- replace(rep(1, G), h, 1 + 1e-5)
    wm <- replace(rep(1, G), h, 1 - 1e-5)
    (ncvg_grad(P, ncvg_state(P, w = wp)) -
      ncvg_grad(P, ncvg_state(P, w = wm))) /
      2e-5
  },
  numeric(K)
))
print(rbind(rho_route = hs$s[seq_len(n_fd), , drop = FALSE], w_route = fd_s))
cat(
  "scaling identity sum_h s_h = grad - H 1: max abs diff",
  max(abs(colSums(hs$s) - (g - drop(hs$H %*% rep(1, K))))),
  " (scale",
  max(abs(colSums(hs$s))),
  ")\n"
)

if (identical(P$family$family, "gaussian")) {
  cat("\n(5) Gaussian: against the Gaussian implementation\n")
  Hg <- ncv_hessian(P)
  sg <- ncv_scores(P)
  cat(
    "max rel diff H",
    max(abs(Hg - hs$H)) / max(abs(Hg)),
    " s",
    max(abs(sg - hs$s)) / max(abs(sg)),
    "\n"
  )
}

free <- diag(hs$H) > BOUNDARY_TOL * max(diag(hs$H)) # as in ncvg_covariances()
cat("\nboundary (treated as fixed):", names(P$sp)[!free], "\n")
D_ij <- matrix(0, G, K)
D_ij[, free] <- ncv_ij_influence(
  P,
  hs$H[free, free, drop = FALSE],
  hs$s[, free, drop = FALSE]
)
cat("\nSD of rho-hat: naive (mgcv hess) / naive (mine) / IJ\n")
print(
  rbind(
    naive_mgcv = sqrt(pmax(diag(pinv_sym(P$H_mgcv)), 0)),
    naive = sqrt(pmax(diag(pinv_sym(hs$H)), 0)),
    ij = sqrt(diag(omega_cov(D_ij, "ij")))
  ),
  digits = 3
)

if (exact) {
  cat(
    "\n(6) leave-curve-out rho-hat: IJ vs exact minimisation of NCV with curve h deleted\n"
  )
  t0 <- proc.time()[[3]]
  refits <- lapply(seq_len(G), \(h) ncvg_refit_rho(P, replace(rep(1, G), h, 0)))
  conv <- vapply(refits, \(r) isTRUE(attr(r, "converged")), TRUE)
  jk <- t(vapply(refits, \(r) r - log(P$sp), numeric(K)))
  cat(sprintf(
    "exact refits %.1f s, converged %d of %d\n",
    proc.time()[[3]] - t0,
    sum(conv),
    G
  ))
  ij <- -D_ij # leave-out orientation
  # converged refits only; boundary components are fixed in the IJ and reported separately
  ijc <- ij[conv, , drop = FALSE]
  jkc <- jk[conv, , drop = FALSE]
  print(
    data.frame(
      par = seq_len(K),
      boundary = !free,
      cor = diag(cor(ijc, jkc)),
      spearman = diag(cor(ijc, jkc, method = "spearman")),
      sd_ratio = apply(ijc, 2, sd) / apply(jkc, 2, sd),
      iqr_ratio = apply(ijc, 2, IQR) / apply(jkc, 2, IQR),
      max_abs_jk = apply(abs(jkc), 2, max)
    ),
    digits = 3
  )
  print(
    rbind(
      ij = sqrt(diag(omega_cov(D_ij, "ij"))),
      jk = sqrt(diag(omega_cov(jkc, "jk")))
    ),
    digits = 3
  )
  saveRDS(
    list(ij = ij, jk = jk, converged = conv, H = hs$H, s = hs$s),
    sprintf("lambda-correction/validation/ncvg-%d-%d.rds", cell_id, rep_id)
  )
}
