# Checks of the Gaussian NCV lambda correction on one replicate (local):
#   Rscript lambda-correction/ncv-validate.R <cell> <rep> [n_fd_curves]
source("scripts/setup.R")
source("lambda-correction/lc-lib.R")
source("lambda-correction/ncv-lib.R")

args <- commandArgs(trailingOnly = TRUE)
cell_id <- as.integer(args[1])
rep_id <- as.integer(args[2])
n_fd <- if (length(args) >= 3) as.integer(args[3]) else 3
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
P <- ncv_pieces(fit)
st <- ncv_state(P)

cat("\n(1) value and gradient vs mgcv\n")
print(c(mine = ncv_value(st), mgcv = unname(fit$gcv.ubre)))
print(rbind(mine = ncv_grad(st), mgcv = drop(fit$outer.info$grad)))

cat("\n(2) gradient vs central differences of the NCV value\n")
K <- length(P$sp)
fd_g <- vapply(
  seq_len(K),
  \(k) {
    e <- replace(rep(0, K), k, 1e-5)
    (ncv_value(ncv_state(P, P$sp * exp(e))) -
      ncv_value(ncv_state(P, P$sp * exp(-e)))) /
      2e-5
  },
  0
)
print(rbind(analytic = ncv_grad(st), fd = fd_g))

cat(
  "\n(3) Hessian: finite differences of the analytic gradient vs mgcv's outer.info$hess\n"
)
t0 <- proc.time()[[3]]
H <- ncv_hessian(P)
cat("time", proc.time()[[3]] - t0, "\n")
print(round(H, 4))
print(round(P$H_mgcv, 4))

cat("\n(4) curve scores vs central differences in the curve weight\n")
t0 <- proc.time()[[3]]
s <- ncv_scores(P, st)
t_s <- proc.time()[[3]] - t0
cat("scores time", t_s, "\n")
G <- length(P$groups)
fd_s <- t(vapply(
  seq_len(n_fd),
  \(h) {
    wp <- replace(rep(1, G), h, 1 + 1e-5)
    wm <- replace(rep(1, G), h, 1 - 1e-5)
    (ncv_grad(ncv_state(P, w = wp)) - ncv_grad(ncv_state(P, w = wm))) / 2e-5
  },
  numeric(K)
))
print(rbind(analytic = s[seq_len(n_fd), ], fd = fd_s))

cat(
  "\n(5) leave-curve-out changes of rho-hat: IJ and one-step vs exact (Newton on NCV_(-h))\n"
)
t0 <- proc.time()[[3]]
os <- ncv_onestep(P, H)
t_os <- proc.time()[[3]] - t0
t0 <- proc.time()[[3]]
jk <- t(vapply(
  seq_len(G),
  \(h) ncv_refit_rho(P, replace(rep(1, G), h, 0)) - log(P$sp),
  numeric(K)
))
t_jk <- proc.time()[[3]] - t0
ij <- -ncv_ij_influence(P, H, s) # leave-out orientation (centred)
cat(sprintf("one-step %.1f s, exact %.1f s\n", t_os, t_jk))
cmp <- \(a, b, lab)
  data.frame(
    what = lab,
    par = seq_len(ncol(a)),
    cor = diag(cor(a, b)),
    sd_ratio = apply(a, 2, sd) / apply(b, 2, sd)
  )
print(rbind(cmp(ij, jk, "IJ"), cmp(os, jk, "one-step")), digits = 3)

cat(
  "\n(6) mgcv zero-weight refit for curve 1 vs exact Newton on my criterion\n"
)
b <- refund:::pffr_model_based_gam(fit)
G0 <- tryCatch(
  mgcv::gam(
    b$formula,
    data = b$model,
    family = b$family,
    method = "NCV",
    nei = refund:::pffr_ncv_nei(P$cluster),
    fit = FALSE
  ),
  error = \(e) e
)
if (!inherits(G0, "error")) {
  G0$w <- replace(rep(1, length(P$y)), P$cluster == P$groups[1], 0)
  m <- tryCatch(
    mgcv::gam(
      G = G0,
      method = "NCV",
      nei = refund:::pffr_ncv_nei(P$cluster),
      in.out = list(sp = P$sp, scale = fit$sig2),
      control = mgcv::gam.control(
        epsilon = 1e-12,
        newton = list(conv.tol = 1e-9)
      )
    ),
    error = \(e) e
  )
  if (!inherits(m, "error"))
    print(rbind(mgcv = log(m$sp) - log(P$sp), mine = jk[1, ])) else
    cat("mgcv refit failed:", conditionMessage(m), "\n")
} else cat("gam setup failed:", conditionMessage(G0), "\n")

cat("\nSD of rho-hat: naive (mgcv hess) / IJ / one-step / exact jackknife\n")
print(
  rbind(
    naive = sqrt(pmax(diag(pinv_sym(P$H_mgcv)), 0)),
    mine_hess = sqrt(pmax(diag(pinv_sym(H)), 0)),
    ij = sqrt(diag(omega_cov(ij, "ij"))),
    os = sqrt(diag(omega_cov(os, "jk"))),
    jk = sqrt(diag(omega_cov(jk, "jk")))
  ),
  digits = 3
)
