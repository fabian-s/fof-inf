# Cluster-robust smoothing-parameter correction for curve-blocked NCV fits of canonical-link
# GLMs (Poisson log, binomial logit; Gaussian identity is the special case W' = 0). Sourced
# after lc-lib.R and ncv-lib.R. Notation as in NOTE-ncvg.qmd.
#
# mgcv's NCV for a GLM (gam.fit3 + C routine Rncv), with prior weights w (one per curve):
#   beta = beta(rho, w), the penalised fit;  u = y - mu,  W = dmu/deta (canonical link);
#   H = X' diag(w W) X + S_lambda;
#   leaving out curve g is one Newton step from beta:
#     H_g = H - w_g X_g' diag(W_g) X_g,  delta_g = H_g^-1 X_g' (w_g u_g),
#     eta_(-g) = X_g (beta - delta_g);
#   NCV(rho, w) = sum_g w_g sum_{i in g} dev(y_i, mu(eta_(-g), i)).
# A zero weight deletes curve g from the fit, from every other fold and from the loss.

ncvg_pieces <- function(fit) {
  if (!identical(fit$method, "NCV")) stop("ncvg_pieces needs an NCV fit")
  if (any(fit$prior.weights != 1)) stop("prior weights must be 1")
  if (any(fit$offset != 0)) stop("offsets are not supported")
  canonical_weights(fit$family, 0.5) # stops for unsupported family/link
  fit$method <- "REML" # lc_pieces() only reads the method for its own check
  P <- lc_pieces(fit, scale_in_H = FALSE) # NCV: no scale parameter in the outer optimisation
  P$H_mgcv <- P$H
  P$idx <- lapply(P$groups, \(g) which(P$cluster == g))
  P$row_curve <- match(P$cluster, P$groups)
  P
}

# Penalised fit at (Slam, row weights wi) by Newton (= Fisher for canonical links), from beta.
# Converged at relative step tol, or after 5 iterations at relative step floor (very large
# smoothing parameters make H ill-conditioned; steps then stall at rounding level ~1e-9).
ncvg_beta <- function(
  P,
  Slam,
  wi,
  beta,
  tol = 1e-12,
  floor = 1e-6,
  maxit = 30
) {
  X <- P$X
  for (it in seq_len(maxit)) {
    mu <- P$family$linkinv(drop(X %*% beta))
    W <- canonical_weights(P$family, mu)$W
    grad <- drop(crossprod(X, wi * (P$y - mu))) - drop(Slam %*% beta)
    step <- drop(solve_pd(crossprod(X * sqrt(wi * W)) + Slam, grad))
    beta <- beta + step
    rel <- max(abs(step)) / (1 + max(abs(beta)))
    if (rel < tol || (it >= 5 && rel < floor)) return(beta)
  }
  stop("ncvg_beta: Newton iterations did not converge")
}

# Everything at (sp, w): the fit, its Hessian, and per curve the fold inverse Hessian,
# the Newton step and the cross-validated predictions.
ncvg_state <- function(
  P,
  sp = P$sp,
  w = rep(1, length(P$groups)),
  beta = P$theta
) {
  X <- P$X
  wi <- w[P$row_curve]
  Sk <- Map(`*`, P$Sk, sp)
  Slam <- Reduce(`+`, Sk)
  beta <- ncvg_beta(P, Slam, wi, beta)
  mu <- P$family$linkinv(drop(X %*% beta))
  cw <- canonical_weights(P$family, mu)
  u <- P$y - mu
  H <- crossprod(X * sqrt(wi * cw$W)) + Slam
  cur <- lapply(seq_along(P$idx), \(j) {
    ix <- P$idx[[j]]
    Xg <- X[ix, , drop = FALSE]
    Hg <- H - w[j] * crossprod(Xg * sqrt(cw$W[ix]))
    Hg_inv <- chol2inv(chol(Hg))
    delta <- drop(Hg_inv %*% crossprod(Xg, w[j] * u[ix]))
    eta_cv <- drop(Xg %*% (beta - delta))
    mu_cv <- P$family$linkinv(eta_cv)
    list(
      Hinv = Hg_inv,
      delta = delta,
      Xdelta = drop(Xg %*% delta),
      mu = mu_cv,
      dev = sum(P$family$dev.resids(P$y[ix], mu_cv, 1))
    )
  })
  list(
    sp = sp,
    w = w,
    wi = wi,
    Sk = Sk,
    beta = beta,
    mu = mu,
    u = u,
    W = cw$W,
    W1 = cw$W1,
    Hinv = chol2inv(chol(H)),
    cur = cur
  )
}

ncvg_value <- function(st) sum(st$w * vapply(st$cur, \(c) c$dev, 0))

# Directional derivatives of NCV at st along the columns of Drho (K x m, in rho = log sp) and
# Dw (G x m, curve weights), jointly: one value per column. Exact (implicit differentiation
# of the penalised fit; for a canonical link d(w u)/d eta = -w W and d(w W)/d eta = w W').
ncvg_dderiv <- function(P, st, Drho, Dw) {
  X <- P$X
  m <- ncol(Drho)
  dwi <- Dw[P$row_curve, , drop = FALSE] # n x m
  Sb <- sapply(st$Sk, \(S) drop(S %*% st$beta)) # p x K
  dbeta <- st$Hinv %*% (crossprod(X, dwi * st$u) - Sb %*% Drho) # p x m
  deta <- X %*% dbeta
  dwW <- dwi * st$W + st$wi * st$W1 * deta # d(w W), n x m
  dwu <- dwi * st$u - st$wi * st$W * deta # d(w u), n x m
  Delta <- sapply(st$cur, \(c) c$delta) # p x G
  # dH_a Delta for every direction a: p x G x m
  dHD <- array(0, c(nrow(Delta), ncol(Delta), m))
  for (a in seq_len(m)) {
    dS <- Reduce(`+`, Map(`*`, st$Sk, Drho[, a]))
    dHD[,, a] <- (crossprod(X, dwW[, a] * X) + dS) %*% Delta
  }
  out <- numeric(m)
  for (j in seq_along(P$idx)) {
    ix <- P$idx[[j]]
    Xg <- X[ix, , drop = FALSE]
    c <- st$cur[[j]]
    rhs <- crossprod(Xg, dwu[ix, , drop = FALSE]) -
      matrix(dHD[, j, ], ncol = m) +
      crossprod(Xg, dwW[ix, , drop = FALSE] * c$Xdelta)
    ddelta <- c$Hinv %*% rhs
    deta_cv <- Xg %*% (dbeta - ddelta)
    out <- out +
      st$w[j] * colSums(-2 * (P$y[ix] - c$mu) * deta_cv) +
      Dw[j, ] * c$dev
  }
  out
}

ncvg_grad <- function(P, st) {
  K <- length(P$sp)
  ncvg_dderiv(P, st, diag(K), matrix(0, length(P$groups), K))
}

# Hessian (K x K) and curve scores s_hk = d^2 NCV / d rho_k d w_h (G x K) at (rho-hat, w = 1),
# by central differences in rho of the exact directional derivatives (2K refits of beta).
# h = 1e-3: the refits are exact only to ~1e-9 when some lambda is very large; truncation
# error O(h^2) is far below that noise.
ncvg_hessian_scores <- function(P, h = 1e-3) {
  K <- length(P$sp)
  G <- length(P$groups)
  Drho <- cbind(diag(K), matrix(0, K, G))
  Dw <- cbind(matrix(0, G, K), diag(G))
  d <- lapply(seq_len(K), \(k) {
    e <- replace(rep(0, K), k, h)
    (ncvg_dderiv(P, ncvg_state(P, P$sp * exp(e)), Drho, Dw) -
      ncvg_dderiv(P, ncvg_state(P, P$sp * exp(-e)), Drho, Dw)) /
      (2 * h)
  })
  H <- sapply(d, \(x) x[seq_len(K)])
  s <- sapply(d, \(x) x[K + seq_len(G)])
  list(H = (H + t(H)) / 2, s = matrix(s, G, K))
}

# Arms of a GLM NCV fit: CL2, mgcv-style naive, IJ and IJ with the cross term.
BOUNDARY_TOL <- 1e-3

ncvg_covariances <- function(fit) {
  tm <- c()
  tick <- \(label, expr) {
    t0 <- proc.time()[["elapsed"]]
    val <- force(expr)
    tm[[label]] <<- proc.time()[["elapsed"]] - t0
    val
  }
  P <- tick("pieces", ncvg_pieces(fit))
  core <- tick(
    "cl2",
    refund:::pffr_influence(fit, sandwich = "cl2", cl2_adjustment = "exact")
  )
  if (!identical(as.character(core$groups), as.character(P$groups))) {
    stop("CL2 influence columns are not in the order of P$groups")
  }
  V_cl2 <- core$correction * tcrossprod(core$K) + core$B2
  st <- tick("state", ncvg_state(P))
  grad <- ncvg_grad(P, st)
  hs <- tick("hessian_scores", ncvg_hessian_scores(P))
  # smoothing parameters on a flat boundary of the criterion (relative curvature < 1e-3; e.g.
  # lambda ~ 1e6, where the fit no longer depends on rho_k) are treated as fixed, as mgcv's
  # Newton step excludes them; their score columns are rounding noise
  free <- diag(hs$H) > BOUNDARY_TOL * max(diag(hs$H))
  D_ij <- matrix(0, length(P$groups), length(P$sp))
  D_ij[, free] <- ncv_ij_influence(
    P,
    hs$H[free, free, drop = FALSE],
    hs$s[, free, drop = FALSE]
  )
  R <- rho_part_linear(P, D_ij)
  list(
    V = list(
      cl2 = V_cl2,
      naive_mgcv = V_cl2 + naive_correction(P),
      naive = V_cl2 + naive_correction(replace(P, "H", list(hs$H))),
      ij_lin = add_rho_part(core, R, "ij"),
      ij_linx = add_rho_part(core, R, "ij", cross = TRUE)
    ),
    omega = list(
      hat = log(P$sp),
      naive_mgcv = sqrt(pmax(diag(pinv_sym(P$H_mgcv)), 0)),
      naive = sqrt(pmax(diag(pinv_sym(hs$H)), 0)),
      ij = sqrt(diag(omega_cov(D_ij, "ij"))),
      boundary = !free
    ),
    grad = grad,
    check = c(
      P$check,
      value = abs(ncvg_value(st) - fit$gcv.ubre) / abs(fit$gcv.ubre),
      grad = max(abs(grad - drop(fit$outer.info$grad))) /
        max(1, abs(fit$gcv.ubre))
    ),
    time = tm
  )
}

# Exact leave-curve-out NCV smoothing parameters (validation reference): minimise
# NCV(rho; w) over rho by Newton steps with the exact gradient and a central-difference
# Hessian, from rho-hat (as ncv_refit_rho; mgcv's own zero-weight NCV refits stall in BFGS).
ncvg_refit_rho <- function(P, w, tol = 1e-8, maxit = 50) {
  rho <- log(P$sp)
  K <- length(rho)
  beta <- P$theta
  at <- \(r) ncvg_state(P, exp(r), w, beta)
  for (it in seq_len(maxit)) {
    st <- at(rho)
    beta <- st$beta
    g <- ncvg_grad(P, st)
    f0 <- ncvg_value(st)
    if (max(abs(g)) < tol * (1 + f0))
      return(structure(
        rho,
        converged = TRUE,
        reason = "gradient",
        iterations = it
      ))
    H <- sapply(seq_len(K), \(k) {
      e <- replace(rep(0, K), k, 1e-4)
      (ncvg_grad(P, at(rho + e)) - ncvg_grad(P, at(rho - e))) / 2e-4
    })
    step <- -drop(pinv_sym((H + t(H)) / 2, 1e-6) %*% g)
    step <- step * min(1, 2 / max(abs(step))) # trust region on log scale
    accepted <- FALSE
    for (half in 0:20) {
      if (ncvg_value(at(rho + step / 2^half)) <= f0) {
        accepted <- TRUE
        break
      }
    }
    if (!accepted)
      return(structure(
        rho,
        converged = FALSE,
        reason = "line search failed",
        iterations = it
      ))
    rho <- rho + step / 2^half
    if (max(abs(step / 2^half)) < 1e-7) {
      st <- at(rho)
      g <- ncvg_grad(P, st)
      return(structure(
        rho,
        converged = max(abs(g)) < 1e-4 * (1 + ncvg_value(st)),
        reason = "small step",
        iterations = it
      ))
    }
  }
  structure(rho, converged = FALSE, reason = "maxit", iterations = maxit)
}
