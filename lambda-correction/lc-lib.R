# Cluster-robust smoothing-parameter-uncertainty correction for CL2 (prototype).
# Sourced after scripts/setup.R. Takes a sandwich = "none" REML pffr fit with a canonical
# link (gaussian identity, poisson log, binomial logit). Notation as in NOTE.qmd:
# outer parameters omega = (rho, tau) with rho = log(sp) and, for the Gaussian family,
# tau = log(scale); V(omega) = mgcv's REML criterion; H = its Hessian (outer.info$hess);
# J = d theta-hat / d rho (fit$db.drho); A = X'WX + S_lambda; curve weights w_g.

# Fit pieces ----------------------------------------------------------------------------

# Canonical-link working weight W(eta) and its first two eta-derivatives.
canonical_weights <- function(family, mu) {
  switch(
    paste(family$family, family$link),
    "gaussian identity" = list(W = 1 + 0 * mu, W1 = 0 * mu, W2 = 0 * mu),
    "poisson log" = list(W = mu, W1 = mu, W2 = mu),
    "binomial logit" = {
      v <- mu * (1 - mu)
      list(W = v, W1 = v * (1 - 2 * mu), W2 = v * (1 - 6 * mu + 6 * mu^2))
    },
    stop(
      "lambda correction: unsupported family/link ",
      family$family,
      "/",
      family$link
    )
  )
}

# Everything the correction needs from a sandwich = "none" REML fit, with consistency checks
# against mgcv (Vp = phi A^-1, db.drho = -A^-1 lambda_k S_k theta).
lc_pieces <- function(fit, check = TRUE, scale_in_H = NULL) {
  if (!identical(fit$method, "REML")) stop("lambda correction needs a REML fit")
  b <- refund:::pffr_model_based_gam(fit)
  X <- model.matrix(b)
  theta <- unname(fit$coefficients)
  p <- length(theta)
  sp <- fit$sp
  # embedded penalty matrices, in fit$sp order
  Sk <- vector("list", length(sp))
  for (sm in fit$smooth) {
    if (!is.null(fit$L) || length(sm$S) != sm$last.sp - sm$first.sp + 1) {
      stop(
        "lambda correction: linked or fixed smoothing parameters not supported"
      )
    }
    for (j in seq_along(sm$S)) {
      M <- matrix(0, p, p)
      idx <- sm$first.para:sm$last.para
      M[idx, idx] <- sm$S[[j]]
      Sk[[sm$first.sp + j - 1]] <- M
    }
  }
  Slam <- Reduce(`+`, Map(`*`, Sk, sp))
  gaussian <- identical(fit$family$family, "gaussian")
  phi <- if (gaussian) fit$sig2 else 1
  if (!gaussian && !isFALSE(fit$scale.estimated))
    stop("scale must be known for GLMs")
  mu <- as.vector(fit$fitted.values)
  ww <- canonical_weights(fit$family, mu)
  A <- crossprod(X * sqrt(ww$W)) + Slam
  Ainv <- chol2inv(chol(A))
  J <- -Ainv %*% sapply(seq_along(sp), \(k) sp[k] * (Sk[[k]] %*% theta))
  H <- fit$outer.info$hess
  if (nrow(H) != length(sp) + (scale_in_H %||% gaussian))
    stop("unexpected outer Hessian dimension")
  if (check) {
    rel <- \(a, b) max(abs(a - b)) / max(abs(b))
    chk <- c(Vp = rel(phi * Ainv, fit$Vp), J = rel(J, fit$db.drho))
    # J is exact for the fitted coefficients; mgcv's db.drho comes from its own fit (for NCV:
    # BFGS path) and agrees only to ~1e-3, so it is a gross-error check (recorded in P$check)
    if (chk[["Vp"]] > 1e-5 || chk[["J"]] > 0.05) {
      stop(
        "lc_pieces consistency check failed: ",
        paste(names(chk), signif(chk, 3), collapse = ", ")
      )
    }
  }
  cid <- refund:::build_cluster_id(fit$pffr)
  list(
    check = if (check) chk,
    X = X,
    y = as.vector(fit$y),
    mu = mu,
    W = ww$W,
    W1 = ww$W1,
    W2 = ww$W2,
    theta = theta,
    sp = sp,
    Sk = Sk,
    Slam = Slam,
    A = A,
    Ainv = Ainv,
    J = J,
    H = H,
    phi = phi,
    gaussian = gaussian,
    family = fit$family,
    cluster = cid,
    groups = unique(cid),
    b = b
  )
}

REG_PRIOR <- 1 / 10 # prior precision on log smoothing parameters (as mgcv's V.sp)

# Pseudo-inverse of a symmetric matrix, non-positive eigenvalues dropped (as mgcv's
# first-order Vc term). prior > 0: inverse of H + prior I after zeroing non-positive
# eigenvalues, exactly mgcv's V.sp for its second-order term (prior = 1/10). This bounds the
# influence along near-flat directions (smoothing parameters drifting to a boundary), but
# gives directions with non-positive eigenvalues weight 1/prior where prior = 0 drops them,
# so it is not a pure shrinkage of the unregularised version.
pinv_sym <- function(H, prior = 0) {
  e <- eigen((H + t(H)) / 2, symmetric = TRUE)
  v <- pmax(e$values, 0) + prior
  d <- ifelse(v > 0, 1 / v, 0)
  e$vectors %*% (d * t(e$vectors))
}

# (a) Infinitesimal jackknife: analytic curve contributions -----------------------------

# s_g = d grad_omega V / d w_g at w = 1 (G x dim(omega)); see NOTE.qmd, eq. (s-glm).
# n_term = FALSE drops the Gaussian n_g / 2 term (mgcv's prior weights leave n fixed; used
# only to compare with finite differences in prior weights).
ij_scores <- function(P, n_term = TRUE) {
  X <- P$X
  K <- length(P$sp)
  r_eta <- (P$y - P$mu) # canonical: W (z - eta) = y - mu
  h <- rowSums((X %*% P$Ainv) * X)
  Xth <- X %*% P$J # x_i' theta_k, n x K
  Ak <- lapply(
    seq_len(K),
    \(k) P$sp[k] * P$Sk[[k]] + crossprod(X, (P$W1 * Xth[, k]) * X)
  )
  q <- sapply(
    seq_len(K),
    \(k) rowSums((X %*% (P$Ainv %*% Ak[[k]] %*% P$Ainv)) * X)
  )
  a_til <- P$Ainv %*% crossprod(X, P$W1 * h)
  Xa <- drop(X %*% a_til)
  SkTh <- sapply(seq_len(K), \(k) P$sp[k] * drop(P$Sk[[k]] %*% P$theta)) # p x K
  Ska <- sapply(seq_len(K), \(k) P$sp[k] * drop(P$Sk[[k]] %*% a_til))
  out <- t(vapply(
    P$groups,
    \(g) {
      ig <- which(P$cluster == g)
      delta <- drop(P$Ainv %*% crossprod(X[ig, , drop = FALSE], r_eta[ig]))
      Xd <- drop(X %*% delta)
      bvec <- P$W1 * Xd
      bvec[ig] <- bvec[ig] + P$W[ig]
      s_rho <- vapply(
        seq_len(K),
        \(k) {
          sum(SkTh[, k] * delta) /
            P$phi +
            0.5 *
              (-sum(bvec * q[, k]) +
                sum(P$W1[ig] * Xth[ig, k] * h[ig]) +
                sum(P$W2 * Xd * Xth[, k] * h) -
                sum(bvec * Xa * Xth[, k]) -
                sum(Ska[, k] * delta))
        },
        numeric(1)
      )
      if (P$gaussian) {
        c(s_rho, -sum(r_eta[ig]^2) / (2 * P$phi) + n_term * length(ig) / 2)
      } else s_rho
    },
    numeric(K + P$gaussian)
  ))
  if (K + P$gaussian == 1) out <- t(out)
  out
}

# Curve influences on omega-hat (rows: curves): d_g = d omega-hat / d w_g = -H^+ s_g,
# with s_g centred over curves (jackknife convention, see NOTE.qmd).
ij_influence <- function(P, s = ij_scores(P), center = TRUE, prior = 0) {
  if (center) s <- sweep(s, 2, colMeans(s))
  -s %*% pinv_sym(P$H, prior)
}

# (c) One-step leave-curve-out --------------------------------------------------------

# w-dependent part of the REML gradient in omega at coefficients theta and curve weights
# w (one per observation): the log|S|_+ term and other w-free terms are omitted, they
# cancel in differences.
reml_grad_w <- function(P, theta, w) {
  X <- P$X
  eta <- drop(X %*% theta)
  mu <- P$family$linkinv(eta)
  ww <- canonical_weights(P$family, mu)
  A <- crossprod(X * sqrt(w * ww$W)) + P$Slam
  Ainv <- chol2inv(chol(A))
  K <- length(P$sp)
  thk <- -Ainv %*% sapply(seq_len(K), \(k) P$sp[k] * (P$Sk[[k]] %*% theta))
  h <- rowSums((X %*% Ainv) * X)
  Xth <- X %*% thk
  g_rho <- vapply(
    seq_len(K),
    \(k) {
      P$sp[k] *
        sum(theta * (P$Sk[[k]] %*% theta)) /
        (2 * P$phi) +
        0.5 * (P$sp[k] * sum(Ainv * P$Sk[[k]]) + sum(w * ww$W1 * Xth[, k] * h))
    },
    numeric(1)
  )
  if (P$gaussian) {
    Dp <- sum(w * (P$y - mu)^2) + sum(theta * (P$Slam %*% theta))
    c(g_rho, -Dp / (2 * P$phi) + sum(w) / 2)
  } else g_rho
}

# One-step (Newton) leave-curve-out changes: theta at fixed rho by one Newton step (exact
# for the Gaussian family), then omega by one Newton step from omega-hat with the full-data
# Hessian. Returns list(d_theta = p x G (fixed rho), d_omega = G x dim(omega)).
# theta_out (p x G): leave-out coefficients at rho-hat to evaluate the leave-out gradient at
# (GLMs: the converged PIRLS solution; default the one Newton step). prior: see pinv_sym().
onestep_loo <- function(P, theta_out = NULL, prior = 0) {
  g0 <- reml_grad_w(P, P$theta, rep(1, length(P$y)))
  Hinv <- pinv_sym(P$H, prior)
  r <- P$y - P$mu
  res <- lapply(P$groups, \(g) {
    ig <- which(P$cluster == g)
    Xg <- P$X[ig, , drop = FALSE]
    Ag <- P$A - crossprod(Xg * sqrt(P$W[ig]))
    dth <- -solve_pd(Ag, drop(crossprod(Xg, r[ig])))
    w <- rep(1, length(P$y))
    w[ig] <- 0
    th <- if (is.null(theta_out)) P$theta + dth else
      theta_out[, match(g, P$groups)]
    grad <- reml_grad_w(P, th, w) - g0
    list(dth = dth, dom = -drop(Hinv %*% grad))
  })
  list(
    d_theta = sapply(res, `[[`, "dth"),
    d_omega = t(sapply(res, `[[`, "dom"))
  )
}

# (b) Exact curve jackknife (refits) ---------------------------------------------------

# mgcv setup of the fit (same parameterisation and constraints for every refit).
gam_setup <- function(fit) {
  b <- refund:::pffr_model_based_gam(fit)
  mf <- b$model
  G0 <- mgcv::gam(
    b$formula,
    data = mf,
    family = b$family,
    method = "REML",
    fit = FALSE
  )
  if (max(abs(G0$X - model.matrix(b))) > 1e-10) {
    stop("gam setup does not reproduce the fit's model matrix")
  }
  G0
}

# Refit with prior weights w (observation level), warm-started at the full fit.
refit_weighted <- function(fit, G0, w, tol = 1e-9) {
  G0$w <- w
  start <- c(
    log(fit$sp),
    if (identical(fit$family$family, "gaussian")) log(fit$sig2)
  )
  ctrl <- mgcv::gam.control(
    epsilon = 1e-12,
    newton = list(conv.tol = tol, maxNstep = 5)
  )
  m <- mgcv::gam(
    G = G0,
    method = "REML",
    control = ctrl,
    in.out = list(
      sp = fit$sp,
      scale = if (identical(fit$family$family, "gaussian")) fit$sig2 else 1
    )
  )
  list(
    theta = unname(m$coefficients),
    omega = c(
      log(m$sp),
      if (identical(fit$family$family, "gaussian")) log(m$sig2)
    ),
    converged = isTRUE(m$converged) && m$outer.info$conv == "full convergence",
    Vp = m$Vp,
    Ve = m$Ve
  )
}

jackknife_refits <- function(fit, P, G0 = gam_setup(fit), groups = P$groups) {
  lapply(groups, \(g) {
    w <- rep(1, length(P$y))
    w[P$cluster == g] <- 0
    refit_weighted(fit, G0, w)
  })
}

# Coefficients at fixed smoothing parameters --------------------------------------------

# Solve M x = b for symmetric M: Cholesky, with a vanishing ridge only if M is singular
# (penalized normal equations at extreme smoothing parameters, e.g. a one-step leave-out
# rho far outside the data's range).
solve_pd <- function(M, b) {
  for (eps in c(0, 10^seq(-12, -6))) {
    R <- tryCatch(
      chol(M + eps * mean(diag(M)) * diag(nrow(M))),
      error = \(e) NULL
    )
    if (!is.null(R)) return(drop(backsolve(R, forwardsolve(t(R), b))))
  }
  stop("solve_pd: matrix not positive definite even with a ridge")
}

# theta-hat(sp, w): penalized IRLS at fixed smoothing parameters sp with observation weights
# w, from `start` (exact solve for the Gaussian family; step-halving Newton otherwise).
theta_fixed_sp <- function(
  P,
  sp,
  w,
  start = P$theta,
  tol = 1e-10,
  maxit = 100
) {
  Slam <- Reduce(`+`, Map(`*`, P$Sk, sp))
  X <- P$X
  if (P$gaussian) {
    return(solve_pd(crossprod(X * sqrt(w)) + Slam, crossprod(X, w * P$y)))
  }
  pdev <- \(th) {
    mu <- P$family$linkinv(drop(X %*% th))
    sum(P$family$dev.resids(P$y, mu, w)) + sum(th * (Slam %*% th))
  }
  th <- start
  obj <- pdev(th)
  for (it in seq_len(maxit)) {
    mu <- P$family$linkinv(drop(X %*% th))
    W <- canonical_weights(P$family, mu)$W
    step <- solve_pd(
      crossprod(X * sqrt(w * W)) + Slam,
      drop(crossprod(X, w * (P$y - mu)) - Slam %*% th)
    )
    accepted <- FALSE
    for (half in 0:30) {
      new <- th + step / 2^half
      obj_new <- pdev(new)
      if (is.finite(obj_new) && obj_new <= obj + 1e-12 * abs(obj)) {
        accepted <- TRUE
        break
      }
    }
    if (!accepted) {
      warning("theta_fixed_sp: no step decreases the penalized deviance")
      return(th)
    }
    th <- new
    done <- abs(obj - obj_new) < tol * (1 + abs(obj)) &&
      max(abs(step)) < 1e-6 * (1 + max(abs(th)))
    obj <- obj_new
    if (done) return(th)
  }
  warning("theta_fixed_sp: no convergence")
  th
}

# Covariances --------------------------------------------------------------------------

# Pseudo-inverse-based omega covariances and the rho part of curve influences ------------

# mgcv's first-order term J [H^+]_rho,rho J'.
naive_correction <- function(P) {
  K <- length(P$sp)
  P$J %*% pinv_sym(P$H)[seq_len(K), seq_len(K), drop = FALSE] %*% t(P$J)
}

# Robust covariance of omega-hat from per-curve changes (G x dim(omega)): IJ form
# (G/(G-1)) or jackknife form ((G-1)/G), centred.
omega_cov <- function(D, form = c("ij", "jk")) {
  form <- match.arg(form)
  G <- nrow(D)
  Dc <- sweep(D, 2, colMeans(D))
  if (form == "ij") G / (G - 1) * crossprod(Dc) else (G - 1) / G * crossprod(Dc)
}

# Linear rho part of the curve influences: J d_g (p x G) from omega influences D (G x .).
rho_part_linear <- function(P, D) {
  P$J %*% t(D[, seq_along(P$sp), drop = FALSE])
}

# Add the rho part R (p x G, sign = influence of including curve g) to CL2:
# cross = FALSE: V_CL2 + f sum_g Rc_g Rc_g'; cross = TRUE: c sum_g (K_g + Rc_g)(...)' + B2,
# with Rc the centred R, f = G/(G-1) (IJ) or (G-1)/G (leave-out based), c = CL2's G/(G-1).
# In the cross form Rc is rescaled by sqrt(f / c), so the rho-rho part is f sum Rc Rc' in both
# forms (the cross part then carries sqrt(f c)).
add_rho_part <- function(core, R, form = c("ij", "jk"), cross = FALSE) {
  form <- match.arg(form)
  G <- ncol(R)
  Rc <- R - rowMeans(R)
  f <- if (form == "ij") G / (G - 1) else (G - 1) / G
  if (cross) {
    Rc <- sqrt(f / core$correction) * Rc
    return(core$correction * tcrossprod(core$K + Rc) + core$B2)
  }
  core$correction * tcrossprod(core$K) + core$B2 + f * tcrossprod(Rc)
}

# All covariance arms of one fit --------------------------------------------------------

# Returns list(V = named list of p x p covariances, omega = per-arm SDs of omega-hat,
# time = seconds per component, n_jk_nonconv). jackknife = FALSE skips the refits.
lc_covariances <- function(fit, jackknife = TRUE) {
  tm <- c()
  tick <- \(label, expr) {
    t0 <- proc.time()[["elapsed"]]
    val <- force(expr)
    tm[[label]] <<- proc.time()[["elapsed"]] - t0
    val
  }
  P <- tick("pieces", lc_pieces(fit))
  if (any(fit$offset != 0)) stop("offsets are not supported")
  core <- tick(
    "cl2",
    refund:::pffr_influence(fit, sandwich = "cl2", cl2_adjustment = "exact")
  )
  # the cross term pairs CL2's influence columns with P$groups by position
  if (!identical(as.character(core$groups), as.character(P$groups))) {
    stop("CL2 influence columns are not in the order of P$groups")
  }
  V_cl2 <- core$correction * tcrossprod(core$K) + core$B2
  K <- length(P$sp)
  n_om <- K + P$gaussian
  omega_hat <- c(log(P$sp), if (P$gaussian) c(log_scale = log(fit$sig2)))
  Hinv <- pinv_sym(P$H)

  D_ij <- tick("ij", ij_influence(P)) # + orientation, centred
  D_ij_reg <- ij_influence(P, prior = REG_PRIOR)
  os0 <- tick("onestep", onestep_loo(P)) # leave-out changes (theta: one Newton step)
  # theta at rho-hat without curve g (exact for Gaussian, PIRLS otherwise)
  w_out <- \(g) {
    w <- rep(1, length(P$y))
    w[P$cluster == g] <- 0
    w
  }
  th_out <- tick(
    "loo_theta",
    sapply(seq_along(P$groups), \(j) {
      start <- P$theta + os0$d_theta[, j]
      if (P$gaussian) start else
        theta_fixed_sp(P, P$sp, w_out(P$groups[j]), start)
    })
  )
  # GLMs: leave-out gradient at the converged leave-out theta (Gaussian: identical)
  os <- if (P$gaussian) os0 else tick("onestep2", onestep_loo(P, th_out))
  os_reg <- onestep_loo(P, if (P$gaussian) NULL else th_out, prior = REG_PRIOR)
  R_os_nl <- tick(
    "onestep_nl",
    -sapply(seq_along(P$groups), \(j) {
      # capped: one-step changes of |log sp| > 10 are outside any meaningful range
      sp_j <- exp(log(P$sp) + pmax(pmin(os$d_omega[j, seq_len(K)], 10), -10))
      theta_fixed_sp(P, sp_j, w_out(P$groups[j]), th_out[, j]) - th_out[, j]
    })
  )
  V <- list(
    cl2 = V_cl2,
    naive = V_cl2 + naive_correction(P),
    vc = V_cl2 + fit$Vc - fit$Vp,
    ij_lin = add_rho_part(core, rho_part_linear(P, D_ij), "ij"),
    ij_linx = add_rho_part(core, rho_part_linear(P, D_ij), "ij", cross = TRUE),
    os_lin = add_rho_part(core, rho_part_linear(P, -os$d_omega), "jk"),
    os_linx = add_rho_part(
      core,
      rho_part_linear(P, -os$d_omega),
      "jk",
      cross = TRUE
    ),
    os_nl = add_rho_part(core, R_os_nl, "jk"),
    os_nlx = add_rho_part(core, R_os_nl, "jk", cross = TRUE),
    ij_reg = add_rho_part(core, rho_part_linear(P, D_ij_reg), "ij"),
    ij_regx = add_rho_part(
      core,
      rho_part_linear(P, D_ij_reg),
      "ij",
      cross = TRUE
    ),
    os_reg = add_rho_part(core, rho_part_linear(P, -os_reg$d_omega), "jk"),
    os_regx = add_rho_part(
      core,
      rho_part_linear(P, -os_reg$d_omega),
      "jk",
      cross = TRUE
    )
  )
  sd_om <- list(
    naive = sqrt(pmax(diag(Hinv), 0)),
    ij = sqrt(diag(omega_cov(D_ij, "ij"))),
    os = sqrt(diag(omega_cov(os$d_omega, "jk"))),
    ij_reg = sqrt(diag(omega_cov(D_ij_reg, "ij"))),
    os_reg = sqrt(diag(omega_cov(os_reg$d_omega, "jk")))
  )
  n_nonconv <- NA_integer_
  if (jackknife) {
    G0 <- tick("jk_setup", gam_setup(fit))
    jk <- tick("jackknife", jackknife_refits(fit, P, G0))
    n_nonconv <- sum(!vapply(jk, `[[`, logical(1), "converged"))
    d_om_jk <- t(vapply(jk, `[[`, numeric(n_om), "omega")) -
      matrix(omega_hat, length(jk), n_om, byrow = TRUE)
    th_jk <- vapply(jk, `[[`, numeric(length(P$theta)), "theta")
    R_jk_nl <- -(th_jk - th_out)
    d_th <- th_jk - P$theta
    V <- c(
      V,
      list(
        jk_lin = add_rho_part(core, rho_part_linear(P, -d_om_jk), "jk"),
        jk_linx = add_rho_part(
          core,
          rho_part_linear(P, -d_om_jk),
          "jk",
          cross = TRUE
        ),
        jk_nl = add_rho_part(core, R_jk_nl, "jk"),
        jk_nlx = add_rho_part(core, R_jk_nl, "jk", cross = TRUE),
        jk_full = (ncol(d_th) - 1) /
          ncol(d_th) *
          tcrossprod(d_th - rowMeans(d_th)) +
          core$B2
      )
    )
    sd_om$jk <- sqrt(diag(omega_cov(d_om_jk, "jk")))
  }
  list(
    V = V,
    omega = c(list(hat = omega_hat), sd_om),
    time = tm,
    n_jk_nonconv = n_nonconv
  )
}
