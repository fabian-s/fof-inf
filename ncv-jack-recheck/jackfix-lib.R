# mgcv's neighbourhood-aware NCV covariance for curve-blocked NCV fits, computed from the
# block fit (no refit, no mgcv change).
#
# mgcv 1.9-5 (gam.outer, mgcv.r ~l.1729-1803) computes its covariance only for the
# per-point neighbourhood form (each point i predicted on its own, d = 1:n, md = 1:n,
# neighbourhood a(i) = the points dropped when predicting i). For the block form refund
# uses (one fold per curve: a(g) = d(g) = the curve's points) it returns Vp unless the
# user sets jackknife = TRUE, and then it passes the n x p matrix of single-point
# perturbations to neicov(), which reads it as if it had one row per neighbourhood (G).
#
# The NCV criterion of the block form and of the per-point form with a(i) = curve(i) is
# the same sum, so both select the same smoothing parameters; the per-point form costs
# D times more per criterion evaluation. jackfix_vcov() returns, from the block fit,
# what mgcv returns for the per-point form:
#   dd_i  = jkw (H - w_i x_i x_i')^-1 x_i u_i          single-point Newton perturbation,
#           H = X'WX + S_lambda, u_i = d loglik_i / d eta_i, jkw = ((n - 1)/n)^0.4
#   rat_i = r_i^cv / r_i                               deviance residuals at the
#           leave-curve-out fit (one Newton step) and at the full fit
#   Vraw  = sum_g (sum_{i in g} dd_i)(...)' n/(n - edf)
#   Vcv   = sum_g (sum_{i in g} rat_i dd_i)(...)' n/(n - edf)
#   alpha = max(tr(Vraw+) / tr(Ve), 1)
#   Vj    = Vcv + alpha (Vp - Ve)  (projected onto the psd cone), Vc_j = Vc - Vp + Vj
# and, as `raw`, mgcv's jackknife = TRUE estimate for the per-point form:
#   dd_g  = (H - X_g'W_g X_g)^-1 X_g'u_g               leave-curve-out Newton step
#   Vjk   = sum_g ((n - D_g)/n)^0.8 dd_g dd_g'  (each curve's fold appears D_g times with
#           weight ((n - D_g)/n)^0.4 / sqrt(D_g)),  alpha_raw = tr(Vjk)/tr(Ve),
#   Vraw_j = Vjk + alpha_raw (Vp - Ve)

# Per-point neighbourhoods: point i is predicted with all points of its curve dropped.
nei_pointwise <- function(cluster_id) {
  groups <- split(seq_along(cluster_id), factor(cluster_id, unique(cluster_id)))
  members <- groups[match(cluster_id, unique(cluster_id))]
  n <- length(cluster_id)
  list(
    a = as.integer(unlist(members, use.names = FALSE)),
    ma = as.integer(cumsum(lengths(members))),
    d = seq_len(n),
    md = seq_len(n)
  )
}

# Penalty matrix sum_j sp_j S_j in the fit's coefficient order (gam/pffr objects,
# single linear predictor, no linked smoothing parameters).
penalty_matrix <- function(fit) {
  p <- length(fit$coefficients)
  S <- matrix(0, p, p)
  sp <- fit$full.sp %||% fit$sp
  for (sm in fit$smooth) {
    if (!length(sm$S)) next
    idx <- sm$first.para:sm$last.para
    for (j in seq_along(sm$S)) {
      S[idx, idx] <- S[idx, idx] + sp[sm$first.sp + j - 1] * sm$S[[j]]
    }
  }
  if (length(fit$paraPen)) stop("paraPen not supported")
  S
}

# Model matrix of the fitting data (pffr and gam objects).
fit_model_matrix <- function(fit) {
  X <- mgcv::predict.gam(fit, type = "lpmatrix")
  off <- fit$offset %||% 0
  stopifnot(
    nrow(X) == length(fit$y),
    max(abs(drop(X %*% fit$coefficients) + off - fit$linear.predictors)) <
      1e-6 * (1 + max(abs(fit$linear.predictors)))
  )
  X
}

# cluster: one id per row of the fitting data (the curves). Returns Vp/Vc of mgcv's
# neighbourhood-aware estimator ("nei") and of its jackknife = TRUE variant ("raw"),
# plus components for checks.
jackfix_vcov <- function(fit, cluster, X = fit_model_matrix(fit)) {
  fam <- fit$family
  if (inherits(fam, "general.family") || inherits(fam, "extended.family")) {
    stop("jackfix_vcov: exponential families only")
  }
  n <- nrow(X)
  stopifnot(length(cluster) == n)
  y <- fit$y
  mu <- fit$fitted.values
  eta <- fit$linear.predictors
  pw <- fit$prior.weights
  w <- fit$weights # working weights at convergence (Fisher = Newton: canonical links)
  # d loglik_i / d eta_i (mgcv's ww): prior weight (y - mu) mu.eta / V(mu)
  u <- pw * (y - mu) * fam$mu.eta(eta) / fam$variance(mu)
  H <- crossprod(X, w * X) + penalty_matrix(fit)
  R <- chol(H)
  HiX <- backsolve(R, forwardsolve(t(R), t(X))) # H^-1 X', p x n
  lev <- w * colSums(t(X) * HiX) # w_i x_i' H^-1 x_i
  jkw <- ((n - 1) / n)^0.4
  dd <- t(HiX) * (jkw * u / (1 - lev)) # n x p single-point perturbations
  # leave-curve-out Newton steps and CV linear predictor
  groups <- split(seq_len(n), factor(cluster, unique(cluster)))
  ddg <- t(vapply(
    groups,
    \(ii) {
      Xg <- X[ii, , drop = FALSE]
      drop(solve(H - crossprod(Xg, w[ii] * Xg), crossprod(Xg, u[ii])))
    },
    numeric(ncol(X))
  ))
  eta_cv <- numeric(n)
  for (k in seq_along(groups)) {
    ii <- groups[[k]]
    eta_cv[ii] <- eta[ii] - drop(X[ii, , drop = FALSE] %*% ddg[k, ])
  }
  dev_res <- \(m) sign(y - m) * sqrt(pmax(fam$dev.resids(y, m, pw), 0))
  rsd0 <- dev_res(mu)
  rsd1 <- dev_res(fam$linkinv(eta_cv))
  rsd1[!is.finite(rsd1)] <- rsd0[!is.finite(rsd1)]
  rat <- rsd1 / rsd0
  rat[!is.finite(rat)] <- 1
  sum_by <- \(M) rowsum(M, factor(cluster, unique(cluster)), reorder = FALSE)
  fac <- n / (n - sum(fit$edf))
  Vraw <- crossprod(sum_by(dd)) * fac
  Vcv <- crossprod(sum_by(dd * rat)) * fac
  ev <- eigen(Vraw, symmetric = TRUE, only.values = TRUE)$values
  alpha <- max(sum(ev[ev > 0]) / sum(diag(fit$Ve)), 1)
  B <- fit$Vp - fit$Ve
  Vj <- psd_project(Vcv + B * alpha)
  # mgcv jackknife = TRUE for the per-point form
  Dg <- lengths(groups)
  Vjk <- crossprod(ddg * ((n - Dg) / n)^0.4)
  alpha_raw <- sum(diag(Vjk)) / sum(diag(fit$Ve))
  Vjr <- Vjk + B * alpha_raw
  Vc0 <- fit$Vc %||% fit$Vp
  list(
    nei = list(Vp = Vj, Vc = Vc0 - fit$Vp + Vj, alpha = alpha),
    raw = list(Vp = Vjr, Vc = Vc0 - fit$Vp + Vjr, alpha = alpha_raw),
    parts = list(
      H = H,
      dd = dd,
      ddg = ddg,
      eta_cv = eta_cv,
      rat = rat,
      Vraw = Vraw,
      Vcv = Vcv,
      Vjk = Vjk
    )
  )
}

psd_project <- function(V) {
  ev <- eigen(V, symmetric = TRUE)
  if (all(ev$values >= 0)) return(V)
  ev$values[ev$values < 0] <- 0
  crossprod(sqrt(ev$values) * t(ev$vectors))
}

# Exact leave-one-curve-out refits at fixed smoothing parameters (penalized IRLS with the
# fit's H structure): returns the G x p matrix of beta_(-g) - beta.
exact_loco <- function(fit, cluster, X = fit_model_matrix(fit), maxit = 100) {
  fam <- fit$family
  S <- penalty_matrix(fit)
  beta <- fit$coefficients
  groups <- split(seq_len(nrow(X)), factor(cluster, unique(cluster)))
  t(vapply(
    groups,
    \(drop) {
      keep <- setdiff(seq_len(nrow(X)), drop)
      Xk <- X[keep, , drop = FALSE]
      yk <- fit$y[keep]
      pwk <- fit$prior.weights[keep]
      b <- beta
      for (it in seq_len(maxit)) {
        e <- drop(Xk %*% b)
        m <- fam$linkinv(e)
        me <- fam$mu.eta(e)
        wk <- pwk * me^2 / fam$variance(m)
        z <- e + (yk - m) / me
        b1 <- drop(solve(crossprod(Xk, wk * Xk) + S, crossprod(Xk, wk * z)))
        if (max(abs(b1 - b)) < 1e-12 * (1 + max(abs(b)))) {
          b <- b1
          break
        }
        b <- b1
      }
      b - beta
    },
    numeric(length(beta))
  ))
}

# Study tasks -------------------------------------------------------------------------

# Curve id per fitted row, from the fit's NCV blocks (refund's block neighbourhoods).
ncv_cluster <- function(fit) {
  nei <- fit$pffr$ncv$nei
  stopifnot(!is.null(nei), identical(sort(nei$a), seq_along(nei$a)))
  id <- integer(length(nei$a))
  id[nei$a] <- rep(seq_along(nei$ma), diff(c(0L, nei$ma)))
  id
}

# One (cell, rep): NCV fit as in the study, its covariances, REML + CL2 on the same data.
# Returns link-scale estimates and SEs per evaluation point and arm.
#   ncv_vp      NCV, Bayesian Vp (the study's model-based NCV interval)
#   ncv_cl2     NCV + CL2 (exact, Bayesian form)
#   jk_nei      mgcv's neighbourhood-aware covariance (per-point form, jackfix_vcov)
#   jk_nei_vc   the same plus mgcv's smoothing-parameter correction (Vc - Vp + Vj)
#   jk_cv       its sampling part alone (Vcv, without alpha (Vp - Ve))
#   jk_raw      mgcv's jackknife = TRUE for the per-point form
#   reml_cl2    REML + CL2 (the study's recipe)
jackfix_task <- function(cell, rep, simulate_fun = simulate) {
  sim <- simulate_fun(cell, rep)
  fit <- fit_model(sim, cell, "NCV")
  t0 <- proc.time()[["elapsed"]]
  jf <- jackfix_vcov(fit, ncv_cluster(fit))
  t_jack <- proc.time()[["elapsed"]] - t0
  t0 <- proc.time()[["elapsed"]]
  Vcl2 <- refund:::pffr_vcov(
    fit,
    sandwich = "cl2",
    cl2_adjustment = "exact",
    freq = FALSE
  )
  t_cl2 <- proc.time()[["elapsed"]] - t0
  fitR <- fit_model(sim, cell, "REML")
  VR <- refund:::pffr_vcov(
    fitR,
    sandwich = "cl2",
    cl2_adjustment = "exact",
    freq = FALSE
  )
  covs <- list(
    ncv_vp = fit$Vp,
    ncv_cl2 = Vcl2,
    jk_nei = jf$nei$Vp,
    jk_nei_vc = jf$nei$Vc,
    jk_cv = jf$parts$Vcv,
    jk_raw = jf$raw$Vp
  )
  rows <- list()
  for (e in model_estimands(cell$model)) {
    L <- evaluation_matrix(fit, sim, e)
    LR <- evaluation_matrix(fitR, sim, e)
    est <- drop(L %*% fit$coefficients)
    se <- lapply(covs, \(V) sqrt(pmax(rowSums((L %*% V) * L), 0)))
    se$reml_cl2 <- sqrt(rowSums((LR %*% VR) * LR))
    truth <- as.vector(sim$truth[[e]])
    for (a in names(se)) {
      rows[[length(rows) + 1]] <- data.frame(
        cell = cell$cell,
        rep = rep,
        estimand = e,
        arm = a,
        point = seq_along(truth),
        truth = truth,
        est = if (a == "reml_cl2") drop(LR %*% fitR$coefficients) else est,
        se = se[[a]]
      )
    }
  }
  list(
    points = do.call(rbind, rows),
    info = data.frame(
      cell = cell$cell,
      rep = rep,
      time_ncv = fit$restart$time,
      time_reml = fitR$restart$time,
      time_jack = t_jack,
      time_cl2 = t_cl2,
      alpha_nei = jf$nei$alpha,
      alpha_raw = jf$raw$alpha,
      tr_ratio_nei = sum(diag(jf$nei$Vp)) / sum(diag(fit$Vp)),
      tr_ratio_raw = sum(diag(jf$raw$Vp)) / sum(diag(fit$Vp)),
      edf = sum(fit$edf),
      p = length(fit$coefficients)
    ),
    time = fit$restart$time + fitR$restart$time + t_jack + t_cl2
  )
}
