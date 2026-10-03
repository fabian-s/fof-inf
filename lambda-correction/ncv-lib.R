# Cluster-robust smoothing-parameter correction for Gaussian curve-blocked NCV fits
# (prototype). Sourced after lc-lib.R. Notation as in NOTE-ncv.qmd: A = X'X + S_lambda,
# b = X'y, for curve g: B_g = A - X_g'X_g, t_g = B_g^-1 (b - X_g'y_g) (leave-curve-out fit),
# r_g = y_g - X_g t_g; NCV(rho) = sum_g w_g ||r_g||^2 (mgcv's Gaussian NCV with refund's
# curve-blocked neighbourhoods; weights w_g on curves, w_g = 0 deletes curve g).

# Pieces of a Gaussian NCV fit (lc_pieces() without the REML-only check).
ncv_pieces <- function(fit) {
  if (!identical(fit$method, "NCV")) stop("ncv_pieces needs an NCV fit")
  if (
    !identical(fit$family$family, "gaussian") ||
      !identical(fit$family$link, "identity")
  ) {
    stop("ncv_pieces: Gaussian identity fits only")
  }
  fit$method <- "REML" # lc_pieces() only reads the method for its own check
  P <- lc_pieces(fit, scale_in_H = FALSE) # NCV: no scale parameter in the outer optimisation
  P$H_mgcv <- P$H
  P
}

# Curve-level quantities at smoothing parameters sp with curve weights w (one per curve,
# in the order of P$groups).
ncv_state <- function(P, sp = P$sp, w = rep(1, length(P$groups))) {
  X <- P$X
  idx <- lapply(P$groups, \(g) which(P$cluster == g))
  wi <- rep(0, nrow(X))
  for (j in seq_along(idx)) wi[idx[[j]]] <- w[j]
  Slam <- Reduce(`+`, Map(`*`, P$Sk, sp))
  A <- crossprod(X * sqrt(wi)) + Slam
  b <- drop(crossprod(X, wi * P$y))
  Sk <- Map(`*`, P$Sk, sp)
  cur <- lapply(seq_along(idx), \(j) {
    Xg <- X[idx[[j]], , drop = FALSE]
    yg <- P$y[idx[[j]]]
    Bg <- A - w[j] * crossprod(Xg)
    Binv <- chol2inv(chol(Bg))
    tg <- drop(Binv %*% (b - w[j] * crossprod(Xg, yg)))
    rg <- yg - drop(Xg %*% tg)
    ug <- drop(Binv %*% crossprod(Xg, rg))
    vg <- sapply(Sk, \(S) drop(Binv %*% (S %*% tg))) # p x K
    list(Binv = Binv, t = tg, r = rg, u = ug, v = vg)
  })
  list(idx = idx, w = w, Sk = Sk, cur = cur)
}

ncv_value <- function(st) sum(st$w * vapply(st$cur, \(c) sum(c$r^2), 0))

# Gradient of NCV in rho: d t_g / d rho_k = -B_g^-1 S_k t_g, so
# dNCV/drho_k = 2 sum_g w_g u_g' S_k t_g with u_g = B_g^-1 X_g' r_g.
ncv_grad <- function(st) {
  K <- length(st$Sk)
  vapply(
    seq_len(K),
    \(k) {
      2 * sum(st$w * vapply(st$cur, \(c) sum(c$u * (st$Sk[[k]] %*% c$t)), 0))
    },
    0
  )
}

# Hessian of NCV in rho by central differences of the analytic gradient.
ncv_hessian <- function(P, h = 1e-4) {
  K <- length(P$sp)
  H <- sapply(seq_len(K), \(k) {
    e <- replace(rep(0, K), k, h)
    (ncv_grad(ncv_state(P, P$sp * exp(e))) -
      ncv_grad(ncv_state(P, P$sp * exp(-e)))) /
      (2 * h)
  })
  (H + t(H)) / 2
}

# Curve scores s_h = d grad NCV / d w_h at w = 1 (G x K), analytic (NOTE-ncv.qmd):
# s_hk = 2 q_hk + 2 sum_{g != h} [ e_{h|g}' X_h (m2_gk - m1_gk) - (X_h u_g)'(X_h v_gk) ],
# q_hk = u_h' S_k t_h, e_{h|g} = y_h - X_h t_g, m1_gk = B_g^-1 X_g'X_g v_gk,
# m2_gk = B_g^-1 S_k u_g.
ncv_scores <- function(P, st = ncv_state(P)) {
  X <- P$X
  G <- length(st$idx)
  K <- length(st$Sk)
  cl <- factor(P$cluster, levels = P$groups)
  s <- matrix(0, G, K)
  for (g in seq_len(G)) {
    c <- st$cur[[g]]
    Xg <- X[st$idx[[g]], , drop = FALSE]
    XtXg <- crossprod(Xg)
    e <- P$y - drop(X %*% c$t)
    Xu <- drop(X %*% c$u)
    for (k in seq_len(K)) {
      v <- c$v[, k]
      m1 <- drop(c$Binv %*% (XtXg %*% v))
      m2 <- drop(c$Binv %*% (st$Sk[[k]] %*% c$u))
      contrib <- rowsum(
        e * drop(X %*% (m2 - m1)) - Xu * drop(X %*% v),
        cl,
        reorder = TRUE
      )[, 1]
      contrib[g] <- 0
      s[, k] <- s[, k] + 2 * contrib
      s[g, k] <- s[g, k] + 2 * sum(c$u * (st$Sk[[k]] %*% c$t))
    }
  }
  s
}

# IJ influences on rho-hat (G x K, + orientation, centred) and one-step leave-curve-out
# changes (G x K): Delta_h = -H^-1 [grad NCV_(-h)(rho-hat) - grad NCV(rho-hat)].
ncv_ij_influence <- function(P, H, s = ncv_scores(P), prior = 0) {
  s <- sweep(s, 2, colMeans(s))
  -s %*% pinv_sym(H, prior)
}
ncv_onestep <- function(P, H, prior = 0) {
  g0 <- ncv_grad(ncv_state(P))
  Hinv <- pinv_sym(H, prior)
  G <- length(P$groups)
  K <- length(P$sp)
  out <- matrix(0, G, K) # explicit shape (vapply/t() would give 1 x G for K = 1)
  for (h in seq_len(G)) {
    w <- replace(rep(1, G), h, 0)
    out[h, ] <- -drop(Hinv %*% (ncv_grad(ncv_state(P, w = w)) - g0))
  }
  out
}

# Exact leave-curve-out NCV smoothing parameters (reference): minimise NCV_(-h) over rho by
# Newton steps with the analytic gradient and a finite-difference Hessian, from rho-hat.
# Returns rho with attributes converged, reason, iterations, value, grad (local minimum
# only; no multiple starts).
ncv_refit_rho <- function(P, w, tol = 1e-8, maxit = 50) {
  rho <- log(P$sp)
  K <- length(rho)
  Pw <- P
  for (it in seq_len(maxit)) {
    Pw$sp <- exp(rho)
    g <- ncv_grad(ncv_state(Pw, w = w))
    H <- sapply(seq_len(K), \(k) {
      e <- replace(rep(0, K), k, 1e-4)
      (ncv_grad(ncv_state(Pw, exp(rho + e), w)) -
        ncv_grad(ncv_state(Pw, exp(rho - e), w))) /
        2e-4
    })
    H <- (H + t(H)) / 2
    step <- -drop(pinv_sym(H, 1e-6) %*% g)
    step <- step * min(1, 2 / max(abs(step))) # trust region on log scale
    f0 <- ncv_value(ncv_state(Pw, exp(rho), w))
    if (max(abs(g)) < tol * (1 + f0))
      return(structure(
        rho,
        converged = TRUE,
        reason = "gradient",
        iterations = it,
        value = f0,
        grad = g
      ))
    accepted <- FALSE
    for (half in 0:20) {
      if (ncv_value(ncv_state(Pw, exp(rho + step / 2^half), w)) <= f0) {
        accepted <- TRUE
        break
      }
    }
    if (!accepted)
      return(structure(
        rho,
        converged = FALSE,
        reason = "line search failed",
        iterations = it,
        value = f0,
        grad = g
      ))
    rho <- rho + step / 2^half
    if (max(abs(step / 2^half)) < 1e-7) {
      Pw$sp <- exp(rho)
      g <- ncv_grad(ncv_state(Pw, w = w))
      f <- ncv_value(ncv_state(Pw, w = w))
      return(structure(
        rho,
        converged = max(abs(g)) < 1e-4 * (1 + f),
        reason = "small step",
        iterations = it,
        value = f,
        grad = g
      ))
    }
  }
  structure(
    rho,
    converged = FALSE,
    reason = "maxit",
    iterations = maxit,
    value = NA,
    grad = g
  )
}

# All covariance arms of a Gaussian NCV fit ---------------------------------------------

ncv_covariances <- function(fit, onestep = TRUE) {
  tm <- c()
  tick <- \(label, expr) {
    t0 <- proc.time()[["elapsed"]]
    val <- force(expr)
    tm[[label]] <<- proc.time()[["elapsed"]] - t0
    val
  }
  P <- tick("pieces", ncv_pieces(fit))
  if (any(fit$offset != 0)) stop("offsets are not supported")
  core <- tick(
    "cl2",
    refund:::pffr_influence(fit, sandwich = "cl2", cl2_adjustment = "exact")
  )
  if (!identical(as.character(core$groups), as.character(P$groups))) {
    stop("CL2 influence columns are not in the order of P$groups")
  }
  V_cl2 <- core$correction * tcrossprod(core$K) + core$B2
  H <- tick("hessian", ncv_hessian(P))
  st <- ncv_state(P)
  s <- tick("scores", ncv_scores(P, st))
  D_ij <- ncv_ij_influence(P, H, s)
  D_ij_reg <- ncv_ij_influence(P, H, s, prior = REG_PRIOR)
  if (onestep) {
    os <- tick("onestep", ncv_onestep(P, H))
    os_reg <- ncv_onestep(P, H, prior = REG_PRIOR)
  }
  P_mine <- replace(P, "H", list(H))
  V <- list(
    cl2 = V_cl2,
    naive_mgcv = V_cl2 + naive_correction(P),
    naive = V_cl2 + naive_correction(P_mine),
    ij_lin = add_rho_part(core, rho_part_linear(P, D_ij), "ij"),
    ij_linx = add_rho_part(core, rho_part_linear(P, D_ij), "ij", cross = TRUE),
    ij_reg = add_rho_part(core, rho_part_linear(P, D_ij_reg), "ij"),
    ij_regx = add_rho_part(
      core,
      rho_part_linear(P, D_ij_reg),
      "ij",
      cross = TRUE
    )
  )
  if (onestep) {
    V$os_lin <- add_rho_part(core, rho_part_linear(P, -os), "jk")
    V$os_linx <- add_rho_part(core, rho_part_linear(P, -os), "jk", cross = TRUE)
    V$os_reg <- add_rho_part(core, rho_part_linear(P, -os_reg), "jk")
  }
  list(
    V = V,
    omega = list(
      hat = log(P$sp),
      naive_mgcv = sqrt(pmax(diag(pinv_sym(P$H)), 0)),
      naive = sqrt(pmax(diag(pinv_sym(H)), 0)),
      ij = sqrt(diag(omega_cov(D_ij, "ij"))),
      os = if (onestep) sqrt(diag(omega_cov(os, "jk"))),
      ij_reg = sqrt(diag(omega_cov(D_ij_reg, "ij"))),
      os_reg = if (onestep) sqrt(diag(omega_cov(os_reg, "jk")))
    ),
    grad = ncv_grad(st),
    time = tm
  )
}

# One task: NCV fit, its arms, and the fixed-lambda oracle (sp fixed at the cell's
# geometric-mean NCV values; CL2, Bayesian, exact). CRN against the study's NCV columns.
ncv_task <- function(
  task,
  cell,
  ref_path,
  sp_table = NULL,
  simulate_fun = simulate,
  onestep = TRUE
) {
  stage <- "simulate"
  tryCatch(
    {
      sim <- simulate_fun(cell, task$rep)
      stage <- "fit"
      fit <- fit_model(sim, cell, "NCV")
      stage <- "covariances"
      cov <- collect_warnings(
        if (identical(fit$family$family, "gaussian"))
          ncv_covariances(fit, onestep) else ncvg_covariances(fit)
      )
      stage <- "fixed"
      fix <- if (is.null(sp_table)) NULL else
        tryCatch(
          {
            rows <- sp_table[sp_table$cell == cell$cell, ]
            sp <- setNames(exp(rows$mean_log_sp), rows$penalty)[names(fit$sp)]
            if (anyNA(sp)) stop("fixed sp do not match fit$sp")
            t0 <- proc.time()[["elapsed"]]
            f <- collect_warnings(fit_pffr(sim, cell, "REML", sp = sp))$value
            V <- refund:::pffr_vcov(
              f,
              sandwich = "cl2",
              freq = FALSE,
              cl2_adjustment = "exact"
            )
            list(fit = f, V = V, time = proc.time()[["elapsed"]] - t0)
          },
          error = \(e) NULL
        )
      stage <- "assemble"
      se <- \(L, V) sqrt(pmax(rowSums((L %*% V) * L), 0))
      records <- lapply(
        setNames(nm = model_estimands(cell$model)),
        \(estimand) {
          L <- evaluation_matrix(fit, sim, estimand)
          cols <- c(
            list(
              truth = sim$truth[[estimand]],
              est_NCV = drop(L %*% fit$coefficients)
            ),
            setNames(
              lapply(cov$value$V, \(V) se(L, V)),
              paste0("se_", names(cov$value$V))
            )
          )
          if (!is.null(fix)) {
            Lf <- evaluation_matrix(fix$fit, sim, estimand)
            cols$est_FIX <- drop(Lf %*% fix$fit$coefficients)
            cols$se_fix <- se(Lf, fix$V)
          }
          do.call(cbind, cols)
        }
      )
      crn <- NULL
      if (file.exists(ref_path)) {
        ref <- readRDS(ref_path)
        if (ref$cell != task$cell || ref$rep != task$rep)
          stop("reference record mismatch")
        crn <- do.call(
          rbind,
          lapply(names(records), \(e) {
            map <- c(
              truth = "truth",
              est_NCV = "est_NCV",
              se_cl2 = "se_cl2_NCV"
            )
            do.call(
              rbind,
              lapply(names(map), \(col) {
                x <- records[[e]][, col]
                y <- ref$records[[e]][, map[[col]]]
                data.frame(
                  estimand = e,
                  column = col,
                  max_abs = max(abs(x - y)),
                  max_rel = max(abs(x - y) / pmax(abs(y), 1e-8))
                )
              })
            )
          })
        )
      }
      list(
        status = "success",
        records = records,
        crn = crn,
        omega = cov$value$omega,
        grad = cov$value$grad,
        check = cov$value$check,
        time_parts = c(cov$value$time, fixed = fix$time %||% NA),
        fit_time = fit$restart$time,
        fixed_ok = !is.null(fix),
        converged = isTRUE(fit$converged) &&
          fit$outer.info$conv == "full convergence",
        warnings = c(fit$restart$warnings, cov$warnings)
      )
    },
    error = \(e)
      list(status = "failure", stage = stage, message = conditionMessage(e))
  )
}
