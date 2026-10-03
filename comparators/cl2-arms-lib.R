# CL2 on the comparator fits: GLS (whitened) + CL2 and pcre + CL2.
# Sourced after comparators-lib.R.
#
# GLS + CL2: the GLS fit of fit_gls() is an ordinary penalized fit of the whitened data
# (each curve's rows of y and X premultiplied by hatSigma^{-1/2}); the transform acts within
# curves, so the curves stay the clusters. CL2 (exact adjustment, Bayesian form B2 = Vp - Ve)
# is computed with refund's fixed-fit core on the whitened design, response and bread,
# exactly as refund's build_cl2_working_standard() builds them for a Gaussian identity fit
# (Xw = X / sigma, z = (y - X b) / sigma). The coefficients are in the REML fit's
# parametrisation, so the REML evaluation matrices map them to the estimand grids.
# pcre + CL2: refund:::pffr_vcov(sandwich = "cl2", exact, Bayesian) on the pcre fit,
# clustered by curve; evaluation matrices embedded as for the pcre arm (curve effect 0).

# GLS set-up and fit, as fit_gls() (comparators-lib.R), also returning the whitened set-up.
fit_gls_setup <- function(fit, Sigma, G, D) {
  G0 <- mgcv::gam(fit$formula, data = fit$model, fit = FALSE, method = "REML")
  lp <- gam_lpmatrix(fit)
  if (!identical(dim(G0$X), dim(lp)) || max(abs(G0$X - lp)) > 1e-10) {
    stop("rebuilt gam set-up does not reproduce the REML model matrix")
  }
  ev <- eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values
  W <- sqrt_sigma_inv(Sigma, GLS_COND_CUTOFF)
  G1 <- refund:::decorrelate_gam_matrices(G0, W, G, D)
  m <- mgcv::gam(G = G1, method = "REML")
  if (!identical(names(m$coefficients), names(fit$coefficients))) {
    stop("GLS coefficients are not in the REML fit's parametrisation")
  }
  m$gls_cond <- if (min(ev) <= 0) Inf else max(ev) / min(ev)
  list(fit = m, X = G1$X, y = G1$y)
}

# Exact CL2 of a Gaussian identity fit from its (working) design X, response y, Bayesian
# (freq = FALSE: + Vp - Ve) or frequentist, clustered by `cluster_id` (one per row).
cl2_from_design <- function(
  X,
  y,
  beta,
  Vp,
  Ve,
  sig2,
  cluster_id,
  freq = FALSE
) {
  s <- 1 / sqrt(sig2)
  Xw <- X * s
  z <- drop(y - X %*% beta) * s
  core <- refund:::pffr_influence_core(Xw, Vp, cluster_id, z, "exact")
  core$B2 <- (Vp + t(Vp)) / 2 - Ve
  V <- refund:::pffr_influence_vcov(core, freq = freq)
  attr(V, "diag") <- core$diagnostics
  V
}

# Leverage monitors of a CL2 covariance (whitened / pcre hat): extreme eigenvalues of
# the curve blocks H_gg, extreme h_ii, clusters with a floored residual-block eigenvalue.
leverage_info <- function(V) {
  d <- attr(V, "diag") %||% NULL
  if (is.null(d)) {
    return(sprintf(
      "hat_eig=[%.3g,%.3g];h_ii=[%.3g,%.3g];n_adjusted=%d",
      attr(V, "min_hat_eig"),
      attr(V, "max_leverage"),
      attr(V, "min_obs_leverage"),
      attr(V, "max_obs_leverage"),
      as.integer(attr(V, "n_adjusted") %||% NA)
    ))
  }
  sprintf(
    "hat_eig=[%.3g,%.3g];h_ii=[%.3g,%.3g];n_adjusted=%d",
    min(d$min_hat_eig),
    max(d$max_leverage),
    min(d$min_obs_leverage),
    max(d$max_obs_leverage),
    sum(d$n_floored > 0)
  )
}

# One task: REML (study path), GLS raw / fpc with Vc and CL2, pcre with Vc and CL2
# (only when `pcre` is TRUE). Arm errors become failure rows.
run_cl2_arms <- function(sim, cell, pcre = TRUE) {
  failures <- new.env()
  failures$rows <- list()
  fit <- fit_model(sim, cell, "REML")
  reml <- method_inputs(fit, "REML", sim, cell)
  Ls <- reml_matrices(fit, sim, cell)
  G <- nrow(sim$data)
  cols <- lapply(reml$estimands, \(x) {
    x[c("est", "se_model", "se_cl2")] |>
      setNames(c("est_REML", "se_model_REML", "se_cl2_REML"))
  })
  add <- \(label, values) {
    for (e in names(cols)) {
      cols[[e]][paste0(names(values[[e]]), "_", label)] <<- values[[e]]
    }
  }
  rows <- list(arm_row("REML", fit, fit$restart$time, fit$restart$warnings))
  E <- residual_curves(fit, G)
  D <- ncol(E)
  cluster_id <- rep(seq_len(G), each = D)

  # Check of cl2_from_design() on the REML fit against refund's pffr_vcov (study path)
  check <- try_arm(
    "CL2check",
    {
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
        cluster_id
      )
      max(abs(V_own - V_ref)) / max(abs(V_ref))
    },
    failures
  )
  rows[[1]]$info <- sprintf("cl2_check=%.2e", check %||% NA_real_)

  start <- elapsed()
  fpc <- try_arm("FPCA", residual_fpca(E, sim$yind), failures)
  fpca_time <- elapsed() - start

  sigmas <- list(GLSraw = stats::cov(E))
  if (!is.null(fpc)) {
    sigmas$GLSfpc <- fpc$efunctions %*%
      (fpc$evalues * t(fpc$efunctions)) +
      fpc$sigma2 * diag(D)
  }
  for (arm in names(sigmas)) {
    gls <- try_arm(
      arm,
      {
        start <- elapsed()
        res <- collect_warnings(fit_gls_setup(fit, sigmas[[arm]], G, D))
        m <- res$value$fit
        fit_time <- elapsed() - start
        fitted_check <- max(abs(
          drop(res$value$X %*% m$coefficients) - m$fitted.values
        ))
        if (fitted_check > 1e-8 * max(abs(m$fitted.values))) {
          stop("whitened design does not reproduce the GLS fitted values")
        }
        start <- elapsed()
        V_cl2 <- cl2_from_design(
          res$value$X,
          res$value$y,
          m$coefficients,
          m$Vp,
          m$Ve,
          m$sig2,
          cluster_id
        )
        cl2_time <- elapsed() - start
        V <- model_cov(m)
        values <- lapply(Ls, \(L) {
          list(
            est = drop(L %*% m$coefficients),
            se = se_rows(L, V),
            se_cl2 = se_rows(L, V_cl2)
          )
        })
        repaired <- any(grepl("condition number", res$warnings))
        list(
          values = values,
          row = arm_row(
            arm,
            m,
            fit_time,
            res$warnings,
            sprintf(
              "cond=%.3g;nearPD=%s;sig2=%.4g;cl2_time=%.2f;%s",
              m$gls_cond,
              repaired,
              m$sig2,
              cl2_time,
              leverage_info(V_cl2)
            )
          )
        )
      },
      failures
    )
    if (!is.null(gls)) {
      add(arm, gls$values)
      rows[[length(rows) + 1]] <- gls$row
    }
  }

  if (pcre && !is.null(fpc)) {
    pc <- try_arm(
      "PCRE",
      {
        start <- elapsed()
        res <- collect_warnings(fit_pcre(sim, cell, fpc))
        fit_time <- elapsed() - start + fpca_time
        m <- res$value
        start <- elapsed()
        cl2 <- collect_warnings(refund:::pffr_vcov(
          m,
          sandwich = "cl2",
          cl2_adjustment = "exact",
          freq = FALSE
        ))
        cl2_time <- elapsed() - start
        V_cl2 <- cl2$value
        if (!identical(attr(V_cl2, "cl2_adjustment"), "exact")) {
          stop("refund did not use exact CL2 for pcre")
        }
        V <- model_cov(m)
        values <- lapply(Ls, \(L) {
          Lp <- embed_matrix(L, names(fit$coefficients), names(m$coefficients))
          list(
            est = drop(Lp %*% m$coefficients),
            se = se_rows(Lp, V),
            se_cl2 = se_rows(Lp, V_cl2)
          )
        })
        list(
          values = values,
          row = arm_row(
            "PCRE",
            m,
            fit_time,
            c(res$warnings, cl2$warnings),
            sprintf(
              "npc=%d;p=%d;cl2_time=%.2f;%s",
              fpc$npc,
              length(m$coefficients),
              cl2_time,
              leverage_info(V_cl2)
            )
          )
        )
      },
      failures
    )
    if (!is.null(pc)) {
      add("PCRE", pc$values)
      rows[[length(rows) + 1]] <- pc$row
    }
  }

  records <- lapply(setNames(nm = names(cols)), \(e) {
    do.call(cbind, c(list(truth = as.vector(sim$truth[[e]])), cols[[e]]))
  })
  list(
    records = records,
    fits = do.call(rbind, rows),
    failures = do.call(rbind, failures$rows)
  )
}

# CRN check against a stored record (study or comparator task file): max differences of
# the shared columns.
crn_against <- function(records, ref_path, cell, rep, columns, source_label) {
  if (!file.exists(ref_path)) return(NULL)
  ref <- readRDS(ref_path)
  if (ref$cell != cell$cell || ref$rep != rep) {
    stop(
      "reference record ",
      ref_path,
      " is not cell/rep ",
      cell$cell,
      "/",
      rep
    )
  }
  do.call(
    rbind,
    lapply(names(records), \(e) {
      shared <- intersect(
        intersect(columns, colnames(records[[e]])),
        colnames(ref$records[[e]])
      )
      do.call(
        rbind,
        lapply(shared, \(col) {
          x <- records[[e]][, col]
          y <- ref$records[[e]][, col]
          data.frame(
            reference = source_label,
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
