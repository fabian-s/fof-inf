# Fits -------------------------------------------------------------------------

# Evaluate expr, collecting (and muffling) its warnings: list(value, warnings).
collect_warnings <- function(expr) {
  warnings <- character()
  value <- withCallingHandlers(expr, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = warnings)
}

fit_model <- function(sim, cell, method) {
  method <- match.arg(
    method,
    c("REML", "NCV", "NCV_point", "AR1", "NCV_block", "NCV_alt")
  )
  if (method == "AR1" && cell$fit_family != "gaussian") {
    stop("AR1 requires a Gaussian fitted family.")
  }
  start <- proc.time()[["elapsed"]]
  result <- collect_warnings(
    if (method == "AR1") profile_ar1(sim, cell) else fit_pffr(sim, cell, method)
  )
  fit <- result$value
  fit$restart <- list(
    time = proc.time()[["elapsed"]] - start,
    warnings = result$warnings
  )
  fit
}

# `sp`: fixed smoothing parameters (all of them, in fit$sp order), passed to gam.
# NCV_block / NCV_alt: NCV leaving out the curve groups sim$data$block / block_alt
# (refund's NCV takes the fit-time `cluster` as its leave-out groups).
fit_pffr <- function(sim, cell, method, rho = 0, sp = NULL) {
  ff <- refund::ff
  s <- mgcv::s
  xind <- sim$xind
  yind <- sim$yind
  bases <- cell_bases(cell)
  ff_basis <- bases$ff
  formula <- switch(
    cell$model,
    M1 = Y ~ ff(X, xind = xind, splinepars = ff_basis) + z,
    # pffr interprets a length-two m vector as one scalar per tensor margin.
    M2 = Y ~
      s(x, bs = BASIS_X$bs, k = BASIS_X$k, m = list(BASIS_X$m, BASIS_T$m)) + z,
    stop("Unknown model: ", cell$model)
  )
  cluster_ids <- switch(
    method,
    NCV_block = sim$data$block,
    NCV_alt = sim$data$block_alt,
    NULL
  )
  if (method %in% c("NCV_block", "NCV_alt") && is.null(cluster_ids)) {
    stop(method, " needs block ids in sim$data")
  }
  refund::pffr(
    formula,
    data = sim$data,
    yind = yind,
    bs.yindex = bases$t,
    bs.int = bases$int %||% bases$t,
    family = FAMILIES[[cell$fit_family]]$fit(),
    sandwich = "none",
    algorithm = if (method == "AR1") "bam" else "gam",
    method = switch(
      method,
      AR1 = "fREML",
      NCV_point = "NCV",
      NCV_block = "NCV",
      NCV_alt = "NCV",
      method
    ),
    ncv_blocks = if (method == "NCV_point") "point" else "cluster",
    knots = if (cell$model == "M2") list(x = x_knots()) else NULL,
    rho = rho,
    sp = sp,
    cluster = cluster_ids
  )
}

# The cell's fitted bases (cells without a `basis` column, e.g. plasmode: the default).
# Larger bases are defined for M1 only (M2's x-margin stays BASIS_X).
cell_bases <- function(cell) {
  name <- cell$basis %||% "default"
  if (!name %in% names(BASES)) stop("unknown basis: ", name)
  if (cell$model == "M2" && name != "default") {
    stop("larger bases are defined for M1 only")
  }
  BASES[[name]]
}

profile_ar1 <- function(sim, cell) {
  candidates <- lapply(AR_RHO_GRID, \(r) fit_pffr(sim, cell, "AR1", r))
  criteria <- vapply(candidates, \(fit) unname(fit$gcv.ubre), numeric(1))
  best <- which.min(criteria)
  neighbours <- AR_RHO_GRID[c(
    max(1, best - 1),
    min(length(criteria), best + 1)
  )]
  # bam.fit subtracts (n - n_curves) * log(1 / sqrt(1 - rho^2))
  # from fREML: gcv.ubre already includes the whitening log-determinant.
  refined <- optimize(
    \(r) fit_pffr(sim, cell, "AR1", r)$gcv.ubre,
    neighbours
  )
  fit <- candidates[[best]]
  if (refined$objective < criteria[best]) {
    fit <- fit_pffr(sim, cell, "AR1", refined$minimum)
  }
  fit$rho_boundary <- fit$AR1.rho %in% range(AR_RHO_GRID)
  fit
}

# Evaluation matrices ----------------------------------------------------------

estimand_smooth <- function(fit, estimand) {
  matches <- vapply(
    fit$smooth,
    \(sm)
      switch(
        estimand,
        beta = "X.smat" %in% sm$term,
        f = "x" %in% sm$term,
        alpha = length(sm$term) == 1 && sm$by == "NA",
        gamma = sm$by == "z",
        FALSE
      ),
    logical(1)
  )
  if (sum(matches) != 1) stop("Cannot identify smooth for ", estimand)
  fit$smooth[[which(matches)]]
}

estimand_grid <- function(sm, estimand) {
  points <- eval_points(estimand)
  grid <- points[names(points) != "curve"]
  names(grid) <- sm$term
  # A unit by-variable evaluates the coefficient, without ff quadrature weights.
  if (sm$by != "NA") grid[[sm$by]] <- 1
  grid
}

# The mean is evaluated at sim$eval_t (default EVAL_T; plasmode apps on native grids set
# their own), which must be a subset of sim$yind.
evaluation_matrix <- function(fit, sim, estimand) {
  points <- eval_points(estimand, sim$eval_t %||% EVAL_T)
  if (estimand == "mean") {
    times <- vapply(points$t, \(t) which.min(abs(sim$yind - t)), integer(1))
    if (any(abs(sim$yind[times] - points$t) > sqrt(.Machine$double.eps))) {
      stop("The mean's evaluation times must be a subset of sim$yind.")
    }
    if (nrow(sim$test) != N_TEST) stop("Test cohort must have N_TEST rows.")
    lp <- predict(fit, sim$test, type = "lpmatrix", reformat = FALSE)
    return(lp[(points$curve - 1) * length(sim$yind) + times, , drop = FALSE])
  }
  sm <- estimand_smooth(fit, estimand)
  L <- matrix(0, nrow(points), length(fit$coefficients))
  L[, sm$first.para:sm$last.para] <- mgcv::PredictMat(
    sm,
    estimand_grid(sm, estimand)
  )
  if (estimand == "alpha") L[, "(Intercept)" == names(fit$coefficients)] <- 1
  L
}

# Covariances, estimates and diagnostics of one fit -------------------------------

# Per-term EDF and log smoothing parameters (one row per penalty).
fit_terms <- function(fit, method) {
  rows <- lapply(fit$smooth, \(sm) {
    sp <- sm$first.sp:sm$last.sp
    data.frame(
      fit = method,
      term = sm$label,
      edf = sum(fit$edf[sm$first.para:sm$last.para]),
      penalty = names(fit$sp)[sp],
      log_sp = log(fit$sp[sp])
    )
  })
  do.call(rbind, rows)
}

# Everything scoring needs from one fit: estimate and SEs (model-based; CL2 for REML
# and NCV) for every estimand, a one-row fit summary and the per-term table.
method_inputs <- function(fit, method, sim, cell) {
  # Model-based covariance: Vc (with smoothing-parameter uncertainty) where valid, else
  # Vp. For NCV fits mgcv's Vc is not a valid covariance (its SEs are ~0.2 x those of Vp),
  # so NCV fits use Vp. refund's pffr_vcov(sandwich = "none") would
  # return that Vc.
  ncv <- identical(fit$method, "NCV")
  model_covariance <- if (!is.null(fit$Vc) && !ncv) "Vc" else "Vp"
  V_model <- fit[[model_covariance]]
  # CL2 (exact adjustment) for REML and NCV fits (and the block-NCV fits of the plasmode
  # block arm), in refund's default Bayesian form (includes Vp - Ve) and in frequentist
  # form (freq = TRUE, without it). Clustered by curve, and additionally by the curve
  # groups sim$data$block ("b") and block_alt ("a") where the cell attaches them.
  cl2_fits <- c("REML", "NCV", "NCV_block", "NCV_alt")
  V_cl2 <- list()
  cl2_time <- NA_real_
  cov_time <- 0
  if (method %in% cl2_fits) {
    cl2 <- \(freq, cluster) {
      out <- collect_warnings(refund:::pffr_vcov(
        fit,
        sandwich = "cl2",
        cl2_adjustment = "exact",
        freq = freq,
        cluster = cluster
      ))
      fit$restart$warnings <<- c(fit$restart$warnings, out$warnings)
      if (!identical(attr(out$value, "cl2_adjustment"), "exact")) {
        stop("Refund did not use exact CL2 for ", method)
      }
      out$value
    }
    # NULL = the fit-time grouping, i.e. curves, except for the block-NCV fits, whose
    # fit-time grouping is the blocks: the identity grouping forces curves there.
    curve <- if (method %in% c("NCV_block", "NCV_alt")) {
      seq_len(nrow(sim$data))
    }
    groupings <- list(
      cl2 = curve,
      cl2b = sim$data$block,
      cl2a = sim$data$block_alt
    )
    groupings <- groupings[c(
      TRUE,
      !is.null(sim$data$block),
      !is.null(sim$data$block_alt)
    )]
    for (name in names(groupings)) {
      start <- proc.time()[["elapsed"]]
      V_cl2[[name]] <- cl2(FALSE, groupings[[name]])
      if (name == "cl2") cl2_time <- proc.time()[["elapsed"]] - start # one matrix
      V_cl2[[paste0(name, "f")]] <- cl2(TRUE, groupings[[name]])
      cov_time <- cov_time + proc.time()[["elapsed"]] - start
    }
  }
  se <- \(L, V) sqrt(rowSums((L %*% V) * L))
  estimands <- model_estimands(cell$model)
  values <- lapply(estimands, \(estimand) {
    L <- evaluation_matrix(fit, sim, estimand)
    c(
      list(est = drop(L %*% fit$coefficients), se_model = se(L, V_model)),
      if (length(V_cl2)) {
        setNames(lapply(V_cl2, \(V) se(L, V)), paste0("se_", names(V_cl2)))
      }
    )
  })
  list(
    estimands = setNames(values, estimands),
    fit_row = data.frame(
      fit = method,
      time = fit$restart$time,
      cl2_time = cl2_time,
      cov_time = cov_time, # all CL2 matrices of this fit
      converged = isTRUE(fit$converged) &&
        (is.null(fit$outer.info$conv) ||
          fit$outer.info$conv == "full convergence"),
      n_warnings = length(fit$restart$warnings),
      criterion = as.numeric(fit$gcv.ubre),
      rho = fit$AR1.rho %||% NA_real_,
      rho_boundary = isTRUE(fit$rho_boundary),
      model_covariance = model_covariance
    ),
    terms = fit_terms(fit, method),
    warnings = fit$restart$warnings
  )
}

# Concurvity ------------------------------------------------------------------------

# Worst-case concurvity (mgcv, full = TRUE: each smooth against all other terms), one row
# per smooth. The "para" column is dropped: it is only pffr's scalar intercept.
fit_concurvity <- function(fit) {
  # model.matrix() through predict.gam, which predict.pffr would reach with a warning
  class(fit) <- setdiff(class(fit), "pffr")
  worst <- mgcv::concurvity(fit, full = TRUE)["worst", ]
  worst <- worst[names(worst) != "para"]
  data.frame(term = names(worst), concurvity = unname(worst))
}

# Concurvity of a fit as a diagnostic: never fails (an error gives NA) and warns, naming
# cell and term, above `threshold`.
check_concurvity <- function(fit, cell, threshold = CONCURVITY_MAX) {
  out <- tryCatch(
    fit_concurvity(fit),
    error = \(e) data.frame(term = NA_character_, concurvity = NA_real_)
  )
  for (i in which(out$concurvity > threshold)) {
    warning(
      sprintf(
        "cell %s: worst-case concurvity %.2f of %s exceeds %.2f",
        cell$cell,
        out$concurvity[i],
        out$term[i],
        threshold
      ),
      call. = FALSE
    )
  }
  out
}

interval_inputs <- function(fits, sim, cell) {
  setNames(
    lapply(names(fits), \(m) method_inputs(fits[[m]], m, sim, cell)),
    names(fits)
  )
}
