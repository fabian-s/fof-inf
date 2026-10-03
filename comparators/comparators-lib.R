# Comparator study: refund's own remedies for dependent residual curves against REML + CL2.
# Sourced after scripts/setup.R, plasmode/sources.R and plasmode/plasmode.R.
#
# Arms on one simulated or plasmode data set (all Gaussian, model M1):
#   REML    the study's REML fit: model-based (Vc) and CL2 (exact, Bayesian) intervals
#   PCRE    REML fit + pcre() functional random effect for the curves, eigenfunctions
#           from fpca.sc() of the REML residual curves (pve = PCRE_PVE); Vc intervals
#   GLSraw  GLS as in refund's (pre-deprecation) pffr_gls(): hatSigma = cov() of the REML
#           residual curves (its documented example); Vc intervals
#   GLSfpc  the same GLS with hatSigma = the fpca.sc() covariance + white noise
#   BOOT    curve bootstrap of the REML fit as in pffr_coefboot(method = "resample"):
#           curves resampled, model refitted with REML smoothing-parameter reselection;
#           percentile (quantile type 8, as refund) and bootstrap-SE Wald intervals

PCRE_PVE <- 0.95 # share of residual variance the pcre eigenfunctions explain
GLS_COND_CUTOFF <- 500 # refund's cond.cutoff default for hatSigma
BOOT_B <- 199
SEED_BOOT_BASE <- 20261001L # bootstrap seed: SEED_BOOT_BASE + 1000 * cell + rep
STUDY_N_CELLS <- 95 # task table of the study run in the synthetic reference directory
COMP_APPS <- c("ecg_gls4b", "af", "running", "dti", "gait_vc", "ecg_roahd")

# Cells and tasks --------------------------------------------------------------------

# Synthetic: the six Gaussian core cells at the default signal (M1, D = 61).
comparator_syn_cells <- function(cells = design_table()) {
  core <- cells[
    cells$block == "core" & cells$signal == "mid" & cells$family == "gaussian",
  ]
  stopifnot(
    nrow(core) == 6,
    all(core$model == "M1"),
    all(core$D == 61),
    all(core$basis == "default")
  )
  core$source <- "synthetic"
  core$key <- sprintf("s%03d", core$cell)
  core
}

# Plasmode: MID truth, REML residual curves attached to their own subjects, curve flips,
# G = all and 40, on the six datasets with independent curves (the cells fGEE used).
comparator_plas_cells <- function(cells = plasmode_cells()) {
  pc <- cells[
    cells$app %in%
      COMP_APPS &
      cells$truth == "MID" &
      cells$residual == "REML" &
      cells$flip == "curve",
  ]
  stopifnot(nrow(pc) == 12, all(pc$model == "M1"))
  pc$source <- "plasmode"
  pc$key <- sprintf("p%03d", pc$cell)
  pc
}

comparator_cells <- function() {
  keep <- \(x) x[, c("key", "source", "cell", "G", "D", "error", "app")]
  syn <- comparator_syn_cells()
  syn$app <- NA_character_
  dplyr::bind_rows(keep(syn), keep(comparator_plas_cells()))
}

# Rep-major task table; tasks are looked up by id.
comparator_tasks <- function(keys, reps = seq_len(N_REP)) {
  tasks <- tidyr::expand_grid(rep = reps, key = keys)
  tasks$task <- seq_len(nrow(tasks))
  dplyr::relocate(tasks, task)
}

# Full cell specification (study design row) of a key.
cell_spec <- function(key) {
  table <- if (startsWith(key, "s")) {
    comparator_syn_cells()
  } else {
    comparator_plas_cells()
  }
  table[table$key == key, ]
}

# Helpers ------------------------------------------------------------------------------

model_cov <- function(fit) if (!is.null(fit$Vc)) fit$Vc else fit$Vp
se_rows <- function(L, V) sqrt(rowSums((L %*% V) * L))
elapsed <- function() proc.time()[["elapsed"]]

fit_ok <- function(fit) {
  isTRUE(fit$converged) &&
    (is.null(fit$outer.info$conv) || fit$outer.info$conv == "full convergence")
}

# Evaluation matrices of every estimand from the REML fit (columns = its coefficients).
reml_matrices <- function(fit, sim, cell) {
  lapply(setNames(nm = model_estimands(cell$model)), \(e) {
    evaluation_matrix(fit, sim, e)
  })
}

# Evaluation matrices for a fit whose coefficients contain the REML fit's (by name) plus
# extra ones (pcre): REML columns mapped by name, zero for the extra coefficients (a new
# curve's functional random effect is zero; the mean is the population mean).
embed_matrix <- function(L, reml_names, target_names) {
  idx <- match(reml_names, target_names)
  if (anyNA(idx))
    stop("fit lacks REML coefficients: ", reml_names[is.na(idx)][1])
  out <- matrix(0, nrow(L), length(target_names))
  out[, idx] <- L
  out
}

# Residual curves (G x D) of a pffr fit.
residual_curves <- function(fit, G) {
  E <- residuals(fit)
  if (!is.matrix(E) || nrow(E) != G)
    stop("residuals(fit) is not a G x D matrix")
  E
}

# FPCA of the REML residual curves (refund's pcre documentation: fpca.sc of resid()).
residual_fpca <- function(E, yind, pve = PCRE_PVE) {
  # var = TRUE: only then does fpca.sc return the white-noise variance sigma2
  fpc <- refund::fpca.sc(Y = E, argvals = yind, pve = pve, var = TRUE)
  if (!is.numeric(fpc$sigma2) || length(fpc$sigma2) != 1) {
    stop("fpca.sc returned no white-noise variance")
  }
  keep <- fpc$evalues > 0
  list(
    efunctions = fpc$efunctions[, keep, drop = FALSE],
    evalues = fpc$evalues[keep],
    sigma2 = fpc$sigma2,
    npc = sum(keep)
  )
}

# Arms ------------------------------------------------------------------------------

# REML + pcre(): same formula, bases and settings as fit_pffr(), plus the curves'
# functional random effect.
fit_pcre <- function(sim, cell, fpc) {
  ff <- refund::ff
  pcre <- refund::pcre
  xind <- sim$xind
  yind <- sim$yind
  bases <- cell_bases(cell)
  ff_basis <- bases$ff
  efunctions <- fpc$efunctions
  evalues <- fpc$evalues
  data <- sim$data
  data$id <- factor(seq_len(nrow(data)))
  refund::pffr(
    Y ~
      ff(X, xind = xind, splinepars = ff_basis) +
        z +
        pcre(id = id, efunctions = efunctions, evalues = evalues, yind = yind),
    data = data,
    yind = yind,
    bs.yindex = bases$t,
    bs.int = bases$int %||% bases$t,
    family = gaussian(),
    sandwich = "none",
    algorithm = "gam",
    method = "REML"
  )
}

# GLS as refund's pffr_gls() did before its deprecation (refund ead7e391^): the gam
# set-up of the model is decorrelated curve by curve with hatSigma^{-1/2} (refund's own
# compute_sqrt_sigma_inv(), incl. its nearPD repair above cond.cutoff, and
# decorrelate_gam_matrices()), then gam(G = G, method = "REML"). The set-up is rebuilt
# from the REML fit's model frame and must reproduce its model matrix exactly.
fit_gls <- function(fit, Sigma, G, D) {
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
  m
}

# refund's compute_sqrt_sigma_inv() (79a346fb), with one fix: a singular hatSigma (e.g.
# the raw residual covariance when G < D) has a smallest eigenvalue of ~ -1e-15, so refund's
# condition number max/min is negative, its nearPD repair is skipped and sqrt() of the
# negative eigenvalue gives NaN. Here a non-positive smallest eigenvalue counts as an
# infinite condition number, i.e. the documented repair is applied.
sqrt_sigma_inv <- function(hat_sigma, cond_cutoff) {
  e_sigma <- eigen(hat_sigma, symmetric = TRUE)
  lo <- min(e_sigma$values)
  cond <- if (lo <= 0) Inf else max(e_sigma$values) / lo
  if (cond <= cond_cutoff) {
    return(refund:::compute_sqrt_sigma_inv(hat_sigma, cond_cutoff))
  }
  hat_sigma_pd <- Matrix::nearPD(
    hat_sigma,
    keepDiag = TRUE,
    ensureSymmetry = TRUE,
    do2eigen = TRUE,
    posd.tol = 1 / cond_cutoff
  )
  warning(
    "hatSigma condition number ",
    signif(cond, 3),
    " > ",
    cond_cutoff,
    ": nearPD repair"
  )
  e_sigma <- eigen(as.matrix(hat_sigma_pd$mat), symmetric = TRUE)
  e_sigma$vectors %*% diag(1 / sqrt(e_sigma$values)) %*% t(e_sigma$vectors)
}

# mgcv's lpmatrix of a pffr fit (bypassing predict.pffr).
gam_lpmatrix <- function(object) {
  class(object) <- setdiff(class(object), "pffr")
  mgcv::predict.gam(object, type = "lpmatrix")
}

# Curve bootstrap (pffr_coefboot(method = "resample") resampling and refit, via boot::boot
# with ordinary resampling of the data rows = curves), returning every estimand's
# replicate values. The REML evaluation matrices are reused for the refits (the bases
# and constraints do not depend on which curves are drawn); the first refit's own
# matrices are compared with them (`L_check`).
curve_bootstrap <- function(sim, cell, Ls, B, seed) {
  n_out <- sum(vapply(Ls, nrow, 1L))
  check <- NA_real_
  stat <- function(data, idx) {
    sim_b <- sim
    sim_b$data <- data[idx, , drop = FALSE]
    fit_b <- tryCatch(
      collect_warnings(fit_pffr(sim_b, cell, "REML"))$value,
      error = \(e) NULL
    )
    if (is.null(fit_b)) return(rep(NA_real_, n_out))
    if (is.na(check)) {
      own <- reml_matrices(fit_b, sim, cell)
      check <<- max(unlist(Map(\(a, b) max(abs(a - b)), own, Ls)))
    }
    unlist(lapply(Ls, \(L) drop(L %*% fit_b$coefficients)))
  }
  start <- elapsed()
  out <- withr::with_seed(seed, boot::boot(sim$data, stat, R = B))
  time <- elapsed() - start
  draws <- out$t
  failed <- rowSums(is.na(draws)) == ncol(draws)
  draws <- draws[!failed, , drop = FALSE]
  q <- \(p) apply(draws, 2, stats::quantile, p, type = 8, names = FALSE)
  alpha <- 1 - LEVEL
  split_by <- rep(names(Ls), vapply(Ls, nrow, 1L))
  pieces <- list(
    lo = q(alpha / 2),
    hi = q(1 - alpha / 2),
    se = apply(draws, 2, sd),
    mean = colMeans(draws)
  )
  list(
    values = lapply(setNames(nm = names(Ls)), \(e) {
      lapply(pieces, \(x) x[split_by == e])
    }),
    row = data.frame(
      arm = "BOOT",
      time = time,
      converged = TRUE,
      n_warnings = 0L,
      info = sprintf("B=%d;failed=%d;L_check=%.1e", B, sum(failed), check)
    ),
    n_failed = sum(failed)
  )
}

# One task ------------------------------------------------------------------------------

arm_row <- function(arm, fit, time, warnings, info = "") {
  data.frame(
    arm = arm,
    time = time,
    converged = fit_ok(fit),
    n_warnings = length(warnings),
    info = info
  )
}

# Run one comparator arm; an error becomes a failure row and NULL (never drops the task).
try_arm <- function(arm, expr, failures) {
  tryCatch(expr, error = \(e) {
    failures$rows[[length(failures$rows) + 1]] <- data.frame(
      stage = arm,
      message = conditionMessage(e)
    )
    NULL
  })
}

run_comparators <- function(sim, cell, rep, boot = TRUE, B = BOOT_B) {
  failures <- new.env()
  failures$rows <- list()
  # REML: estimates, model-based and CL2 SEs exactly as the study computes them
  fit <- fit_model(sim, cell, "REML")
  reml <- method_inputs(fit, "REML", sim, cell)
  Ls <- reml_matrices(fit, sim, cell)
  rows <- list(arm_row(
    "REML",
    fit,
    fit$restart$time,
    fit$restart$warnings,
    sprintf("cl2_time=%.2f", reml$fit_row$cl2_time)
  ))
  cols <- lapply(reml$estimands, \(x) {
    x[c("est", "se_model", "se_cl2", "se_cl2f")] |>
      setNames(c("est_REML", "se_model_REML", "se_cl2_REML", "se_cl2f_REML"))
  })
  add <- \(label, values) {
    for (e in names(cols)) {
      cols[[e]][paste0(names(values[[e]]), "_", label)] <<- values[[e]]
    }
  }
  E <- residual_curves(fit, nrow(sim$data))

  # FPCA of the residual curves (shared by PCRE and GLSfpc)
  start <- elapsed()
  fpc <- try_arm("FPCA", residual_fpca(E, sim$yind), failures)
  fpca_time <- elapsed() - start

  # PCRE
  if (!is.null(fpc)) {
    pcre <- try_arm(
      "PCRE",
      {
        start <- elapsed()
        res <- collect_warnings(fit_pcre(sim, cell, fpc))
        time <- elapsed() - start
        m <- res$value
        V <- model_cov(m)
        values <- lapply(Ls, \(L) {
          Lp <- embed_matrix(L, names(fit$coefficients), names(m$coefficients))
          list(est = drop(Lp %*% m$coefficients), se = se_rows(Lp, V))
        })
        list(
          values = values,
          row = arm_row(
            "PCRE",
            m,
            time + fpca_time,
            res$warnings,
            sprintf(
              "npc=%d;fpca_time=%.2f;cov=%s",
              fpc$npc,
              fpca_time,
              if (is.null(m$Vc)) "Vp" else "Vc"
            )
          )
        )
      },
      failures
    )
    if (!is.null(pcre)) {
      add("PCRE", pcre$values)
      rows[[length(rows) + 1]] <- pcre$row
    }
  }

  # GLS with the two covariance estimates
  sigmas <- list(GLSraw = stats::cov(E))
  if (!is.null(fpc)) {
    sigmas$GLSfpc <- fpc$efunctions %*%
      (fpc$evalues * t(fpc$efunctions)) +
      fpc$sigma2 * diag(ncol(E))
  }
  for (arm in names(sigmas)) {
    gls <- try_arm(
      arm,
      {
        start <- elapsed()
        res <- collect_warnings(fit_gls(
          fit,
          sigmas[[arm]],
          nrow(E),
          ncol(E)
        ))
        m <- res$value
        V <- model_cov(m)
        values <- lapply(Ls, \(L) {
          list(est = drop(L %*% m$coefficients), se = se_rows(L, V))
        })
        repaired <- any(grepl("condition number", res$warnings)) # refund's or ours
        list(
          values = values,
          row = arm_row(
            arm,
            m,
            elapsed() - start,
            res$warnings,
            sprintf("cond=%.3g;nearPD=%s", m$gls_cond, repaired)
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

  # Curve bootstrap
  if (boot) {
    bs <- try_arm(
      "BOOT",
      curve_bootstrap(
        sim,
        cell,
        Ls,
        B,
        SEED_BOOT_BASE + 1000L * cell$cell + rep
      ),
      failures
    )
    if (!is.null(bs)) {
      add("BOOT", bs$values)
      rows[[length(rows) + 1]] <- bs$row
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

# CRN check -----------------------------------------------------------------------------

CRN_COLUMNS <- c(
  "truth",
  "est_REML",
  "se_model_REML",
  "se_cl2_REML",
  "se_cl2f_REML"
)

reference_path <- function(task, cell, ref_syn, ref_plas) {
  if (cell$source == "synthetic") {
    return(file.path(
      ref_syn,
      sprintf("task-%d.rds", (task$rep - 1) * STUDY_N_CELLS + cell$cell)
    ))
  }
  file.path(
    ref_plas,
    cell$app,
    sprintf("task-c%03d-r%03d.rds", cell$cell, task$rep)
  )
}

# Max absolute and relative differences between our REML columns and the study record.
crn_check <- function(records, ref_path, cell, rep) {
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
      do.call(
        rbind,
        lapply(CRN_COLUMNS, \(col) {
          x <- records[[e]][, col]
          y <- ref$records[[e]][, col]
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
