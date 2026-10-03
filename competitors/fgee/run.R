# fastFGEE comparison, both experiments. Run from the repository root with the BLAS pinned:
#   Rscript competitors/fgee/run.R <A|B> <worker> <n_workers> <out_dir> [max_rep]
# Workers claim tasks dynamically (atomic lock directories), rep-major order.
#
# A: fastFGEE's own simulation design (their DGP, their pffr call, cluster = subject).
#    Every interval is scored on the same data; the pffr-based ones on the same fit.
# B: the study's data sets (M1, one curve per cluster), fastFGEE one-step
#    started from the study's REML fit, against REML + model-based and REML + CL2.
# Every estimand is L theta for one evaluation matrix L in pffr's parametrisation;
# fastFGEE returns coefficients in that parametrisation (checked in smoke-a.R: its
# fitted curves and SEs equal ours to machine precision).

.libPaths(c(
  Sys.getenv("RESTART_LIB", path.expand("~/R/restart-lib")),
  path.expand("~/R/fastfgee-lib"),
  .libPaths()
))
args <- commandArgs(TRUE)
experiment <- args[1]
worker <- as.integer(args[2])
n_workers <- as.integer(args[3])
out_dir <- args[4]
max_rep <- as.integer(if (length(args) >= 5) args[5] else 200)
stopifnot(experiment %in% c("A", "B", "P"))

source("scripts/setup.R") # study code (R/)
suppressPackageStartupMessages(library(fastFGEE))
source("competitors/fgee/fgee-dgp.R")
if (experiment == "P") {
  source("plasmode/sources.R")
  source("plasmode/plasmode.R")
}

Z <- qnorm(0.975)
BOOT <- 3000 # fastFGEE default; the package may raise it internally

# Tasks --------------------------------------------------------------------------

cells_A <- function() {
  cells <- expand.grid(
    N = c(25L, 100L, 500L),
    family = c("gaussian", "poisson", "binomial"),
    stringsAsFactors = FALSE
  )
  # n_i = 5: their design (Tables 1-4); n_i = 1 (cells 10-18): one curve per subject,
  # the bridge to our setting (the longitudinal correlation then plays no role).
  cells <- rbind(transform(cells, n_i = 5L), transform(cells, n_i = 1L))
  cells$long_corr <- "exchangeable"
  cells$cell <- seq_len(nrow(cells))
  cells
}

cells_B <- function() {
  cells <- design_table()
  cells[cells$block == "core" & cells$signal == "mid", ]
}

# P: plasmode cells of the public datasets with independent curves, neutral MID truth
# with REML residual curves attached to their own subjects, G = all and G = 40.
# FGEE_P_APPS overrides (e.g. the de-identified cardiology apps "ecg_gls4b,af", which run on
# LRZ only with PLASMODE_CARDIO=deid; their per-task files never leave LRZ).
P_APPS <- strsplit(
  Sys.getenv("FGEE_P_APPS", "running,dti,gait_vc,ecg_roahd"),
  ","
)[[1]]
if (
  experiment == "P" &&
    any(P_APPS %in% PLASMODE_LOCAL_ONLY) &&
    !identical(PLASMODE_CARDIO, "deid")
) {
  stop(
    "cardiology apps only from the de-identified input (PLASMODE_CARDIO=deid)"
  )
}
cells_P <- function() {
  cells <- plasmode_cells()
  cells[
    cells$app %in%
      P_APPS &
      cells$truth == "MID" &
      cells$residual == "REML" &
      cells$flip == "curve",
  ]
}

tasks_for <- function(cells, reps) {
  tasks <- expand.grid(cell = cells$cell, rep = reps)
  tasks[order(tasks$rep, tasks$cell), ]
}

claim <- function(path) dir.create(paste0(path, ".lock"), showWarnings = FALSE)

save_atomic <- function(x, path) {
  tmp <- paste0(path, ".tmp", Sys.getpid())
  saveRDS(x, tmp)
  if (!file.rename(tmp, path)) stop("could not write ", path)
}

# Interval records -----------------------------------------------------------------

# One record per (method, estimand): estimates, SEs and critical values per point.
record <- function(
  method,
  estimand,
  est,
  se,
  crit,
  truth,
  lower = NULL,
  upper = NULL
) {
  data.frame(
    method = method,
    estimand = estimand,
    point = seq_along(est),
    est = est,
    se = se,
    lower = lower %||% (est - crit * se),
    upper = upper %||% (est + crit * se),
    truth = truth
  )
}

# fastFGEE's shipped pointwise intervals for the terms its CI object covers: its own
# centre and bounds (wild-bootstrap critical value x t_{N - edf} adjustment).
their_pointwise <- function(obj, term) {
  ci <- obj$crit$ci %||% obj$model$crit$ci
  if (is.null(ci) || !term %in% names(ci$df)) {
    return(NULL)
  }
  ci$df[[term]]
}

# Experiment A -----------------------------------------------------------------------

run_A <- function(cell, rep, GenerateData) {
  family <- switch(
    cell$family,
    gaussian = gaussian(),
    poisson = poisson(),
    binomial = binomial()
  )
  seed <- 20261001L + 1000L * cell$cell + rep
  sim <- simulate_fgee(
    GenerateData,
    cell$family,
    cell$N,
    cell$n_i,
    cell$long_corr,
    seed
  )
  timing <- list()
  tic <- \() proc.time()[[3]]
  t0 <- tic()
  fit <- fit_initial(sim, family)
  timing$pffr <- tic() - t0
  Ls <- coef_matrices(fit, sim$s)
  ID <- sim$data$ID

  t0 <- tic()
  base <- fgee_baseline(fit, sim, family, joint_ci = "wild", boot = BOOT)
  timing$baseline <- tic() - t0
  w <- base$working0
  meat <- w$dd_sum - w$N * tcrossprod(w$d_bar)
  Ai <- solve(w$W_bar + base$pen.mat / (w$N * penalty_scale(fit)))
  t0 <- tic()
  V <- list(
    model = fit$Vp,
    fgee_base = base$vb,
    fgee_base_fixed = Ai %*% meat %*% Ai / (w$N * (w$N - 1)),
    refund_cr = refund:::pffr_vcov(fit, sandwich = "cluster", cluster = ID),
    # Convention decomposition: frequentist form (no Vp - Ve term), and
    # additionally centred scores, i.e. the fastFGEE meat on refund's bread.
    refund_cr_freq = refund:::pffr_vcov(
      fit,
      sandwich = "cluster",
      cluster = ID,
      freq = TRUE
    ),
    refund_cr_freq_centred = refund:::pffr_vcov(
      fit,
      sandwich = "cluster",
      cluster = ID,
      freq = TRUE,
      center_scores = TRUE
    ),
    refund_cl2 = refund:::pffr_vcov(
      fit,
      sandwich = "cl2",
      cl2_adjustment = "exact",
      cluster = ID
    )
  )
  timing$refund_cov <- tic() - t0

  # The same model fitted as in the study (gam, REML, same bases): separates
  # the fitting pipeline (bam/fREML/discrete) from the covariance question.
  t0 <- tic()
  reml <- tryCatch(fit_initial(sim, family, reml = TRUE), error = \(e) e)
  timing$reml <- tic() - t0
  V_reml <- if (!inherits(reml, "error")) {
    list(
      reml_model = reml$Vc %||% reml$Vp,
      reml_cr = refund:::pffr_vcov(reml, sandwich = "cluster", cluster = ID),
      reml_cl2 = refund:::pffr_vcov(
        reml,
        sandwich = "cl2",
        cl2_adjustment = "exact",
        cluster = ID
      )
    )
  }

  onestep <- list()
  for (spec in list(
    c("independent", "independent"),
    c("ar1", "exchangeable")
  )) {
    name <- paste0("fgee_", spec[1], "_", spec[2])
    t0 <- tic()
    onestep[[name]] <- tryCatch(
      fgee_onestep(fit, sim, family, spec[1], spec[2], boot = BOOT),
      error = \(e) e
    )
    timing[[name]] <- tic() - t0
  }

  terms <- c(beta0 = "fn_intercept", beta1 = "X1", beta2 = "X2")
  recs <- list()
  for (k in seq_along(Ls)) {
    r <- names(Ls)[k]
    L <- Ls[[k]]
    truth <- sim$beta[k, ]
    est <- drop(L %*% fit$coefficients)
    for (m in names(V)) {
      recs[[length(recs) + 1]] <- record(m, r, est, se_of(L, V[[m]]), Z, truth)
    }
    if (!is.null(V_reml)) {
      est_r <- drop(L %*% reml$coefficients)
      for (m in names(V_reml)) {
        recs[[length(recs) + 1]] <- record(
          m,
          r,
          est_r,
          se_of(L, V_reml[[m]]),
          Z,
          truth
        )
      }
    }
    tp <- their_pointwise(base, terms[k])
    if (!is.null(tp)) {
      recs[[length(recs) + 1]] <- record(
        "fgee_base_wild",
        r,
        tp$fit,
        tp$se,
        NA,
        truth,
        tp$lower_pt,
        tp$upper_pt
      )
    }
    for (m in names(onestep)) {
      obj <- onestep[[m]]
      if (inherits(obj, "error")) next
      est_m <- drop(L %*% obj$beta)
      recs[[length(recs) + 1]] <- record(
        m,
        r,
        est_m,
        se_of(L, obj$vb),
        Z,
        truth
      )
      tp <- their_pointwise(obj, terms[k])
      if (!is.null(tp)) {
        recs[[length(recs) + 1]] <- record(
          paste0(m, "_wild"),
          r,
          tp$fit,
          tp$se,
          NA,
          truth,
          tp$lower_pt,
          tp$upper_pt
        )
      }
    }
  }
  list(
    records = do.call(rbind, recs),
    timing = unlist(timing),
    errors = vapply(
      onestep,
      \(o) if (inherits(o, "error")) conditionMessage(o) else "",
      ""
    ),
    sp = fit$sp,
    sig2 = fit$sig2,
    converged = c(
      pffr = isTRUE(fit$converged),
      reml = !inherits(reml, "error") && isTRUE(reml$converged),
      vapply(onestep, \(o) !inherits(o, "error") && isTRUE(o$converged), TRUE)
    ),
    n_boot = vapply(
      c(list(base = base), onestep),
      \(o) if (inherits(o, "error")) NA_real_ else NROW(o$qn %||% NA),
      1
    )
  )
}

# Experiment B -----------------------------------------------------------------------

run_B <- function(cell, rep, simulate_fun = simulate) {
  sim <- simulate_fun(cell, rep)
  set.seed(20261001L + 10000L * cell$cell + rep) # simulate() restores the RNG state
  fam <- FAMILIES[[cell$fit_family]]$fit()
  tic <- \() proc.time()[[3]]
  t0 <- tic()
  fit <- fit_model(sim, cell, "REML")
  timing <- list(reml = tic() - t0)
  t0 <- tic()
  inputs <- method_inputs(fit, "REML", sim, cell) # model-based (Vc) and CL2 SEs
  timing$cl2 <- tic() - t0

  d <- sim$data
  d$ID <- seq_len(nrow(d))
  xind <- sim$xind
  ff <- refund::ff
  ff_basis <- cell_bases(cell)$ff # the cell's own bases (plasmode apps may differ)
  form <- Y ~ ff(X, xind = xind, splinepars = ff_basis) + z
  onestep <- list()
  for (cf in c("independent", "ar1", "fpca")) {
    name <- paste0("fgee_", cf)
    t0 <- tic()
    onestep[[name]] <- tryCatch(
      suppressMessages(fastFGEE::fgee(
        formula = form,
        data = d,
        cluster = "ID",
        family = fam,
        corr_fn = cf,
        corr_long = "independent",
        pffr.mod = fit,
        joint.CI = "wild",
        boot.samps = BOOT,
        verbose.tuning = FALSE
      )),
      error = \(e) e
    )
    timing[[name]] <- tic() - t0
  }

  terms <- c(alpha = "fn_intercept", gamma = "z")
  recs <- list()
  for (estimand in model_estimands(cell$model)) {
    L <- evaluation_matrix(fit, sim, estimand)
    truth <- sim$truth[[estimand]]
    if (estimand == "mean") truth <- fam$linkfun(truth) # scored on the link scale
    x <- inputs$estimands[[estimand]]
    recs[[length(recs) + 1]] <- record(
      "reml_model",
      estimand,
      x$est,
      x$se_model,
      Z,
      truth
    )
    recs[[length(recs) + 1]] <- record(
      "reml_cl2",
      estimand,
      x$est,
      x$se_cl2,
      Z,
      truth
    )
    for (m in names(onestep)) {
      obj <- onestep[[m]]
      if (inherits(obj, "error")) next
      est <- drop(L %*% obj$beta)
      se <- se_of(L, obj$vb)
      recs[[length(recs) + 1]] <- record(m, estimand, est, se, Z, truth)
      # Their calibrated critical value (wild bootstrap x t adjustment) where their CI
      # object covers the term, applied to the same estimate and SE.
      tp <- if (estimand %in% names(terms))
        their_pointwise(obj, terms[[estimand]])
      if (!is.null(tp)) {
        # pointwise_crit already includes the t_{N - edf} adjustment (their
        # lower_pt = fit - pointwise_crit * se); transferred to our evaluation points.
        crit <- tp$pointwise_crit[1]
        recs[[length(recs) + 1]] <- record(
          paste0(m, "_wild"),
          estimand,
          est,
          se,
          crit,
          truth
        )
      }
    }
  }
  list(
    records = do.call(rbind, recs),
    timing = unlist(timing),
    errors = vapply(
      onestep,
      \(o) if (inherits(o, "error")) conditionMessage(o) else "",
      ""
    ),
    lambda = lapply(onestep, \(o) if (!inherits(o, "error")) o$lambda),
    sp = fit$sp,
    converged = c(
      reml = isTRUE(fit$converged),
      vapply(onestep, \(o) !inherits(o, "error") && isTRUE(o$converged), TRUE)
    ),
    warnings = fit$restart$warnings
  )
}

# Main loop ---------------------------------------------------------------------------

cells <- switch(experiment, A = cells_A(), B = cells_B(), P = cells_P())
tasks <- tasks_for(cells, seq_len(max_rep))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
GenerateData <- if (experiment == "A") source_upstream("competitors/fgee")
versions <- list(
  refund = as.character(packageVersion("refund")),
  refund_sha = Sys.getenv("REFUND_SHA"),
  mgcv = as.character(packageVersion("mgcv")),
  mgcv_patched = packageDescription("mgcv")$Patched %||% "",
  fastFGEE = as.character(packageVersion("fastFGEE"))
)
for (i in seq_len(nrow(tasks))) {
  cell <- cells[cells$cell == tasks$cell[i], ]
  rep <- tasks$rep[i]
  path <- file.path(
    out_dir,
    sprintf("%s-c%03d-r%03d.rds", experiment, cell$cell, rep)
  )
  if (file.exists(path) || !claim(path)) next
  t0 <- proc.time()[[3]]
  out <- tryCatch(
    switch(
      experiment,
      A = run_A(cell, rep, GenerateData),
      B = run_B(cell, rep),
      P = run_B(cell, rep, plasmode_simulate)
    ),
    error = \(e) list(error = conditionMessage(e))
  )
  out$cell <- as.list(cell)
  out$rep <- rep
  out$elapsed <- proc.time()[[3]] - t0
  out$versions <- versions
  out$worker <- worker
  save_atomic(out, path)
  cat(sprintf(
    "[%s] worker %d: %s cell %d rep %d in %.0f s%s\n",
    format(Sys.time(), "%H:%M:%S"),
    worker,
    experiment,
    cell$cell,
    rep,
    out$elapsed,
    if (!is.null(out[["error"]])) paste(" ERROR:", out[["error"]]) else ""
  ))
}
