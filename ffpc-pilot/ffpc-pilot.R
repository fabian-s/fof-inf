# Pilot (2026-10-03): inference for the coefficient surface from pffr's ffpc() term (regression on
# the leading FPC scores of X, beta(s,t) = sum_k psi_k(s) b_k(t)) vs the study's ff fits, on the same
# simulated data. Run from the repository root:
#   Rscript ffpc-pilot/ffpc-pilot.R <cell> <rep_from> <rep_to> <out_dir> [pve]
# pve: share of X's variance the retained FPCs explain (ffpc's decomppars; default 0.99, at most 15 PCs).
# Per replicate: the ffpc REML fit's beta on the evaluation grid, model-based and exact-CL2 (Bayesian)
# SEs, the true beta and its projection onto the estimated FPC span (the target ffpc identifies),
# plus the study's stored ff REML/NCV estimates and CL2 SEs for the same cell and replicate.
source("scripts/setup.R")
suppressMessages(library(refund))
args <- commandArgs(trailingOnly = TRUE)
cell <- design_table()
cell <- cell[cell$cell == as.integer(args[1]), ]
reps <- seq(as.integer(args[2]), as.integer(args[3]))
out_dir <- args[4]
pve <- if (length(args) >= 5) as.numeric(args[5]) else 0.99
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
stopifnot(cell$model == "M1")

study_record <- function(cell, rep) {
  # the first study run (cells 1-95) names files by its 95-cell task table
  f <- c(
    if (cell$cell <= 95) {
      sprintf("results/373f6fd65cac/task-%d.rds", (rep - 1) * 95 + cell$cell)
    },
    file.path(
      list.dirs("results", recursive = FALSE),
      sprintf("task-c%03d-r%03d.rds", cell$cell, rep)
    )
  )
  f <- f[file.exists(f)]
  if (!length(f)) return(NULL)
  r <- readRDS(f[1])
  stopifnot(r$cell == cell$cell, r$rep == rep)
  r$records$beta
}

ffpc_task <- function(cell, rep) {
  sim <- simulate(cell, rep)
  ffpc <- refund::ffpc
  xind <- sim$xind
  t0 <- proc.time()[["elapsed"]]
  fit <- refund::pffr(
    Y ~
      ffpc(
        X,
        xind = xind,
        splinepars = BASIS_T,
        decomppars = list(pve = pve, useSymm = TRUE)
      ) +
        z,
    data = sim$data,
    yind = sim$yind,
    bs.yindex = BASIS_T,
    bs.int = BASIS_T,
    family = FAMILIES[[cell$fit_family]]$fit(),
    sandwich = "none",
    method = "REML"
  )
  fit_time <- proc.time()[["elapsed"]] - t0
  # scores and eigenfunctions exactly as ffpc() computes them for this X
  pc <- refund::ffpc(
    sim$data$X,
    xind = xind,
    decomppars = list(pve = pve, useSymm = TRUE)
  )
  xi <- pc$data
  stopifnot(max(abs(pc$PCMat - fit$pffr$ffpc[[1]]$PCMat)) < 1e-10)
  K <- ncol(xi)
  Xc <- sweep(sim$data$X, 2, pc$meanX)
  w <- simpson_weights(xind)
  # psi: functions in span(PCMat) with int X_c psi_k = xi_k (Simpson, as the truth's predictor)
  A <- qr.solve(Xc %*% (w * pc$PCMat), xi)
  psi <- pc$PCMat %*% A
  # Phi: X_c ~ xi Phi'; the identified target is b_k(t) = int Phi_k(s) beta(s,t) ds
  Phi <- t(qr.solve(xi, Xc))
  pts <- eval_points("beta")
  s_idx <- match(round(pts$s, 10), round(xind, 10))
  t_idx <- match(pts$t, EVAL_T)
  L <- matrix(0, nrow(pts), length(fit$coefficients))
  for (k in seq_len(K)) {
    by <- paste0("X.PC", k)
    j <- which(vapply(fit$smooth, \(sm) sm$by == by, TRUE))
    stopifnot(length(j) == 1)
    sm <- fit$smooth[[j]]
    nd <- data.frame(EVAL_T, 1)
    names(nd) <- c(sm$term, by)
    B <- mgcv::PredictMat(sm, nd)
    L[, sm$first.para:sm$last.para] <- psi[s_idx, k] * B[t_idx, ]
  }
  truth <- sim$truth$beta
  beta_fun <- matrix(truth, length(EVAL_S)) # only for checks on the EVAL grid
  truth_full <- make_truth(
    cell$model,
    cell$truth,
    covariate_sd(cell$covariate)
  )$beta(xind, EVAL_T) *
    signal_scale(cell)$scale
  b_target <- t(Phi) %*% (w * truth_full) # K x EVAL_T
  proj <- as.vector((psi %*% b_target)[s_idx[seq_len(length(EVAL_S))], ])
  V_m <- fit$Vp
  V_c <- refund:::pffr_vcov(
    fit,
    sandwich = "cl2",
    cl2_adjustment = "exact",
    freq = FALSE
  )
  se <- \(V) sqrt(rowSums((L %*% V) * L))
  rec <- cbind(
    truth = truth,
    truth_proj = proj,
    est_ffpc = drop(L %*% fit$coefficients),
    se_model_ffpc = se(V_m),
    se_cl2_ffpc = se(V_c)
  )
  st <- study_record(cell, rep)
  if (!is.null(st)) {
    stopifnot(max(abs(st[, "truth"] - truth)) < 1e-10)
    rec <- cbind(
      rec,
      st[, c("est_REML", "se_cl2_REML", "est_NCV", "se_cl2_NCV")]
    )
  }
  list(
    cell = cell$cell,
    rep = rep,
    K = K,
    pve = sum(pc$eigenvalues[seq_len(K)]) / sum(pc$eigenvalues),
    fit_time = fit_time,
    records = rec
  )
}

for (rep in reps) {
  f <- file.path(
    out_dir,
    sprintf("ffpc-p%03d-c%03d-r%03d.rds", round(100 * pve), cell$cell, rep)
  )
  if (file.exists(f)) next
  res <- tryCatch(
    ffpc_task(cell, rep),
    error = \(e) list(cell = cell$cell, rep = rep, error = conditionMessage(e))
  )
  saveRDS(res, f)
  cat(sprintf(
    "cell %d rep %d: %s\n",
    cell$cell,
    rep,
    res$error %||% sprintf("K = %d, %.1f s", res$K, res$fit_time)
  ))
}
