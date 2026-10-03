# ffpc() worker (2026-10-03): inference for the coefficient surface from pffr's ffpc() term
# (regression on the leading FPC scores of X, beta(s,t) = sum_k psi_k(s) b_k(t)) on the
# study's simulated data (common random numbers), next to the study's stored ff fits.
# Run from the repository root (locally or from the frozen LRZ code copy):
#   Rscript ffpc-pilot/run-ffpc.R <out_dir> <cells, R expression> <rules> <methods> <rep_from> <rep_to>
#   rules:   comma-separated; p<100*pve> (ffpc's pve rule, at most 15 PCs, e.g. p095) or
#            k<K> (fixed number of PCs, e.g. k3)
#   methods: comma-separated smoothing-parameter selection for the b_k: REML, NCV
#            (NCV leaves out whole curves, as the study's ff NCV fits)
# Walks the rep-major task table and claims tasks through atomic lock directories, so any
# number of workers can run or be restarted. Per task: the ffpc fit's beta on the study's
# evaluation grid, model-based (Vp) and exact Bayesian CL2 SEs, the true beta and its
# projection onto the estimated FPC span, plus the study's stored ff REML/NCV estimates and
# CL2 SEs for the same cell and replicate (searched in the colon-separated directories of
# FFPC_REF_DIRS, default results/, and their subdirectories).
# ffpc-pilot.R is the single-cell pilot version of the same computation (identical records).
source("scripts/setup.R")
suppressMessages(library(refund))

parse_rule <- function(rule) {
  stopifnot(grepl("^[pk][0-9]+$", rule))
  value <- as.numeric(sub("^[pk]", "", rule))
  if (startsWith(rule, "p")) list(pve = value / 100, useSymm = TRUE) else
    list(npc = as.integer(value), useSymm = TRUE)
}

# FFPC_REF_DIRS (default results/) and their immediate subdirectories, searched once.
ref_dirs <- local({
  dirs <- NULL
  function() {
    if (is.null(dirs)) {
      env <- Sys.getenv("FFPC_REF_DIRS", "results")
      top <- strsplit(env, ":")[[1]]
      dirs <<- unique(c(top, unlist(lapply(top, list.dirs, recursive = FALSE))))
    }
    dirs
  }
})

# The study's ff REML/NCV columns for this cell and replicate (NULL if not found). The first
# study run (cells 1-95) names files by its 95-cell task table; other tables reuse those
# names for other cells, so cell and rep are checked rather than assumed.
study_record <- function(cell, rep) {
  name <- sprintf("task-c%03d-r%03d.rds", cell$cell, rep)
  if (cell$cell <= 95) {
    name <- c(name, sprintf("task-%d.rds", (rep - 1) * 95 + cell$cell))
  }
  for (path in Filter(
    file.exists,
    as.vector(outer(ref_dirs(), name, file.path))
  )) {
    r <- readRDS(path)
    if (
      isTRUE(r$cell == cell$cell) &&
        isTRUE(r$rep == rep) &&
        !is.null(r$records$beta)
    ) {
      return(r$records$beta)
    }
  }
  NULL
}

ffpc_task <- function(cell, rep, decomp, method) {
  sim <- simulate(cell, rep)
  ffpc <- refund::ffpc
  xind <- sim$xind
  t0 <- proc.time()[["elapsed"]]
  fit <- refund::pffr(
    Y ~ ffpc(X, xind = xind, splinepars = BASIS_T, decomppars = decomp) + z,
    data = sim$data,
    yind = sim$yind,
    bs.yindex = BASIS_T,
    bs.int = BASIS_T,
    family = FAMILIES[[cell$fit_family]]$fit(),
    sandwich = "none",
    method = method,
    ncv_blocks = "cluster"
  )
  fit_time <- proc.time()[["elapsed"]] - t0
  # scores and eigenfunctions exactly as ffpc() computes them for this X
  pc <- refund::ffpc(sim$data$X, xind = xind, decomppars = decomp)
  xi <- pc$data
  stopifnot(max(abs(pc$PCMat - fit$pffr$ffpc[[1]]$PCMat)) < 1e-10)
  K <- ncol(xi)
  Xc <- sweep(sim$data$X, 2, pc$meanX)
  w <- simpson_weights(xind)
  # psi: functions in span(PCMat) with int X_c psi_k = xi_k (Simpson, as the truth's
  # predictor). fpca.sc's scores are shrunken BLUPs, so this holds in least squares only.
  A <- qr.solve(Xc %*% (w * pc$PCMat), xi)
  psi <- pc$PCMat %*% A
  psi_resid <- norm(Xc %*% (w * psi) - xi, "F") / norm(xi, "F")
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
  truth_full <- make_truth(
    cell$model,
    cell$truth,
    covariate_sd(cell$covariate)
  )$beta(xind, EVAL_T) *
    signal_scale(cell)$scale
  b_target <- t(Phi) %*% (w * truth_full) # K x EVAL_T
  proj <- as.vector((psi %*% b_target)[s_idx[seq_len(length(EVAL_S))], ])
  t1 <- proc.time()[["elapsed"]]
  V_c <- refund:::pffr_vcov(
    fit,
    sandwich = "cl2",
    cl2_adjustment = "exact",
    freq = FALSE
  )
  cl2_time <- proc.time()[["elapsed"]] - t1
  se <- \(V) sqrt(rowSums((L %*% V) * L))
  rec <- cbind(
    truth = truth,
    truth_proj = proj,
    est_ffpc = drop(L %*% fit$coefficients),
    se_model_ffpc = se(fit$Vp),
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
    K = K,
    pve = sum(pc$eigenvalues[seq_len(K)]) / sum(pc$eigenvalues),
    psi_resid = psi_resid,
    edf_b = vapply(
      fit$smooth[grepl("PC[0-9]+$", vapply(fit$smooth, `[[`, "", "by"))],
      \(sm) sum(fit$edf[sm$first.para:sm$last.para]),
      numeric(1)
    ),
    fit_time = fit_time,
    cl2_time = cl2_time,
    records = rec
  )
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 6) {
  stop("usage: run-ffpc.R out_dir cells rules methods rep_from rep_to")
}
out_dir <- args[1]
cells <- design_table()
cells <- cells[cells$cell %in% eval(parse(text = paste0("c(", args[2], ")"))), ]
if (nrow(cells) == 0) stop("no cells selected")
stopifnot(all(cells$model == "M1"))
rules <- strsplit(args[3], ",")[[1]]
methods <- strsplit(args[4], ",")[[1]]
stopifnot(all(methods %in% c("REML", "NCV")))
invisible(lapply(rules, parse_rule))
reps <- seq(as.integer(args[5]), as.integer(args[6]))
software <- software_versions()
if (nzchar(Sys.getenv("SLURM_JOB_ID"))) check_software(software)
dir.create(file.path(out_dir, "locks"), recursive = TRUE, showWarnings = FALSE)
tasks <- expand.grid(
  cell = cells$cell,
  rule = rules,
  method = methods,
  rep = reps,
  stringsAsFactors = FALSE
)
tasks <- tasks[order(tasks$rep, tasks$cell, tasks$rule, tasks$method), ]
for (i in seq_len(nrow(tasks))) {
  tk <- tasks[i, ]
  # file names: ffpc-<rule>-c<cell>-r<rep>.rds for REML (as the pilot), -ncv- for NCV
  path <- file.path(
    out_dir,
    sprintf(
      "ffpc-%s%s-c%03d-r%03d.rds",
      tk$rule,
      if (tk$method == "NCV") "-ncv" else "",
      tk$cell,
      tk$rep
    )
  )
  if (file.exists(path)) next
  lock <- file.path(out_dir, "locks", sub("\\.rds$", ".lock", basename(path)))
  if (!dir.create(lock, showWarnings = FALSE)) next
  start <- proc.time()[["elapsed"]]
  res <- tryCatch(
    ffpc_task(
      cells[cells$cell == tk$cell, ],
      tk$rep,
      parse_rule(tk$rule),
      tk$method
    ),
    error = \(e) list(error = conditionMessage(e))
  )
  res <- c(
    res,
    list(
      cell = tk$cell,
      rep = tk$rep,
      rule = tk$rule,
      method = tk$method,
      software = software,
      time = proc.time()[["elapsed"]] - start
    )
  )
  atomic_save(res, path)
  cat(sprintf(
    "cell %d rep %d %s %s: %s\n",
    tk$cell,
    tk$rep,
    tk$rule,
    tk$method,
    res$error %||% sprintf("K = %d, %.1f s", res$K, res$time)
  ))
}
