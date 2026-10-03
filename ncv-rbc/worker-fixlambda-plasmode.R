# Oracle fixed-lambda check for NCV-centred intervals on the plasmode cells:
# refit with every smoothing parameter fixed at the cell's median NCV lambda-hat (from the study's
# plasmode task records, replicates 1-100) and compute CL2. If NCV + CL2's short SEs come from
# lambda-hat variability, these intervals recover; if CL2 itself is biased at NCV's EDF, they do not.
#   Rscript ncv-rbc/worker-fixlambda-plasmode.R <out_dir> <records_dir> <app,app,...> <n_rep>
source("scripts/setup.R")
source("plasmode/sources.R")
source("plasmode/plasmode.R")
args <- commandArgs(trailingOnly = TRUE)
out <- args[1]
records <- args[2]
apps <- strsplit(args[3], ",")[[1]]
n_rep <- as.integer(args[4])
if (any(apps %in% PLASMODE_DEID_APPS) && !identical(PLASMODE_CARDIO, "deid"))
  stop("cardiology apps only de-identified")
cells <- plasmode_cells()
cells <- cells[
  cells$app %in%
    apps &
    cells$flip == "curve" &
    cells$G > 40 &
    cells$error == "attached",
]
dir.create(file.path(out, "locks"), recursive = TRUE, showWarnings = FALSE)

# median NCV log sp per penalty, in fit$sp order (rows of the stored terms table follow it)
cell_sp <- function(cell) {
  f <- file.path(out, sprintf("sp-c%03d.rds", cell$cell))
  if (file.exists(f)) return(readRDS(f))
  files <- file.path(
    records,
    cell$app,
    sprintf("task-c%03d-r%03d.rds", cell$cell, seq_len(n_rep))
  )
  files <- files[file.exists(files)]
  terms <- lapply(files, \(p) {
    r <- readRDS(p)
    stopifnot(r$cell == cell$cell)
    r$terms[r$terms$fit == "NCV", ]
  })
  pen <- terms[[1]]$penalty
  ls <- sapply(terms, \(t) {
    stopifnot(identical(t$penalty, pen))
    t$log_sp
  })
  sp <- setNames(exp(apply(matrix(ls, nrow = length(pen)), 1, median)), pen)
  tmp <- tempfile(tmpdir = out)
  saveRDS(sp, tmp)
  file.rename(tmp, f)
  sp
}

task <- function(cell, rep, sp) {
  sim <- plasmode_simulate(cell, rep)
  start <- proc.time()[["elapsed"]]
  fit <- fit_pffr(sim, cell, "REML", sp = unname(sp))
  if (!isTRUE(all.equal(unname(fit$full.sp), unname(sp), tolerance = 1e-8)))
    stop("sp not kept fixed")
  V <- list(
    ncvfix_cl2 = refund:::pffr_vcov(
      fit,
      sandwich = "cl2",
      cl2_adjustment = "exact",
      freq = FALSE
    ),
    ncvfix_cl2f = refund:::pffr_vcov(
      fit,
      sandwich = "cl2",
      cl2_adjustment = "exact",
      freq = TRUE
    )
  )
  linkinv <- FAMILIES[[cell$fit_family]]$linkinv
  z <- qnorm(1 - (1 - LEVEL) / 2)
  rows <- list()
  for (e in model_estimands(cell$model)) {
    L <- evaluation_matrix(fit, sim, e)
    truth <- as.vector(sim$truth[[e]])
    tf <- if (e == "mean") linkinv else identity
    est <- drop(L %*% fit$coefficients)
    for (nm in names(V)) {
      se <- sqrt(rowSums((L %*% V[[nm]]) * L))
      rows[[length(rows) + 1]] <- data.frame(
        cell = cell$cell,
        rep = rep,
        estimand = e,
        arm = nm,
        point = seq_along(truth),
        truth = truth,
        est = tf(est),
        lower = tf(est - z * se),
        upper = tf(est + z * se)
      )
    }
  }
  list(points = do.call(rbind, rows), time = proc.time()[["elapsed"]] - start)
}

tasks <- expand.grid(cell = cells$cell, rep = seq_len(n_rep))
tasks <- tasks[order(tasks$rep, tasks$cell), ]
for (i in seq_len(nrow(tasks))) {
  cell <- cells[cells$cell == tasks$cell[i], ]
  rep <- tasks$rep[i]
  f <- file.path(out, sprintf("c%03d-r%03d.rds", cell$cell, rep))
  if (file.exists(f)) next
  if (!dir.create(file.path(out, "locks", basename(f)), showWarnings = FALSE))
    next
  res <- tryCatch(
    task(cell, rep, cell_sp(cell)),
    error = \(e) list(error = conditionMessage(e))
  )
  tmp <- tempfile(tmpdir = out)
  saveRDS(res, tmp)
  file.rename(tmp, f)
  cat(sprintf(
    "cell %d rep %d %s\n",
    cell$cell,
    rep,
    if (is.null(res$error)) sprintf("%.1fs", res$time) else res$error
  ))
}
