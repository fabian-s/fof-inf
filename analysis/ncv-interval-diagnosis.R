# Why do intervals centred at the NCV fit miss? Per cell, from the stored task records:
# pointwise bias and SD over replicates of the REML and NCV estimates, CL2 SE / SD,
# bias-eliminated coverage, how well delta = NCV - REML tracks NCV's bias, and candidate
# NCV-centred intervals that need no refit. Aggregates only. Run from the repository root:
#   Rscript analysis/ncv-interval-diagnosis.R synthetic <out_prefix> [cell ...]
#   Rscript analysis/ncv-interval-diagnosis.R plasmode <out_prefix> <app> [<app> ...]
# Writes <out_prefix>-cells.csv (one row per cell x estimand x interval) and
# <out_prefix>-points.csv (per cell x estimand: grid summaries of bias/SD diagnostics).
source("scripts/setup.R")
suppressPackageStartupMessages(library(dplyr))
args <- commandArgs(trailingOnly = TRUE)
mode <- args[1]
prefix <- args[2]
cores <- as.integer(Sys.getenv("RESTART_CORES", "1"))
z <- qnorm(1 - (1 - LEVEL) / 2)

if (mode == "synthetic") {
  cells <- design_table()
  dirs <- list.dirs("results", recursive = FALSE)
  if (nzchar(Sys.getenv("DIAG_DIRS")))
    dirs <- strsplit(Sys.getenv("DIAG_DIRS"), ",")[[1]]
  keep <- as.integer(args[-(1:2)])
} else {
  source("plasmode/sources.R")
  source("plasmode/plasmode.R")
  cells <- plasmode_cells()
  results <- Sys.getenv(
    "PLASMODE_RESULTS",
    file.path(PLASMODE_CACHE, "results")
  )
  dirs <- file.path(results, args[-(1:2)])
  keep <- integer(0)
}
files <- unlist(lapply(
  dirs,
  list.files,
  pattern = "^task-.*\\.rds$",
  full.names = TRUE
))
# file -> cell: synthetic names carry the task id only, so read the cell from each file
file_cell <- unlist(parallel::mclapply(files, mc.cores = cores, \(f) {
  r <- readRDS(f)
  if (identical(r$status, "success")) r$cell else NA_integer_
}))
ok <- !is.na(file_cell) & (length(keep) == 0 | file_cell %in% keep)
by_cell <- split(files[ok], file_cell[ok])

# Interval constructions on the link scale (est, half-width); `rec` is a list of
# replicate x point matrices.
intervals <- function(rec) {
  d <- rec$est_NCV - rec$est_REML
  vN <- rec$se_cl2_NCV^2
  vR <- rec$se_cl2_REML^2
  # noise allowance for delta: Var(N - R) ~ Var(R) - Var(N) if N is a smoother of R
  vd <- pmax(vR - vN, 0)
  list(
    reml_cl2 = list(est = rec$est_REML, half = z * rec$se_cl2_REML),
    ncv_cl2 = list(est = rec$est_NCV, half = z * rec$se_cl2_NCV),
    ncv_bias = list(est = rec$est_NCV, half = z * sqrt(vN + d^2)),
    # delta^2 minus its own noise variance, floored at 0
    ncv_bias_denoised = list(
      est = rec$est_NCV,
      half = z * sqrt(vN + pmax(d^2 - vd, 0))
    ),
    # NCV centre, larger of the two CL2 SEs
    ncv_maxse = list(
      est = rec$est_NCV,
      half = z * pmax(rec$se_cl2_NCV, rec$se_cl2_REML)
    ),
    # NCV centre, REML's SE (the paper's hybrid)
    hybrid = list(est = rec$est_NCV, half = z * rec$se_cl2_REML),
    # NCV centre, SE inflated by the replicate-free plug-in |delta| only where it exceeds
    # its noise level: bias-aware where delta is detectable, NCV + CL2 elsewhere
    ncv_bias_screened = list(
      est = rec$est_NCV,
      half = z * sqrt(vN + d^2 * (d^2 > z^2 * vd))
    )
  )
}

cell_summary <- function(fs) {
  rs <- lapply(fs, readRDS)
  cell <- cells[cells$cell == rs[[1]]$cell, ]
  linkinv <- FAMILIES[[cell$fit_family]]$linkinv
  cols <- c("truth", "est_REML", "se_cl2_REML", "est_NCV", "se_cl2_NCV")
  out_c <- list()
  out_p <- list()
  for (e in intersect(
    c("beta", "mean", "alpha", "gamma"),
    scored_estimands(cell, names(rs[[1]]$records))
  )) {
    if (!all(cols %in% colnames(rs[[1]]$records[[e]]))) next
    rec <- lapply(
      setNames(cols, cols),
      \(k) do.call(rbind, lapply(rs, \(r) r$records[[e]][, k]))
    )
    truth <- rec$truth
    resp <- e == "mean"
    tf <- if (resp) linkinv else identity
    # pointwise diagnostics on the estimand's own scale
    errN <- tf(rec$est_NCV) - truth
    errR <- tf(rec$est_REML) - truth
    biasN <- colMeans(errN)
    biasR <- colMeans(errR)
    sdN <- apply(tf(rec$est_NCV), 2, sd)
    sdR <- apply(tf(rec$est_REML), 2, sd)
    d <- tf(rec$est_NCV) - tf(rec$est_REML)
    n <- nrow(truth)
    # squared biases with the Monte Carlo noise of the replicate mean removed (unbiased)
    b2 <- \(bias, sd) mean(bias^2 - sd^2 / n)
    # per replicate: REML's relative error against the REML-NCV disagreement, a flag
    # available on real data (both standardised by the grid mean of the squared truth)
    den <- if (e %in% c("mean", "alpha"))
      mean((truth[1, ] - mean(truth[1, ]))^2) else mean(truth[1, ]^2)
    relR <- sqrt(rowMeans(errR^2) / den)
    relN <- sqrt(rowMeans(errN^2) / den)
    disagree <- sqrt(rowMeans(d^2) / den)
    top <- relR >= quantile(relR, 0.8)
    auc <- \(score, label) {
      r <- rank(score)
      (sum(r[label]) - sum(label) * (sum(label) + 1) / 2) /
        (sum(label) * sum(!label))
    }
    out_p[[e]] <- tibble(
      cell = cell$cell,
      estimand = e,
      n_rep = n,
      absbias_sd_ncv = mean(abs(biasN) / sdN),
      absbias_sd_reml = mean(abs(biasR) / sdR),
      rmsbias_ncv = sqrt(max(b2(biasN, sdN), 0)),
      rmsbias_reml = sqrt(max(b2(biasR, sdR), 0)),
      sd_ncv = sqrt(mean(sdN^2)),
      sd_reml = sqrt(mean(sdR^2)),
      # delta = NCV - REML as a bias plug-in: its mean is biasN - biasR, its noise SD
      delta_signal = sqrt(max(b2(colMeans(d), apply(d, 2, sd)), 0)),
      delta_noise_sd = sqrt(mean(apply(d, 2, var))),
      # share of NCV's squared bias that REML shares (projection; noise-corrected numerator)
      shared_bias = sum(
        biasR *
          biasN -
          vapply(seq_len(ncol(errR)), \(j) cov(errR[, j], errN[, j]), 0) / n
      ) /
        max(sum(biasN^2 - sdN^2 / n), 1e-12),
      # flag: Spearman correlation and AUC (top 20% of REML rel. error) of the disagreement
      flag_spearman = cor(disagree, relR, method = "spearman"),
      flag_auc = auc(disagree, top),
      rel_reml_median = median(relR),
      rel_ncv_median = median(relN)
    )
    for (nm in names(ints <- intervals(rec))) {
      it <- ints[[nm]]
      lo <- tf(it$est - it$half)
      hi <- tf(it$est + it$half)
      cen <- tf(it$est)
      bias_it <- colMeans(cen - truth)
      cov_pt <- colMeans(lo <= truth & truth <= hi)
      # bias-eliminated: shift each interval by the arm's pointwise bias
      be_pt <- colMeans(
        lo - rep(bias_it, each = nrow(lo)) <= truth &
          truth <= hi - rep(bias_it, each = nrow(lo))
      )
      out_c[[paste(e, nm)]] <- tibble(
        cell = cell$cell,
        estimand = e,
        interval = nm,
        n_rep = nrow(truth),
        coverage = mean(cov_pt),
        coverage_q05 = unname(quantile(cov_pt, 0.05)),
        coverage_bias_elim = mean(be_pt),
        halfwidth = mean((hi - lo) / 2),
        # interval score, alpha = 0.05
        score = mean(
          (hi - lo) +
            2 /
              (1 - LEVEL) *
              ((lo - truth) * (truth < lo) + (truth - hi) * (truth > hi))
        )
      )
    }
  }
  list(cells = bind_rows(out_c), points = bind_rows(out_p))
}

res <- parallel::mclapply(by_cell, mc.cores = cores, cell_summary)
keys <- intersect(
  c("cell", "app", "truth", "residual", "flip", "error", "G"),
  names(cells)
)
add_keys <- \(x)
  if (mode == "plasmode") left_join(x, cells[keys], by = "cell") else x
dir.create(dirname(prefix), recursive = TRUE, showWarnings = FALSE)
readr::write_csv(
  add_keys(bind_rows(lapply(res, `[[`, "cells"))),
  paste0(prefix, "-cells.csv")
)
readr::write_csv(
  add_keys(bind_rows(lapply(res, `[[`, "points"))),
  paste0(prefix, "-points.csv")
)
cat(sprintf("%d cells -> %s-*.csv\n", length(res), prefix))
