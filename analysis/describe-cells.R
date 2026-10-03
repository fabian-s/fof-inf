# Descriptive per-cell quantities of the data-generating process (no model fits):
# realized response-scale design effect (1'S1 / tr S of the residual covariance), and
# for warp cells the phase share of the residual variance and the gap between the
# observed-clock (smeared) and registered truths.
# Run from the repository root: Rscript analysis/describe-cells.R <out.csv> [n_rep]
source("scripts/setup.R")
source("analysis/dgp-helpers.R")
args <- commandArgs(trailingOnly = TRUE)
n_rep <- if (length(args) >= 2) as.integer(args[2]) else 20L

design_effect <- function(residuals) {
  S <- cov(residuals)
  sum(S) / sum(diag(S))
}

# Relative sup-norm gap between observed-clock and registered targets (warp cells).
target_gap <- function(cell) {
  iid <- modifyList(cell, list(error = "iid"))
  a <- simulate(cell, 1)$truth
  b <- simulate(iid, 1)$truth
  sapply(
    c("mean", "alpha", "beta"),
    \(e) max(abs(a[[e]] - b[[e]])) / max(abs(b[[e]]))
  )
}

describe_cell <- function(cell) {
  residuals <- do.call(
    rbind,
    lapply(seq_len(n_rep), \(r) residual_curves(cell, r))
  )
  out <- tibble::tibble(
    cell = cell$cell,
    design_effect = design_effect(residuals)
  )
  # The phase share subtracts the white-noise variance; this decomposition is exact only
  # for Gaussian responses (other families' noise variance depends on the warped mean).
  if (cell$error != "warp" || cell$family != "gaussian") return(out)
  noise <- modifyList(cell, list(error = "iid"))
  noise_var <- mean(
    do.call(rbind, lapply(seq_len(n_rep), \(r) residual_curves(noise, r)))^2
  )
  gap <- target_gap(cell)
  dplyr::mutate(
    out,
    phase_share = 1 - noise_var / mean(residuals^2),
    gap_mean = gap[["mean"]],
    gap_alpha = gap[["alpha"]],
    gap_beta = gap[["beta"]]
  )
}

cells <- design_table()
cells <- cells[cells$block != "families" | cells$fit_family == cells$family, ]
description <- dplyr::bind_rows(lapply(
  seq_len(nrow(cells)),
  \(i) describe_cell(as.list(cells[i, ]))
))
readr::write_csv(description, args[1])
cat("described", nrow(description), "cells\n")
