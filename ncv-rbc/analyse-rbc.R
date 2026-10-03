# Summarise the bias-corrected NCV results per cell x estimand x arm (processed cell by
# cell): grid-averaged coverage (MC SE over replicates), 5% quantile of pointwise coverage,
# mean half-width, relative error and half-width (as analysis/relative-error.R),
# detection rate. REML + CL2 for the same cells and replicates comes from the study
# summaries (summaries/final/cells.csv, relative-error.csv).
#   Rscript ncv-rbc/analyse-rbc.R <results_dir> <out.csv> [rep_from] [rep_to]
# Also: interval score (width + 2/alpha x miss distance, alpha = 0.05; lower is better), a derived
# arm `hybrid` (NCV estimate with REML + CL2's half-width, where both arms are in the files; not for
# the mean, whose intervals are asymmetric on the response scale), and an oracle matched-coverage
# view: each arm's half-widths scaled by the factor c95 that makes the cell's pooled coverage 0.95,
# with detection and half-width at that common coverage (link scale; not for the mean).
suppressPackageStartupMessages(library(dplyr))
args <- commandArgs(trailingOnly = TRUE)
dir <- args[1] # task files c<cell>-r<rep>.rds; raw outputs stay on LRZ (rbc-results, rbc-fresh-results, rbc-pm2-results, fixlam-pm-results, combo*-results)
out <- if (length(args) >= 2) args[2] else "ncv-rbc/summaries/rbc-cells.csv"
files <- list.files(dir, "^c[0-9]+-r[0-9]+\\.rds$", full.names = TRUE)
rep_of <- as.integer(sub(".*-r([0-9]+)\\.rds$", "\\1", basename(files)))
if (length(args) >= 4)
  files <- files[rep_of >= as.integer(args[3]) & rep_of <= as.integer(args[4])]
cell_of <- as.integer(sub("^c([0-9]+)-.*", "\\1", basename(files)))
ALPHA <- 0.05
CENTRED <- c("mean", "alpha")

summarise_cell <- function(fs) {
  pts <- bind_rows(lapply(fs, \(f) {
    r <- readRDS(f)
    if (!is.null(r$error)) return(NULL)
    r$points
  }))
  # hybrid: NCV estimate, REML + CL2 half-width (same replicate and point)
  if (all(c("ncv_cl2", "reml_cl2") %in% pts$arm)) {
    h <- pts |>
      filter(arm == "ncv_cl2", estimand != "mean") |>
      inner_join(
        pts |>
          filter(arm == "reml_cl2") |>
          transmute(estimand, rep, point, hR = (upper - lower) / 2),
        by = c("estimand", "rep", "point")
      ) |>
      mutate(arm = "hybrid", lower = est - hR, upper = est + hR) |>
      select(-hR)
    pts <- bind_rows(pts, h)
  }
  # denominators and the detection threshold per replicate, as in analysis/relative-error.R
  pts |>
    group_by(estimand, rep, arm) |>
    mutate(
      denom = if (first(estimand) %in% CENTRED)
        mean((truth - mean(truth))^2) else mean(truth^2),
      big = abs(truth) > 0.25 * max(abs(truth))
    ) |>
    ungroup() |>
    mutate(
      covered = lower <= truth & truth <= upper,
      half = (upper - lower) / 2,
      hit = (lower > 0 & truth > 0) | (upper < 0 & truth < 0),
      wrong = (lower > 0 & truth <= 0) | (upper < 0 & truth >= 0),
      score = (upper - lower) +
        2 /
          ALPHA *
          ((lower - truth) *
            (truth < lower) +
            (truth - upper) * (truth > upper))
    ) |>
    group_by(cell, estimand, arm) |>
    group_modify(\(d, k) {
      per_rep <- d |>
        group_by(rep) |>
        summarise(
          cov = mean(covered),
          half = mean(half),
          rel = sqrt(mean((est - truth)^2) / first(denom)),
          relhw = sqrt(mean(half^2) / first(denom)),
          det = if (k$estimand %in% CENTRED) NA_real_ else mean(hit[big]),
          wrong = if (k$estimand %in% CENTRED) NA_real_ else mean(wrong),
          score = mean(score),
          .groups = "drop"
        )
      # oracle matched coverage (symmetric intervals only)
      matched <- if (k$estimand == "mean")
        list(c95 = NA_real_, det = NA_real_, half = NA_real_) else {
        c95 <- unname(quantile(abs(d$est - d$truth) / d$half, 0.95))
        hitm <- (d$est - c95 * d$half > 0 & d$truth > 0) |
          (d$est + c95 * d$half < 0 & d$truth < 0)
        list(
          c95 = c95,
          det = if (k$estimand %in% CENTRED) NA_real_ else
            mean(tapply(hitm[d$big], d$rep[d$big], mean)),
          half = c95 * mean(d$half)
        )
      }
      z975 <- qnorm(0.975)
      pc <- d |>
        group_by(point) |>
        mutate(
          se = half / z975,
          zc = (est - mean(est)) / se,
          zt = (est - truth) / se
        ) |>
        summarise(
          pc = mean(covered),
          sd_est = sd(est),
          rms_se = sqrt(mean(se^2)),
          zc = list(zc),
          zt = list(zt),
          .groups = "drop"
        )
      tibble(
        n_rep = nrow(per_rep),
        coverage = mean(per_rep$cov),
        coverage_se = sd(per_rep$cov) / sqrt(nrow(per_rep)),
        coverage_q05 = unname(quantile(pc$pc, 0.05)),
        halfwidth = mean(per_rep$half),
        rel_rmse_median = median(per_rep$rel),
        rel_halfwidth_median = median(per_rep$relhw),
        detect = mean(per_rep$det),
        wrong_sign = mean(per_rep$wrong),
        score = mean(per_rep$score),
        c95 = matched$c95,
        detect_matched = matched$det,
        halfwidth_matched = matched$half,
        # Variance calibration (link scale; the mean's intervals are not symmetric on the
        # response scale, so it is left out): median over grid points of RMS SE / replicate SD
        # of the estimate, and the pooled SD of the studentized error with (sd_z) and without
        # (sd_zc, centred at the pointwise replicate mean) the bias. sd_zc ~ 1: SE right.
        se_sd = if (k$estimand == "mean") NA_real_ else
          median(pc$rms_se / pc$sd_est),
        sd_zc = if (k$estimand == "mean") NA_real_ else sd(unlist(pc$zc)),
        sd_z = if (k$estimand == "mean") NA_real_ else
          sqrt(mean(unlist(pc$zt)^2))
      )
    }) |>
    ungroup()
}

res <- bind_rows(lapply(split(files, cell_of), summarise_cell))
n_fail <- sum(vapply(files, \(f) !is.null(readRDS(f)$error), TRUE))
dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
readr::write_csv(res, out)
cat(sprintf(
  "%d files (%d failed), %d cells -> %s\n",
  length(files),
  n_fail,
  n_distinct(res$cell),
  out
))
