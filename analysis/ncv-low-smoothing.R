# Low-smoothing NCV fits of the coefficient surface in the first 95 settings with a function-on-function term: how often the surface's smoothing
# parameter lies far below the setting's median, and what these fits do to error and coverage, with the REML fit
# of the same replicates as control. Reads existing records only; no fits. Aggregates only. From the repository root:
#   Rscript analysis/ncv-low-smoothing.R <results_dir> <study.rds> summaries/final/ncv-low-smoothing.csv [per_rep.rds]
source("scripts/setup.R")
suppressPackageStartupMessages(library(dplyr))
args <- commandArgs(trailingOnly = TRUE)
results <- args[1]
st <- readRDS(args[2])
cells <- design_table()
# The term-type block has no function-on-function term.
cells <- cells[cells$cell <= 95 & cells$block != "term_type", ]
stopifnot(nrow(cells) == 86)
Z_GAP <- c(2, 3, 4)

# First penalty of the surface in both fits, and its EDF. `log_sp1` is computed after `log_sp2`
# because summarise() evaluates sequentially (no column is overwritten).
surface_sp <- st$terms |>
  filter(fit %in% c("NCV", "REML"), grepl("X", term), cell %in% cells$cell) |>
  group_by(fit, cell, rep) |>
  summarise(
    log_sp2 = last(log_sp),
    log_sp1 = first(log_sp),
    edf = first(edf),
    .groups = "drop"
  ) |>
  tidyr::pivot_wider(
    id_cols = c(cell, rep),
    names_from = fit,
    values_from = c(log_sp1, log_sp2, edf)
  )
stopifnot(all(count(surface_sp, cell)$n == N_REP), !anyNA(surface_sp))

rel_err <- function(rec, fit) {
  tr <- rec[, "truth"]
  sqrt(mean((rec[, paste0("est_", fit)] - tr)^2) / mean(tr^2))
}

task_row <- function(cell, rep) {
  r <- readRDS(file.path(
    results,
    sprintf("task-%d.rds", (rep - 1) * 95 + cell$cell)
  ))
  stopifnot(identical(r$status, "success"), r$cell == cell$cell, r$rep == rep)
  rec <- r$records$beta
  linkinv <- FAMILIES[[cell$fit_family]]$linkinv
  tibble(
    cell = cell$cell,
    rep = rep,
    flagged = !r$fits$converged[r$fits$fit == "NCV"],
    rel_ncv = rel_err(rec, "NCV"),
    rel_reml = rel_err(rec, "REML"),
    cov_ncv = mean(point_scores(rec, "ncv_cl2", "beta", linkinv)$covered),
    cov_reml = mean(point_scores(rec, "reml_cl2", "beta", linkinv)$covered)
  )
}

per_rep <- bind_rows(parallel::mclapply(
  seq_len(nrow(cells)),
  \(i) bind_rows(lapply(seq_len(N_REP), \(rep) task_row(cells[i, ], rep))),
  mc.cores = 4
)) |>
  left_join(surface_sp, by = c("cell", "rep")) |>
  left_join(
    as_tibble(cells)[c(
      "cell",
      "block",
      "family",
      "fit_family",
      "error",
      "G",
      "signal",
      "truth"
    )],
    by = "cell"
  ) |>
  group_by(cell) |>
  mutate(
    gap = median(log_sp1_NCV) - log_sp1_NCV,
    gap2 = median(log_sp2_NCV) - log_sp2_NCV,
    gap_reml = median(log_sp1_REML) - log_sp1_REML,
    sd_ncv = sd(log_sp1_NCV),
    sd2_ncv = sd(log_sp2_NCV),
    sd2_reml = sd(log_sp2_REML),
    within_edf = edf_NCV / median(edf_NCV),
    # shares of the setting's squared error, with each setting weighted equally
    wmse_ncv = rel_ncv^2 / mean(rel_ncv^2),
    wmse_reml = rel_reml^2 / mean(rel_reml^2),
    sd_reml = sd(log_sp1_REML),
    spread_ncv = unname(quantile(rel_ncv, 0.9)) / median(rel_ncv),
    spread_reml = unname(quantile(rel_reml, 0.9)) / median(rel_reml),
    # relative to the setting's median, so that settings with more low-smoothing fits do not drive the contrast
    within_ncv = rel_ncv / median(rel_ncv),
    within_reml = rel_reml / median(rel_reml),
    cov_ncv_c = cov_ncv - mean(cov_ncv),
    cov_reml_c = cov_reml - mean(cov_reml)
  ) |>
  ungroup() |>
  mutate(dependence = if_else(error == "iid", "independent", "dependent"))
stopifnot(!anyNA(per_rep$rel_ncv), !anyNA(per_rep$gap))
# As in scored_estimands(): beta is not a valid target under misregistration with a non-identity link, so only
# the smoothing parameters and EDFs of these settings are summarised.
UNSCORED <- c(
  "rel_ncv",
  "rel_reml",
  "cov_ncv",
  "cov_reml",
  "spread_ncv",
  "spread_reml",
  "within_ncv",
  "within_reml",
  "wmse_ncv",
  "wmse_reml",
  "cov_ncv_c",
  "cov_reml_c"
)
per_rep <- per_rep |>
  mutate(across(
    all_of(UNSCORED),
    \(x) if_else(error == "warp" & fit_family != "gaussian", NA_real_, x)
  ))
if (length(args) >= 4) saveRDS(per_rep, args[4])

summarise_low <- function(d, thr) {
  d |>
    mutate(low = gap > thr) |>
    summarise(
      threshold = thr,
      n_rep = n(),
      share_low = mean(low),
      share_low_either = mean(low | gap2 > thr),
      share_low_reml = mean(gap_reml > thr),
      reml_low_given_ncv_low = mean(gap_reml[low] > thr),
      reml_low_given_ncv_other = mean(gap_reml[!low] > thr),
      sd_ncv = median(sd_ncv[!duplicated(cell)]),
      sd_reml = median(sd_reml[!duplicated(cell)]),
      sd2_ncv = median(sd2_ncv[!duplicated(cell)]),
      sd2_reml = median(sd2_reml[!duplicated(cell)]),
      within_edf_low = median(within_edf[low]),
      gap2_low = median(gap2[low]),
      mse_share_ncv_low = sum(wmse_ncv[low]) / sum(wmse_ncv),
      mse_share_reml_low = sum(wmse_reml[low]) / sum(wmse_reml),
      # geometric mean over settings; REML's MSE over NCV's, on all replicates and without NCV's low-smoothing fits
      mse_ratio_all = exp(mean(log(
        tapply(rel_reml^2, cell, mean) / tapply(rel_ncv^2, cell, mean)
      ))),
      mse_ratio_other = exp(mean(log(
        tapply(rel_reml[!low]^2, cell[!low], mean) /
          tapply(rel_ncv[!low]^2, cell[!low], mean)
      ))),
      # the same contrasts, selected on REML's own tail
      within_ncv_remltail = median(within_ncv[gap_reml > thr]),
      within_reml_remltail = median(within_reml[gap_reml > thr]),
      cov_diff_ncv_remltail = mean(cov_ncv_c[gap_reml > thr]) -
        mean(cov_ncv_c[gap_reml <= thr]),
      cov_diff_reml_remltail = mean(cov_reml_c[gap_reml > thr]) -
        mean(cov_reml_c[gap_reml <= thr]),
      spread_ncv = median(spread_ncv[!duplicated(cell)]),
      spread_reml = median(spread_reml[!duplicated(cell)]),
      edf_ncv_low = median(edf_NCV[low]),
      edf_reml_low = median(edf_REML[low]),
      smoother_than_reml_low = mean(edf_NCV[low] < edf_REML[low]),
      n_cells = n_distinct(cell),
      within_ncv_low = median(within_ncv[low]),
      within_ncv_other = median(within_ncv[!low]),
      within_reml_low = median(within_reml[low]),
      within_reml_other = median(within_reml[!low]),
      ratio_ncv_reml_low = median(rel_ncv[low] / rel_reml[low]),
      ratio_ncv_reml_other = median(rel_ncv[!low] / rel_reml[!low]),
      cov_diff_ncv = mean(cov_ncv_c[low]) - mean(cov_ncv_c[!low]),
      cov_diff_reml = mean(cov_reml_c[low]) - mean(cov_reml_c[!low]),
      flagged_low = mean(flagged[low]),
      rel_ncv_low = median(rel_ncv[low]),
      rel_ncv_other = median(rel_ncv[!low]),
      rel_ncv_all = median(rel_ncv),
      rel_reml_low = median(rel_reml[low]),
      rel_reml_other = median(rel_reml[!low]),
      rel_mean_ncv_all = mean(rel_ncv),
      rel_mean_ncv_other = mean(rel_ncv[!low]),
      rel_mean_reml_all = mean(rel_reml),
      ncv_beats_reml_low = mean(rel_ncv[low] < rel_reml[low]),
      ncv_beats_reml_other = mean(rel_ncv[!low] < rel_reml[!low]),
      cov_ncv_low = mean(cov_ncv[low]),
      cov_ncv_other = mean(cov_ncv[!low]),
      cov_ncv_all = mean(cov_ncv),
      cov_reml_low = mean(cov_reml[low]),
      cov_reml_other = mean(cov_reml[!low]),
      .groups = "drop"
    )
}

out <- bind_rows(lapply(Z_GAP, \(thr) {
  bind_rows(
    per_rep |>
      group_by(dependence) |>
      summarise_low(thr) |>
      mutate(family = "all", error = "all", G = NA),
    per_rep |>
      group_by(dependence, family, error) |>
      summarise_low(thr) |>
      mutate(G = NA),
    per_rep |>
      filter(block != "warp") |>
      group_by(dependence) |>
      summarise_low(thr) |>
      mutate(family = "no warp", error = "all", G = NA),
    per_rep |>
      group_by(dependence, G) |>
      summarise_low(thr) |>
      mutate(family = "all", error = "all")
  )
}))
readr::write_csv(out, args[3])
