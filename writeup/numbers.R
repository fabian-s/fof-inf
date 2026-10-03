# Headline objects shared by writeup/RESULTS.qmd and paper/paper.qmd (sourced after
# writeup/helpers.R, the repository root as working directory).
# Core dependent cells at G = 100 (both signal levels of GLMs, three of Gaussian).
core_dep <- S |>
  filter(
    block == "core",
    G == 100,
    error != "iid",
    estimand %in% c("mean", "beta")
  )
core_iid <- S |>
  filter(
    block == "core",
    G == 100,
    error == "iid",
    estimand %in% c("mean", "beta")
  )
core_dep40 <- S |>
  filter(
    block == "core",
    G == 40,
    error != "iid",
    estimand %in% c("mean", "beta")
  )
core_iid40 <- S |>
  filter(
    block == "core",
    G == 40,
    error == "iid",
    estimand %in% c("mean", "beta")
  )
cov_of <- \(d, a) d$coverage_mean[d$arm == a]
q05_of <- \(d, a) d$pointwise_q05[d$arm == a]
se_of <- \(d, a) d$coverage_se[d$arm == a]
dense <- MC |> filter(block == "dense_grid", estimand == "beta")
# plasmode pool means; AR(1) principal cell = (MID truth, REML residuals)
POOL_MEANS_MAIN <- POOL_MEANS |>
  filter(arm != "ar1_model") |>
  bind_rows(
    P_POOL |>
      filter(arm == "ar1_model", truth == "MID", residual == "REML") |>
      transmute(
        app,
        role,
        estimand,
        arm,
        n_cells = 1L,
        coverage = coverage_mean,
        coverage_min = coverage_mean,
        coverage_max = coverage_mean,
        q05_min = pointwise_q05,
        q05_mean = pointwise_q05
      )
  )
pm_counted <- POOL_MEANS_MAIN |>
  filter(role == "counted", estimand %in% c("mean", "beta"))
pm_stress <- POOL_MEANS_MAIN |>
  filter(role == "stress test", estimand %in% c("mean", "beta"))
pm_all <- POOL_MEANS_MAIN |> filter(estimand %in% c("mean", "beta"))
pmc <- \(a, e = c("mean", "beta"))
  pm_counted$coverage[pm_counted$arm == a & pm_counted$estimand %in% e]
pms <- \(a, e = c("mean", "beta"))
  pm_stress$coverage[pm_stress$arm == a & pm_stress$estimand %in% e]
pma <- \(a, e = c("mean", "beta"))
  pm_all$coverage[pm_all$arm == a & pm_all$estimand %in% e]
# G = 40 vs all: paired per (app, truth, residual, estimand, arm)
g40 <- P |>
  filter(
    error == "attached",
    flip == "curve",
    estimand %in% c("mean", "beta"),
    arm %in% c(MAIN_ARMS, "ar1_model")
  ) |>
  select(app, role, truth, residual, estimand, arm, G_label, coverage_mean) |>
  pivot_wider(names_from = G_label, values_from = coverage_mean) |>
  filter(!is.na(all), !is.na(`40`)) |>
  mutate(loss = all - `40`)
g40_app <- g40 |>
  group_by(app, role, estimand, arm) |>
  summarise(loss = mean(loss), .groups = "drop")
# NCV / REML MSE ratios, synthetic dependent core cells at G = 100
r_dep <- MSE_RATIO |> filter(block == "core", G == 100, error != "iid")
r_iid <- MSE_RATIO |> filter(block == "core", G == 100, error == "iid")
ratio_of <- \(d, e) d$ratio[d$estimand == e]
rt <- \(tr, e, role_sel = c("counted", "stress test")) {
  RATIO_BY_TRUTH$ratio[
    RATIO_BY_TRUTH$truth == tr &
      RATIO_BY_TRUTH$estimand == e &
      RATIO_BY_TRUTH$role %in% role_sel
  ]
}
# per truth source, beta, G = all, model residuals: coverage and 5% quantile per arm
by_truth <- P_POOL |>
  filter(estimand %in% c("mean", "beta"), residual %in% c("REML", "NCV")) |>
  group_by(app, role, truth, estimand, arm) |>
  summarise(
    coverage = mean(coverage_mean),
    q05 = mean(pointwise_q05),
    width = mean(width_mean),
    .groups = "drop"
  )
bt <- \(
  a,
  tr,
  e = "beta",
  col = "coverage",
  role_sel = c("counted", "stress test")
) {
  by_truth[[col]][
    by_truth$arm == a &
      by_truth$truth == tr &
      by_truth$estimand == e &
      by_truth$role %in% role_sel
  ]
}
# AR(1): principal cell and sensitivity cells
ar1_main <- P_POOL |>
  filter(
    arm == "ar1_model",
    truth == "MID",
    residual == "REML",
    estimand %in% c("mean", "beta")
  )
ar1_sens <- P_POOL |>
  filter(
    arm == "ar1_model",
    !(truth == "MID" & residual == "REML"),
    estimand %in% c("mean", "beta")
  )
ar1_s <- S |> filter(arm == "ar1_model", estimand %in% c("mean", "beta"))
# hybrid vs REML + CL2 under the REML-fitted truth (beta, G = all, model residuals)
hyb_reml <- HSUM_P |>
  filter(
    flip == "curve",
    error == "attached",
    G > 40,
    residual %in% c("REML", "NCV"),
    estimand == "beta",
    truth == "REML"
  )
hyb_q05_reml <- mean(hyb_reml$pointwise_q05[hyb_reml$arm == "hybrid"])
cl2_q05_reml <- mean(hyb_reml$pointwise_q05[hyb_reml$arm == "reml_cl2"])
# Interval-score ratios everywhere in this document: geometric mean over cells of the
# per-cell ratio of mean interval scores (the definition of the recommendation rule).
hyb_by_truth <- HSUM_P |>
  filter(
    flip == "curve",
    error == "attached",
    G > 40,
    residual %in% c("REML", "NCV"),
    estimand == "beta",
    arm %in% c("reml_cl2", "hybrid")
  ) |>
  select(cell, truth, arm, coverage_mean, score_mean) |>
  pivot_wider(names_from = arm, values_from = c(coverage_mean, score_mean)) |>
  group_by(truth) |>
  summarise(
    cov_hybrid = mean(coverage_mean_hybrid),
    cov_reml = mean(coverage_mean_reml_cl2),
    is_ratio = exp(mean(log(score_mean_hybrid / score_mean_reml_cl2))),
    .groups = "drop"
  )
hv2s <- \(what, tr) {
  h <- hyb_by_truth
  switch(
    what,
    "cov hybrid" = h$cov_hybrid[h$truth == tr],
    "cov REML+CL2" = h$cov_reml[h$truth == tr],
    "IS hybrid/REML+CL2" = h$is_ratio[h$truth == tr]
  )
}
za_bin_core <- if (!is.null(ZS))
  ZS |>
    filter(
      fit == "REML",
      block == "core",
      family == "binomial",
      error != "iid",
      estimand == "alpha"
    ) else NULL
# MC precision of stored grid-average coverages
se_dep_cl2 <- se_of(core_dep, "reml_cl2")
se_dep_model <- se_of(core_dep, "reml_model")
se_plas_cl2 <- P_POOL$coverage_se[
  P_POOL$arm == "reml_cl2" & P_POOL$estimand %in% c("mean", "beta")
]
# exceptions for summary item 1
bin_ag <- S |>
  filter(
    block == "core",
    G == 100,
    family == "binomial",
    estimand %in% c("alpha", "gamma"),
    arm %in% c("reml_cl2", "reml_model")
  )
warp_pois <- S |>
  filter(
    block == "warp",
    family == "poisson",
    estimand == "mean",
    arm %in% c("reml_cl2", "ncv_cl2_bias", "reml_model")
  )
iid40_cost <- S |>
  filter(
    block == "core",
    G == 40,
    error == "iid",
    estimand %in% c("mean", "beta"),
    arm %in% c("reml_model", "reml_cl2")
  ) |>
  select(family, signal, estimand, arm, coverage_mean) |>
  pivot_wider(names_from = arm, values_from = coverage_mean) |>
  mutate(cost = reml_model - reml_cl2)
# bias-aware and hybrid widths relative to REML + CL2 (paired cell means)
width_ratio <- S |>
  filter(
    estimand %in% c("mean", "beta"),
    arm %in% c("reml_cl2", "ncv_cl2_bias", "ncv_cl2")
  ) |>
  select(
    cell,
    block,
    family,
    error,
    G,
    signal,
    truth,
    basis,
    estimand,
    arm,
    width_mean
  ) |>
  pivot_wider(names_from = arm, values_from = width_mean) |>
  mutate(bias_aware = ncv_cl2_bias / reml_cl2, ncv_cl2 = ncv_cl2 / reml_cl2)
wr_dep <- width_ratio |> filter(block == "core", G == 100, error != "iid")
