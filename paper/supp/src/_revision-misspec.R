MS <- readr::read_csv(
  "misspec/summaries/misspec-cells.csv",
  show_col_types = FALSE
)
MSC <- readr::read_csv(
  "misspec/summaries/misspec-contrasts.csv",
  show_col_types = FALSE
)
MST <- readr::read_csv(
  "misspec/summaries/misspec-truths.csv",
  show_col_types = FALSE
)
MSCRN <- readr::read_csv(
  "misspec/summaries/misspec-crn.csv",
  show_col_types = FALSE
)
MS_ARMS <- c(
  reml_model = "REML model-based",
  reml_cl2 = "REML + CL2",
  reml_cl2f = "REML + CL2 (freq.)",
  ncv_cl2 = "NCV + CL2"
)
MS_TRUTH <- c(T0 = "T0 smooth", T1 = "T1 bump", T2 = "T2 weak")
ms_cell <- \(error, G)
  paste0(if_else(error == "iid", "iid", "smooth"), ", G = ", G)
MS <- MS |>
  mutate(
    truth = MS_TRUTH[truth_label],
    setting = ms_cell(error, G),
    arm_l = factor(MS_ARMS[arm], MS_ARMS)
  )
msv <- \(tr, err, g, b, est, a, col = "coverage") {
  v <- MS[[col]][
    MS$truth_label == tr &
      MS$error == err &
      MS$G == g &
      MS$basis == b &
      MS$estimand == est &
      MS$arm == a
  ]
  stopifnot(length(v) == 1)
  v
}

d <- MS |>
  filter(estimand %in% c("beta", "mean"), arm != "reml_cl2f") |>
  mutate(
    estimand = label_estimand(estimand),
    arm_l = factor(arm_l, rev(levels(droplevels(arm_l))))
  )
db <- d |> filter(truth_label == "T1", estimand == label_estimand("beta"))
fig <- ggplot(d, aes(coverage, arm_l, colour = basis)) +
  geom_vline(xintercept = LEVEL, linetype = 2, colour = "grey40") +
  geom_pointrange(
    aes(xmin = coverage - 2 * coverage_se, xmax = coverage + 2 * coverage_se),
    position = position_dodge(width = 0.6),
    size = 0.15
  ) +
  geom_point(
    data = db,
    aes(x = coverage_bump),
    shape = 1,
    size = 1.5,
    position = position_dodge(width = 0.6)
  ) +
  facet_grid(truth + estimand ~ setting) +
  scale_colour_manual(
    values = c(default = "#D55E00", xlarge = "#0072B2"),
    name = "basis"
  ) +
  labs(x = "coverage of the exact truth (grid average)", y = NULL) +
  theme_doc
save_fig(fig, "misspec-coverage", 7, 6)
fig

ms_truths <- MST |>
  transmute(
    truth = MS_TRUTH[truth],
    amplitude = feature_amplitude,
    `share ∫β²` = feature_l2_share,
    `rel. L2 err. default` = rel_l2_default,
    `resid. signal default` = resid_signal_share_default,
    `rel. L2 err. xlarge` = rel_l2_xlarge,
    `resid. signal xlarge` = resid_signal_share_xlarge
  )
save_tab(ms_truths, "misspec-truths")
ms_wide <- MS |>
  mutate(cov = sprintf("%.3f (%.3f)", coverage, coverage_se)) |>
  select(estimand, truth, setting, basis, arm_l, cov) |>
  pivot_wider(names_from = arm_l, values_from = cov) |>
  arrange(estimand, truth, setting, basis)
save_tab(ms_wide, "misspec-coverage-wide")
ms_bump <- MS |>
  filter(truth_label == "T1", estimand == "beta") |>
  select(
    setting,
    basis,
    arm = arm_l,
    bump_region = coverage_bump,
    elsewhere = coverage_off,
    bump_region_bias_elim = coverage_bump_bias_elim
  ) |>
  arrange(setting, basis, arm)
save_tab(ms_bump, "misspec-bump-region")
ms_bias <- MS |>
  filter(
    estimand %in% c("beta", "mean"),
    arm %in% c("reml_cl2", "ncv_cl2"),
    truth_label != "T0"
  ) |>
  select(
    estimand,
    truth,
    setting,
    basis,
    arm = arm_l,
    coverage,
    coverage_bias_elim,
    pointwise_q05,
    mean_absbias_sd,
    se_sd_ratio
  ) |>
  arrange(estimand, truth, setting, basis, arm)
save_tab(ms_bias, "misspec-bias")
ms_c <- MSC |>
  filter(
    contrast == "xlarge - default",
    estimand %in% c("beta", "mean"),
    arm %in% c("reml_cl2", "ncv_cl2")
  ) |>
  mutate(
    truth = MS_TRUTH[truth_label],
    setting = ms_cell(error, G),
    arm = MS_ARMS[arm]
  ) |>
  select(estimand, truth, setting, arm, n_pairs, diff, diff_se) |>
  arrange(estimand, truth, setting, arm)
save_tab(ms_c, "misspec-basis-contrasts")
ms_mse <- MS |>
  filter(estimand %in% c("beta", "mean"), arm %in% c("reml_cl2", "ncv_cl2")) |>
  mutate(
    fit = if_else(arm == "reml_cl2", "REML", "NCV"),
    bias_share = bias2 / mse
  ) |>
  select(estimand, truth, setting, basis, fit, mse, bias_share) |>
  pivot_wider(names_from = fit, values_from = c(mse, bias_share)) |>
  arrange(estimand, truth, setting, basis)
save_tab(ms_mse, "misspec-mse")

# the contrast file must hold the basis contrast used below
stopifnot(nrow(ms_c) > 0)

show(
  ms_truths |> mutate(across(where(is.numeric), \(x) signif(x, 2))),
  digits = 8
)

show(ms_wide)

show(ms_bump)

show(ms_bias)

show(ms_c)

show(ms_mse, digits = 4)
