REL <- read_summary("final", "relative-error") |>
  left_join(
    S |> distinct(cell, block, family, signal, error, G, truth, basis, extra),
    by = "cell"
  ) |>
  mutate(dependence = if_else(error == "iid", "independent", "dependent"))
LSP <- read_summary("final", "log-sp-spread") |>
  left_join(
    S |> distinct(cell, block, family, signal, error, G, truth, basis),
    by = "cell"
  )
REL_ARMS <- c(
  reml_model = "REML, model-based",
  reml_cl2 = "REML + CL2",
  ncv_cl2 = "NCV + CL2",
  ncv_cl2_bias = "NCV + CL2, bias-aware",
  ar1_model = "AR(1) working model"
)
rel_sum <- REL |>
  filter(estimand %in% c("beta", "mean", "gamma", "alpha")) |>
  group_by(estimand, dependence, arm) |>
  summarise(
    cells = n(),
    coverage = median(coverage),
    rel_error_max = max(rel_rmse_median),
    cells_rel_gt_0.5 = sum(rel_rmse_median > 0.5),
    rel_error = median(rel_rmse_median),
    reps_rel_gt_1 = mean(share_rel_gt_1),
    halfwidth = median(rel_halfwidth_median),
    detect_min = min(detect),
    wrong_sign_max = max(wrong_sign),
    detect = median(detect),
    .groups = "drop"
  ) |>
  mutate(
    arm = factor(REL_ARMS[arm], REL_ARMS),
    estimand = label_estimand(estimand)
  ) |>
  select(
    estimand,
    dependence,
    arm,
    cells,
    coverage,
    rel_error,
    rel_error_max,
    cells_rel_gt_0.5,
    reps_rel_gt_1,
    halfwidth,
    detect,
    detect_min,
    wrong_sign_max
  ) |>
  arrange(estimand, dependence, arm)
save_tab(rel_sum, "relative-error-summary")
rs <- \(e, dep, a, col)
  rel_sum[[col]][
    rel_sum$estimand == label_estimand(e) &
      rel_sum$dependence == dep &
      rel_sum$arm == REL_ARMS[[a]]
  ]
rel_core <- REL |>
  filter(
    block == "core",
    estimand == "beta",
    arm %in% c("reml_cl2", "ncv_cl2")
  ) |>
  mutate(fit = if_else(arm == "reml_cl2", "REML", "NCV")) |>
  select(family, signal, error, G, fit, rel_rmse_median) |>
  pivot_wider(names_from = fit, values_from = rel_rmse_median) |>
  mutate(ratio = REML / NCV) |>
  arrange(family, error, G, signal)
save_tab(rel_core, "relative-error-core-beta")
rr <- \(dep)
  rel_core$ratio[
    if (dep == "iid") rel_core$error == "iid" else rel_core$error != "iid"
  ]
lsp_core <- LSP |>
  filter(
    block == "core",
    signal == "mid",
    fit %in% c("REML", "NCV"),
    grepl("^te\\(", term)
  ) |>
  group_by(family, error, G, fit) |>
  summarise(log_sp_sd = mean(log_sp_sd), .groups = "drop") |>
  pivot_wider(names_from = fit, values_from = log_sp_sd) |>
  arrange(family, error, G)
save_tab(lsp_core, "log-sp-spread-core")
ls <- \(f, dep)
  lsp_core[[f]][
    if (dep == "iid") lsp_core$error == "iid" else lsp_core$error != "iid"
  ]

d <- REL |>
  filter(estimand == "beta", arm %in% c("reml_cl2", "ncv_cl2")) |>
  mutate(fit = if_else(arm == "reml_cl2", "REML", "NCV")) |>
  select(cell, family, dependence, fit, rel_rmse_median) |>
  pivot_wider(names_from = fit, values_from = rel_rmse_median) |>
  mutate(family = label_family(family))
fig <- ggplot(d, aes(NCV, REML, colour = family)) +
  geom_abline(linetype = 2, colour = "grey40") +
  geom_hline(yintercept = 1, linetype = 3, colour = "grey40") +
  geom_vline(xintercept = 1, linetype = 3, colour = "grey40") +
  geom_point(size = 1.2) +
  facet_wrap(~dependence) +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    x = "NCV: relative error of beta",
    y = "REML: relative error of beta",
    colour = NULL
  ) +
  theme_doc
save_fig(fig, "relative-error-beta", 6.5, 3.6)
fig

show(
  rel_sum |>
    rename(
      dep. = dependence,
      cov. = coverage,
      rel = rel_error,
      `rel max` = rel_error_max,
      `cells >0.5` = cells_rel_gt_0.5,
      `reps >1` = reps_rel_gt_1,
      hw = halfwidth,
      det = detect,
      `det min` = detect_min,
      `wrong max` = wrong_sign_max
    ) |>
    mutate(dep. = substr(dep., 1, 3)),
  digits = 2
)

show(rel_core, digits = 2)

show(lsp_core, digits = 2)
