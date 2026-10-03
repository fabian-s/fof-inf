PREL <- read_summary("plasmode", "relative-error") |>
  mutate(
    role = if_else(app %in% APP_COUNTED, "counted", "stress test"),
    G_label = if_else(G > 40, "all", "40")
  )
prel_app <- PREL |>
  filter(
    flip == "curve",
    G_label == "all",
    error == "attached",
    arm %in% c("reml_cl2", "ncv_cl2", "ncv_cl2_bias")
  ) |>
  group_by(role, app, estimand, arm) |>
  summarise(
    cells = n(),
    coverage = mean(coverage),
    rel_error_max = max(rel_rmse_median),
    rel_error = median(rel_rmse_median),
    halfwidth = median(rel_halfwidth_median),
    detect = median(detect),
    .groups = "drop"
  )
save_tab(prel_app, "plasmode-relative-error")
pv <- \(e, a, col, role_sel = c("counted", "stress test"))
  prel_app[[col]][
    prel_app$estimand == e & prel_app$arm == a & prel_app$role %in% role_sel
  ]
prel_beta <- prel_app |>
  filter(estimand == "beta") |>
  mutate(
    arm = recode(
      arm,
      reml_cl2 = "REML+CL2",
      ncv_cl2 = "NCV+CL2",
      ncv_cl2_bias = "bias-aware"
    ),
    app = label_app(app)
  ) |>
  select(
    role,
    app,
    arm,
    cov = coverage,
    rel = rel_error,
    hw = halfwidth,
    det = detect
  ) |>
  pivot_wider(
    names_from = arm,
    values_from = c(cov, rel, hw, det),
    names_sep = " "
  ) |>
  select(-`rel bias-aware`) |>
  arrange(role, app)

show(prel_beta, digits = 2)
