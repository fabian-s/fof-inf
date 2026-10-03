NR <- \(f)
  readr::read_csv(file.path("ncv-rbc/summaries", f), show_col_types = FALSE)
NCVI_ARMS <- c(
  reml_cl2 = "REML + CL2",
  ncv_cl2 = "NCV + CL2",
  rbc_p = "bias-corrected (1 step)",
  rbc2_f = "bias-corrected (2 steps, freq.)",
  ncv_ij = "NCV + λ term",
  ncv_ijx = "NCV + λ term (cross)",
  rbc_p_ij = "1 step + λ term",
  rbc_p_ijx = "1 step + λ term (cross)",
  rbc2_p_ijx = "2 steps + λ term (cross)",
  ncvfix_cl2 = "NCV, λ fixed at median (oracle)",
  hybrid = "NCV estimate, REML + CL2 width"
)
pm_keys <- read_summary("plasmode", "ncv-diag-cells") |>
  distinct(cell, app, truth)
syn_keys <- S |> distinct(cell, block, family, error)
ncvi_fmt <- \(d) d |> mutate(arm = factor(NCVI_ARMS[arm], NCVI_ARMS))
# synthetic, single-fit corrections (all families), fresh replicates
nv_syn <- NR("rbc-fresh-cells.csv") |>
  left_join(syn_keys, by = "cell") |>
  filter(
    estimand == "beta",
    error != "iid",
    arm %in% c("reml_cl2", "ncv_cl2", "rbc_p", "rbc2_f")
  ) |>
  group_by(cell) |>
  mutate(width = halfwidth / halfwidth[arm == "reml_cl2"]) |>
  ungroup() |>
  mutate(
    group = paste(
      if_else(block == "core", "smooth truth", "rough truth"),
      if_else(family == "binomial", "binary", "Gaussian/Poisson"),
      sep = ", "
    )
  ) |>
  group_by(group, arm) |>
  summarise(
    cells = n(),
    cov_min = min(coverage),
    coverage = mean(coverage),
    q05 = mean(coverage_q05),
    width = mean(width),
    se_sd = mean(se_sd),
    detect = mean(detect),
    .groups = "drop"
  ) |>
  ncvi_fmt() |>
  arrange(group, arm)
save_tab(nv_syn, "ncv-intervals-synthetic")
# synthetic combinations (Gaussian), reps 201-300; REML + CL2 from the same cells (reps 201-400)
nv_comb <- bind_rows(
  NR("combo2-syn-cells.csv"),
  NR("rbc-fresh-cells.csv") |>
    filter(arm == "reml_cl2", cell %in% NR("combo2-syn-cells.csv")$cell)
) |>
  left_join(syn_keys, by = "cell") |>
  filter(estimand == "beta") |>
  group_by(cell) |>
  mutate(width = halfwidth / halfwidth[arm == "reml_cl2"]) |>
  ungroup() |>
  mutate(group = if_else(block == "core", "smooth truth", "rough truth")) |>
  group_by(group, arm) |>
  summarise(
    coverage = mean(coverage),
    q05 = mean(coverage_q05),
    width = mean(width),
    se_sd = mean(se_sd),
    .groups = "drop"
  ) |>
  filter(arm %in% names(NCVI_ARMS)) |>
  ncvi_fmt() |>
  arrange(group, arm)
save_tab(nv_comb, "ncv-intervals-combination-synthetic")
# plasmode per dataset
pm_all <- bind_rows(
  NR("combo2-pm-cells.csv"),
  NR("fixlam-plasmode-cells.csv"),
  NR("rbc-plasmode-cells.csv") |> filter(arm == "reml_cl2")
) |>
  left_join(pm_keys, by = "cell") |>
  filter(estimand == "beta", !is.na(app))
nv_pm <- pm_all |>
  filter(
    arm %in%
      c(
        "reml_cl2",
        "ncv_cl2",
        "ncvfix_cl2",
        "rbc_p",
        "ncv_ij",
        "rbc_p_ij",
        "rbc_p_ijx"
      )
  ) |>
  group_by(app, arm) |>
  summarise(
    v = sprintf(
      "%.2f (%.2f) %.2f",
      mean(coverage),
      mean(coverage_q05),
      mean(se_sd)
    ),
    .groups = "drop"
  ) |>
  mutate(arm = NCVI_ARMS[arm], app = label_app(app)) |>
  pivot_wider(names_from = arm, values_from = v) |>
  arrange(app)
save_tab(nv_pm, "ncv-intervals-plasmode")
pmv <- \(a, col = "coverage")
  pm_all |>
    filter(arm == a) |>
    group_by(app) |>
    summarise(v = mean(.data[[col]]), .groups = "drop") |>
    pull(v)
# misspecified truths
ms_keys <- readr::read_csv(
  "misspec/summaries/misspec-cells.csv",
  show_col_types = FALSE
) |>
  distinct(cell, truth_label, error, G, basis)
nv_ms <- bind_rows(
  NR("combo2-ms-cells.csv") |>
    filter(arm %in% c("ncv_cl2", "rbc_p", "rbc_p_ij", "rbc2_p_ijx")),
  readr::read_csv(
    "misspec/summaries/misspec-cells.csv",
    show_col_types = FALSE
  ) |>
    filter(
      arm == "reml_cl2",
      basis == "default",
      cell %in% NR("combo2-ms-cells.csv")$cell
    ) |>
    transmute(cell, estimand, arm, coverage, coverage_q05 = pointwise_q05)
) |>
  left_join(ms_keys, by = "cell") |>
  filter(estimand == "beta") |>
  mutate(
    setting = paste(
      truth_label,
      if_else(error == "iid", "independent", "smooth errors"),
      paste0("G = ", G)
    ),
    v = sprintf("%.2f (%.2f)", coverage, coverage_q05),
    arm = NCVI_ARMS[arm]
  ) |>
  select(setting, arm, v) |>
  pivot_wider(names_from = arm, values_from = v) |>
  arrange(setting)
save_tab(nv_ms, "ncv-intervals-misspec")
# interval scores (paired replicates)
sc <- bind_rows(
  bind_rows(
    NR("score-fresh-cells.csv") |>
      filter(arm %in% c("reml_cl2", "ncv_cl2", "hybrid")),
    NR("score-combo-syn-cells.csv") |> filter(arm == "rbc_p_ijx")
  ) |>
    filter(cell %in% NR("score-combo-syn-cells.csv")$cell) |>
    left_join(syn_keys, by = "cell") |>
    mutate(
      setting = if_else(
        block == "core",
        "synthetic, smooth truth",
        "synthetic, rough truth"
      )
    ),
  bind_rows(
    NR("score-pm2-cells.csv") |>
      filter(arm %in% c("reml_cl2", "ncv_cl2", "hybrid")),
    NR("score-combo-pm-cells.csv") |> filter(arm == "rbc_p_ijx")
  ) |>
    left_join(pm_keys, by = "cell") |>
    mutate(setting = paste0("plasmode, ", truth, "-fitted truth"))
) |>
  filter(estimand == "beta") |>
  group_by(cell, setting) |>
  mutate(
    score_rel = score / score[arm == "reml_cl2"],
    det_rel = detect - detect[arm == "reml_cl2"]
  ) |>
  ungroup() |>
  group_by(setting, arm) |>
  summarise(
    coverage = mean(coverage),
    q05 = mean(coverage_q05),
    score_vs_REML = exp(mean(log(score_rel))),
    detection = mean(detect),
    .groups = "drop"
  ) |>
  ncvi_fmt() |>
  arrange(setting, arm)
save_tab(sc, "ncv-intervals-scores")

show(nv_syn, digits = 3)

show(nv_comb, digits = 3)

show(nv_pm)

show(nv_ms)

show(sc, digits = 3)
