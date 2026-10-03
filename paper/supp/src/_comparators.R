CMP <- readr::read_csv(
  "comparators/summaries/comparators-cells.csv",
  show_col_types = FALSE
)
CMPF <- readr::read_csv(
  "comparators/summaries/comparators-fits.csv",
  show_col_types = FALSE
)
CMP_ARMS <- c(
  reml_model = "REML model-based",
  reml_cl2 = "REML + CL2",
  pcre = "pcre",
  gls_fpc = "GLS (FPCA cov.)",
  gls_raw = "GLS (raw cov.)",
  boot_pct = "bootstrap percentile",
  boot_wald = "bootstrap Wald"
)
# NCV estimates for the same cells (study records; same replicates by common random numbers)
ncv_syn <- S |>
  filter(
    cell %in% CMP$cell[CMP$source == "synthetic"],
    arm %in% c("ncv_cl2", "reml_cl2")
  ) |>
  transmute(source = "synthetic", cell, estimand, arm, mse)
ncv_pm <- read_summary("plasmode", "coverage") |>
  filter(
    truth == "MID",
    residual == "REML",
    flip == "curve",
    error == "attached",
    arm %in% c("ncv_cl2", "reml_cl2")
  ) |>
  select(app, G, estimand, arm, mse) |>
  inner_join(
    CMP |> filter(source == "plasmode") |> distinct(cell, app, G),
    by = c("app", "G")
  ) |>
  transmute(source = "plasmode", cell, estimand, arm, mse)
# check: REML MSE identical in both sources
cmp_chk <- bind_rows(ncv_syn, ncv_pm) |>
  filter(arm == "reml_cl2") |>
  inner_join(
    CMP |>
      filter(arm == "reml_cl2") |>
      select(source, cell, estimand, mse_c = mse),
    by = c("source", "cell", "estimand")
  )
stopifnot(nrow(cmp_chk) > 0, max(abs(cmp_chk$mse / cmp_chk$mse_c - 1)) < 1e-6)
cmp_den <- bind_rows(
  read_summary("final", "relative-error") |> mutate(source = "synthetic"),
  read_summary("plasmode", "relative-error") |> mutate(source = "plasmode")
) |>
  distinct(source, cell, estimand, denom)
cmp_est_arms <- c(
  REML = "reml_cl2",
  NCV = "ncv_cl2",
  pcre = "pcre",
  `GLS (FPCA cov.)` = "gls_fpc",
  `GLS (raw cov.)` = "gls_raw"
)
CMP_EST <- bind_rows(
  CMP |>
    filter(arm %in% c("reml_cl2", "pcre", "gls_fpc", "gls_raw")) |>
    select(source, cell, estimand, arm, mse),
  bind_rows(ncv_syn, ncv_pm) |> filter(arm == "ncv_cl2")
) |>
  left_join(
    CMP |> distinct(source, cell, label, error, app, G),
    by = c("source", "cell")
  ) |>
  left_join(cmp_den, by = c("source", "cell", "estimand")) |>
  mutate(
    rel = sqrt(mse / denom),
    fit = factor(
      names(cmp_est_arms)[match(arm, cmp_est_arms)],
      names(cmp_est_arms)
    ),
    setting = case_when(
      source == "plasmode" ~ "plasmode",
      error == "iid" ~ "synthetic, independent",
      TRUE ~ "synthetic, dependent"
    )
  )
stopifnot(!anyNA(CMP_EST$rel))
cmp_ratio <- CMP_EST |>
  group_by(source, cell, estimand) |>
  mutate(
    vs_REML = mse / mse[fit == "REML"],
    vs_NCV = mse / mse[fit == "NCV"]
  ) |>
  ungroup()
cmp_er <- \(f, e, set, col = "vs_NCV")
  cmp_ratio[[col]][
    cmp_ratio$fit == f & cmp_ratio$estimand == e & cmp_ratio$setting == set
  ]
cmp_cv <- \(a, e, set) {
  d <- CMP |>
    mutate(
      setting = case_when(
        source == "plasmode" ~ "plasmode",
        error == "iid" ~ "synthetic, independent",
        TRUE ~ "synthetic, dependent"
      )
    )
  d$coverage[d$arm == a & d$estimand == e & d$setting == set]
}

d <- CMP |>
  mutate(
    setting = case_when(
      source == "plasmode" ~ "plasmode",
      error == "iid" ~ "synthetic, independent",
      TRUE ~ "synthetic, dependent"
    ),
    arm = factor(CMP_ARMS[arm], rev(CMP_ARMS)),
    estimand = label_estimand(estimand)
  )
fig <- ggplot(d, aes(coverage, arm)) +
  geom_vline(xintercept = LEVEL, linetype = 2, colour = "grey40") +
  geom_point(
    position = position_jitter(height = 0.15, width = 0, seed = 1),
    size = 1,
    alpha = 0.7
  ) +
  facet_grid(estimand ~ setting) +
  labs(x = "coverage (grid average)", y = NULL) +
  theme_doc
save_fig(fig, "comparators-coverage", 7, 5.5)
fig

d <- cmp_ratio |>
  filter(fit != "NCV") |>
  mutate(estimand = label_estimand(estimand))
fig <- ggplot(d, aes(vs_NCV, fit)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey40") +
  geom_point(
    position = position_jitter(height = 0.15, width = 0, seed = 1),
    size = 1,
    alpha = 0.7
  ) +
  scale_x_log10() +
  facet_grid(estimand ~ setting) +
  labs(x = "MSE / MSE of the NCV fit (log scale)", y = NULL) +
  theme_doc
save_fig(fig, "comparators-estimation", 7, 5)
fig

comp_cov <- CMP |>
  mutate(
    setting = case_when(
      source == "plasmode" ~ "plasmode",
      error == "iid" ~ "synthetic, independent",
      TRUE ~ "synthetic, dependent"
    )
  ) |>
  group_by(setting, estimand, cell) |>
  mutate(
    width_rel = width / width[arm == "reml_cl2"],
    score_rel = score / score[arm == "reml_cl2"]
  ) |>
  group_by(setting, estimand, arm) |>
  summarise(
    cells = n(),
    coverage_min = min(coverage),
    coverage = mean(coverage),
    q05 = mean(pointwise_q05),
    width_rel = exp(mean(log(width_rel))),
    score_rel = exp(mean(log(score_rel))),
    .groups = "drop"
  ) |>
  mutate(
    arm = factor(CMP_ARMS[arm], CMP_ARMS),
    estimand = label_estimand(estimand)
  ) |>
  select(
    setting,
    estimand,
    arm,
    cells,
    coverage,
    coverage_min,
    q05,
    width_rel,
    score_rel
  ) |>
  arrange(setting, estimand, arm)
save_tab(comp_cov, "comparators-coverage")
comp_est <- cmp_ratio |>
  group_by(setting, estimand, fit) |>
  summarise(
    cells = n(),
    rel_error = median(rel),
    mse_vs_REML = exp(mean(log(vs_REML))),
    vs_NCV_min = min(vs_NCV),
    vs_NCV_max = max(vs_NCV),
    mse_vs_NCV = exp(mean(log(vs_NCV))),
    .groups = "drop"
  ) |>
  mutate(estimand = label_estimand(estimand)) |>
  arrange(setting, estimand, fit)
save_tab(comp_est, "comparators-estimation")
comp_est_pm <- cmp_ratio |>
  filter(source == "plasmode", estimand == "beta") |>
  select(app, G, fit, rel) |>
  pivot_wider(names_from = fit, values_from = rel) |>
  mutate(app = label_app(app)) |>
  arrange(app, desc(G))
save_tab(comp_est_pm, "comparators-estimation-plasmode-beta")
comp_time <- CMPF |>
  group_by(arm) |>
  summarise(
    median_s = median(time),
    q90_s = unname(quantile(time, 0.9)),
    .groups = "drop"
  )
save_tab(comp_time, "comparators-time")
cmp_gm <- \(x) exp(mean(log(x)))

show(
  comp_cov |>
    mutate(
      setting = recode(
        setting,
        `synthetic, dependent` = "syn. dep.",
        `synthetic, independent` = "syn. iid"
      )
    ),
  digits = 2
)

show(
  comp_est |>
    mutate(
      setting = recode(
        setting,
        `synthetic, dependent` = "syn. dep.",
        `synthetic, independent` = "syn. iid"
      )
    ),
  digits = 2
)

show(comp_est_pm, digits = 2)

CA <- readr::read_csv(
  "comparators/summaries/cl2-arms-cells.csv",
  show_col_types = FALSE
) |>
  filter(window == "1-100") |>
  mutate(
    setting = case_when(
      source == "plasmode" ~ "plasmode",
      error == "iid" ~ "synthetic, independent",
      TRUE ~ "synthetic, dependent"
    )
  )
CA_ARMS <- c(
  reml_cl2 = "REML + CL2",
  ncv_cl2 = "NCV + CL2",
  gls_fpc = "GLS (FPCA cov.)",
  gls_fpc_cl2 = "GLS (FPCA cov.) + CL2",
  gls_raw_cl2 = "GLS (raw cov.) + CL2",
  pcre = "pcre",
  pcre_cl2 = "pcre + CL2"
)
ca_tab <- CA |>
  filter(arm %in% names(CA_ARMS)) |>
  group_by(setting, estimand, cell) |>
  mutate(
    width_rel = width / width[arm == "reml_cl2"],
    score_rel = score / score[arm == "reml_cl2"],
    mse_vs_NCV = mse / mse[arm == "ncv_cl2"]
  ) |>
  group_by(setting, estimand, arm) |>
  summarise(
    cov_min = min(coverage),
    coverage = mean(coverage),
    q05 = mean(pointwise_q05),
    width_rel = exp(mean(log(width_rel))),
    score_rel = exp(mean(log(score_rel))),
    mse_vs_NCV = exp(mean(log(mse_vs_NCV))),
    detection = mean(detection),
    .groups = "drop"
  ) |>
  mutate(
    arm = factor(CA_ARMS[arm], CA_ARMS),
    estimand = label_estimand(estimand)
  ) |>
  select(
    setting,
    estimand,
    arm,
    coverage,
    cov_min,
    q05,
    width_rel,
    score_rel,
    mse_vs_NCV,
    detection
  ) |>
  arrange(setting, estimand, arm)
save_tab(ca_tab, "comparators-cl2-arms")
cav <- \(a, set, col = "coverage", e = "beta")
  ca_tab[[col]][
    ca_tab$arm == CA_ARMS[[a]] &
      ca_tab$setting == set &
      ca_tab$estimand == label_estimand(e)
  ]
ca_pm_app <- CA |>
  filter(source == "plasmode", estimand == "beta", arm == "gls_fpc_cl2") |>
  group_by(app) |>
  summarise(v = mean(coverage), .groups = "drop")

show(
  ca_tab |>
    mutate(
      setting = recode(
        setting,
        `synthetic, dependent` = "syn. dep.",
        `synthetic, independent` = "syn. iid"
      )
    ),
  digits = 2
)
