ABL <- readr::read_csv(
  "ablation/summaries/ablation-cells.csv",
  show_col_types = FALSE
)
ABLC <- readr::read_csv(
  "ablation/summaries/ablation-contrasts.csv",
  show_col_types = FALSE
)
SAT <- readr::read_csv(
  "satterthwaite/summaries/satt-cells.csv",
  show_col_types = FALSE
)
SATC <- readr::read_csv(
  "satterthwaite/summaries/satt-contrasts.csv",
  show_col_types = FALSE
)
SATD <- readr::read_csv(
  "satterthwaite/summaries/satt-df.csv",
  show_col_types = FALSE
)
JACK <- readr::read_csv(
  "ncv-jack-recheck/jack-reps.csv",
  show_col_types = FALSE
)
ABL_ARMS <- c(
  model = "model-based (Vc)",
  cr1_b = "CR1",
  cr1_f = "CR1 (freq.)",
  cl2s_b = "CL2, I - H_gg approx.",
  cl2s_f = "CL2, I - H_gg approx. (freq.)",
  cl2e_b = "CL2 (full block)",
  cl2e_f = "CL2 (full block, freq.)",
  fix_model = "fixed lambda, model-based",
  fix_cl2e_b = "fixed lambda, CL2"
)
dep_label <- \(e) if_else(e == "iid", "independent", "dependent")

d <- ABL |>
  filter(crit == "z", estimand %in% c("mean", "beta")) |>
  mutate(
    arm = factor(ABL_ARMS[arm], rev(ABL_ARMS)),
    estimand = label_estimand(estimand),
    family = label_family(family),
    dependence = dep_label(error),
    G = paste0("G = ", G)
  )
fig <- ggplot(d, aes(coverage, arm, colour = dependence)) +
  geom_vline(xintercept = LEVEL, linetype = 2, colour = "grey40") +
  geom_pointrange(
    aes(xmin = coverage - 2 * coverage_se, xmax = coverage + 2 * coverage_se),
    position = position_dodge(width = 0.6),
    size = 0.15
  ) +
  facet_grid(estimand + G ~ family) +
  scale_colour_manual(
    values = c(dependent = "#D55E00", independent = "#0072B2"),
    name = NULL
  ) +
  labs(x = "coverage (grid average)", y = NULL) +
  theme_doc
save_fig(fig, "ablation-coverage", 7, 6)
fig

abl_wide <- ABL |>
  filter(estimand %in% c("mean", "beta", "alpha", "gamma"), crit == "z") |>
  mutate(arm = ABL_ARMS[arm]) |>
  select(family, error, G, estimand, arm, coverage) |>
  pivot_wider(names_from = arm, values_from = coverage) |>
  arrange(estimand, family, error, G)
save_tab(abl_wide, "ablation-coverage-wide")
abl_c <- ABLC |>
  mutate(dependence = dep_label(error)) |>
  group_by(contrast, estimand, dependence, G) |>
  summarise(
    diff_min = min(diff),
    diff_max = max(diff),
    diff = mean(diff),
    .groups = "drop"
  ) |>
  arrange(contrast, estimand, dependence, G)
save_tab(abl_c, "ablation-contrasts-summary")
abl_se <- ABL |>
  filter(crit == "z", estimand %in% c("mean", "beta")) |>
  group_by(arm, dependence = dep_label(error), estimand) |>
  summarise(se_sd = mean(se_sd_ratio), width = mean(width), .groups = "drop") |>
  mutate(arm = ABL_ARMS[arm])
save_tab(abl_se, "ablation-se-sd")
ac <- \(ct, dep = "dependent", e = c("mean", "beta"), g = c(40, 100)) {
  x <- ABLC$diff[
    ABLC$contrast == ct &
      dep_label(ABLC$error) == dep &
      ABLC$estimand %in% e &
      ABLC$G %in% g
  ]
  x
}
ab_bin_a <- \(a, err = "fpc", g = 40, e = "alpha")
  ABL$coverage[
    ABL$arm == a &
      ABL$crit == "z" &
      ABL$family == "binomial" &
      ABL$error == err &
      ABL$G == g &
      ABL$estimand == e
  ]

show(abl_c, digits = 4)

show(abl_se)

d <- SAT |>
  filter(crit %in% c("z", "t", "satt")) |>
  mutate(
    crit = factor(
      c(z = "z", t = "t(G−1)", satt = "Satterthwaite")[crit],
      c("z", "t(G−1)", "Satterthwaite")
    ),
    estimand = label_estimand(estimand),
    family = label_family(family),
    dependence = dep_label(error),
    G = paste0("G = ", G)
  )
fig <- ggplot(
  d,
  aes(coverage, interaction(error, G, sep = ", "), colour = crit)
) +
  geom_vline(xintercept = LEVEL, linetype = 2, colour = "grey40") +
  geom_pointrange(
    aes(xmin = coverage - 2 * coverage_se, xmax = coverage + 2 * coverage_se),
    position = position_dodge(width = 0.6),
    size = 0.15
  ) +
  facet_grid(estimand ~ family) +
  scale_colour_manual(
    values = c(
      "z" = "#999999",
      "t(G−1)" = "#56B4E9",
      "Satterthwaite" = "#0072B2"
    ),
    name = NULL
  ) +
  labs(x = "coverage (grid average)", y = "error process, G") +
  theme_doc
save_fig(fig, "satterthwaite-coverage", 7, 6)
fig

satt_c <- SATC |>
  mutate(dependence = dep_label(error)) |>
  group_by(contrast, estimand, dependence, G) |>
  summarise(
    diff_min = min(diff),
    diff_max = max(diff),
    diff = mean(diff),
    width_ratio = mean(width_ratio),
    .groups = "drop"
  )
save_tab(satt_c, "satterthwaite-contrasts-summary")
satt_df <- SATD |>
  group_by(estimand, G) |>
  summarise(
    df_q05 = mean(df_q05),
    df_median = mean(df_median),
    df_q95 = mean(df_q95),
    df_rel = mean(df_median_rel),
    df_b2 = mean(df_b2_median),
    .groups = "drop"
  )
save_tab(satt_df, "satterthwaite-df")
sc <- \(ct, e = c("mean", "beta", "gamma"), g = c(40, 100))
  SATC$diff[SATC$contrast == ct & SATC$estimand %in% e & SATC$G %in% g]
sworst <- SAT |> filter(crit == "satt") |> slice_min(coverage, n = 1)

show(satt_c, digits = 4)

show(satt_df, digits = 2)

jack_tab <- JACK |>
  group_by(family, estimand, arm) |>
  summarise(
    coverage = mean(coverage),
    coverage_se = sd(coverage) / sqrt(n()),
    se_ratio_to_bayes = mean(se_ratio_to_bayes),
    width = mean(width),
    reps = n(),
    .groups = "drop"
  ) |>
  mutate(
    arm = c(
      bayes_Vp = "Bayesian Vp",
      jack_Vp = "jackknife Vp (jackknife = TRUE)",
      jack_Vc = "jackknife Vc (patched mgcv)",
      cl2 = "CL2"
    )[arm]
  )
save_tab(jack_tab, "ncv-jackknife-recheck")

show(jack_tab |> mutate(family = label_family(family)))
