# (helpers: paper/_common.R)
options(knitr.kable.NA = "", scipen = 100)

# (numbers: paper/_common.R)

# relative error and interval informativeness (details: @sec-relerr and the plasmode subsection)
rel_s <- read_summary("final", "relative-error") |>
  left_join(S |> distinct(cell, error), by = "cell") |>
  filter(estimand == "beta", error != "iid")
rs_med <- \(a, col) median(rel_s[[col]][rel_s$arm == a])
rel_p <- read_summary("plasmode", "relative-error") |>
  filter(estimand == "beta", flip == "curve", G > 40, error == "attached") |>
  group_by(app, arm) |>
  summarise(
    rel = median(rel_rmse_median),
    det = median(detect),
    .groups = "drop"
  )
rp <- \(a, col) rel_p[[col]][rel_p$arm == a]

blocks <- DESIGN |>
  mutate(
    fitted = if_else(
      family == fit_family,
      FAMILY_LABELS[family],
      paste0(
        FAMILY_LABELS[family],
        " (fitted as ",
        FAMILY_LABELS[fit_family],
        ")"
      )
    )
  ) |>
  group_by(block) |>
  summarise(
    cells = n(),
    families = paste(unique(fitted), collapse = ", "),
    errors = paste(unique(SHORT_ERROR[error]), collapse = ", "),
    G = paste(sort(unique(G)), collapse = ", "),
    signal = paste(unique(signal), collapse = ", "),
    truth = paste(unique(truth), collapse = ", "),
    D = paste(unique(D), collapse = ", "),
    model = paste(
      unique(c(M1 = "ff(X)", M2 = "f(x,t)")[model]),
      collapse = ", "
    ),
    covariate = paste(unique(covariate), collapse = ", "),
    basis = paste(unique(basis), collapse = ", "),
    .groups = "drop"
  ) |>
  arrange(match(
    block,
    c(
      "core",
      "warp",
      "dense_grid",
      "rough_truth",
      "term_type",
      "families",
      "ar1_home",
      "oscillating",
      "warp_ar1",
      "heteroskedastic",
      "lowrank_covariate",
      "basis_size"
    )
  ))
save_tab(blocks, "synthetic-design-blocks")
fit_summary <- FITS |>
  group_by(fit_family, fit) |>
  summarise(
    converged = weighted.mean(converged, n),
    median_time_s = median(median_time),
    .groups = "drop"
  )
save_tab(fit_summary, "synthetic-fits")

show(
  blocks |>
    select(
      block,
      cells,
      families,
      errors,
      G,
      signal,
      truth,
      D,
      model,
      covariate,
      basis
    )
)

d <- S |>
  filter(
    block == "core",
    signal == "mid",
    estimand %in% c("mean", "beta"),
    arm %in% MAIN_ARMS
  ) |>
  mutate(
    arm = label_arm(arm),
    estimand = label_estimand(estimand),
    error = label_error(error),
    family = label_family(family),
    G = paste0("G = ", G)
  )
fig <- coverage_plot(d, "error", vars(G, estimand), vars(family))
save_fig(fig, "coverage-core", 6.5, 5.5)
fig

q1 <- Q1 |>
  mutate(value = ci(estimate, lower, upper)) |>
  select(family, estimand, name, value) |>
  pivot_wider(names_from = name, values_from = value) |>
  rename(
    `CE model-based (REML)` = ce_reml_model,
    `CE REML + CL2` = ce_reml_cl2,
    `CE model-based (NCV)` = ce_ncv_model,
    `CE NCV + CL2` = ce_ncv_cl2,
    `gain REML` = ce_gain_reml,
    `gain NCV` = ce_gain_ncv,
    `iid cost REML` = iid_cost_reml,
    `iid cost NCV` = iid_cost_ncv
  ) |>
  pretty()
save_tab(Q1, "synthetic-question1")
by_fam <- core_dep |>
  group_by(family, arm) |>
  summarise(lo = min(coverage_mean), hi = max(coverage_mean), .groups = "drop")
bf <- \(fam, a)
  c(
    by_fam$lo[by_fam$family == fam & by_fam$arm == a],
    by_fam$hi[by_fam$family == fam & by_fam$arm == a]
  )

show(q1)

d <- S |>
  filter(
    block == "core",
    signal == "mid",
    estimand %in% c("mean", "beta"),
    arm %in% MAIN_ARMS
  ) |>
  mutate(
    arm = label_arm(arm),
    estimand = label_estimand(estimand),
    error = label_error(error),
    family = label_family(family),
    G = paste0("G = ", G)
  )
fig <- coverage_plot(
  d,
  "error",
  vars(G, estimand),
  vars(family),
  x = "pointwise_q05",
  se = "pointwise_q05_se",
  xlab = "5% quantile of pointwise coverage"
)
save_fig(fig, "lower-tail-core", 6.5, 5.5)
fig
tail_tab <- S |>
  filter(
    block == "core",
    estimand %in% c("mean", "beta"),
    arm %in% MAIN_ARMS
  ) |>
  group_by(
    dependence = if_else(error == "iid", "independent", "dependent"),
    G,
    estimand,
    arm
  ) |>
  summarise(
    q05_min = min(pointwise_q05),
    q05_median = median(pointwise_q05),
    coverage_min = min(coverage_mean),
    .groups = "drop"
  ) |>
  mutate(arm = label_arm(arm), estimand = label_estimand(estimand))
save_tab(tail_tab, "synthetic-lower-tail-core")
tq <- \(a, e, g, dep = "dependent")
  tail_tab$q05_min[
    tail_tab$arm == ARM_LABELS[a] &
      tail_tab$estimand == ESTIMAND_LABELS[e] &
      tail_tab$G == g &
      tail_tab$dependence == dep
  ]

zs_avail <- !is.null(ZS)
# the mean's stored records hold link-scale estimates and SEs next to the response-scale
# truth, so z is defined for the mean only in Gaussian cells
zs_tab <- if (zs_avail)
  ZS |>
    filter(fit == "REML", !(estimand == "mean" & fit_family != "gaussian")) |>
    group_by(block, family, error, G, basis, covariate, estimand) |>
    summarise(
      across(
        c(
          z_mean_cl2,
          z_sd_cl2,
          share_big_cl2,
          z_mean_model,
          z_sd_model,
          share_big_model
        ),
        mean
      ),
      .groups = "drop"
    ) |>
    arrange(block, family, error, G, basis, estimand) else NULL
if (zs_avail) save_tab(zs_tab, "synthetic-zstats-reml")
zs_blocks <- if (zs_avail) sort(unique(ZS$block)) else character()
zs_dep <- if (zs_avail) zs_tab |> filter(error != "iid") else NULL
zs_iid <- if (zs_avail) zs_tab |> filter(error == "iid") else NULL

zc <- if (zs_avail) zs_tab |> filter(block == "core") else NULL
zcv <- \(col, dep = TRUE, g = 100, e = "beta")
  zc[[col]][(zc$error != "iid") == dep & zc$G == g & zc$estimand == e]
zx <- if (zs_avail)
  zs_tab |>
    filter(
      block %in%
        c(
          "warp",
          "rough_truth",
          "dense_grid",
          "oscillating",
          "families",
          "heteroskedastic",
          "lowrank_covariate",
          "term_type"
        ),
      error != "iid"
    ) else NULL
zb <- if (zs_avail)
  zs_tab |>
    filter(block == "basis_size", error != "iid", estimand == "beta") else NULL
za <- if (zs_avail)
  ZS |> filter(fit == "REML", estimand == "alpha", error != "iid") else NULL
zav <- \(fam, blk = "core", g = 100, col = "z_mean_cl2", b = "default")
  za[[col]][za$family == fam & za$block == blk & za$G == g & za$basis == b]
zg <- if (zs_avail)
  ZS |>
    filter(
      fit == "REML",
      estimand == "gamma",
      error != "iid",
      block == "core"
    ) else NULL

d <- S |>
  filter(
    (block == "core" & G == 100 & signal == "mid") | block == "dense_grid",
    estimand %in% c("beta", "mean"),
    arm %in% c("reml_model", "reml_cl2", "ncv_cl2_bias")
  ) |>
  mutate(
    arm = label_arm(arm),
    estimand = label_estimand(estimand),
    error = label_error_short(error),
    family = label_family(family),
    D = factor(D)
  )
fig <- ggplot(d, aes(D, coverage_mean, colour = arm, group = arm)) +
  geom_hline(yintercept = LEVEL, linetype = 2, colour = "grey40") +
  geom_line() +
  geom_pointrange(
    aes(
      ymin = coverage_mean - 2 * coverage_se,
      ymax = coverage_mean + 2 * coverage_se
    ),
    size = 0.2
  ) +
  facet_grid(estimand ~ family + error) +
  scale_colour_manual(values = ARM_COLOURS, name = NULL) +
  labs(x = "grid points per curve D", y = "coverage") +
  theme_doc +
  theme(strip.text = element_text(size = 7))
save_fig(fig, "coverage-grid", 6.5, 3.5)
fig
grid_tab <- MC |>
  filter(
    block == "dense_grid",
    estimand == "beta",
    value %in%
      c(
        "coverage_reml_model",
        "coverage_reml_cl2",
        "coverage_ncv_cl2",
        "coverage_ncv_cl2_bias",
        "log_mse_ratio_ncv_reml"
      )
  ) |>
  mutate(diff = ci(estimate, lower, upper)) |>
  select(family, error, value, diff) |>
  pivot_wider(names_from = value, values_from = diff) |>
  rename(
    `model-based` = coverage_reml_model,
    `REML + CL2` = coverage_reml_cl2,
    `NCV + CL2` = coverage_ncv_cl2,
    `bias-aware` = coverage_ncv_cl2_bias,
    `log MSE ratio` = log_mse_ratio_ncv_reml
  ) |>
  pretty()
save_tab(MC |> filter(block == "dense_grid"), "synthetic-grid-contrasts")
gd <- \(v, err = c("ou", "fpc"))
  dense$estimate[dense$value == v & dense$error %in% err]

show(grid_tab)

g_tab <- S |>
  filter(
    block == "core",
    estimand %in% c("mean", "beta"),
    arm %in% c(MAIN_ARMS, "ncv_cl2f_bias")
  ) |>
  select(family, error, signal, estimand, arm, G, coverage_mean) |>
  pivot_wider(
    names_from = G,
    values_from = coverage_mean,
    names_prefix = "G"
  ) |>
  mutate(diff = G100 - G40)
g_sum <- g_tab |>
  group_by(
    dependence = if_else(error == "iid", "independent", "dependent"),
    family,
    estimand,
    arm
  ) |>
  summarise(
    `G = 40` = mean(G40),
    `G = 100` = mean(G100),
    `loss at G = 40` = mean(diff),
    .groups = "drop"
  ) |>
  mutate(
    arm = label_arm(arm),
    estimand = label_estimand(estimand),
    family = label_family(family)
  ) |>
  arrange(dependence, family, estimand, arm)
save_tab(g_sum, "synthetic-G40-vs-100")
gl <- \(a, dep = "dependent")
  g_sum$`loss at G = 40`[g_sum$arm == ARM_LABELS[a] & g_sum$dependence == dep]
g40_bin <- g_sum |> filter(dependence == "dependent", family == "binary")
g40_bin_ncv_beta <- g40_bin[["G = 40"]][
  g40_bin$arm == ARM_LABELS["ncv_cl2"] &
    g40_bin$estimand == ESTIMAND_LABELS["beta"]
]

d <- S |>
  filter(
    block %in% c("core", "rough_truth"),
    G == 100,
    error %in% c("iid", "fpc"),
    estimand == "beta",
    arm %in% c("reml_cl2", "ncv_cl2", "ncv_cl2_bias", "ncv_cl2f_bias")
  ) |>
  mutate(
    arm = label_arm(arm),
    family = label_family(family),
    error = label_error(error),
    signal = label_signal(signal),
    truth = paste(if_else(truth == "wiggly", "rough", "smooth"), "truth")
  )
fig <- coverage_plot(d, "signal", vars(truth, error), vars(family))
save_fig(fig, "coverage-signal-truth", 6.5, 5)
fig
st <- S |>
  filter(
    block %in% c("core", "rough_truth"),
    G == 100,
    error %in% c("iid", "fpc"),
    estimand == "beta"
  )
stv <- \(a, tr, err = "fpc")
  st$coverage_mean[st$arm == a & st$truth == tr & st$error == err]
st_tab <- cov_wide(
  st$cell,
  "beta",
  c("reml_cl2", "ncv_cl2", "ncv_cl2_bias", "ncv_cl2f_bias")
) |>
  select(-D, -basis, -G) |>
  pretty()
save_tab(st_tab, "synthetic-signal-truth-beta")

est_tab <- S |>
  filter(block == "core", G == 100, signal == "mid", arm %in% MAIN_ARMS) |>
  group_by(
    dependence = if_else(error == "iid", "independent", "dependent"),
    family,
    estimand,
    arm
  ) |>
  summarise(
    coverage = mean(coverage_mean),
    q05 = mean(pointwise_q05),
    .groups = "drop"
  ) |>
  mutate(
    arm = label_arm(arm),
    estimand = label_estimand(estimand),
    family = label_family(family)
  )
est_wide <- est_tab |>
  select(-q05) |>
  pivot_wider(names_from = arm, values_from = coverage) |>
  arrange(dependence, family, estimand)
save_tab(est_tab, "synthetic-coverage-by-estimand")
ev <- \(a, e, dep = "dependent")
  est_tab$coverage[
    est_tab$arm == ARM_LABELS[a] &
      est_tab$estimand == ESTIMAND_LABELS[e] &
      est_tab$dependence == dep
  ]
m2 <- MC |> filter(block == "term_type", estimand == "main")
m2_tab <- m2 |>
  filter(
    value %in%
      c(
        "coverage_reml_model",
        "coverage_reml_cl2",
        "coverage_ncv_cl2_bias",
        "log_mse_ratio_ncv_reml"
      )
  ) |>
  mutate(diff = ci(estimate, lower, upper)) |>
  select(family, error, value, diff) |>
  pivot_wider(names_from = value, values_from = diff) |>
  rename(
    `model-based` = coverage_reml_model,
    `REML + CL2` = coverage_reml_cl2,
    `bias-aware` = coverage_ncv_cl2_bias,
    `log MSE ratio` = log_mse_ratio_ncv_reml
  ) |>
  pretty()
save_tab(MC |> filter(block == "term_type"), "synthetic-term-type-contrasts")
m2v <- \(v, err = c("ou", "fpc")) m2$estimate[m2$value == v & m2$error %in% err]

show(est_wide)

show(m2_tab)

het_tab <- cov_wide(
  DESIGN$cell[
    DESIGN$block == "heteroskedastic" |
      (DESIGN$block == "core" &
        DESIGN$error == "fpc" &
        DESIGN$family == "gaussian" &
        DESIGN$signal == "mid")
  ],
  c("mean", "beta", "gamma"),
  c(MAIN_ARMS, "ar1_model")
) |>
  select(-D, -basis, -signal, -truth, -family)
hv <- \(a, e, err) het_tab[[a]][het_tab$estimand == e & het_tab$error %in% err]
het_show <- pretty(het_tab)
save_tab(het_show, "synthetic-heteroskedastic-coverage")
# misregistration: the AR(1)-on-warp cell (same DGP as the Gaussian rough-truth warp cell at G = 100) supplies the AR(1) column of that row
warp_raw <- cov_wide(
  DESIGN$cell[DESIGN$block %in% c("warp", "warp_ar1")],
  c("mean", "beta"),
  c(MAIN_ARMS, "ncv_cl2f_bias", "ar1_model")
)
warp_tab <- warp_raw |>
  group_by(family, G, signal, truth, estimand) |>
  summarise(
    across(
      c(all_of(MAIN_ARMS), ncv_cl2f_bias, ar1_model),
      \(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
    ),
    .groups = "drop"
  )
wv <- \(a, fam, e = "mean")
  warp_tab[[a]][warp_tab$family == fam & warp_tab$estimand == e]
warp_show <- pretty(warp_tab)
save_tab(warp_show, "synthetic-misregistration-coverage")
osc_tab <- cov_wide(
  DESIGN$cell[DESIGN$block == "oscillating"],
  c("mean", "beta"),
  c(MAIN_ARMS, "ar1_model")
) |>
  select(-D, -basis, -signal, -truth, -G)
osc_show <- pretty(osc_tab)
save_tab(osc_show, "synthetic-oscillating-coverage")

show(warp_show)

fam_cov <- cov_wide(
  DESIGN$cell[DESIGN$block == "families"],
  c("mean", "beta"),
  MAIN_ARMS
) |>
  left_join(distinct(DESIGN, cell, fit_family), by = "cell") |>
  select(family, fit_family, error, estimand, all_of(MAIN_ARMS))
fam_mse <- MSE_RATIO |>
  filter(block == "families", estimand %in% c("mean", "beta")) |>
  mutate(`MSE ratio NCV/REML` = ci(ratio, ratio_lower, ratio_upper, 2)) |>
  select(family, fit_family, error, estimand, `MSE ratio NCV/REML`)
fam_tab <- left_join(
  fam_cov,
  fam_mse,
  by = c("family", "fit_family", "error", "estimand")
)
fv <- \(a, err = "fpc")
  fam_tab[[a]][fam_tab$error == err & fam_tab$family == fam_tab$fit_family]
fam_show <- pretty(fam_tab)
save_tab(fam_show, "synthetic-extended-families")

show(fam_show)

basis_levels <- c("default", "large", "xlarge")
bcells <- DESIGN |>
  filter(
    block == "basis_size" |
      (block == "core" &
        G == 100 &
        signal == "mid" &
        error %in% c("iid", "fpc") &
        family %in% c("gaussian", "poisson", "binomial"))
  )
edf_b <- EDF |>
  filter(
    cell %in% bcells$cell,
    term == "te(X.smat,X.tmat):L.X",
    fit %in% c("REML", "NCV")
  ) |>
  mutate(
    basis = factor(basis, basis_levels),
    family = label_family(family),
    error = label_error(error)
  )
p1 <- ggplot(edf_b, aes(basis, edf_mean, colour = fit, group = fit)) +
  geom_line() +
  geom_point() +
  facet_grid(error ~ family) +
  scale_colour_manual(
    values = c(REML = "#0072B2", NCV = "#D55E00"),
    name = NULL
  ) +
  labs(x = NULL, y = "EDF of the ff term") +
  theme_doc
mse_b <- MSE_RATIO |>
  filter(cell %in% bcells$cell, estimand %in% c("beta", "mean")) |>
  mutate(
    basis = factor(basis, basis_levels),
    family = label_family(family),
    error = label_error(error),
    estimand = label_estimand(estimand)
  )
p2 <- ggplot(mse_b, aes(basis, ratio, colour = estimand, group = estimand)) +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey40") +
  geom_line() +
  geom_pointrange(aes(ymin = ratio_lower, ymax = ratio_upper), size = 0.2) +
  facet_grid(error ~ family) +
  scale_y_log10(breaks = c(0.01, 0.03, 0.1, 0.3, 1)) +
  scale_colour_manual(values = c("#000000", "#56B4E9"), name = NULL) +
  labs(x = NULL, y = "MSE NCV / MSE REML") +
  theme_doc
cov_b <- S |>
  filter(cell %in% bcells$cell, estimand == "beta", arm %in% MAIN_ARMS) |>
  mutate(
    basis = factor(basis, basis_levels),
    family = label_family(family),
    error = label_error(error),
    arm = label_arm(arm)
  )
p3 <- ggplot(cov_b, aes(basis, coverage_mean, colour = arm, group = arm)) +
  geom_hline(yintercept = LEVEL, linetype = 2, colour = "grey40") +
  geom_line() +
  geom_pointrange(
    aes(
      ymin = coverage_mean - 2 * coverage_se,
      ymax = coverage_mean + 2 * coverage_se
    ),
    size = 0.2
  ) +
  facet_grid(error ~ family) +
  scale_colour_manual(values = ARM_COLOURS, name = NULL) +
  guides(colour = guide_legend(nrow = 2)) +
  labs(x = NULL, y = "coverage of beta") +
  theme_doc
fig <- (p1 | p2 | p3) &
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 30, hjust = 1),
    strip.text = element_text(size = 6)
  )
save_fig(fig, "basis-size", 9, 4.5)
fig
basis_tab <- cov_b |>
  select(family, error, basis, arm, coverage_mean, score_mean, width_mean) |>
  pivot_wider(
    names_from = arm,
    values_from = c(coverage_mean, score_mean, width_mean)
  ) |>
  left_join(
    mse_b |>
      filter(estimand == ESTIMAND_LABELS["beta"]) |>
      select(family, error, basis, mse_ratio_beta = ratio),
    by = c("family", "error", "basis")
  ) |>
  arrange(family, error, basis)
names(basis_tab) <- sub(
  "coverage_mean_",
  "cov: ",
  sub("score_mean_", "IS: ", sub("width_mean_", "width: ", names(basis_tab)))
)
save_tab(basis_tab, "synthetic-basis-size")
ed <- \(fit, b, err = "fpc")
  edf_b$edf_mean[
    edf_b$fit == fit & edf_b$basis == b & edf_b$error == ERROR_LABELS[err]
  ]
bs <- \(a, b, err = "fpc")
  cov_b$coverage_mean[
    cov_b$arm == ARM_LABELS[a] &
      cov_b$basis == b &
      cov_b$error == ERROR_LABELS[err]
  ]
bis <- MC |>
  filter(
    block == "basis_size",
    estimand == "beta",
    value == "log_score_ratio_proposal_fallback"
  )
wr_b <- width_ratio |>
  filter(cell %in% bcells$cell, estimand == "beta", error == "fpc")

show(pretty(basis_tab), digits = 2)

d <- MSE_RATIO |>
  filter(
    block %in% c("core", "rough_truth", "term_type", "lowrank_covariate"),
    G == 100,
    estimand %in% c("beta", "f", "mean", "gamma")
  ) |>
  mutate(
    estimand = factor(
      recode(
        estimand,
        mean = "E(Y | X)",
        beta = "surface: beta(s,t) or f(x,t)",
        f = "surface: beta(s,t) or f(x,t)",
        gamma = "gamma(t)"
      ),
      c("E(Y | X)", "surface: beta(s,t) or f(x,t)", "gamma(t)")
    ),
    setting = factor(
      recode(
        block,
        core = "smooth truth",
        rough_truth = "rough truth",
        term_type = "f(x,t) model",
        lowrank_covariate = "low-rank X"
      ),
      c("smooth truth", "rough truth", "f(x,t) model", "low-rank X")
    ),
    error = label_error(error),
    signal = label_signal(signal),
    family = label_family(family)
  )
fig <- ggplot(d, aes(ratio, error, colour = family, shape = signal)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey40") +
  geom_pointrange(
    aes(xmin = ratio_lower, xmax = ratio_upper),
    position = position_dodge(width = 0.7),
    size = 0.2
  ) +
  facet_grid(setting ~ estimand, scales = "free_y", space = "free_y") +
  scale_shape_manual(
    values = c(low = 1, mid = 16, high = 17),
    name = "signal"
  ) +
  scale_x_log10(breaks = c(0.05, 0.1, 0.2, 0.5, 1, 2)) +
  scale_colour_manual(
    values = c(Gaussian = "#000000", Poisson = "#56B4E9", binary = "#CC79A7"),
    name = NULL
  ) +
  labs(x = "MSE of NCV / MSE of REML (log scale; G = 100)", y = NULL) +
  theme_doc
save_fig(fig, "mse-ratio", 6.5, 6)
fig
prof_err <- PROF |>
  filter(value == "log_mse_ratio_ncv_reml", factor == "error") |>
  mutate(ratio = ci(exp(estimate), exp(lower), exp(upper), 2)) |>
  select(family, estimand, error = name, ratio) |>
  pivot_wider(names_from = error, values_from = ratio) |>
  pretty()
save_tab(prof_err, "synthetic-mse-ratio-by-error")
prof_G <- PROF |>
  filter(value == "log_mse_ratio_ncv_reml", factor == "G") |>
  mutate(ratio = exp(estimate)) |>
  select(family, estimand, G = name, ratio) |>
  pivot_wider(names_from = G, values_from = ratio, names_prefix = "G")
pe <- \(e, err)
  exp(PROF$estimate[
    PROF$value == "log_mse_ratio_ncv_reml" &
      PROF$factor == "error" &
      PROF$name %in% err &
      PROF$estimand == e
  ])
rough <- MSE_RATIO |> filter(block == "rough_truth", error == "fpc")
lowr <- MSE_RATIO |> filter(block == "lowrank_covariate")

show(prof_err)

bv <- S |>
  filter(
    (block == "core" | block == "rough_truth") & G == 100 & signal == "mid",
    estimand %in% c("beta", "mean"),
    arm %in% c("reml_cl2", "ncv_cl2")
  ) |>
  transmute(
    family = label_family(family),
    error = label_error_short(error),
    truth = if_else(truth == "wiggly", "rough", "smooth"),
    estimand = label_estimand(estimand),
    fit = if_else(arm == "reml_cl2", "REML", "NCV"),
    `squared bias` = bias2,
    variance = variance
  ) |>
  mutate(
    x = factor(
      paste(fit, truth),
      c("REML smooth", "NCV smooth", "REML rough", "NCV rough")
    )
  ) |>
  pivot_longer(
    c(`squared bias`, variance),
    names_to = "component",
    values_to = "value"
  )
fig <- ggplot(bv, aes(x, value, fill = component)) +
  geom_col(position = "fill") +
  facet_grid(estimand ~ family + error) +
  scale_fill_manual(
    values = c(`squared bias` = "#D55E00", variance = "#0072B2"),
    name = NULL
  ) +
  labs(x = "fit and truth", y = "share of the integrated squared error") +
  theme_doc +
  theme(
    strip.text = element_text(size = 6),
    axis.text.x = element_text(angle = 60, hjust = 1, size = 6)
  )
save_fig(fig, "bias-variance", 7.5, 4.5)
fig
bv_tab <- S |>
  filter(
    block %in% c("core", "rough_truth"),
    estimand %in% c("beta", "mean"),
    arm %in% c("reml_cl2", "ncv_cl2")
  ) |>
  transmute(
    family,
    error,
    G,
    signal,
    truth,
    estimand,
    fit = if_else(arm == "reml_cl2", "REML", "NCV"),
    mse,
    mse_se,
    bias2,
    bias2_se,
    variance,
    bias_share = pmax(bias2, 0) / mse
  ) |>
  arrange(family, error, G, signal, truth, estimand, fit)
save_tab(bv_tab, "synthetic-bias-variance")
bshare <- \(fit_sel, tr, e = "beta", err = c("ou", "fpc"))
  bv_tab$bias_share[
    bv_tab$fit == fit_sel &
      bv_tab$truth == tr &
      bv_tab$estimand == e &
      bv_tab$error %in% err &
      bv_tab$G == 100
  ]

edf_core <- EDF |>
  filter(
    block == "core",
    G == 100,
    signal == "mid",
    fit %in% c("REML", "NCV")
  ) |>
  mutate(
    term = recode(
      term,
      `s(yindex.vec)` = "alpha(t)",
      `te(X.smat,X.tmat):L.X` = "beta(s,t)",
      `s(yindex.vec):z` = "gamma(t)"
    )
  ) |>
  select(family, error, term, fit, edf_mean) |>
  pivot_wider(names_from = fit, values_from = edf_mean) |>
  mutate(ratio = REML / NCV) |>
  arrange(term, family, error)
save_tab(
  EDF |>
    select(
      cell,
      block,
      family,
      error,
      G,
      signal,
      truth,
      basis,
      fit,
      term,
      edf_mean,
      edf_se,
      mean_log_sp
    ),
  "synthetic-term-edf"
)
er <- \(term_sel, err)
  edf_core$ratio[edf_core$term == term_sel & edf_core$error %in% err]

show(pretty(edf_core))

pw <- W |>
  filter(
    contrast == "pointwise ncv vs reml: squared error",
    estimand %in% c("beta", "mean")
  ) |>
  transmute(
    error,
    estimand,
    `MSE pointwise NCV / REML` = ci(ratio, ratio_lower, ratio_upper, 2)
  ) |>
  pretty()
save_tab(pw, "synthetic-pointwise-ncv")
pwv <- \(e, err)
  W$ratio[
    W$contrast == "pointwise ncv vs reml: squared error" &
      W$estimand == e &
      W$error %in% err
  ]
lowrank_tab <- MSE_RATIO |>
  filter(
    block == "lowrank_covariate" |
      (block == "core" &
        G == 100 &
        signal == "mid" &
        error %in% c("iid", "fpc") &
        family %in% c("gaussian", "poisson", "binomial"))
  ) |>
  mutate(
    covariate = if_else(block == "core", "rich", "low-rank"),
    ratio_ci = ci(ratio, ratio_lower, ratio_upper, 2)
  ) |>
  select(family, error, estimand, covariate, ratio_ci) |>
  pivot_wider(names_from = covariate, values_from = ratio_ci) |>
  arrange(family, error, estimand) |>
  pretty()
save_tab(lowrank_tab, "synthetic-lowrank-mse")

show(pw)

rv_syn <- RV |>
  filter(data == "synthetic") |>
  select(
    family = group,
    estimand,
    criterion,
    ce_proposal,
    ce_fallback,
    score_ratio,
    decision
  ) |>
  arrange(family, estimand, criterion) |>
  pretty()
save_tab(rv_syn, "synthetic-rule-variants")
save_tab(REC, "synthetic-recommendation")
save_tab(RECSIG, "synthetic-recommendation-by-signal")
rsym <- rv_syn |> filter(criterion == "symmetric")
runder <- rv_syn |> filter(criterion == "under only")

show(rv_syn)

mc_rough <- HMC_S |> filter(truth == "wiggly")
mcv <- \(d, m, col) d[[col]][d$method == m]
fbias <- S |>
  filter(
    block %in% c("core", "rough_truth"),
    G == 100,
    estimand == "beta",
    arm == "ncv_cl2f_bias"
  )
# hybrid in the synthetic study: dependent core cells at G = 100 (from the per-cell hybrid summary)
hs_cells <- HSUM_S |>
  filter(
    block == "core",
    G == 100,
    error != "iid",
    estimand %in% c("beta", "mean"),
    arm %in% c("reml_cl2", "ncv_cl2_bias", "hybrid")
  ) |>
  select(
    cell,
    family,
    error,
    signal,
    estimand,
    arm,
    coverage_mean,
    pointwise_q05,
    width_mean,
    score_mean
  )
hs_wide <- hs_cells |>
  pivot_wider(
    names_from = arm,
    values_from = c(coverage_mean, pointwise_q05, width_mean, score_mean)
  )
hsv <- \(a, col, e = "beta")
  hs_wide[[paste0(col, "_", a)]][hs_wide$estimand == e]
hs_rough <- HSUM_S |>
  filter(
    block == "rough_truth",
    error == "fpc",
    estimand == "beta",
    arm %in% c("reml_cl2", "ncv_cl2_bias", "hybrid")
  ) |>
  select(cell, arm, coverage_mean, pointwise_q05, width_mean, score_mean) |>
  pivot_wider(
    names_from = arm,
    values_from = c(coverage_mean, pointwise_q05, width_mean, score_mean)
  )
hyb_tab <- bind_rows(
  hs_wide |> mutate(setting = "dependent core, G = 100"),
  HSUM_S |>
    filter(
      block == "rough_truth",
      error == "fpc",
      estimand %in% c("beta", "mean"),
      arm %in% c("reml_cl2", "ncv_cl2_bias", "hybrid")
    ) |>
    select(
      cell,
      family,
      error,
      signal,
      estimand,
      arm,
      coverage_mean,
      pointwise_q05,
      width_mean,
      score_mean
    ) |>
    pivot_wider(
      names_from = arm,
      values_from = c(coverage_mean, pointwise_q05, width_mean, score_mean)
    ) |>
    mutate(setting = "rough truth, smooth errors")
) |>
  group_by(setting, estimand) |>
  summarise(
    n_cells = n(),
    `cov REML+CL2` = mean(coverage_mean_reml_cl2),
    `cov bias-aware` = mean(coverage_mean_ncv_cl2_bias),
    `cov hybrid` = mean(coverage_mean_hybrid),
    `q05 REML+CL2` = mean(pointwise_q05_reml_cl2),
    `q05 bias-aware` = mean(pointwise_q05_ncv_cl2_bias),
    `q05 hybrid` = mean(pointwise_q05_hybrid),
    `width bias-aware/REML+CL2` = exp(mean(log(
      width_mean_ncv_cl2_bias / width_mean_reml_cl2
    ))),
    `width hybrid/REML+CL2` = exp(mean(log(
      width_mean_hybrid / width_mean_reml_cl2
    ))),
    `IS bias-aware/REML+CL2` = exp(mean(log(
      score_mean_ncv_cl2_bias / score_mean_reml_cl2
    ))),
    `IS hybrid/REML+CL2` = exp(mean(log(
      score_mean_hybrid / score_mean_reml_cl2
    ))),
    .groups = "drop"
  ) |>
  mutate(estimand = ESTIMAND_LABELS[estimand])
save_tab(hyb_tab, "synthetic-hybrid-bias-aware")
ht <- \(col, setting_sel = "dependent core, G = 100", e = "beta(s,t)")
  hyb_tab[[col]][hyb_tab$setting == setting_sel & hyb_tab$estimand == e]

show(hyb_tab, digits = 3)

maps_s <- read_maps("synthetic")
map_cells <- maps_s |>
  filter(
    family == "gaussian",
    signal == "mid",
    error == "fpc",
    basis == "default"
  ) |>
  distinct(cell, truth)
mp <- maps_s |>
  filter(cell %in% map_cells$cell) |>
  mutate(
    method = label_arm(method),
    truth = paste(if_else(truth == "wiggly", "rough", "smooth"), "truth")
  )
curv <- mp |>
  filter(method == ARM_LABELS["reml_cl2"]) |>
  mutate(method = "curvature of truth")
p_cov <- ggplot(mp, aes(s, t)) +
  geom_raster(aes(fill = coverage)) +
  geom_contour(
    aes(z = coverage),
    breaks = 0.9,
    colour = "white",
    linewidth = 0.3
  ) +
  facet_grid(truth ~ method) +
  scale_fill_viridis_c(limits = c(0, 1), name = "coverage") +
  coord_equal() +
  theme_doc +
  theme(strip.text = element_text(size = 6))
p_curv <- ggplot(curv, aes(s, t, fill = curvature)) +
  geom_raster() +
  facet_grid(truth ~ method) +
  scale_fill_viridis_c(option = "magma", name = "curvature") +
  coord_equal() +
  theme_doc +
  theme(strip.text = element_text(size = 6), axis.title.y = element_blank())
fig <- p_cov + p_curv + plot_layout(widths = c(4, 1.15))
save_fig(fig, "coverage-maps-synthetic", 8, 4.2)
fig

ar1_cells <- DESIGN$cell[
  !is.na(DESIGN$extra) &
    grepl("ar1", DESIGN$extra) |
    DESIGN$block %in% c("ar1_home", "warp_ar1")
]
ar1_tab <- cov_wide(
  ar1_cells,
  c("mean", "beta"),
  c("reml_model", "reml_cl2", "ncv_cl2_bias", "ar1_model")
) |>
  filter(!is.na(ar1_model)) |>
  select(-D, -basis, -family, -signal)
ar1_is <- W |>
  filter(
    contrast == "ar1 vs proposal: interval score",
    estimand %in% c("mean", "beta")
  ) |>
  transmute(
    error,
    truth,
    estimand,
    `IS AR(1) / bias-aware NCV` = ci(ratio, ratio_lower, ratio_upper, 2)
  )
ar1_mse <- W |>
  filter(
    contrast == "ar1 vs ncv: squared error",
    estimand %in% c("mean", "beta")
  ) |>
  transmute(
    error,
    truth,
    estimand,
    `MSE AR(1) / NCV` = ci(ratio, ratio_lower, ratio_upper, 2)
  )
ar1_full <- ar1_tab |>
  left_join(ar1_is, by = c("error", "truth", "estimand")) |>
  left_join(ar1_mse, by = c("error", "truth", "estimand")) |>
  pretty()
save_tab(ar1_full, "synthetic-ar1")
rho <- FITS |>
  filter(fit == "AR1") |>
  select(cell, error, truth, rho_mean, rho_boundary)
av <- \(a, err) ar1_tab[[a]][ar1_tab$error %in% err]

show(ar1_full)

n_cells_p <- P |> distinct(app, truth, residual, flip, error, G) |> nrow()
n_tasks_p <- PFIT |> filter(fit == "REML") |> summarise(n = sum(n)) |> pull(n)
conv_p <- PFIT |>
  filter(fit == "NCV") |>
  summarise(c = weighted.mean(converged, n)) |>
  pull(c)
cohort <- PD_MODEL |>
  transmute(
    dataset = label_app(app),
    G,
    `cohort` = if_else(
      test_cohort_first50 == 1,
      "first 50 rows",
      "seeded sample of 50"
    ),
    `one-sided at 5%` = test_cohort_one_sided == 1,
    `permutation p` = test_cohort_p,
    `SMD z` = test_smd_z,
    `SMD PC1` = test_smd_PC1,
    `SMD PC2` = test_smd_PC2,
    `SMD PC3` = test_smd_PC3
  )
save_tab(cohort, "plasmode-test-cohort")
mid_above <- DATASETS |> filter(edf_ff_MID > pmax(edf_ff_REML, edf_ff_NCV))
save_tab(DATASETS, "plasmode-datasets")

# Dataset descriptions (source, access, variables, preprocessing): from the dataset screening,
# plasmode/sources.R and plasmode/export-public.R. Every number is taken from the
# diagnostics tables below; this table holds the textual facts only.
DS_META <- tibble::tribble(
  ~app,
  ~source,
  ~access,
  ~response,
  ~native,
  ~grid_note,
  ~covariate,
  ~z,
  ~prep,
  ~interest,
  "ecg_gls4b",
  "ECG-to-strain study (consulting project with N. Bouchahda's group, 2025–26); first beat per subject",
  "local only (patient and volunteer data; de-identified arrays on LRZ, aggregates here)",
  "global longitudinal strain, 4-chamber view (%), one cardiac cycle normalised to [0, 1]",
  "about 52 frames per beat (29–83)",
  "interpolated to 61 points",
  "simultaneously recorded ECG of the same beat (digitised pixel intensity, shape usable, amplitude up to scale), 51 points",
  "age (standardised)",
  "cycle time normalised to the R–R interval; subjects with at least two beats (the second beat supplies the beat-to-beat difference curves)",
  "real cardiac misregistration (largest phase share), strongest dependence among the cardiology sets with sign-changing correlations and near-zero variance at both ends of the cycle; the only model-free error source",
  "af",
  "atrial-fibrillation trial (N. Bouchahda's group): patients randomised to digoxin or beta-blocker",
  "local only (patient data; de-identified arrays on LRZ, aggregates here)",
  "left-atrial strain, 4-chamber view, at day 30 (%), one cycle on [0, 1]",
  "30 points per cycle",
  "used as is",
  "the same strain curve at day 0, 51 points",
  "randomised treatment (digoxin = 1)",
  "none beyond the cycle normalisation; a baseline-adjusted functional ANCOVA",
  "a randomised binary covariate, so γ(t) is a model-defined target informed by a trial; the coarsest grid; large phase variability with a broad, noisy peak",
  "running",
  "Fukuchi, Fukuchi & Duarte 2017 (PeerJ 5:e3298) via tidyfundata::running; one row per person, first condition",
  "public (tidyfundata, MIT; original data public)",
  "sagittal knee moment over normalised stance",
  "101 points",
  "interpolated to 61 points (as in the original 16 cells)",
  "sagittal knee angle over the same stance, 51 points",
  "body mass (standardised)",
  "one row per person (post-intervention rows dropped); speed and backpack conditions differ between persons and stay in the residuals",
  "the lowest-rank functional covariate together with gait, very smooth residuals (strongest dependence of the counted datasets)",
  "dti",
  "refund::DTI, baseline visit (Goldsmith et al. 2011): multiple-sclerosis patients and controls",
  "public (refund, GPL ≥ 2)",
  "fractional anisotropy along the right corticospinal tract",
  "55 points",
  "used as is",
  "fractional anisotropy along the corpus callosum, 51 points",
  "case status (MS = 1)",
  "complete cases of both tracts",
  "the study's application and the source of the synthetic smooth error process; no phase variation (registration control), richest covariate rank",
  "gait_vc",
  "Van Criekinge et al. 2023 (Sci Data 10:852; figshare 24192489), able-bodied adults walking barefoot at preferred speed",
  "public (CC0)",
  "sagittal knee moment over one stride, one mean stride per subject",
  "1001 points (0.1% of the stride)",
  "every 10th point (101)",
  "sagittal knee angle over the stride, 51 points",
  "age (standardised; from the subject-characteristics file)",
  "the authors' stride normalisation and averaging over strides (within-subject stride variability is averaged out)",
  "the largest G; the running data's model on an independent cohort with a walking task; a third of the residual correlations negative; no phase variation",
  "ecg_roahd",
  "roahd::mfD_healthy and mfD_LBBB (Ieva et al.), 50 healthy and 50 left-bundle-branch-block subjects",
  "public (roahd, GPL-3)",
  "ECG lead V5 (µV) over one beat",
  "1024 points at 1 kHz",
  "interpolated to 101 points",
  "ECG lead I (presumed from the documented lead order), 51 points",
  "LBBB indicator",
  "registered and smoothed by the package authors; lead identities resolved from the R-wave progression",
  "the strongest sign-changing residual correlation and the most extreme variance profile along the beat; residual phase (QRS width, LBBB morphology) survives the registration; larger intercept and γ bases",
  "ocean",
  "Hawaii Ocean Time-series profiles (FRegSigCom::ocean), monthly casts",
  "public (FRegSigCom, GPL ≥ 2; HOT data public)",
  "dissolved oxygen over depth",
  "101 points",
  "used as is",
  "temperature over depth, 51 points",
  "cast order (standardised)",
  "none; casts kept in time order",
  "stress test: consecutive monthly casts (no residual-score autocorrelation found, so the block arm is a pure cost check); very smooth residuals",
  "aemet",
  "AEMET open data via fda.usc::aemet, 73 Spanish stations, daily means 1980–2009",
  "public (fda.usc, GPL-2)",
  "log precipitation over the year",
  "365 days",
  "every 5th day (73)",
  "daily mean temperature over the year, 51 points",
  "altitude (standardised)",
  "none beyond the day thinning",
  "stress test: spatially correlated stations (residual scores correlated up to about 250 km); the lowest covariate rank of all datasets; rough residuals with almost no negative correlations; NCV rougher than REML",
  "adelaide",
  "Adelaide electricity demand (fds::mondaydemand, fds::mondaytempairport; Magnano, Boland & Hyndman 2008), every 5th Monday 1997–2007",
  "public (fds, GPL-3)",
  "half-hourly electricity demand over the day",
  "48 half-hours",
  "used as is",
  "airport temperature over the same day, 51 points",
  "week index (standardised; absorbs the linear trend)",
  "every 5th Monday (thins, not removes, the week-to-week dependence; the annual cycle stays)",
  "stress test: a weekly time series with dependence to about five weeks; strongly positively correlated residuals (2 FPCs for 90%), the largest REML–NCV EDF gap"
)
DS <- DATASETS |>
  mutate(app = as.character(app)) |>
  left_join(DS_META, by = "app") |>
  left_join(
    PD_RESID |> filter(residual == "REML") |> select(app, lag10),
    by = "app"
  ) |>
  left_join(
    PD_TRUTH |>
      select(
        app,
        truth,
        edf_alpha = `edf:s(yindex.vec)`,
        edf_gamma = `edf:s(yindex.vec):z`
      ) |>
      pivot_wider(names_from = truth, values_from = c(edf_alpha, edf_gamma)),
    by = "app"
  ) |>
  left_join(
    PD_BLOCKS |>
      select(
        app,
        n_blocks,
        block_size_min,
        block_size_max,
        score_corr_within,
        score_corr_between,
        n_pairs_between,
        n_blocks_alt
      ),
    by = "app"
  ) |>
  left_join(
    PD_MODEL |>
      select(app, test_cohort_first50, test_cohort_one_sided, test_cohort_p),
    by = "app"
  ) |>
  mutate(dataset = label_app(app), public = !startsWith(access, "local"))
# (b) data table
data_tab <- DS |>
  transmute(
    dataset,
    role,
    access = if_else(public, "public", "local only"),
    G,
    `native grid` = native,
    `fitted D` = D,
    `grid used` = grid_note,
    response,
    covariate,
    z
  )
save_tab(data_tab, "plasmode-data-table")
# (c) model table: bases per dataset (R/constants.R BASES$default and BASES$ecg8)
model_tab <- DS |>
  transmute(
    dataset,
    model = "Y(t) ~ α(t) + ∫X(s)β(s,t)ds + z γ(t), Gaussian, identity link",
    `intercept basis k` = if_else(app == "ecg_roahd", 40L, 12L),
    `gamma(t) basis k` = if_else(app == "ecg_roahd", 30L, 12L),
    `ff basis (s × t)` = "8 × 10",
    `R²` = R2,
    `ff EDF REML / NCV / MID` = sprintf(
      "%s / %s / %s",
      f(edf_ff_REML, 1),
      f(edf_ff_NCV, 1),
      f(edf_ff_MID, 1)
    ),
    `total EDF REML / NCV / MID` = sprintf(
      "%s / %s / %s",
      f(edf_total_REML, 1),
      f(edf_total_NCV, 1),
      f(edf_total_MID, 1)
    ),
    `MID interior` = mid_interior
  )
save_tab(model_tab, "plasmode-model-table")
edf_tab <- DS |>
  transmute(
    dataset,
    `alpha(t): REML / NCV / MID` = sprintf(
      "%s / %s / %s",
      f(edf_alpha_REML, 1),
      f(edf_alpha_NCV, 1),
      f(edf_alpha_MID, 1)
    ),
    `ff REML / NCV / MID` = sprintf(
      "%s / %s / %s",
      f(edf_ff_REML, 1),
      f(edf_ff_NCV, 1),
      f(edf_ff_MID, 1)
    ),
    `gamma(t): REML / NCV / MID` = sprintf(
      "%s / %s / %s",
      f(edf_gamma_REML, 1),
      f(edf_gamma_NCV, 1),
      f(edf_gamma_MID, 1)
    ),
    `total REML / NCV / MID` = sprintf(
      "%s / %s / %s",
      f(edf_total_REML, 1),
      f(edf_total_NCV, 1),
      f(edf_total_MID, 1)
    )
  )
save_tab(edf_tab, "plasmode-edf-table")
# (d) residual / dependence table
dep_tab <- DS |>
  transmute(
    dataset,
    role,
    DE,
    `DE/D` = DE_per_D,
    `lag 1 step` = lag1,
    `lag 10%` = lag10,
    `lag 25%` = lag25,
    `neg. corr. share` = neg_share,
    `min r` = min_r,
    `sd ratio` = sd_ratio,
    FPC90 = fpc90,
    `rank X` = rank_X,
    `phase share` = phase_share
  )
save_tab(dep_tab, "plasmode-dependence-table")
block_tab <- DS |>
  filter(role == "stress test") |>
  transmute(
    dataset,
    `block rule` = c(
      ocean = "3 consecutive monthly casts (cost check; the autocorrelation rule gives single casts)",
      aemet = "complete-linkage clusters of station distance cut at 250 km (mismatched partition: 150 km)",
      adelaide = "2 consecutive frame curves (10 weeks)"
    )[app],
    `blocks` = n_blocks,
    `block size` = if_else(
      block_size_min == block_size_max,
      as.character(block_size_min),
      sprintf("%d–%d", block_size_min, block_size_max)
    ),
    `alt. blocks` = n_blocks_alt,
    `score corr. within` = score_corr_within,
    `score corr. between` = score_corr_between,
    `between pairs` = n_pairs_between
  )
save_tab(block_tab, "plasmode-block-table")
ecg8_note <- "intercept k = 40 and γ k = 30 for ECG 8-lead because with k = 12 neither the functional intercept nor the LBBB effect can follow the QRS complex (22% of the REML residual sum of squares was mean-curve misfit; 0.2% with the larger bases), so that mean structure does not leak into the frozen truth's residuals"

for (i in seq_len(nrow(DS))) {
  d <- DS[i, ]
  cat(sprintf(
    "**%s** (%s; %s). %s: %d curves; response: %s, native grid %s (fitted grid: %s); functional covariate: %s; scalar covariate z: %s. Preprocessing: %s. Fit: R² %s, ff EDF %s (REML) vs %s (NCV); residual design effect %s (DE/D %s), lag-1 correlation %s, %s%% negative correlations (minimum r %s), SD ratio along t %s, %d FPCs for 90%% of the residual variance, covariate rank %d, phase share %s%%%s. Why it is here: %s.\n\n",
    d$dataset,
    d$source,
    d$access,
    if (d$role == "counted") "Counted dataset" else "Stress test",
    d$G,
    d$response,
    d$native,
    d$grid_note,
    d$covariate,
    d$z,
    d$prep,
    f(d$R2, 2),
    f(d$edf_ff_REML, 1),
    f(d$edf_ff_NCV, 1),
    f(d$DE, 1),
    f(d$DE_per_D, 2),
    f(d$lag1, 2),
    f(100 * d$neg_share, 0),
    f(d$min_r, 2),
    f(d$sd_ratio, 0),
    as.integer(d$fpc90),
    as.integer(d$rank_X),
    f(100 * d$phase_share, 0),
    if (d$role == "stress test")
      sprintf(
        "; %d blocks of %s curves, residual-score correlation %s within and %s between blocks",
        d$n_blocks,
        if (d$block_size_min == d$block_size_max)
          as.character(d$block_size_min) else
          sprintf("%d–%d", d$block_size_min, d$block_size_max),
        f(d$score_corr_within, 2),
        f(d$score_corr_between, 2)
      ) else "",
    d$interest
  ))
}

show(
  data_tab |>
    select(-covariate, -response) |>
    mutate(
      role = if_else(role == "stress test", "stress", role),
      z = sub(" \\(standardised.*\\)", " (std.)", z)
    ) |>
    rename(`native` = `native grid`)
)

show(model_tab |> select(-model), digits = 2)

show(dep_tab |> select(-role), digits = 2)

show(block_tab, digits = 2)

show(cohort)

d <- P_POOL |>
  filter(
    estimand %in% c("mean", "beta"),
    arm %in% c(MAIN_ARMS, "ar1_model"),
    residual %in% c("REML", "NCV")
  ) |>
  mutate(
    arm = label_arm(arm),
    estimand = label_estimand(estimand),
    app = label_app(app),
    truth = factor(truth, TRUTH_LEVELS),
    residual = RESIDUAL_LABELS[residual]
  )
fig <- ggplot(d, aes(coverage_mean, truth, colour = arm, shape = residual)) +
  geom_vline(xintercept = LEVEL, linetype = 2, colour = "grey40") +
  geom_pointrange(
    aes(
      xmin = coverage_mean - 2 * coverage_se,
      xmax = coverage_mean + 2 * coverage_se
    ),
    position = position_dodge(width = 0.7),
    size = 0.2
  ) +
  facet_grid(app ~ estimand) +
  scale_colour_manual(values = ARM_COLOURS, name = NULL) +
  scale_shape_manual(values = c(16, 1), name = NULL) +
  guides(colour = guide_legend(nrow = 2), shape = guide_legend(nrow = 2)) +
  labs(x = "coverage (grid average, nominal 0.95)", y = "truth source") +
  theme_doc
save_fig(fig, "plasmode-coverage", 6.5, 8.4)
fig
pool_wide <- POOL_MEANS_MAIN |>
  filter(arm %in% c(MAIN_ARMS, "ar1_model"), estimand %in% c("mean", "beta")) |>
  select(app, role, estimand, arm, coverage) |>
  pivot_wider(names_from = arm, values_from = coverage) |>
  mutate(app = label_app(app)) |>
  arrange(role, app, estimand) |>
  pretty()
save_tab(POOL_MEANS_MAIN, "plasmode-pool-means")
pv <- \(a, e, r = "counted")
  pool_wide[[ARM_LABELS[a]]][
    pool_wide$estimand == ESTIMAND_LABELS[e] & pool_wide$role == r
  ]

show(pool_wide)

tails_p <- POOL_MEANS |>
  filter(arm %in% MAIN_ARMS, estimand %in% c("mean", "beta")) |>
  select(app, role, estimand, arm, q05_mean, q05_min) |>
  pivot_wider(names_from = arm, values_from = c(q05_mean, q05_min)) |>
  mutate(app = label_app(app)) |>
  arrange(role, app, estimand)
names(tails_p) <- sub(
  "q05_mean_",
  "q05 mean: ",
  sub("q05_min_", "q05 min: ", names(tails_p))
)
save_tab(pretty(tails_p), "plasmode-lower-tails")
tv <- \(a, e, col = "q05 mean", r = "counted")
  tails_p[[paste0(col, ": ", a)]][tails_p$estimand == e & tails_p$role == r]
# REML + CL2 per truth source (beta, G = all, model residuals)
cl2_truth <- by_truth |>
  filter(arm == "reml_cl2", estimand == "beta") |>
  select(app, role, truth, coverage, q05) |>
  pivot_wider(names_from = truth, values_from = c(coverage, q05)) |>
  mutate(app = label_app(app)) |>
  select(
    app,
    role,
    `cov REML truth` = coverage_REML,
    `cov MID` = coverage_MID,
    `cov NCV truth` = coverage_NCV,
    `q05 REML truth` = q05_REML,
    `q05 MID` = q05_MID,
    `q05 NCV truth` = q05_NCV
  ) |>
  arrange(role, app)
share90 <- HMC_P |>
  filter(method == "reml_cl2") |>
  mutate(dataset = label_app(app)) |>
  select(dataset, truth, share_below_90, coverage_q05)
save_tab(cl2_truth, "plasmode-remlcl2-by-truth")
ct <- \(col, r = "counted") cl2_truth[[col]][cl2_truth$role == r]

show(cl2_truth)

zp_avail <- !is.null(ZP)
zp_tab <- if (zp_avail)
  ZP |>
    filter(
      fit == "REML",
      flip == "curve",
      error == "attached",
      G > 40,
      residual %in% c("REML", "NCV"),
      estimand %in% c("mean", "beta")
    ) |>
    group_by(app, truth, estimand) |>
    summarise(
      across(
        c(
          z_mean_cl2,
          z_sd_cl2,
          share_big_cl2,
          z_mean_model,
          z_sd_model,
          share_big_model
        ),
        mean
      ),
      .groups = "drop"
    ) |>
    mutate(
      role = if_else(app %in% APP_COUNTED, "counted", "stress test"),
      dataset = label_app(app),
      truth = factor(truth, TRUTH_LEVELS)
    ) |>
    select(dataset, role, truth, estimand, everything(), -app) |>
    arrange(estimand, role, dataset, truth) else NULL
if (zp_avail)
  save_tab(
    zp_tab |> mutate(estimand = ESTIMAND_LABELS[estimand]),
    "plasmode-zstats-reml"
  )
zpb <- if (zp_avail) zp_tab |> filter(estimand == "beta") else NULL
zpv <- \(col, r = "counted", tr = TRUTH_LEVELS)
  zpb[[col]][zpb$role == r & zpb$truth %in% tr]
zp_worst <- if (zp_avail) zpb |> arrange(desc(z_sd_cl2)) |> slice(1) else NULL
zp_cell_max <- if (zp_avail)
  ZP |>
    filter(
      fit == "REML",
      flip == "curve",
      error == "attached",
      G > 40,
      estimand == "beta",
      app %in% APP_STRESS
    ) |>
    summarise(sd = max(z_sd_cl2), big = share_big_cl2[which.max(z_sd_cl2)]) else
  NULL
zp_norm_share <- \(sd) 2 * (1 - pnorm(1.96 / sd))

g40_tab <- g40_app |>
  filter(arm %in% MAIN_ARMS) |>
  mutate(arm = label_arm(arm), app = label_app(app)) |>
  pivot_wider(names_from = arm, values_from = loss) |>
  arrange(role, app, estimand)
save_tab(g40 |> mutate(app = label_app(app)), "plasmode-G40-paired")
g40v <- \(a, e = c("mean", "beta"), r = "counted")
  g40_app$loss[g40_app$arm == a & g40_app$estimand %in% e & g40_app$role == r]
# the rule on the G = 40 cells (point estimates, no bootstrap flags)
rule40 <- P_G40 |>
  filter(
    estimand %in% c("mean", "beta"),
    arm %in% c("ncv_cl2_bias", "reml_cl2")
  ) |>
  group_by(app, role, estimand) |>
  group_modify(\(d, k) rule_both(d)) |>
  ungroup() |>
  mutate(app = label_app(app)) |>
  arrange(estimand, role, app)
save_tab(rule40, "plasmode-rule-G40-point")

d <- PM |>
  filter(
    error == "attached",
    flip == "curve",
    G_label == "all",
    residual %in% c("REML", "NCV")
  ) |>
  mutate(
    estimand = label_estimand(estimand),
    app = label_app(app),
    truth = factor(truth, TRUTH_LEVELS),
    residual = RESIDUAL_LABELS[residual]
  )
fig <- ggplot(
  d,
  aes(mse_ratio_ncv_reml, app, colour = truth, shape = residual)
) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey40") +
  geom_pointrange(
    aes(xmin = ratio_lower, xmax = ratio_upper),
    position = position_dodge(width = 0.7),
    size = 0.2
  ) +
  facet_wrap(~estimand, nrow = 1) +
  scale_x_log10(breaks = c(0.03, 0.1, 0.3, 1, 3)) +
  scale_colour_manual(
    values = c(REML = "#0072B2", MID = "#009E73", NCV = "#D55E00"),
    name = "truth"
  ) +
  scale_shape_manual(values = c(16, 1), name = NULL) +
  labs(x = "MSE of NCV / MSE of REML", y = NULL) +
  theme_doc
save_fig(fig, "plasmode-mse-ratio", 7.5, 5)
fig
ratio_wide <- RATIO_BY_TRUTH |>
  select(app, role, estimand, truth, ratio) |>
  pivot_wider(names_from = truth, values_from = ratio) |>
  mutate(app = label_app(app)) |>
  select(role, app, estimand, REML, MID, NCV) |>
  arrange(estimand, role, app) |>
  pretty()
save_tab(RATIO_BY_TRUTH, "plasmode-mse-ratio-by-truth")
order_ok <- RATIO_BY_TRUTH |>
  filter(estimand == "beta") |>
  select(app, truth, ratio) |>
  pivot_wider(names_from = truth, values_from = ratio) |>
  mutate(ok = NCV < MID & MID < REML)
losers <- RATIO_BY_TRUTH |>
  filter(truth == "REML", estimand == "beta", ratio > 1)
# univariate gains with intervals (NCV truth, REML residuals)
uni <- PM |>
  filter(
    error == "attached",
    flip == "curve",
    G_label == "all",
    truth == "NCV",
    residual == "REML",
    estimand %in% c("alpha", "gamma"),
    ratio_upper < 0.9
  ) |>
  mutate(
    dataset = label_app(app),
    value = ci(mse_ratio_ncv_reml, ratio_lower, ratio_upper, 2)
  ) |>
  select(dataset, estimand, value)
wr <- DATASETS |> filter(app %in% c("aemet", "running"))

show(ratio_wide)

rule_all <- PR |>
  filter(pool == "all") |>
  transmute(
    dataset = label_app(app),
    role,
    estimand,
    n_cells,
    `CE proposal` = ci(ce_proposal, ce_proposal_lower, ce_proposal_upper),
    `CE fallback` = ci(ce_fallback, ce_fallback_lower, ce_fallback_upper),
    `score ratio` = ci(score_ratio, score_ratio_lower, score_ratio_upper, 2),
    decision,
    class,
    `Q1 class (REML)` = q1_class_reml,
    `Q1 class (NCV)` = q1_class_ncv
  ) |>
  arrange(role, estimand, dataset)
save_tab(PR, "plasmode-rule")
save_tab(PC, "plasmode-consistency")
prc <- PR |> filter(pool == "all", role == "counted")
prs <- PR |> filter(pool == "all", role == "stress test")
n_labels <- prc |> filter(class == "N", estimand %in% c("mean", "beta"))
q1_gain <- prc |> filter(estimand %in% c("mean", "beta"))
q1_gain_s <- prs |> filter(estimand %in% c("mean", "beta"))
# the symmetric criterion on the nine-dataset pools (point estimates)
rule_sym <- P_POOL |>
  filter(
    estimand %in% c("mean", "beta"),
    arm %in% c("ncv_cl2_bias", "reml_cl2")
  ) |>
  group_by(app, role, estimand) |>
  group_modify(\(d, k) rule_both(d)) |>
  ungroup() |>
  mutate(app = label_app(app)) |>
  arrange(estimand, role, app)
save_tab(rule_sym, "plasmode-rule-both-criteria-point")
sym_diff <- rule_sym |> filter(decision_under != decision_sym)
rv_pl <- RV |>
  filter(data == "plasmode") |>
  mutate(
    group = recode(
      group,
      ecg_gls4 = "ECG strain (79-subject frame)",
      running = "running (first-50 cohort)"
    )
  )

show(rule_all)

show(
  rule_sym |>
    select(
      app,
      role,
      estimand,
      n_cells,
      ce_proposal_under,
      ce_fallback_under,
      ce_proposal_sym,
      ce_fallback_sym,
      score_ratio,
      decision_under,
      decision_sym
    )
)

att <- PA |>
  filter(
    !grepl("beat", effect),
    quantity %in%
      c(
        "log_mse_ratio",
        "coverage:reml_cl2",
        "coverage:ncv_cl2_bias",
        "coverage:reml_model"
      )
  ) |>
  mutate(dataset = label_app(app), value = ci(estimate, lower, upper)) |>
  select(dataset, estimand, quantity, effect, value, verdict)
att_wide <- att |>
  select(-verdict) |>
  pivot_wider(names_from = effect, values_from = value) |>
  left_join(
    att |>
      filter(effect == "residual") |>
      select(dataset, estimand, quantity, `residual verdict` = verdict),
    by = c("dataset", "estimand", "quantity")
  ) |>
  arrange(quantity, estimand, dataset)
save_tab(PA, "plasmode-attribution")
lm_truth <- PA |>
  filter(quantity == "log_mse_ratio", effect == "truth", estimand == "beta")
lm_resid <- PA |>
  filter(quantity == "log_mse_ratio", effect == "residual", estimand == "beta")
lm_int <- PA |>
  filter(
    quantity == "log_mse_ratio",
    effect == "interaction",
    estimand == "beta"
  )
cov_resid <- PA |>
  filter(startsWith(quantity, "coverage:"), effect == "residual")
cov_exc <- cov_resid |>
  filter(verdict != "small") |>
  mutate(
    dataset = label_app(app),
    arm = ARM_LABELS[sub("coverage:", "", quantity)]
  )
flips <- PF |> mutate(dataset = label_app(app))
save_tab(flips, "plasmode-flips")
n_pf <- sum(
  flips$flips &
    ((flips$REML == "P" & flips$NCV == "F") |
      (flips$REML == "F" & flips$NCV == "P"))
)

beat_cov <- P_POOL |>
  filter(
    app == "ecg_gls4b",
    estimand %in% c("mean", "beta"),
    arm %in% c(MAIN_ARMS, "ar1_model")
  ) |>
  select(truth, residual, estimand, arm, coverage_mean) |>
  pivot_wider(names_from = arm, values_from = coverage_mean) |>
  mutate(
    residual = RESIDUAL_LABELS[residual],
    truth = factor(truth, TRUTH_LEVELS)
  ) |>
  arrange(estimand, truth, residual)
beat_mse <- PM |>
  filter(
    app == "ecg_gls4b",
    error == "attached",
    G_label == "all",
    estimand %in% c("mean", "beta")
  ) |>
  select(truth, residual, estimand, mse_ratio_ncv_reml) |>
  mutate(residual = RESIDUAL_LABELS[residual])
beat_tab <- left_join(
  beat_cov,
  beat_mse,
  by = c("truth", "residual", "estimand")
) |>
  pretty()
save_tab(beat_tab, "plasmode-ecg-beats")
bc <- PA |>
  filter(
    grepl("beat", effect),
    quantity %in%
      c("coverage:reml_cl2", "coverage:ncv_cl2_bias", "coverage:reml_model")
  )
bcv <- \(q, which) bc$estimate[bc$quantity == q & grepl(which, bc$effect)]
bm <- \(res, e = "beta")
  beat_mse$mse_ratio_ncv_reml[
    beat_mse$residual == RESIDUAL_LABELS[res] & beat_mse$estimand == e
  ]
resid_de <- PD_RESID |>
  filter(app == "ecg_gls4b") |>
  select(residual, DE, residual_variance, sd_ratio)

show(beat_tab)

ar1_pt <- P_POOL |>
  filter(arm == "ar1_model", estimand %in% c("mean", "beta", "gamma")) |>
  select(app, role, truth, residual, estimand, coverage_mean, pointwise_q05) |>
  left_join(
    PFIT |>
      filter(fit == "AR1") |>
      select(app, truth, residual, rho, median_time),
    by = c("app", "truth", "residual")
  ) |>
  left_join(
    P_POOL |>
      filter(arm == "reml_cl2") |>
      select(app, truth, residual, estimand, reml_cl2 = coverage_mean),
    by = c("app", "truth", "residual", "estimand")
  ) |>
  mutate(
    dataset = label_app(app),
    cell_role = if_else(
      truth == "MID" & residual == "REML",
      "principal",
      "sensitivity"
    )
  ) |>
  select(
    dataset,
    role,
    cell_role,
    truth,
    residual,
    estimand,
    `AR(1)` = coverage_mean,
    `AR(1) q05` = pointwise_q05,
    `REML + CL2` = reml_cl2,
    rho,
    `fit time s` = median_time
  ) |>
  arrange(role, dataset, cell_role, estimand)
save_tab(pretty(ar1_pt), "plasmode-ar1")
ar1_v <- \(e, r = c("counted", "stress test"))
  ar1_pt[["AR(1)"]][
    ar1_pt$estimand == e & ar1_pt$role %in% r & ar1_pt$cell_role == "principal"
  ]
ar1_n_ok <- sum(
  (ar1_pt |>
    filter(cell_role == "principal", estimand %in% c("mean", "beta")) |>
    group_by(dataset) |>
    summarise(ok = all(.data[["AR(1)"]] >= 0.93)))$ok
)
ar1_sens_tab <- ar1_pt |>
  filter(cell_role == "sensitivity", estimand %in% c("mean", "beta"))

blk <- PB |>
  filter(
    name %in%
      c(
        "mse_ratio_block",
        "mse_ratio_mismatch",
        "coverage_loss",
        "coverage_loss_recipe"
      ),
    cells %in% c("block_flips", "curve_flips")
  ) |>
  mutate(dataset = label_app(app), value = ci(estimate, lower, upper)) |>
  select(dataset, estimand, verdict, cells, name, value) |>
  pivot_wider(names_from = name, values_from = value) |>
  arrange(dataset, estimand, cells)
blk_cov <- PB |>
  filter(
    startsWith(name, "coverage:"),
    cells %in% c("block_flips", "curve_flips")
  ) |>
  mutate(dataset = label_app(app), recipe = sub("coverage:", "", name)) |>
  select(dataset, estimand, cells, recipe, estimate) |>
  pivot_wider(names_from = recipe, values_from = estimate) |>
  arrange(dataset, estimand, cells)
save_tab(PB, "plasmode-block-arm")
save_tab(blk_cov, "plasmode-block-arm-coverage")
bgain <- PB |> filter(name == "mse_ratio_block", cells == "block_flips")
bcost <- PB |> filter(name == "coverage_loss", cells == "curve_flips")
blocks_d <- PD_BLOCKS |>
  select(
    app,
    n_blocks,
    block_size_min,
    block_size_max,
    score_corr_within,
    score_corr_between
  )
# between-curve dependence effect: block-flip minus curve-flip coverage of the standard recipes
bf_effect <- PB |>
  filter(
    name %in% c("coverage:standard", "coverage:reml_curvecl2"),
    cells %in% c("block_flips", "curve_flips")
  ) |>
  select(app, estimand, name, cells, estimate) |>
  pivot_wider(names_from = cells, values_from = estimate) |>
  mutate(
    diff = block_flips - curve_flips,
    dataset = label_app(app),
    recipe = recode(
      name,
      `coverage:standard` = "bias-aware NCV (curve)",
      `coverage:reml_curvecl2` = "REML + curve CL2"
    )
  )

show(blk)

hyb_p <- HSUM_P |>
  filter(
    flip == "curve",
    error == "attached",
    G > 40,
    residual %in% c("REML", "NCV"),
    estimand %in% c("beta", "mean"),
    arm %in% c("reml_cl2", "ncv_cl2_bias", "hybrid")
  ) |>
  select(
    cell,
    app,
    truth,
    estimand,
    arm,
    coverage_mean,
    pointwise_q05,
    width_mean,
    score_mean
  ) |>
  pivot_wider(
    names_from = arm,
    values_from = c(coverage_mean, pointwise_q05, width_mean, score_mean)
  )
hyb_truth <- hyb_p |>
  group_by(truth, estimand) |>
  summarise(
    n_cells = n(),
    `cov REML+CL2` = mean(coverage_mean_reml_cl2),
    `cov bias-aware` = mean(coverage_mean_ncv_cl2_bias),
    `cov hybrid` = mean(coverage_mean_hybrid),
    `q05 REML+CL2` = mean(pointwise_q05_reml_cl2),
    `q05 bias-aware` = mean(pointwise_q05_ncv_cl2_bias),
    `q05 hybrid` = mean(pointwise_q05_hybrid),
    `width bias-aware/REML+CL2` = exp(mean(log(
      width_mean_ncv_cl2_bias / width_mean_reml_cl2
    ))),
    `width hybrid/REML+CL2` = exp(mean(log(
      width_mean_hybrid / width_mean_reml_cl2
    ))),
    `IS bias-aware/REML+CL2` = exp(mean(log(
      score_mean_ncv_cl2_bias / score_mean_reml_cl2
    ))),
    `IS hybrid/REML+CL2` = exp(mean(log(
      score_mean_hybrid / score_mean_reml_cl2
    ))),
    .groups = "drop"
  ) |>
  mutate(
    truth = factor(truth, TRUTH_LEVELS),
    estimand = ESTIMAND_LABELS[estimand]
  ) |>
  arrange(estimand, truth)
save_tab(hyb_truth, "plasmode-hybrid-by-truth")
hyb_app <- HP |>
  filter(estimand %in% c("mean", "beta")) |>
  mutate(dataset = label_app(app)) |>
  select(
    dataset,
    estimand,
    n_cells,
    incl,
    cov_reml_cl2,
    cov_ncv_cl2_bias,
    cov_hybrid,
    q05_reml_cl2,
    q05_ncv_cl2_bias,
    q05_hybrid,
    `IS_ncv_cl2_bias/reml_cl2`,
    `IS_hybrid/reml_cl2`
  ) |>
  arrange(estimand, dataset)
save_tab(HP, "plasmode-hybrid-digest")
hv2 <- \(col, tr, e = "beta(s,t)")
  hyb_truth[[col]][hyb_truth$truth == tr & hyb_truth$estimand == e]
mcp <- HMC_P |> mutate(dataset = label_app(app))
mcpv <- \(m, tr, col, a = "running")
  mcp[[col]][mcp$method == m & mcp$truth == tr & mcp$app == a]

show(hyb_truth, digits = 3)

maps_p <- read_maps("plasmode") |>
  filter(truth %in% c("REML", "NCV")) |>
  mutate(
    method = label_arm(method),
    dataset = label_app(app),
    panel = paste(dataset, truth, "truth")
  )
fig <- ggplot(maps_p, aes(s, t)) +
  geom_raster(aes(fill = coverage)) +
  geom_contour(
    aes(z = coverage),
    breaks = 0.9,
    colour = "white",
    linewidth = 0.3
  ) +
  facet_grid(panel ~ method) +
  scale_fill_viridis_c(limits = c(0, 1), name = "coverage") +
  coord_equal() +
  theme_doc +
  theme(strip.text = element_text(size = 6))
save_fig(fig, "coverage-maps-plasmode", 7, 6.5)
fig

syn_core100 <- S |>
  filter(block == "core", G == 100, estimand %in% c("mean", "beta"))
syn_cov <- \(a, err, e = c("mean", "beta"))
  syn_core100$coverage_mean[
    syn_core100$arm == a &
      syn_core100$error %in% err &
      syn_core100$estimand %in% e
  ]
syn_ratio <- \(err, e = "beta", blk = "core")
  MSE_RATIO$ratio[
    MSE_RATIO$block == blk &
      MSE_RATIO$G == 100 &
      MSE_RATIO$error %in% err &
      MSE_RATIO$estimand == e
  ]

# Per-dataset outcomes: pool means over the attached curve-flip cells at G = all (all truth x
# residual combinations), MSE ratios by truth source (geometric mean over residual sources),
# the AR(1) principal cell, and the dataset diagnostics (REML residuals).
cross <- DATASETS |>
  mutate(app = as.character(app)) |>
  left_join(
    POOL_MEANS_MAIN |>
      filter(
        arm %in%
          c("reml_model", "reml_cl2", "ncv_cl2", "ncv_cl2_bias", "ar1_model"),
        estimand %in% c("mean", "beta")
      ) |>
      select(app, estimand, arm, coverage) |>
      pivot_wider(names_from = c(arm, estimand), values_from = coverage),
    by = "app"
  ) |>
  left_join(
    RATIO_BY_TRUTH |>
      filter(estimand %in% c("beta", "mean", "gamma", "alpha")) |>
      select(app, truth, estimand, ratio) |>
      pivot_wider(
        names_from = c(truth, estimand),
        values_from = ratio,
        names_prefix = "ratio_"
      ),
    by = "app"
  ) |>
  left_join(
    g40_app |>
      filter(arm == "reml_cl2", estimand == "beta") |>
      select(app, g40_loss_beta = loss),
    by = "app"
  ) |>
  left_join(
    by_truth |>
      filter(arm == "ncv_cl2_bias", truth == "REML", estimand == "beta") |>
      select(app, ba_beta_remltruth = coverage),
    by = "app"
  ) |>
  mutate(
    shortfall_model_beta = LEVEL - reml_model_beta,
    shortfall_model_mean = LEVEL - reml_model_mean,
    shortfall_cl2_beta = LEVEL - reml_cl2_beta,
    log_ratio_MID_beta = log(ratio_MID_beta),
    dataset = label_app(app)
  )
save_tab(cross, "cross-study-dataset-outcomes")
cc <- cross |> filter(role == "counted")
A <- \(x, y, d = cross) rho_text(d[[x]], d[[y]])
Av <- \(x, y, d = cross) rho_perm(d[[x]], d[[y]])$rho
Ap <- \(x, y, d = cross) rho_perm(d[[x]], d[[y]])$p

# Synthetic cell-level quantities for the scatter plots.
syn_cells <- S |>
  filter(
    estimand %in% c("mean", "beta"),
    arm %in% c("reml_model", "reml_cl2", "ncv_cl2_bias")
  ) |>
  select(
    cell,
    block,
    family,
    error,
    G,
    signal,
    truth,
    D,
    covariate,
    basis,
    estimand,
    arm,
    coverage_mean
  ) |>
  pivot_wider(names_from = arm, values_from = coverage_mean) |>
  left_join(
    MSE_RATIO |> select(cell, estimand, ratio),
    by = c("cell", "estimand")
  ) |>
  left_join(
    EDF |>
      filter(term == "te(X.smat,X.tmat):L.X", fit %in% c("REML", "NCV")) |>
      select(cell, fit, edf_mean) |>
      pivot_wider(
        names_from = fit,
        values_from = edf_mean,
        names_prefix = "edf_"
      ) |>
      mutate(edf_ratio = edf_REML / edf_NCV),
    by = "cell"
  )
# Design effect: realized values (analysis/describe-cells.R on LRZ) where available, else the
# design values of the error kernels (iid 1; OU 8; smooth 11.8 at D = 61 and 47 at D = 241; AR(1) 8).
DE_DESIGN <- c(iid = 1, ou = 8, fpc = 11.8, ou0 = 8)
syn_cells <- syn_cells |>
  mutate(
    DE_design = if_else(error == "fpc" & D == 241, 47, unname(DE_DESIGN[error]))
  )
if (!is.null(SDE))
  syn_cells <- left_join(
    syn_cells,
    SDE |> select(cell, DE = design_effect),
    by = "cell"
  )
has_de <- "DE" %in% names(syn_cells)
if (!has_de) syn_cells$DE <- NA_real_
syn_cells <- syn_cells |>
  mutate(
    DE_realized = !is.na(DE),
    DE_plot = coalesce(DE, DE_design),
    DE_per_D = DE_plot / D
  )
de_missing <- syn_cells |>
  filter(!DE_realized) |>
  distinct(cell, block, family, error)
save_tab(syn_cells, "cross-study-synthetic-cells")
de_label <- if (has_de) "design effect 1'S1/tr(S) of the residual curves" else
  "design effect (synthetic: design values, realized values pending)"
syn_de_span <- syn_cells |>
  filter(block == "core", G == 100, error != "iid") |>
  summarise(lo = min(DE_plot), hi = max(DE_plot))

label_pts <- \(d)
  geom_text(
    data = d,
    aes(x, y, label = dataset),
    size = 2.2,
    vjust = -0.8,
    colour = "black"
  )
sy <- syn_cells |>
  filter(
    estimand == "beta",
    G == 100,
    block %in%
      c(
        "core",
        "dense_grid",
        "oscillating",
        "ar1_home",
        "heteroskedastic",
        "lowrank_covariate",
        "rough_truth",
        "basis_size"
      )
  )
fam_cols <- c(gaussian = "#000000", poisson = "#56B4E9", binomial = "#CC79A7")
base <- \(xl, yl, ylim = NULL)
  list(
    labs(x = xl, y = yl),
    theme_doc,
    scale_colour_manual(values = fam_cols, name = "synthetic family"),
    scale_shape_manual(
      values = c(counted = 16, `stress test` = 17),
      name = "plasmode dataset"
    )
  )
sy_de <- sy |>
  mutate(
    de_shape = if_else(
      error %in% c("osc", "warp"),
      "sign-changing / misregistration cell",
      "other synthetic cell"
    )
  )
pa <- ggplot() +
  geom_hline(yintercept = LEVEL, linetype = 2, colour = "grey40") +
  geom_point(
    data = sy_de |> filter(de_shape == "other synthetic cell"),
    aes(DE_plot, reml_model, colour = family),
    size = 0.9,
    alpha = 0.6
  ) +
  geom_point(
    data = sy_de |> filter(de_shape != "other synthetic cell"),
    aes(DE_plot, reml_model, colour = family),
    size = 1.6,
    shape = 4,
    show.legend = FALSE
  ) +
  geom_point(
    data = cross |> mutate(x = DE, y = reml_model_beta),
    aes(x, y, shape = role),
    size = 2.5
  ) +
  label_pts(cross |> mutate(x = DE, y = reml_model_beta)) +
  base(de_label, "model-based coverage of beta") +
  ggtitle("(a)")
pb <- ggplot() +
  geom_hline(yintercept = LEVEL, linetype = 2, colour = "grey40") +
  geom_point(
    data = sy_de |> filter(de_shape == "other synthetic cell"),
    aes(DE_plot, reml_cl2, colour = family),
    size = 0.9,
    alpha = 0.6
  ) +
  geom_point(
    data = sy_de |> filter(de_shape != "other synthetic cell"),
    aes(DE_plot, reml_cl2, colour = family),
    size = 1.6,
    shape = 4,
    show.legend = FALSE
  ) +
  geom_point(
    data = cross |> mutate(x = DE, y = reml_cl2_beta),
    aes(x, y, shape = role),
    size = 2.5
  ) +
  label_pts(cross |> mutate(x = DE, y = reml_cl2_beta)) +
  base(de_label, "REML + CL2 coverage of beta") +
  ggtitle("(b)")
pc <- ggplot() +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey40") +
  geom_point(
    data = sy |> filter(!is.na(edf_ratio), truth == "smooth"),
    aes(edf_ratio, ratio, colour = family),
    size = 0.9,
    alpha = 0.6
  ) +
  geom_point(
    data = cross |> mutate(x = edf_ratio, y = ratio_MID_beta),
    aes(x, y, shape = role),
    size = 2.5
  ) +
  label_pts(cross |> mutate(x = edf_ratio, y = ratio_MID_beta)) +
  scale_y_log10(breaks = c(0.01, 0.03, 0.1, 0.3, 1, 3)) +
  scale_x_log10() +
  scale_colour_manual(values = fam_cols, name = "synthetic family") +
  scale_shape_manual(
    values = c(counted = 16, `stress test` = 17),
    name = "plasmode dataset"
  ) +
  labs(
    x = "EDF ratio REML / NCV of the ff term",
    y = "MSE NCV / MSE REML for beta"
  ) +
  theme_doc +
  ggtitle("(c)")
syG <- syn_cells |> filter(estimand == "beta", block == "core")
pd_ <- ggplot() +
  geom_hline(yintercept = LEVEL, linetype = 2, colour = "grey40") +
  geom_jitter(
    data = syG,
    aes(G, reml_cl2, colour = family),
    size = 0.9,
    alpha = 0.6,
    width = 2,
    height = 0
  ) +
  geom_point(
    data = cross |> mutate(x = G, y = reml_cl2_beta),
    aes(x, y, shape = role),
    size = 2.5
  ) +
  geom_point(
    data = cross |> mutate(x = 40, y = reml_cl2_beta - g40_loss_beta),
    aes(x, y, shape = role),
    size = 2.5,
    colour = "grey50",
    show.legend = FALSE
  ) +
  label_pts(cross |> mutate(x = G, y = reml_cl2_beta)) +
  base("number of curves G", "REML + CL2 coverage of beta") +
  ggtitle("(d)")
fig <- (pa | pb) /
  (pc | pd_) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
save_fig(fig, "cross-study-scatter", 8, 7)
fig

# Rankings and associations computed here; the table cells are short, the evidence follows as a list.
rank_cl2 <- cc |> arrange(reml_cl2_beta) |> mutate(rank = row_number())
rank_gain_all <- cross |> arrange(ratio_MID_beta) |> mutate(rank = row_number())
low_rank <- cc |> arrange(rank_X) |> slice(1:2)
alpha_gain_n <- sum(cross$ratio_NCV_alpha < 0.9)
signal_p <- c(Ap("R2", "shortfall_model_beta"), Ap("R2", "log_ratio_MID_beta"))
n_adequate_beta <- sum(
  PR$ce_fallback[
    PR$pool == "all" & PR$role == "counted" & PR$estimand == "beta"
  ] <=
    ADEQUATE
)
de_osc <- syn_cells$DE_plot[
  syn_cells$error == "osc" & syn_cells$estimand == "beta"
]
de_warp <- syn_cells$DE_plot[
  syn_cells$error == "warp" &
    syn_cells$estimand == "beta" &
    syn_cells$DE_realized
]
cov_oscwarp <- syn_cells$reml_model[
  syn_cells$error %in% c("osc", "warp") & syn_cells$estimand == "beta"
]
cross_rows <- tibble::tribble(
  ~property,
  ~`synthetic: inference`,
  ~`synthetic: estimation`,
  ~`plasmode`,
  ~verdict,
  "Dependence strength (design effect DE)",
  sprintf(
    "model-based %s (iid) → %s (dependent); REML + CL2 %s",
    rng(syn_cov("reml_model", "iid")),
    rng(syn_cov("reml_model", c("ou", "fpc"))),
    rng(syn_cov("reml_cl2", c("iid", "ou", "fpc")))
  ),
  sprintf(
    "β ratio %s (iid) vs %s (dependent)",
    rng(syn_ratio("iid")),
    rng(syn_ratio(c("ou", "fpc")))
  ),
  sprintf(
    "DE %s; model-based shortfall vs DE %s; CL2 improves coverage on every dataset and meets the point adequacy bound for β on %s of 6 counted datasets",
    rng(cross$DE, 0),
    A("DE", "shortfall_model_beta"),
    n_adequate_beta
  ),
  "reinforces the direction; DE is a grid-mean variance inflation that misses sign-changing and misregistration dependence and is confounded with D",
  "Dependence shape (sign changes; variance along t or by covariate)",
  sprintf(
    "sign-changing: model-based %s, REML + CL2 %s; covariate-dependent variance breaks model-based γ (%s) and AR(1) (%s), not CL2 (%s)",
    rng(osc_tab$reml_model),
    rng(osc_tab$reml_cl2),
    rng(hv("reml_model", "gamma", c("fpc_het", "fpc_raw_het"))),
    rng(na.omit(hv("ar1_model", "gamma", c("fpc_het", "fpc_raw_het")))),
    rng(hv("reml_cl2", "gamma", c("fpc_het", "fpc_raw_het")))
  ),
  sprintf(
    "one cell per family; the realized DE of the sign-changing (%s) and misregistration (%s) cells is near 1 although model-based coverage is %s",
    rng(de_osc, 1),
    rng(de_warp, 1),
    rng(cov_oscwarp)
  ),
  sprintf(
    "negative-correlation share %s%%, SD ratio %s; shortfall vs SD ratio %s; CL2 β coverage vs SD ratio %s (counted only %s)",
    rng(100 * cross$neg_share, 0),
    rng(cross$sd_ratio, 0),
    A("sd_ratio", "shortfall_model_beta"),
    A("sd_ratio", "reml_cl2_beta"),
    A("sd_ratio", "reml_cl2_beta", cc)
  ),
  "not separately testable (shape, DE and dataset role are confounded)",
  "Misregistration (phase variability)",
  sprintf(
    "model-based %s, REML + CL2 %s, AR(1) mean %s; Poisson mean undercovers for every recipe (%s)",
    rng(wv("reml_model", "gaussian")),
    rng(wv("reml_cl2", "gaussian")),
    rng(na.omit(warp_tab$ar1_model[warp_tab$estimand == "mean"])),
    rng(wv("reml_cl2", "poisson"))
  ),
  sprintf(
    "β ratio %s",
    rng(MSE_RATIO$ratio[
      MSE_RATIO$block == "warp" & MSE_RATIO$estimand == "beta"
    ])
  ),
  sprintf(
    "phase share %s%%; AR(1) β vs phase %s; CL2 β vs phase %s",
    rng(100 * cross$phase_share, 0),
    A("phase_share", "ar1_model_beta"),
    A("phase_share", "reml_cl2_beta")
  ),
  "partly (CL2 near nominal on the two cardiac datasets with real misregistration; AR(1)'s real-data shortfall has no single attributable cause)",
  "Number of curves G",
  sprintf(
    "REML + CL2 loses %s at G = 40; NCV + CL2 %s",
    rng_pp(gl("reml_cl2")),
    rng_pp(gl("ncv_cl2"))
  ),
  sprintf(
    "β ratio %s (G = 40) vs %s (G = 100)",
    rng(prof_G$G40[prof_G$estimand == "beta"]),
    rng(prof_G$G100[prof_G$estimand == "beta"])
  ),
  sprintf(
    "within-dataset G = all vs 40: REML + CL2 loses %s, NCV + CL2 %s; across datasets (G %s) CL2 β vs G %s",
    rng_pp(g40_app$loss[
      g40_app$arm == "reml_cl2" & g40_app$estimand == "beta"
    ]),
    rng_pp(g40_app$loss[g40_app$arm == "ncv_cl2" & g40_app$estimand == "beta"]),
    rng(cross$G, 0),
    A("G", "reml_cl2_beta")
  ),
  "reinforces (paired within-dataset contrast of the same size)",
  "Grid density D",
  sprintf(
    "D 61 → 241: model-based β −%s; REML + CL2 %s",
    rng_pp(-gd("coverage_reml_model")),
    rng_pp(gd("coverage_reml_cl2"))
  ),
  sprintf("NCV β advantage × %s", rng(exp(-gd("log_mse_ratio_ncv_reml")), 1)),
  sprintf(
    "not varied within a dataset (D %s); DE/D = (1'S1/trS)/D %s; shortfall vs DE/D %s",
    rng(cross$D, 0),
    rng(cross$DE_per_D, 2),
    A("DE_per_D", "shortfall_model_beta")
  ),
  "not testable",
  "Covariate rank (identifiability of β)",
  sprintf(
    "low-rank X: REML + CL2 β %s",
    rng(S$coverage_mean[
      S$block == "lowrank_covariate" &
        S$error == "fpc" &
        S$estimand == "beta" &
        S$arm == "reml_cl2"
    ])
  ),
  sprintf(
    "β ratio %s (low-rank) vs %s (rich), smooth process",
    rng(lowr$ratio[lowr$estimand == "beta" & lowr$error == "fpc"]),
    rng(MSE_RATIO$ratio[
      MSE_RATIO$block == "core" &
        MSE_RATIO$G == 100 &
        MSE_RATIO$error == "fpc" &
        MSE_RATIO$signal == "mid" &
        MSE_RATIO$estimand == "beta" &
        MSE_RATIO$family %in% c("gaussian", "poisson", "binomial")
    ])
  ),
  sprintf(
    "rank %s; MID-truth β gain vs rank %s; CL2 β vs rank %s; the two lowest-rank counted datasets (%s) rank %s of 6 in CL2 β coverage",
    rng(cross$rank_X, 0),
    A("rank_X", "log_ratio_MID_beta"),
    A("rank_X", "reml_cl2_beta"),
    paste(APP_LABELS[as.character(low_rank$app)], collapse = ", "),
    paste(rank_cl2$rank[match(low_rank$app, rank_cl2$app)], collapse = " and ")
  ),
  "partly (estimation direction consistent; coverage association uninformative at n = 9)",
  "Truth roughness / EDF gap between REML and NCV",
  sprintf(
    "rough truth: NCV + CL2 β %s, bias-aware %s, REML + CL2 %s",
    rng(stv("ncv_cl2", "wiggly")),
    rng(stv("ncv_cl2_bias", "wiggly")),
    rng(stv("reml_cl2", "wiggly"))
  ),
  sprintf(
    "β ratio %s (rough) vs %s (smooth); EDF ratio %s",
    rng(rough$ratio[rough$estimand == "beta"]),
    rng(syn_ratio("fpc")),
    rng(er("beta(s,t)", c("ou", "fpc")))
  ),
  sprintf(
    "β ratio ordered NCV < MID < REML truth on %s of 9; NCV loses under REML truth on %s; MID gain vs EDF ratio %s; weather: NCV rougher than REML yet wins under all truths",
    sum(order_ok$ok),
    nrow(losers),
    A("edf_ratio", "log_ratio_MID_beta")
  ),
  "partly (strong truth-source sensitivity; the EDF-gap mechanism is not identified in the plasmode study)",
  "Signal strength",
  sprintf(
    "REML + CL2 within %s of nominal across signal levels",
    pp(max(abs(
      PROF$estimate[
        PROF$factor == "signal" & PROF$value == "coverage_reml_cl2"
      ] -
        LEVEL
    )))
  ),
  sprintf(
    "β ratio %s (low) to %s (high), Gaussian smooth process",
    rng(MSE_RATIO$ratio[
      MSE_RATIO$block == "core" &
        MSE_RATIO$G == 100 &
        MSE_RATIO$error == "fpc" &
        MSE_RATIO$family == "gaussian" &
        MSE_RATIO$signal == "low" &
        MSE_RATIO$estimand == "beta"
    ]),
    rng(MSE_RATIO$ratio[
      MSE_RATIO$block == "core" &
        MSE_RATIO$G == 100 &
        MSE_RATIO$error == "fpc" &
        MSE_RATIO$family == "gaussian" &
        MSE_RATIO$signal == "high" &
        MSE_RATIO$estimand == "beta"
    ])
  ),
  sprintf(
    "R² %s; shortfall vs R² %s; MID gain vs R² %s",
    rng(cross$R2, 2),
    A("R2", "shortfall_model_beta"),
    A("R2", "log_ratio_MID_beta")
  ),
  "uninformative at n = 9 (confounded)",
  "Term type and estimand",
  sprintf(
    "same pattern for α, γ, β, mean, f(x,t); REML + CL2 differs by %s between f(x,t) and ff",
    rng_pp(m2v("coverage_reml_cl2"))
  ),
  sprintf(
    "gain on the surface: γ(t) ratio %s, α(t) %s (dependent, G = 100)",
    rng(ratio_of(r_dep, "gamma")),
    rng(ratio_of(r_dep, "alpha"))
  ),
  sprintf(
    "ff model only; γ ratio %s / %s (NCV / REML truth), α %s / %s; α gains < 0.9 on %s datasets under the NCV truth",
    rng(rt("NCV", "gamma")),
    rng(rt("REML", "gamma")),
    rng(rt("NCV", "alpha")),
    rng(rt("REML", "alpha")),
    alpha_gain_n
  ),
  "partly (gains most consistent for the surface; univariate gains on several datasets; f(x,t) not testable)",
  "Basis size",
  sprintf(
    "bias-aware β %s → %s; REML + CL2 %s → %s (default → xlarge)",
    rng(bs("ncv_cl2_bias", "default")),
    rng(bs("ncv_cl2_bias", "xlarge")),
    rng(bs("reml_cl2", "default")),
    rng(bs("reml_cl2", "xlarge"))
  ),
  sprintf(
    "β ratio %s → %s",
    rng(mse_b$ratio[
      mse_b$basis == "default" &
        mse_b$error == ERROR_LABELS["fpc"] &
        mse_b$estimand == ESTIMAND_LABELS["beta"]
    ]),
    rng(mse_b$ratio[
      mse_b$basis == "xlarge" &
        mse_b$error == ERROR_LABELS["fpc"] &
        mse_b$estimand == ESTIMAND_LABELS["beta"]
    ])
  ),
  sprintf(
    "default bases except ECG 8-lead (k = 40/30 by design), whose MID-truth β gain ranks %s of 9 (1 = largest)",
    rank_gain_all$rank[rank_gain_all$app == "ecg_roahd"]
  ),
  "not testable (one dataset, two changes)",
  "Response family",
  sprintf(
    "six families: model-based %s, REML + CL2 %s (smooth process, G = 100)",
    rng(c(syn_cov("reml_model", "fpc"), fv("reml_model"))),
    rng(c(syn_cov("reml_cl2", "fpc"), fv("reml_cl2")))
  ),
  sprintf(
    "β ratio %s across families",
    rng(c(
      syn_ratio("fpc"),
      MSE_RATIO$ratio[
        MSE_RATIO$block == "families" &
          MSE_RATIO$error == "fpc" &
          MSE_RATIO$estimand == "beta" &
          MSE_RATIO$fit_family == MSE_RATIO$family
      ]
    ))
  ),
  "Gaussian only",
  "not testable",
  "Between-curve dependence",
  "not in the design (independent curves)",
  "—",
  sprintf(
    "block-flip minus curve-flip coverage: REML + curve CL2 %s; block NCV/CL2 MSE ratio %s",
    rng_pp(bf_effect$diff[bf_effect$name == "coverage:reml_curvecl2"]),
    rng(bgain$estimate)
  ),
  "plasmode only (imposed block dependence; curve-flip shortfalls are not evidence about it)"
)
save_tab(cross_rows, "cross-study-assessment")

show(cross_rows)

z_de <- rho_perm(cross$DE, cross$shortfall_model_beta)
z_ded <- rho_perm(cross$DE_per_D, cross$shortfall_model_beta)
z_edf <- rho_perm(cross$edf_ratio, cross$log_ratio_MID_beta)
z_sd <- rho_perm(cross$sd_ratio, cross$reml_cl2_beta)
z_sd_c <- rho_perm(cc$sd_ratio, cc$reml_cl2_beta)
resid_excess <- cross |>
  mutate(
    pred = predict(lm(shortfall_model_beta ~ log(DE), data = cross)),
    excess = shortfall_model_beta - pred
  ) |>
  arrange(desc(excess))
min_de <- cross |> filter(DE == min(DE))
gain_order <- cross |> arrange(desc(ratio_MID_beta))
gain_order_c <- cc |> arrange(desc(ratio_MID_beta))

n_rep <- unique(S$n_rep)

show(tail_tab)

show(g_sum |> filter(dependence == "dependent"))

show(st_tab)

show(het_show)

show(osc_show)

show(lowrank_tab)

show(edf_tab)

show(pretty(tails_p))

show(g40_tab)

show(
  rule40 |>
    select(
      app,
      role,
      estimand,
      n_cells,
      ce_proposal_under,
      ce_fallback_under,
      ce_proposal_sym,
      ce_fallback_sym,
      score_ratio,
      decision_under,
      decision_sym
    )
)

show(
  rv_pl |>
    select(
      dataset = group,
      estimand,
      criterion,
      n_cells,
      ce_proposal,
      ce_fallback,
      score_ratio,
      decision
    ) |>
    pretty()
)

show(att_wide)

show(pretty(ar1_pt))

show(blk_cov)

show(hyb_app)

if (zs_avail) show(pretty(zs_tab)) else cat("pending")

if (zp_avail) show(zp_tab |> mutate(estimand = ESTIMAND_LABELS[estimand])) else
  cat("pending")

supp_cov <- S |>
  filter(
    estimand %in% c("mean", "beta"),
    arm %in% c(MAIN_ARMS, "ncv_model", "ncv_cl2f_bias", "ar1_model")
  ) |>
  select(
    cell,
    block,
    family,
    fit_family,
    error,
    G,
    signal,
    truth,
    D,
    model,
    covariate,
    basis,
    estimand,
    arm,
    coverage_mean,
    coverage_se,
    pointwise_q05,
    width_mean,
    score_mean
  )
save_tab(supp_cov, "S1-synthetic-coverage-per-cell")
supp_cov_wide <- supp_cov |>
  mutate(v = sprintf("%s (%s)", f(coverage_mean, 3), f(coverage_se, 3))) |>
  select(
    cell,
    block,
    family,
    error,
    G,
    signal,
    truth,
    D,
    basis,
    estimand,
    arm,
    v
  ) |>
  pivot_wider(names_from = arm, values_from = v) |>
  pretty()
supp_ag <- S |>
  filter(estimand %in% c("alpha", "gamma", "f"), arm %in% MAIN_ARMS) |>
  select(
    cell,
    block,
    family,
    error,
    G,
    signal,
    truth,
    D,
    basis,
    estimand,
    arm,
    coverage_mean,
    coverage_se,
    pointwise_q05
  )
save_tab(supp_ag, "S2-synthetic-coverage-alpha-gamma-f")
supp_ag_wide <- supp_ag |>
  mutate(v = sprintf("%s (%s)", f(coverage_mean, 3), f(coverage_se, 3))) |>
  select(-coverage_mean, -coverage_se, -pointwise_q05) |>
  pivot_wider(names_from = arm, values_from = v) |>
  pretty()
supp_mse <- S |>
  filter(arm %in% c("reml_cl2", "ncv_cl2")) |>
  transmute(
    cell,
    block,
    family,
    error,
    G,
    signal,
    truth,
    D,
    basis,
    estimand,
    fit = if_else(arm == "reml_cl2", "REML", "NCV"),
    mse,
    mse_se,
    bias2,
    bias2_se,
    variance
  ) |>
  left_join(
    MSE_RATIO |> select(cell, estimand, ratio, ratio_lower, ratio_upper),
    by = c("cell", "estimand")
  )
save_tab(supp_mse, "S3-synthetic-mse-per-cell")
supp_mse_wide <- supp_mse |>
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
    fit,
    mse
  ) |>
  pivot_wider(names_from = fit, values_from = mse, names_prefix = "MSE ") |>
  left_join(
    MSE_RATIO |>
      transmute(
        cell,
        estimand,
        `ratio NCV/REML` = ci(ratio, ratio_lower, ratio_upper, 2)
      ),
    by = c("cell", "estimand")
  ) |>
  pretty()
supp_p <- P |>
  filter(estimand %in% c("mean", "beta", "alpha", "gamma")) |>
  mutate(dataset = label_app(app)) |>
  select(
    dataset,
    role,
    truth,
    residual,
    flip,
    error,
    G,
    estimand,
    arm,
    n_rep,
    coverage_mean,
    coverage_se,
    pointwise_q05,
    width_mean,
    score_mean,
    mse,
    bias2,
    variance
  )
save_tab(supp_p, "S4-plasmode-coverage-per-cell")
supp_p_wide <- supp_p |>
  filter(arm %in% c(MAIN_ARMS, "ar1_model"), estimand %in% c("mean", "beta")) |>
  mutate(v = sprintf("%s (%s)", f(coverage_mean, 3), f(coverage_se, 3))) |>
  select(dataset, truth, residual, flip, error, G, estimand, arm, v) |>
  pivot_wider(names_from = arm, values_from = v) |>
  arrange(dataset, estimand, desc(G), truth, residual) |>
  pretty()
supp_pm <- PM |>
  mutate(dataset = label_app(app)) |>
  select(
    dataset,
    role,
    truth,
    residual,
    flip,
    error,
    G,
    estimand,
    n_pairs,
    mse_ratio_ncv_reml,
    ratio_lower,
    ratio_upper
  )
save_tab(supp_pm, "S5-plasmode-mse-ratio-per-cell")
supp_pm_wide <- supp_pm |>
  mutate(v = ci(mse_ratio_ncv_reml, ratio_lower, ratio_upper, 2)) |>
  select(dataset, truth, residual, flip, error, G, estimand, v) |>
  pivot_wider(names_from = estimand, values_from = v) |>
  arrange(dataset, desc(G), truth, residual) |>
  pretty()

show(supp_cov_wide)

show(supp_ag_wide)

show(supp_mse_wide, digits = 4)

show(supp_p_wide)

show(supp_pm_wide)

s7_full <- REL |>
  select(
    cell,
    block,
    family,
    signal,
    error,
    G,
    truth,
    basis,
    estimand,
    arm,
    coverage,
    rel_rmse_median,
    rel_rmse_q90,
    share_rel_gt_1,
    rel_halfwidth_median,
    detect,
    wrong_sign
  )
save_tab(s7_full, "relative-error-cells")
s7 <- REL |>
  filter(estimand == "beta", arm %in% c("reml_cl2", "ncv_cl2")) |>
  mutate(fit = if_else(arm == "reml_cl2", "R", "N")) |>
  select(
    cell,
    block,
    family,
    signal,
    error,
    G,
    fit,
    rel = rel_rmse_median,
    hw = rel_halfwidth_median,
    det = detect
  ) |>
  pivot_wider(
    names_from = fit,
    values_from = c(rel, hw, det),
    names_sep = "_"
  ) |>
  arrange(cell)
show(s7, digits = 2)

tibble::tibble(file = list.files(TAB_DIR, pattern = "\\.csv$")) |> show()
