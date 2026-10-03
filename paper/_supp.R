# Derivations for paper/supplement.qmd that are not in the frozen write-up code
# (paper/supp/src): Satterthwaite summaries, intervals centred at the NCV fit, the
# smoothing-parameter uncertainty terms, timing. Sourced after paper/_common.R and
# paper/supp/_objects.R, the repository root as working directory. Reads only summary CSVs.

dep_lab <- \(e) if_else(e == "iid", "independent", "dependent")
fam_lab <- \(x) unname(FAMILY_LABELS[x])
app_lab <- \(x) unname(APP_LABELS[x])
est_lab <- \(x) unname(ESTIMAND_LABELS[x])
pct <- \(x, d = 1) f(100 * x, d)

# ---- S1.7 Satterthwaite ------------------------------------------------------------------
SAT <- readr::read_csv(
  "satterthwaite/summaries/satt-cells.csv",
  show_col_types = FALSE
)
SATD <- readr::read_csv(
  "satterthwaite/summaries/satt-df.csv",
  show_col_types = FALSE
)
CRIT_LAB <- c(
  z = "z",
  t = "t(G-1)",
  satt = "Satterthwaite",
  satt_b2 = "Satterthwaite incl. Vp - Ve"
)
sat_w <- SAT |>
  filter(crit %in% names(CRIT_LAB)) |>
  mutate(dependence = dep_lab(error)) |>
  group_by(dependence, G, family, estimand, crit) |>
  summarise(
    cov = mean(coverage),
    q05 = mean(pointwise_q05),
    width = mean(width),
    .groups = "drop"
  ) |>
  group_by(dependence, G, family, estimand) |>
  mutate(wr = width / width[crit == "z"]) |>
  ungroup()
sat_tab <- sat_w |>
  mutate(
    v = sprintf("%s (%s) %s", f(cov, 3), f(q05, 2), f(wr, 2)),
    crit = factor(CRIT_LAB[crit], CRIT_LAB)
  ) |>
  select(dependence, G, family, estimand, crit, v) |>
  arrange(crit) |>
  tidyr::pivot_wider(names_from = crit, values_from = v) |>
  mutate(family = fam_lab(family), estimand = est_lab(estimand)) |>
  arrange(desc(dependence), G, family, estimand)
satv <- \(
  cr,
  g,
  e = c("mean", "beta", "gamma"),
  dep = "dependent",
  col = "cov",
  fam = c("gaussian", "poisson", "binomial")
)
  sat_w[[col]][
    sat_w$crit == cr &
      sat_w$G %in% g &
      sat_w$estimand %in% e &
      sat_w$dependence == dep &
      sat_w$family %in% fam
  ]
sat_df <- SATD |>
  group_by(G) |>
  summarise(
    med = median(df_median),
    lo = min(df_median),
    hi = max(df_median),
    b2 = median(df_b2_median),
    .groups = "drop"
  )

# ---- S5.1 intervals centred at the NCV fit -------------------------------------------------
# (a) diagnosis: bias-eliminated coverage (pointwise bias estimated from the replicates removed)
ND_S <- read_summary("final", "ncv-diag-cells") |>
  left_join(DESIGN, by = "cell")
ND_P <- read_summary("plasmode", "ncv-diag-cells")
diag_s <- ND_S |>
  filter(
    estimand == "beta",
    error != "iid",
    G == 100,
    block %in% c("core", "rough_truth"),
    interval %in% c("ncv_cl2", "reml_cl2")
  ) |>
  group_by(
    setting = if_else(block == "core", "smooth truth", "rough truth"),
    interval
  ) |>
  summarise(
    n = n(),
    cov = mean(coverage),
    q05 = mean(coverage_q05),
    be = mean(coverage_bias_elim),
    .groups = "drop"
  )
dsv <- \(st, iv, col)
  diag_s[[col]][diag_s$setting == st & diag_s$interval == iv]

RBC_P <- read_pcells("ncv-rbc/summaries/rbc-plasmode-cells.csv")
FIX_P <- read_pcells("ncv-rbc/summaries/fixlam-plasmode-cells.csv")
COMBO_P <- read_pcells("ncv-rbc/summaries/combo-plasmode-cells.csv")
diag_p <- bind_rows(
  RBC_P |> filter(estimand == "beta", arm %in% c("ncv_cl2", "reml_cl2")),
  FIX_P |> filter(estimand == "beta", arm == "ncvfix_cl2")
) |>
  group_by(app, role, arm) |>
  summarise(
    cov = mean(coverage),
    q05 = mean(coverage_q05),
    se_sd = median(se_sd),
    .groups = "drop"
  )
be_p <- ND_P |>
  filter(
    estimand == "beta",
    interval == "ncv_cl2",
    flip == "curve",
    error == "attached",
    G > 40,
    residual %in% c("REML", "NCV")
  ) |>
  group_by(app) |>
  summarise(be = mean(coverage_bias_elim), .groups = "drop")
diag_p_tab <- diag_p |>
  mutate(v = sprintf("%s (%s) %s", f(cov, 2), f(q05, 2), f(se_sd, 2))) |>
  select(app, role, arm, v) |>
  tidyr::pivot_wider(names_from = arm, values_from = v) |>
  left_join(be_p, by = "app") |>
  mutate(app = factor(app, names(APP_LABELS))) |>
  arrange(desc(role == "counted"), app) |>
  transmute(
    dataset = app_lab(as.character(app)),
    `NCV + CL2` = ncv_cl2,
    `bias removed` = f(be, 2),
    `NCV, λ fixed` = ncvfix_cl2,
    `REML + CL2` = reml_cl2
  )
dpv <- \(a, col, r = "counted")
  diag_p[[col]][diag_p$arm == a & diag_p$role == r]

# (b), (c) bias allowance and NCV estimate with REML + CL2 width
hyb_s <- bind_rows(
  HSUM_S |>
    filter(block == "core", G == 100, error != "iid") |>
    mutate(setting = "dependent, smooth truth"),
  HSUM_S |>
    filter(block == "rough_truth", error == "fpc") |>
    mutate(setting = "dependent, rough truth"),
  HSUM_S |>
    filter(block == "core", G == 100, error == "iid") |>
    mutate(setting = "independent")
) |>
  filter(
    estimand == "beta",
    arm %in% c("reml_cl2", "ncv_cl2", "ncv_cl2_bias", "hybrid")
  ) |>
  select(
    setting,
    cell,
    arm,
    coverage_mean,
    pointwise_q05,
    width_mean,
    score_mean
  ) |>
  group_by(setting, cell) |>
  mutate(
    wr = width_mean / width_mean[arm == "reml_cl2"],
    sr = score_mean / score_mean[arm == "reml_cl2"]
  ) |>
  group_by(setting, arm) |>
  summarise(
    n = n(),
    cov = mean(coverage_mean),
    q05 = mean(pointwise_q05),
    wr = geo(wr),
    sr = geo(sr),
    .groups = "drop"
  )
hyb_p <- HSUM_P |>
  filter(
    flip == "curve",
    error == "attached",
    G > 40,
    residual %in% c("REML", "NCV"),
    estimand == "beta",
    arm %in% c("reml_cl2", "ncv_cl2", "ncv_cl2_bias", "hybrid")
  ) |>
  select(
    truth,
    cell,
    arm,
    coverage_mean,
    pointwise_q05,
    width_mean,
    score_mean
  ) |>
  group_by(truth, cell) |>
  mutate(
    wr = width_mean / width_mean[arm == "reml_cl2"],
    sr = score_mean / score_mean[arm == "reml_cl2"]
  ) |>
  group_by(
    setting = paste(
      "plasmode,",
      c(REML = "REML-fitted", MID = "midpoint", NCV = "NCV-fitted")[truth],
      "truth"
    ),
    arm
  ) |>
  summarise(
    n = n(),
    cov = mean(coverage_mean),
    q05 = mean(pointwise_q05),
    wr = geo(wr),
    sr = geo(sr),
    .groups = "drop"
  )
hyb_all <- bind_rows(hyb_s, hyb_p)
HYB_ARMS <- c(
  reml_cl2 = "REML + CL2",
  ncv_cl2 = "NCV + CL2",
  ncv_cl2_bias = "NCV–REML widening",
  hybrid = "NCV estimate, REML + CL2 width"
)
hyb_tab_s5 <- hyb_all |>
  mutate(
    v = sprintf("%s (%s) %s / %s", f(cov, 3), f(q05, 2), f(wr, 2), f(sr, 2)),
    arm = factor(HYB_ARMS[arm], HYB_ARMS)
  ) |>
  select(setting, n, arm, v) |>
  arrange(arm) |>
  tidyr::pivot_wider(names_from = arm, values_from = v)
hv5 <- \(st, a, col) hyb_all[[col]][hyb_all$setting == st & hyb_all$arm == a]
# detection of the NCV estimate with REML + CL2 width vs REML + CL2, plasmode (paired reps)
SC_P <- read_pcells("ncv-rbc/summaries/score-pm2-cells.csv")
det_p <- SC_P |>
  filter(estimand == "beta", arm %in% c("hybrid", "reml_cl2")) |>
  group_by(app, arm) |>
  summarise(det = mean(detect), .groups = "drop") |>
  tidyr::pivot_wider(names_from = arm, values_from = det)
det_lower_n <- sum(det_p$hybrid < det_p$reml_cl2)

# (d)-(f) one-step bias correction, smoothing-parameter term, combination
RBC_S <- read_scells("ncv-rbc/summaries/rbc-fresh-cells.csv")
COMBO_S <- read_scells("ncv-rbc/summaries/combo-syn-cells.csv")
syn_grp <- \(d)
  d |>
    mutate(
      group = case_when(
        block == "core" &
          error != "iid" &
          family %in% c("gaussian", "poisson") ~
          "smooth truth, Gaussian/Poisson",
        block == "rough_truth" &
          error != "iid" &
          family %in% c("gaussian", "poisson") ~
          "rough truth, Gaussian/Poisson",
        error != "iid" & family == "binomial" ~ "binary",
        TRUE ~ "independent errors"
      )
    )
rbc_s <- syn_grp(RBC_S) |>
  filter(estimand == "beta", arm %in% c("ncv_cl2", "rbc_p", "reml_cl2")) |>
  group_by(cell) |>
  mutate(wr = halfwidth / halfwidth[arm == "reml_cl2"]) |>
  group_by(group, arm) |>
  summarise(
    n = n(),
    cov = mean(coverage),
    q05 = mean(coverage_q05),
    wr = geo(wr),
    det = mean(detect),
    rel = median(rel_rmse_median),
    .groups = "drop"
  )
rsv <- \(g, a, col) rbc_s[[col]][rbc_s$group == g & rbc_s$arm == a]
# REML + CL2 half-width per synthetic cell from the fresh run (same data as twicing; reps 201-400)
reml_hw_s <- RBC_S |>
  filter(arm == "reml_cl2") |>
  select(cell, estimand, hw_reml = halfwidth)
combo_s <- syn_grp(COMBO_S) |>
  filter(estimand == "beta") |>
  left_join(reml_hw_s, by = c("cell", "estimand")) |>
  mutate(wr = halfwidth / hw_reml) |>
  group_by(group, arm) |>
  summarise(
    n = n(),
    reps = min(n_rep),
    cov = mean(coverage),
    q05 = mean(coverage_q05),
    wr = geo(wr),
    se_sd = median(se_sd),
    det = mean(detect),
    .groups = "drop"
  )
csv_ <- \(g, a, col) combo_s[[col]][combo_s$group == g & combo_s$arm == a]
reml_cov_s <- syn_grp(RBC_S) |>
  filter(estimand == "beta", arm == "reml_cl2", cell %in% COMBO_S$cell) |>
  group_by(group) |>
  summarise(
    cov = mean(coverage),
    q05 = mean(coverage_q05),
    det = mean(detect),
    .groups = "drop"
  )
# the combination was run for Gaussian responses only
combo_s$group <- sub(
  ", Gaussian/Poisson",
  ", Gaussian",
  combo_s$group,
  fixed = TRUE
)
reml_cov_s$group <- sub(
  ", Gaussian/Poisson",
  ", Gaussian",
  reml_cov_s$group,
  fixed = TRUE
)
COMBO_ARMS <- c(
  ncv_cl2 = "NCV + CL2",
  ncv_ijx = "+ λ term",
  rbc_p = "bias correction",
  rbc_p_ijx = "bias correction + λ term"
)
combo_s_tab <- combo_s |>
  filter(arm %in% names(COMBO_ARMS)) |>
  mutate(
    v = sprintf("%s (%s) %s", f(cov, 3), f(q05, 2), f(wr, 2)),
    arm = factor(COMBO_ARMS[arm], COMBO_ARMS)
  ) |>
  select(group, n, arm, v) |>
  arrange(arm) |>
  tidyr::pivot_wider(names_from = arm, values_from = v) |>
  left_join(
    reml_cov_s |>
      transmute(group, `REML + CL2` = sprintf("%s (%s)", f(cov, 3), f(q05, 2))),
    by = "group"
  )
# plasmode: combination (reps 1-50) against REML + CL2 of the bias-correction run (reps 1-100)
reml_p <- RBC_P |>
  filter(estimand == "beta", arm == "reml_cl2") |>
  select(
    cell,
    hw_reml = halfwidth,
    cov_reml = coverage,
    q05_reml = coverage_q05
  )
combo_p <- COMBO_P |>
  filter(estimand == "beta") |>
  left_join(reml_p, by = "cell") |>
  mutate(wr = halfwidth / hw_reml) |>
  group_by(app, role, arm) |>
  summarise(
    cov = mean(coverage),
    q05 = mean(coverage_q05),
    wr = geo(wr),
    se_sd = median(se_sd),
    .groups = "drop"
  )
reml_p_app <- reml_p |>
  left_join(PCELL, by = "cell") |>
  group_by(app) |>
  summarise(cov = mean(cov_reml), q05 = mean(q05_reml), .groups = "drop")
combo_p_tab <- combo_p |>
  filter(arm %in% names(COMBO_ARMS)) |>
  mutate(
    v = sprintf("%s (%s) %s", f(cov, 2), f(q05, 2), f(wr, 2)),
    arm = factor(COMBO_ARMS[arm], COMBO_ARMS)
  ) |>
  select(app, role, arm, v) |>
  arrange(arm) |>
  tidyr::pivot_wider(names_from = arm, values_from = v) |>
  left_join(
    reml_p_app |>
      transmute(app, `REML + CL2` = sprintf("%s (%s)", f(cov, 2), f(q05, 2))),
    by = "app"
  ) |>
  mutate(app = factor(app, names(APP_LABELS))) |>
  arrange(desc(role == "counted"), app) |>
  mutate(dataset = app_lab(as.character(app))) |>
  select(dataset, all_of(unname(COMBO_ARMS)), `REML + CL2`)
cpv <- \(a, col, r = "counted")
  combo_p[[col]][combo_p$arm == a & combo_p$role %in% r]
combo_missing_p <- COMBO_P |>
  filter(estimand == "beta", arm == "rbc_p_ijx") |>
  summarise(miss = sum(50 - n_rep), tasks = 50 * n())
combo_missing_s <- COMBO_S |>
  filter(estimand == "beta", arm == "rbc_p_ijx") |>
  summarise(miss = sum(100 - n_rep), tasks = 100 * n())
combo_truth <- COMBO_P |>
  filter(estimand == "beta", arm == "rbc_p_ijx") |>
  group_by(truth) |>
  summarise(cov = mean(coverage), .groups = "drop")

# (e) NCV smoothing-parameter term alone (Gaussian synthetic, GLM synthetic, plasmode)
NCVL <- readr::read_csv(
  "lambda-correction/summaries-ncv/ncv-cells.csv",
  show_col_types = FALSE
)
NCVG <- readr::read_csv(
  "lambda-correction/summaries-ncvg/ncv-cells.csv",
  show_col_types = FALSE
)
NCVP <- readr::read_csv(
  "lambda-correction/summaries-ncv-plasmode/ncvp-cells.csv",
  show_col_types = FALSE
) |>
  select(-family, -error, -G) |>
  left_join(PCELL, by = "cell")
lam_syn <- bind_rows(NCVL, NCVG) |>
  filter(estimand == "beta", crit == "z", arm %in% c("cl2", "ij_lin", "fix")) |>
  mutate(dependence = dep_lab(error)) |>
  group_by(cell) |>
  mutate(wr = width / width[arm == "cl2"], sr = score / score[arm == "cl2"]) |>
  group_by(family, dependence, G, arm) |>
  summarise(
    cov = mean(coverage),
    q05 = mean(pointwise_q05),
    wr = mean(wr),
    sr = mean(sr),
    se_sd = mean(se_sd_ratio),
    .groups = "drop"
  )
lsv <- \(fam, g, a, col = "cov", dep = "dependent")
  lam_syn[[col]][
    lam_syn$family == fam &
      lam_syn$G == g &
      lam_syn$arm == a &
      lam_syn$dependence == dep
  ]
lam_p <- NCVP |>
  filter(
    estimand == "beta",
    crit == "z",
    arm %in% c("cl2", "ij_lin"),
    G_label == "all"
  ) |>
  group_by(cell) |>
  mutate(wr = width / width[arm == "cl2"]) |>
  group_by(role, arm) |>
  summarise(
    cov = mean(coverage),
    cmin = min(coverage),
    cmax = max(coverage),
    wr = median(wr),
    .groups = "drop"
  )
lpv <- \(a, col, r = "counted") lam_p[[col]][lam_p$arm == a & lam_p$role == r]
LAM_ARMS <- c(
  cl2 = "NCV + CL2",
  ij_lin = "+ λ term",
  fix = "λ fixed (diagnostic)"
)
lam_syn_tab <- lam_syn |>
  filter(dependence == "dependent") |>
  mutate(
    v = sprintf("%s (%s) %s", f(cov, 3), f(q05, 2), f(wr, 2)),
    arm = factor(LAM_ARMS[arm], LAM_ARMS),
    family = fam_lab(family)
  ) |>
  select(family, G, arm, v) |>
  arrange(arm) |>
  tidyr::pivot_wider(names_from = arm, values_from = v)

# ---- S5.2 mgcv's jackknife covariance for curve-blocked NCV --------------------------------
JACK <- readr::read_csv(
  "ncv-jack-recheck/jack-reps.csv",
  show_col_types = FALSE
)
jack_s5 <- JACK |>
  group_by(family, estimand, arm) |>
  summarise(
    reps = n(),
    cov = mean(coverage),
    se = sd(coverage) / sqrt(n()),
    ser = mean(se_ratio_to_bayes),
    .groups = "drop"
  )
jkv <- \(a, col = "cov") jack_s5[[col]][jack_s5$arm == a]

# corrected neighbourhood covariance from the block fit (ncv-jack-recheck/jackfix-lib.R)
JF_DIR <- "ncv-jack-recheck/summaries/"
JF_ARMS <- c(
  ncv_vp = "Bayesian Vp",
  ncv_cl2 = "NCV + CL2",
  jk_nei = "corrected neighbourhood covariance",
  jk_nei_vc = "... with λ̂ correction",
  jk_cv = "its sampling part",
  jk_raw = "mgcv per-point jackknife = TRUE",
  reml_cl2 = "REML + CL2"
)
JF_S <- readr::read_csv(
  paste0(JF_DIR, "jackfix-syn-cells.csv"),
  show_col_types = FALSE
)
JF_SI <- readr::read_csv(
  paste0(JF_DIR, "jackfix-syn-info.csv"),
  show_col_types = FALSE
)
JF_PI <- readr::read_csv(
  paste0(JF_DIR, "jackfix-pm-info.csv"),
  show_col_types = FALSE
)
JF_P <- readr::read_csv(
  paste0(JF_DIR, "jackfix-pm-cells.csv"),
  show_col_types = FALSE
)
JF_PW <- readr::read_csv(
  paste0(JF_DIR, "jackfix-pointwise-check-c12-r1.csv"),
  show_col_types = FALSE
)
JF_CRN <- readr::read_csv(
  paste0(JF_DIR, "jackfix-syn-crn.csv"),
  show_col_types = FALSE
)
JF_FAIL <- dplyr::bind_rows(
  readr::read_csv(
    paste0(JF_DIR, "jackfix-syn-failures.csv"),
    show_col_types = FALSE
  ),
  readr::read_csv(
    paste0(JF_DIR, "jackfix-pm-failures.csv"),
    show_col_types = FALSE
  )
)
jf_cell <- \(cov, q05, sesd, w1, w2)
  sprintf(
    "%s (%s) %s [%s/%s]",
    f(cov, 3),
    f(q05, 2),
    f(sesd, 2),
    f(w1, 2),
    f(w2, 2)
  )
jf_syn_tab <- JF_S |>
  dplyr::filter(arm %in% names(JF_ARMS)) |>
  mutate(
    v = jf_cell(coverage, q05, se_sd, width_vs_ncv_cl2, width_vs_reml_cl2),
    family = recode(
      family,
      gaussian = "Gaussian",
      poisson = "Poisson",
      binomial = "binary"
    ),
    estimand = factor(ESTIMAND_LABELS[estimand], ESTIMAND_LABELS),
    arm = factor(JF_ARMS[arm], JF_ARMS)
  ) |>
  select(estimand, interval = arm, family, v) |>
  tidyr::pivot_wider(names_from = family, values_from = v) |>
  arrange(estimand, interval) |>
  select(estimand, interval, Gaussian, Poisson, binary) |>
  mutate(across(c(estimand, interval), as.character))
jf_sv <- \(a, e, col = "coverage", fam = NULL)
  JF_S[[col]][
    JF_S$arm == a & JF_S$estimand %in% e & (is.null(fam) | JF_S$family %in% fam)
  ]
jf_p_app <- JF_P |>
  group_by(app, estimand, arm) |>
  summarise(
    cells = n(),
    cov = mean(coverage),
    q05 = mean(q05),
    se_sd = median(se_sd),
    w1 = geo(width_vs_ncv_cl2),
    w2 = geo(width_vs_reml_cl2),
    .groups = "drop"
  )
jf_p_tab <- jf_p_app |>
  dplyr::filter(
    estimand %in% c("beta", "mean"),
    arm %in% c("ncv_vp", "ncv_cl2", "jk_nei", "jk_raw", "reml_cl2")
  ) |>
  mutate(
    v = ifelse(
      arm %in% c("jk_nei", "jk_raw"),
      jf_cell(cov, q05, se_sd, w1, w2),
      sprintf("%s (%s) %s", f(cov, 2), f(q05, 2), f(se_sd, 2))
    ),
    estimand = factor(
      ESTIMAND_LABELS[estimand],
      ESTIMAND_LABELS[c("beta", "mean")]
    ),
    arm = factor(JF_ARMS[arm], JF_ARMS)
  ) |>
  select(estimand, app, cells, arm, v) |>
  tidyr::pivot_wider(names_from = arm, values_from = v) |>
  arrange(estimand, app) |>
  mutate(dataset = app_lab(app), estimand = as.character(estimand)) |>
  select(estimand, dataset, combinations = cells, any_of(unname(JF_ARMS)))
jf_p_all <- JF_P |>
  dplyr::filter(estimand == "beta") |>
  group_by(arm) |>
  summarise(
    cells = n(),
    cov = mean(coverage),
    q05 = mean(q05),
    min = min(coverage),
    w1 = geo(width_vs_ncv_cl2),
    w2 = geo(width_vs_reml_cl2),
    .groups = "drop"
  )
jf_pa <- \(a, col = "cov") jf_p_all[[col]][jf_p_all$arm == a]
jf_pb <- jf_p_app |>
  dplyr::filter(estimand == "beta") |>
  select(app, arm, cov) |>
  tidyr::pivot_wider(names_from = arm, values_from = cov)

# ---- S5.3 smoothing-parameter term for REML --------------------------------------------------
LC <- readr::read_csv(
  "lambda-correction/summaries/lc-cells.csv",
  show_col_types = FALSE
)
LCT <- readr::read_csv(
  "lambda-correction/summaries/lc-targeting.csv",
  show_col_types = FALSE
)
LCP <- readr::read_csv(
  "lambda-correction/summaries-plasmode/lcp-cells.csv",
  show_col_types = FALSE
) |>
  select(-family, -error, -G) |>
  left_join(PCELL, by = "cell")
lc_w <- LC |>
  filter(crit == "z", arm %in% c("cl2", "vc", "ij_lin", "fix")) |>
  mutate(dependence = dep_lab(error)) |>
  group_by(cell, estimand) |>
  mutate(wr = width / width[arm == "cl2"], sr = score / score[arm == "cl2"]) |>
  group_by(estimand, dependence, G, arm) |>
  summarise(
    cov = mean(coverage),
    wr = mean(wr),
    sr = mean(sr),
    .groups = "drop"
  )
lcv <- \(a, col = "cov", e = "beta", g = c(40, 100), dep = "dependent")
  lc_w[[col]][
    lc_w$arm == a & lc_w$estimand == e & lc_w$G %in% g & lc_w$dependence == dep
  ]
LC_ARMS <- c(
  cl2 = "CL2",
  vc = "CL2 + mgcv's term",
  ij_lin = "CL2 + curve-robust λ term",
  fix = "λ fixed (diagnostic)"
)
lc_tab <- lc_w |>
  filter(dependence == "dependent", estimand %in% c("mean", "beta")) |>
  mutate(
    v = sprintf("%s %s / %s", f(cov, 3), f(wr, 2), f(sr, 3)),
    arm = factor(LC_ARMS[arm], LC_ARMS),
    estimand = est_lab(estimand)
  ) |>
  select(estimand, G, arm, v) |>
  arrange(arm) |>
  tidyr::pivot_wider(names_from = arm, values_from = v)
lc_target <- LCT |>
  filter(arm == "ij_lin", error != "iid") |>
  group_by(G) |>
  summarise(
    m = mean(spearman),
    lo = min(spearman),
    hi = max(spearman),
    .groups = "drop"
  )
lcp_w <- LCP |>
  filter(crit == "z", arm %in% c("cl2", "ij_lin"), estimand == "beta") |>
  group_by(cell) |>
  mutate(
    gain = coverage - coverage[arm == "cl2"],
    sr = score / score[arm == "cl2"]
  ) |>
  ungroup() |>
  filter(arm == "ij_lin")
lcp_cl2_beta <- LCP |> filter(crit == "z", arm == "cl2", estimand == "beta")

# ---- S7 computing time -----------------------------------------------------------------------
TIM <- read_summary("final", "timing-benchmark")
tim_tab <- TIM |>
  group_by(method) |>
  summarise(
    fit = median(fit_time),
    cl2 = median(cl2_time, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )
CMPF <- readr::read_csv(
  "comparators/summaries/comparators-fits.csv",
  show_col_types = FALSE
)
cmp_time <- CMPF |>
  group_by(arm) |>
  summarise(
    median = median(time),
    q90 = unname(quantile(time, 0.9)),
    n = n(),
    .groups = "drop"
  )
LCTIME <- readr::read_csv(
  "lambda-correction/summaries/lc-times.csv",
  show_col_types = FALSE
)
NCVTIME_P <- readr::read_csv(
  "lambda-correction/summaries-ncv-plasmode/ncvp-times.csv.gz",
  show_col_types = FALSE
)
NCVTIME <- readr::read_csv(
  "lambda-correction/summaries-ncv/ncv-times.csv",
  show_col_types = FALSE
)

# ---- helpers formerly in paper/_sim.R --------------------------------------------------------
conv_ncv <- FITS |>
  filter(fit == "NCV") |>
  group_by(fit_family) |>
  summarise(conv = weighted.mean(converged, n), .groups = "drop")
ar1_vs_cl2 <- S |>
  filter(
    cell %in% ar1_cells,
    estimand %in% c("mean", "beta"),
    arm %in% c("ar1_model", "reml_cl2")
  ) |>
  select(cell, estimand, arm, score_mean) |>
  tidyr::pivot_wider(names_from = arm, values_from = score_mean) |>
  filter(!is.na(ar1_model)) |>
  mutate(ratio = ar1_model / reml_cl2)
wb <- \(a, b, fam = "gaussian", err = "fpc")
  S$width_mean[
    S$cell %in%
      bcells$cell &
      S$basis == b &
      S$family == fam &
      S$error == err &
      S$estimand == "beta" &
      S$arm == a
  ]

# ---- S5.1 (e)-(g): smoothing-parameter term with and without the cross term; misspecified truths
# (ncv-rbc/summaries/combo2-*.csv; objects pm_all, nv_comb from the frozen write-up code)
c2_pm <- pm_all |>
  filter(
    arm %in% c("ncv_ij", "ncv_ijx", "rbc_p_ij", "rbc_p_ijx", "reml_cl2")
  ) |>
  mutate(role = if_else(app %in% APP_COUNTED, "counted", "stress test")) |>
  group_by(role, app, arm) |>
  summarise(
    cov = mean(coverage),
    q05 = mean(coverage_q05),
    se_sd = mean(se_sd),
    .groups = "drop"
  )
c2v <- \(a, col = "cov", r = c("counted", "stress test"))
  c2_pm[[col]][c2_pm$arm == a & c2_pm$role %in% r]
c2_gap <- c2_pm |>
  select(role, app, arm, cov) |>
  tidyr::pivot_wider(names_from = arm, values_from = cov) |>
  mutate(
    cross_gain = ncv_ijx - ncv_ij,
    cross_gain_bc = rbc_p_ijx - rbc_p_ij,
    short = reml_cl2 - rbc_p_ij
  )
c2_fail <- c2_gap |>
  filter(rbc_p_ij < 0.85) |>
  mutate(app = factor(app, names(APP_LABELS))) |>
  arrange(app)
c2_tab <- c2_pm |>
  mutate(
    v = sprintf("%s (%s) %s", f(cov, 2), f(q05, 2), f(se_sd, 2)),
    arm = factor(
      c(
        ncv_ij = "+ λ term",
        ncv_ijx = "+ λ term, cross",
        rbc_p_ij = "bias corr. + λ term",
        rbc_p_ijx = "bias corr. + λ term, cross",
        reml_cl2 = "REML + CL2"
      )[arm],
      c(
        "+ λ term",
        "+ λ term, cross",
        "bias corr. + λ term",
        "bias corr. + λ term, cross",
        "REML + CL2"
      )
    )
  ) |>
  select(role, app, arm, v) |>
  arrange(arm) |>
  tidyr::pivot_wider(names_from = arm, values_from = v) |>
  mutate(app = factor(app, names(APP_LABELS))) |>
  arrange(desc(role == "counted"), app) |>
  transmute(dataset = app_lab(as.character(app)), across(-c(role, app)))
c2_syn <- readr::read_csv(
  "ncv-rbc/summaries/combo2-syn-cells.csv",
  show_col_types = FALSE
) |>
  left_join(DESIGN, by = "cell") |>
  filter(
    estimand == "beta",
    arm %in% c("ncv_ij", "ncv_ijx", "rbc_p_ij", "rbc_p_ijx")
  ) |>
  group_by(truth, arm) |>
  summarise(cov = mean(coverage), se_sd = mean(se_sd), .groups = "drop")
c2s <- \(a, tr, col = "cov") c2_syn[[col]][c2_syn$arm == a & c2_syn$truth == tr]
# misspecified truths, default basis: every NCV-centred interval
MSK <- readr::read_csv(
  "misspec/summaries/misspec-cells.csv",
  show_col_types = FALSE
) |>
  distinct(cell, truth_label, error, G, basis)
c2_ms <- readr::read_csv(
  "ncv-rbc/summaries/combo2-ms-cells.csv",
  show_col_types = FALSE
) |>
  left_join(MSK, by = "cell") |>
  filter(estimand == "beta")
c2m <- \(tr, err, col = "coverage", arms = unique(c2_ms$arm))
  c2_ms[[col]][
    c2_ms$truth_label == tr & c2_ms$error %in% err & c2_ms$arm %in% arms
  ]
c2_ms_n_arms <- n_distinct(c2_ms$arm)

# ---- round 2: SE/SD of GLS-FPCA + CL2 on the plasmode data (replicates 1-100)
CA_se <- CA |>
  filter(source == "plasmode", estimand == "beta", arm == "gls_fpc_cl2")

# ---- 20 curves, zero surface and Satterthwaite beyond the main factorial ---------------------
# design blocks incl. the G = 20 and beta = 0 blocks (Appendix B.1)
blocks <- bind_rows(DESIGN, DESIGN_X) |>
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
      "small_G",
      "warp",
      "dense_grid",
      "rough_truth",
      "null_beta",
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

# G ladder: coverage (mean over settings), relative error and detection (medians over settings)
# per family, dependence, estimand and interval, main factorial at G = 20, 40 and 100
lad_keys <- SGL |> distinct(cell, family, signal, error, G)
lad_rel <- bind_rows(
  read_summary("final", "relative-error") |> filter(cell %in% lad_keys$cell),
  SGREL |> filter(cell %in% lad_keys$cell)
) |>
  inner_join(lad_keys, by = "cell")
lad_cov <- SGL |>
  filter(arm %in% MAIN_ARMS) |>
  group_by(family, error = dep_lab(error), estimand, arm, G) |>
  summarise(cov = mean(coverage_mean), .groups = "drop")
lad_re <- lad_rel |>
  filter(arm %in% MAIN_ARMS) |>
  group_by(family, error = dep_lab(error), estimand, arm, G) |>
  summarise(
    rel = median(rel_rmse_median),
    det = median(detect),
    .groups = "drop"
  )
ladder_tab <- lad_cov |>
  left_join(lad_re, by = c("family", "error", "estimand", "arm", "G")) |>
  pivot_wider(names_from = G, values_from = c(cov, rel, det)) |>
  mutate(
    family = factor(
      fam_lab(family),
      FAMILY_LABELS[c("gaussian", "poisson", "binomial")]
    ),
    estimand = factor(
      est_lab(estimand),
      ESTIMAND_LABELS[c("mean", "beta", "alpha", "gamma")]
    ),
    arm = factor(ARM_LABELS[arm], ARM_LABELS[MAIN_ARMS])
  ) |>
  arrange(family, desc(error), estimand, arm) |>
  mutate(across(c(family, estimand, arm), as.character)) |>
  select(
    family,
    errors = error,
    estimand,
    interval = arm,
    cov_20,
    cov_40,
    cov_100,
    rel_20,
    rel_40,
    rel_100,
    det_20,
    det_40,
    det_100
  )
# critical values by G (default signal): coverage per family, dependence and estimand
lad_satt_tab <- SGLS |>
  filter(crit %in% c("z", "t", "satt")) |>
  group_by(family, error = dep_lab(error), estimand, crit, G) |>
  summarise(
    cov = mean(coverage),
    q05 = mean(pointwise_q05),
    width = mean(width),
    .groups = "drop"
  ) |>
  group_by(family, error, estimand, G) |>
  mutate(wr = width / width[crit == "z"]) |>
  ungroup()
lst <- \(cr, g, col = "cov") {
  d <- lad_satt_tab |>
    filter(crit == cr, G == g) |>
    select(family, error, estimand, v = all_of(col))
  names(d)[4] <- paste(cr, g, col)
  d
}
ladder_satt_tab <- lst("z", 20) |>
  left_join(lst("t", 20), by = c("family", "error", "estimand")) |>
  left_join(lst("satt", 20), by = c("family", "error", "estimand")) |>
  left_join(lst("satt", 20, "q05"), by = c("family", "error", "estimand")) |>
  left_join(lst("satt", 20, "wr"), by = c("family", "error", "estimand")) |>
  left_join(lst("z", 40), by = c("family", "error", "estimand")) |>
  left_join(lst("satt", 40), by = c("family", "error", "estimand")) |>
  left_join(lst("z", 100), by = c("family", "error", "estimand")) |>
  left_join(lst("satt", 100), by = c("family", "error", "estimand")) |>
  mutate(
    family = factor(
      fam_lab(family),
      FAMILY_LABELS[c("gaussian", "poisson", "binomial")]
    ),
    estimand = factor(
      est_lab(estimand),
      ESTIMAND_LABELS[c("mean", "beta", "alpha", "gamma")]
    )
  ) |>
  arrange(family, desc(error), estimand) |>
  mutate(across(c(family, estimand), as.character)) |>
  rename(errors = error)
lsv <- \(
  cr,
  g,
  dep = "dependent",
  e = c("mean", "beta", "gamma"),
  fam = c("gaussian", "poisson", "binomial"),
  col = "cov"
)
  lad_satt_tab[[col]][
    lad_satt_tab$crit == cr &
      lad_satt_tab$G == g &
      lad_satt_tab$error == dep &
      lad_satt_tab$estimand %in% e &
      lad_satt_tab$family %in% fam
  ]
lcv <- \(
  a,
  g,
  dep = "dependent",
  e = c("mean", "beta", "gamma"),
  fam = c("gaussian", "poisson", "binomial"),
  col = "cov"
) {
  d <- lad_cov |>
    left_join(lad_re, by = c("family", "error", "estimand", "arm", "G"))
  d[[col]][
    d$arm == a &
      d$G == g &
      d$error == dep &
      d$estimand %in% e &
      d$family %in% fam
  ]
}

# Satterthwaite in the extension blocks: coverage with z and Satterthwaite per block, G and
# estimand (mean over the block's settings; f(x,t) with beta), and the mean width ratio
SATT_BLOCKS <- c(
  warp = "misregistration",
  dense_grid = "dense grid",
  rough_truth = "rough truth",
  term_type = "smooth effect",
  families = "extended families",
  ar1_home = "exact AR(1)",
  oscillating = "sign-changing",
  warp_ar1 = "misreg. (AR(1) fit)",
  heteroskedastic = "heteroskedastic",
  lowrank_covariate = "low-rank covariate",
  basis_size = "larger bases"
)
satx_est <- \(e)
  c(mean = "mean", beta = "beta", f = "beta", alpha = "alpha", gamma = "gamma")[
    e
  ]
satx_e <- SATX |>
  filter(!block %in% c("small_G", "null_beta"), crit %in% c("z", "satt")) |>
  mutate(est = satx_est(estimand))
satx_tab <- satx_e |>
  group_by(block, G) |>
  mutate(settings = n_distinct(cell)) |>
  group_by(block, G, settings, est, crit) |>
  summarise(cov = mean(coverage), .groups = "drop") |>
  mutate(col = paste(est, crit)) |>
  select(-est, -crit) |>
  pivot_wider(names_from = col, values_from = cov) |>
  left_join(
    SATXC |>
      filter(!block %in% c("small_G", "null_beta"), contrast == "Satt - z") |>
      group_by(block, G) |>
      summarise(width = mean(width_ratio), .groups = "drop"),
    by = c("block", "G")
  ) |>
  mutate(block = factor(SATT_BLOCKS[block], SATT_BLOCKS)) |>
  arrange(block, G) |>
  mutate(block = as.character(block)) |>
  select(
    block,
    G,
    settings,
    `mean z`,
    `mean satt`,
    `beta z`,
    `beta satt`,
    `alpha z`,
    `alpha satt`,
    `gamma z`,
    `gamma satt`,
    width
  )
satx_low_tab <- satx_e |>
  select(
    cell,
    block,
    family,
    error,
    G,
    signal,
    truth,
    estimand,
    crit,
    coverage,
    pointwise_q05
  ) |>
  pivot_wider(names_from = crit, values_from = c(coverage, pointwise_q05)) |>
  filter(coverage_satt < LEVEL - ADEQUATE) |>
  mutate(block = factor(SATT_BLOCKS[block], SATT_BLOCKS)) |>
  arrange(block, family, G, estimand) |>
  transmute(
    block = as.character(block),
    family = fam_lab(family),
    errors = unname(SHORT_ERROR[error]),
    G,
    signal,
    truth = recode(truth, wiggly = "rough"),
    estimand = est_lab(estimand),
    z = coverage_z,
    Satt. = coverage_satt,
    `q05 Satt.` = pointwise_q05_satt
  )
satx_con <- SATXC |>
  filter(!block %in% c("small_G", "null_beta"), contrast == "Satt - z")
# per-point df relative to G - 1, all synthetic settings and plasmode, by G and estimand
satt_df_tab <- bind_rows(
  synthetic = bind_rows(SATD, SATXD) |> mutate(G = as.character(G)),
  plasmode = SATPD |> mutate(G = G_label),
  .id = "data"
) |>
  mutate(est = satx_est(estimand)) |>
  group_by(data, G, est) |>
  summarise(
    settings = n(),
    df_lo = min(df_median),
    df_hi = max(df_median),
    rel_lo = min(df_median_rel),
    rel_med = median(df_median_rel),
    rel_hi = max(df_median_rel),
    .groups = "drop"
  ) |>
  arrange(
    data,
    factor(G, c("20", "40", "100", "all")),
    factor(est, c("mean", "beta", "alpha", "gamma"))
  ) |>
  mutate(est = est_lab(est))

# plasmode: z and Satterthwaite per dataset and number of subjects (mean over truth x residual)
satp_tab <- satp_app |>
  filter(crit %in% c("z", "satt")) |>
  mutate(col = paste(estimand, crit)) |>
  select(app, G_label, col, cov) |>
  pivot_wider(names_from = col, values_from = cov) |>
  left_join(
    satp_app |>
      filter(crit == "satt", estimand == "beta") |>
      select(app, G_label, width = wr),
    by = c("app", "G_label")
  ) |>
  left_join(SATP |> distinct(app, G_label, G), by = c("app", "G_label")) |>
  mutate(app = factor(APP_LABELS[app], APP_LABELS)) |>
  arrange(app, G) |>
  mutate(app = as.character(app)) |>
  select(
    dataset = app,
    G,
    `mean z`,
    `mean satt`,
    `beta z`,
    `beta satt`,
    `alpha z`,
    `alpha satt`,
    `gamma z`,
    `gamma satt`,
    width
  )

# zero surface: all intervals
NULL_ARMS <- c(
  reml_model = "REML model",
  reml_cl2 = "REML+CL2",
  reml_cl2_t = "REML+CL2, t(G-1)",
  reml_cl2_satt = "REML+CL2, Satt.",
  ncv_model = "NCV model",
  ncv_cl2 = "NCV+CL2",
  ncv_cl2_bias = "NCV+CL2 widened"
)
null_s <- SGN |>
  mutate(
    arm_l = factor(NULL_ARMS[arm], NULL_ARMS),
    fam = factor(fam_lab(family), FAMILY_LABELS[c("gaussian", "binomial")]),
    err = factor(unname(SHORT_ERROR[error]), SHORT_ERROR[c("fpc", "iid")])
  ) |>
  arrange(fam, err, G, arm_l)
null_fpr_tab <- null_s |>
  mutate(v = sprintf("%s (%s)", f(fpr, 3), f(fpr_se, 3))) |>
  select(family = fam, errors = err, G, arm_l, v) |>
  pivot_wider(names_from = arm_l, values_from = v) |>
  mutate(across(c(family, errors), as.character))
null_desc_tab <- null_s |>
  filter(arm %in% c("reml_cl2", "reml_cl2_satt", "ncv_cl2")) |>
  mutate(v = sprintf("%s / %s", f(any_excl, 2), f(share_gt_5pct, 2))) |>
  select(family = fam, errors = err, G, arm_l, v) |>
  pivot_wider(names_from = arm_l, values_from = v) |>
  left_join(
    null_s |>
      filter(arm %in% c("reml_cl2", "ncv_cl2")) |>
      mutate(fit = if_else(arm == "reml_cl2", "REML", "NCV")) |>
      select(family = fam, errors = err, G, fit, rms, rms_rel) |>
      pivot_wider(names_from = fit, values_from = c(rms, rms_rel)),
    by = c("family", "errors", "G")
  ) |>
  mutate(across(c(family, errors), as.character)) |>
  select(
    family,
    errors,
    G,
    `REML+CL2`,
    `REML+CL2, Satt.`,
    `NCV+CL2`,
    `RMS REML` = rms_REML,
    `RMS NCV` = rms_NCV,
    `rel. REML` = rms_rel_REML,
    `rel. NCV` = rms_rel_NCV
  )
nsv <- \(
  a,
  err = "fpc",
  fam = c("gaussian", "binomial"),
  g = c(40, 100),
  col = "fpr"
)
  SGN[[col]][
    SGN$arm == a & SGN$error %in% err & SGN$family %in% fam & SGN$G %in% g
  ]
# the other estimands in the zero-surface settings, and REML + CL2 in the matching non-zero setting
null_keys <- DESIGN_X |> filter(block == "null_beta")
null_core <- S |>
  filter(
    block == "core",
    truth == "smooth",
    arm == "reml_cl2",
    estimand %in% c("alpha", "gamma")
  ) |>
  select(family, error, G, signal, estimand, core = coverage_mean)
null_other <- SG |>
  filter(block == "null_beta", estimand != "beta", arm %in% MAIN_ARMS) |>
  select(family, error, G, signal, estimand, arm, coverage_mean) |>
  pivot_wider(names_from = arm, values_from = coverage_mean) |>
  left_join(null_core, by = c("family", "error", "G", "signal", "estimand"))
null_other_tab <- null_other |>
  mutate(
    fam = factor(fam_lab(family), FAMILY_LABELS[c("gaussian", "binomial")]),
    err = factor(unname(SHORT_ERROR[error]), SHORT_ERROR[c("fpc", "iid")]),
    est = factor(
      est_lab(estimand),
      ESTIMAND_LABELS[c("mean", "alpha", "gamma")]
    )
  ) |>
  arrange(fam, err, G, est) |>
  transmute(
    family = as.character(fam),
    errors = as.character(err),
    G,
    estimand = as.character(est),
    `REML model` = reml_model,
    `REML+CL2` = reml_cl2,
    `NCV+CL2` = ncv_cl2,
    `NCV+CL2 widened` = ncv_cl2_bias,
    `REML+CL2, non-zero β` = core
  )

# smallest Satterthwaite - z difference anywhere (synthetic and plasmode)
satt_min_gain <- min(c(
  SATXC$diff[SATXC$contrast == "Satt - z"],
  SATPC$diff[SATPC$contrast == "Satt - z"],
  readr::read_csv(
    "satterthwaite/summaries/satt-contrasts.csv",
    show_col_types = FALSE
  ) |>
    filter(contrast == "Satt - z") |>
    pull(diff)
))
# plasmode combinations (per dataset and number of subjects) below the adequacy margin
sp_low <- satp_app |>
  filter(crit == "satt", cov < LEVEL - ADEQUATE) |>
  mutate(
    est = c(
      mean = "the mean",
      beta = "$\\beta$",
      alpha = "$\\alpha(t)$",
      gamma = "$\\gamma(t)$"
    )[estimand],
    est = factor(
      est,
      c("$\\beta$", "$\\gamma(t)$", "$\\alpha(t)$", "the mean")
    ),
    app = factor(APP_LABELS[app], APP_LABELS)
  ) |>
  arrange(app, est, G_label) |>
  group_by(app, est) |>
  summarise(
    txt = sprintf(
      "%s of the %s data (%s)",
      first(est),
      first(app),
      paste(
        sprintf(
          "%s with %s",
          f(cov, 3),
          if_else(G_label == "all", "all subjects", "40 subjects")
        ),
        collapse = ", "
      )
    ),
    .groups = "drop"
  )
sp_low_txt <- paste(sp_low$txt, collapse = "; ")

# Satterthwaite on the plasmode data (curve flips): mean over datasets of the per-dataset means
spa <- \(cr, g, e = "beta", col = "cov", apps = names(APP_LABELS))
  satp_app[[col]][
    satp_app$crit == cr &
      satp_app$G_label == g &
      satp_app$estimand == e &
      satp_app$app %in% apps
  ]
spm <- \(cr, g, e = "beta", col = "cov") mean(spa(cr, g, e, col))
