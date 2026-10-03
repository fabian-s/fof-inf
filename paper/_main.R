# Derived objects for paper/paper.qmd (sourced after paper/_common.R; the repository root is the working
# directory). Reads summary CSVs only; writes nothing.

dep_lab <- \(e) if_else(e == "iid", "independent", "dependent")
MAIN3 <- c("reml_model", "reml_cl2", "ncv_cl2")
LAB3 <- c(
  reml_model = "model-based (REML)",
  reml_cl2 = "REML + CL2",
  ncv_cl2 = "NCV + CL2"
)
COL3 <- c(
  "model-based (REML)" = "#999999",
  "REML + CL2" = "#0072B2",
  "NCV + CL2" = "#D55E00"
)
ERR3 <- c(iid = "independent", ou = "rough stationary", fpc = "smooth")

# ---- Section 3: design -----------------------------------------------------------------------
n_cells <- nrow(DESIGN) + nrow(DESIGN_X)
n_rep <- unique(S$n_rep)
n_tasks <- sum(FITS$n[FITS$fit == "REML"]) +
  sum(distinct(SG, cell, n_rep)$n_rep)
BLOCK_LABELS <- c(
  core = "main factorial",
  warp = "misregistration",
  dense_grid = "dense grid (D = 241)",
  rough_truth = "rough truth",
  term_type = "smooth effect f(x,t) instead of ff",
  families = "scaled t, beta, negative binomial",
  ar1_home = "exact AR(1) errors",
  oscillating = "sign-changing correlation",
  warp_ar1 = "misregistration, AR(1) fit",
  heteroskedastic = "heteroskedastic errors",
  lowrank_covariate = "low-rank functional covariate",
  basis_size = "larger bases"
)
# the G = 20 block after the main factorial, the zero surface after the rough truth
BLOCK_LABELS <- append(BLOCK_LABELS, c(small_G = "main factorial, G = 20"), 1)
BLOCK_LABELS <- append(
  BLOCK_LABELS,
  c(null_beta = "zero surface (β ≡ 0)"),
  which(names(BLOCK_LABELS) == "rough_truth")
)
design_tab <- bind_rows(DESIGN, DESIGN_X) |>
  group_by(block) |>
  summarise(
    settings = n(),
    families = paste(unique(FAMILY_LABELS[family]), collapse = ", "),
    errors = paste(unique(SHORT_ERROR[error]), collapse = ", "),
    G = paste(sort(unique(G)), collapse = ", "),
    .groups = "drop"
  ) |>
  mutate(
    errors = gsub("OU", "rough", gsub("iid", "indep.", errors)),
    families = gsub("negative binomial", "neg. binomial", families),
    block = factor(BLOCK_LABELS[block], BLOCK_LABELS)
  ) |>
  arrange(block)

# ---- Section 4.1: intervals ------------------------------------------------------------------
fig1 <- S |>
  filter(
    block == "core",
    signal == "mid",
    estimand %in% c("mean", "beta"),
    arm %in% MAIN3
  ) |>
  mutate(
    arm = factor(LAB3[arm], LAB3),
    estimand = label_estimand(estimand),
    error = factor(ERR3[error], rev(ERR3)),
    family = label_family(family),
    G = paste0("G = ", G)
  )
# dense grid (D = 241 minus D = 61, beta, G = 100), dependent errors
gd <- \(v, err = c("ou", "fpc"))
  dense$estimate[dense$value == v & dense$error %in% err]
grid_cov <- \(a, d, dep = TRUE) {
  S$coverage_mean[
    ((S$block == "core" & S$G == 100 & S$signal == "mid") |
      S$block == "dense_grid") &
      S$estimand == "beta" &
      S$arm == a &
      S$D == d &
      (S$error != "iid") == dep
  ]
}
# iid cost of CL2 relative to model-based (G = 100), core cells
iid_cost100 <- Q1$estimate[Q1$name == "iid_cost_reml"]
# binary alpha/gamma, misregistered counts
ag <- S |>
  filter(block == "core", estimand %in% c("alpha", "gamma"), arm == "reml_cl2")
ag_bin <- \(dep, g = c(40, 100)) {
  ag$coverage_mean[
    ag$family == "binomial" & (ag$error != "iid") == dep & ag$G %in% g
  ]
}
warp_pois <- S |>
  filter(block == "warp", family == "poisson", estimand == "mean")
wp <- \(a) warp_pois$coverage_mean[warp_pois$arm == a]
# lower tails (5% quantile of pointwise coverage), core cells
tail_core <- S |>
  filter(block == "core", estimand %in% c("mean", "beta"), arm %in% MAIN3) |>
  mutate(dependence = dep_lab(error))
tq <- \(a, g, dep = "dependent", e = c("mean", "beta")) {
  tail_core$pointwise_q05[
    tail_core$arm == a &
      tail_core$G == g &
      tail_core$dependence == dep &
      tail_core$estimand %in% e
  ]
}
# G = 40 vs 100 (core, dependent): mean loss per family x estimand
g_loss <- S |>
  filter(
    block == "core",
    estimand %in% c("mean", "beta"),
    arm %in% MAIN3,
    error != "iid"
  ) |>
  select(family, error, signal, estimand, arm, G, coverage_mean) |>
  pivot_wider(
    names_from = G,
    values_from = coverage_mean,
    names_prefix = "G"
  ) |>
  group_by(family, estimand, arm) |>
  summarise(loss = mean(G100 - G40), .groups = "drop")
gl <- \(a) g_loss$loss[g_loss$arm == a]

# Satterthwaite degrees of freedom (satterthwaite/): REML + CL2 (Bayesian, full block) with
# z, t_{G-1} and per-point Satterthwaite critical values, core cells at the default signal
SAT <- readr::read_csv(
  "satterthwaite/summaries/satt-cells.csv",
  show_col_types = FALSE
)
SATD <- readr::read_csv(
  "satterthwaite/summaries/satt-df.csv",
  show_col_types = FALSE
)
# with the G = 20 settings at the default signal
m_sat_w <- bind_rows(SAT, SGLS |> filter(G == 20)) |>
  filter(crit %in% c("z", "satt")) |>
  mutate(dependence = dep_lab(error)) |>
  group_by(dependence, G, family, estimand, crit) |>
  summarise(
    cov = mean(coverage),
    q05 = mean(pointwise_q05),
    width = mean(width),
    .groups = "drop"
  ) |>
  pivot_wider(names_from = crit, values_from = c(cov, q05, width)) |>
  mutate(wr = width_satt / width_z)
m_sat_tab <- m_sat_w |>
  filter(dependence == "dependent") |>
  mutate(
    family = factor(
      FAMILY_LABELS[family],
      FAMILY_LABELS[c("gaussian", "poisson", "binomial")]
    ),
    estimand = factor(
      ESTIMAND_LABELS[estimand],
      ESTIMAND_LABELS[c("mean", "beta", "alpha", "gamma")]
    )
  ) |>
  select(G, family, estimand, cov_z, cov_satt, wr) |>
  pivot_wider(
    names_from = G,
    values_from = c(cov_z, cov_satt, wr)
  ) |>
  arrange(family, estimand) |>
  select(
    family,
    estimand,
    cov_z_20,
    cov_satt_20,
    wr_20,
    cov_z_40,
    cov_satt_40,
    wr_40,
    cov_z_100,
    cov_satt_100,
    wr_100
  )
sv <- \(
  col,
  g,
  dep = "dependent",
  e = c("mean", "beta", "gamma"),
  fam = c("gaussian", "poisson", "binomial")
) {
  m_sat_w[[col]][
    m_sat_w$dependence == dep &
      m_sat_w$G == g &
      m_sat_w$estimand %in% e &
      m_sat_w$family %in% fam
  ]
}
sat_df_med <- SATD |>
  group_by(G) |>
  summarise(lo = min(df_median), hi = max(df_median), .groups = "drop")

# ---- Section 4.2-4.3: estimation and the trade-off -------------------------------------------
rel_sum <- REL |>
  filter(estimand %in% c("beta", "mean", "gamma", "alpha"), arm %in% MAIN3) |>
  group_by(estimand, dependence, arm) |>
  summarise(
    settings = n(),
    coverage = median(coverage),
    rel = median(rel_rmse_median),
    rel_max = max(rel_rmse_median),
    gt1 = mean(share_rel_gt_1),
    hw = median(rel_halfwidth_median),
    det = median(detect),
    wrong = max(wrong_sign),
    .groups = "drop"
  )
# 5% quantile of pointwise coverage per setting, joined from cells.csv
q05_dep <- S |>
  filter(estimand %in% c("beta", "mean", "gamma", "alpha"), arm %in% MAIN3) |>
  mutate(dependence = dep_lab(error)) |>
  group_by(estimand, dependence, arm) |>
  summarise(q05 = median(pointwise_q05), .groups = "drop")
rel_sum <- left_join(rel_sum, q05_dep, by = c("estimand", "dependence", "arm"))
rs <- \(e, dep, a, col)
  rel_sum[[col]][
    rel_sum$estimand == e & rel_sum$dependence == dep & rel_sum$arm == a
  ]
trade_tab <- rel_sum |>
  filter(dependence == "dependent") |>
  mutate(
    estimand = factor(
      ESTIMAND_LABELS[estimand],
      ESTIMAND_LABELS[c("beta", "mean", "alpha", "gamma")]
    ),
    arm = factor(LAB3[arm], LAB3)
  ) |>
  arrange(estimand, arm) |>
  select(estimand, arm, settings, rel, coverage, q05, hw, det)
# REML / NCV ratio of the relative error of beta, main factorial
rel_core <- REL |>
  filter(
    block == "core",
    estimand == "beta",
    arm %in% c("reml_cl2", "ncv_cl2")
  ) |>
  select(family, signal, error, G, arm, rel_rmse_median) |>
  pivot_wider(names_from = arm, values_from = rel_rmse_median) |>
  mutate(ratio = reml_cl2 / ncv_cl2)
rr <- \(dep) rel_core$ratio[(rel_core$error != "iid") == dep]
# REML / NCV ratio of the relative error per estimand under dependence: main factorial for the
# mean, alpha and gamma, term-type block for the smooth effect f(x, t)
rr_est <- function(est, blk = if (est == "f") "term_type" else "core") {
  REL |>
    filter(
      block == blk,
      estimand == est,
      error != "iid",
      arm %in% c("reml_cl2", "ncv_cl2")
    ) |>
    select(cell, arm, rel_rmse_median) |>
    pivot_wider(names_from = arm, values_from = rel_rmse_median) |>
    mutate(ratio = reml_cl2 / ncv_cl2) |>
    pull(ratio)
}
# model-based (REML) and REML + CL2 coverage of every estimand under dependence, G = 100
# (main factorial and term-type block)
cov_all_dep <- S |>
  filter(block %in% c("core", "term_type"), G == 100, error != "iid")
# per-setting trade-off data for the figure (dependent errors)
fig2 <- REL |>
  filter(
    estimand == "beta",
    arm %in% c("reml_cl2", "ncv_cl2"),
    error != "iid"
  ) |>
  select(
    cell,
    block,
    family,
    truth,
    arm,
    rel_rmse_median,
    coverage,
    detect,
    rel_halfwidth_median
  )
fig2a <- fig2 |>
  select(cell, family, truth, arm, rel_rmse_median) |>
  pivot_wider(names_from = arm, values_from = rel_rmse_median) |>
  mutate(
    family = label_family(family),
    truth = if_else(truth == "wiggly", "rough truth", "smooth truth")
  )
fig2b <- fig2 |>
  mutate(
    arm = factor(LAB3[arm], LAB3[c("reml_cl2", "ncv_cl2")]),
    truth = if_else(truth == "wiggly", "rough truth", "smooth truth")
  )
# EDF of the ff term: REML / NCV, core cells at G = 100 and default signal
edf_core <- EDF |>
  filter(
    block == "core",
    G == 100,
    signal == "mid",
    fit %in% c("REML", "NCV"),
    term == "te(X.smat,X.tmat):L.X"
  ) |>
  select(family, error, fit, edf_mean) |>
  pivot_wider(names_from = fit, values_from = edf_mean) |>
  mutate(ratio = REML / NCV)
er <- \(dep) edf_core$ratio[(edf_core$error != "iid") == dep]
# squared-bias share of the MSE of beta, core cells at G = 100, dependent errors
bv <- S |>
  filter(
    block %in% c("core", "rough_truth"),
    G == 100,
    estimand == "beta",
    arm %in% c("reml_cl2", "ncv_cl2"),
    error != "iid"
  ) |>
  mutate(
    fit = if_else(arm == "reml_cl2", "REML", "NCV"),
    share = pmax(bias2, 0) / mse
  )
bshare <- \(f, tr) bv$share[bv$fit == f & bv$truth == tr]
# basis size: EDF of the ff term and REML + CL2 width for beta, Gaussian/Poisson/binary, smooth errors
bcells <- DESIGN |>
  filter(
    block == "basis_size" |
      (block == "core" &
        G == 100 &
        signal == "mid" &
        error %in% c("iid", "fpc"))
  )
edf_b <- EDF |>
  filter(
    cell %in% bcells$cell,
    term == "te(X.smat,X.tmat):L.X",
    fit %in% c("REML", "NCV")
  )
ed <- \(f, b, err = "fpc")
  edf_b$edf_mean[edf_b$fit == f & edf_b$basis == b & edf_b$error == err]
hw_b <- REL |>
  filter(
    cell %in% bcells$cell,
    estimand == "beta",
    arm %in% c("reml_cl2", "ncv_cl2"),
    error == "fpc"
  )
hwb <- \(a, b, col = "rel_halfwidth_median")
  hw_b[[col]][hw_b$arm == a & hw_b$basis == b]
cov_bb <- \(a, b, err = "fpc")
  S$coverage_mean[
    S$cell %in%
      bcells$cell &
      S$basis == b &
      S$error == err &
      S$estimand == "beta" &
      S$arm == a
  ]
# NCV + CL2: coverage with the bias removed (bias estimated from the same replicates)
ND <- read_summary("final", "ncv-diag-cells") |> left_join(DESIGN, by = "cell")
nd <- \(blk, col, iv = "ncv_cl2", g = 100) {
  ND[[col]][
    ND$block == blk &
      ND$interval == iv &
      ND$estimand == "beta" &
      ND$error != "iid" &
      ND$G %in% g
  ]
}
# rough truth: NCV + CL2 coverage and lower tail, smooth errors, G = 100
rough_ncv <- S |>
  filter(block == "rough_truth", error == "fpc", estimand == "beta")
rn <- \(a, col = "coverage_mean") rough_ncv[[col]][rough_ncv$arm == a]
# binary NCV + CL2 at G = 40 (dependent)
bin40_ncv <- S$coverage_mean[
  S$block == "core" &
    S$family == "binomial" &
    S$G == 40 &
    S$error != "iid" &
    S$estimand == "beta" &
    S$arm == "ncv_cl2"
]

# ---- Section 5: plasmode ---------------------------------------------------------------------
n_pcells <- P |> distinct(app, truth, residual, flip, error, G) |> nrow()
DS_TAB <- DATASETS |>
  mutate(app = as.character(app)) |>
  transmute(
    dataset = APP_LABELS[app],
    role = if_else(role == "counted", "independent", "stress test"),
    response = c(
      ecg_gls4b = "myocardial strain",
      af = "atrial strain (day 30)",
      running = "knee moment",
      dti = "FA, corticospinal tract",
      gait_vc = "knee moment",
      ecg_roahd = "ECG lead V5",
      ocean = "dissolved oxygen",
      aemet = "log precipitation",
      adelaide = "electricity demand"
    )[app],
    covariate = c(
      ecg_gls4b = "ECG of the same beat",
      af = "atrial strain (day 0)",
      running = "knee angle",
      dti = "FA, corpus callosum",
      gait_vc = "knee angle",
      ecg_roahd = "ECG lead I",
      ocean = "temperature",
      aemet = "temperature",
      adelaide = "temperature"
    )[app],
    z = c(
      ecg_gls4b = "age",
      af = "treatment",
      running = "body mass",
      dti = "MS status",
      gait_vc = "age",
      ecg_roahd = "LBBB",
      ocean = "cast order",
      aemet = "altitude",
      adelaide = "week"
    )[app],
    G,
    D,
    DE,
    rank_X
  )
pool <- POOL_MEANS |> filter(arm %in% MAIN3, estimand %in% c("mean", "beta"))
pm <- \(a, e, r = "counted", col = "coverage")
  pool[[col]][pool$arm == a & pool$estimand %in% e & pool$role == r]
fig3 <- pool |>
  mutate(
    arm = factor(LAB3[arm], LAB3),
    estimand = label_estimand(estimand),
    app = factor(APP_LABELS[app], rev(APP_LABELS)),
    role = factor(
      if_else(
        role == "counted",
        "indep. source",
        "dep. source"
      ),
      c("indep. source", "dep. source")
    )
  )
g40v <- \(a, e = c("mean", "beta"), r = "counted")
  g40_app$loss[g40_app$arm == a & g40_app$estimand %in% e & g40_app$role == r]
# REML + CL2 lower tail by truth source (beta, G = all, model residuals)
ct <- \(tr, r = "counted", col = "q05")
  bt("reml_cl2", tr, col = col, role_sel = r)
# relative error, half-width and detection per dataset, beta, all subjects, curve flips
prel <- PREL |>
  filter(
    flip == "curve",
    G_label == "all",
    error == "attached",
    arm %in% MAIN3,
    residual %in% c("REML", "NCV")
  )
prel_app <- prel |>
  group_by(role, app, estimand, arm) |>
  summarise(
    coverage = mean(coverage),
    rel = median(rel_rmse_median),
    hw = median(rel_halfwidth_median),
    det = median(detect),
    .groups = "drop"
  )
pv <- \(e, a, col, r = c("counted", "stress test"))
  prel_app[[col]][
    prel_app$estimand == e & prel_app$arm == a & prel_app$role %in% r
  ]
prel_truth <- prel |>
  filter(estimand == "beta", arm %in% c("reml_cl2", "ncv_cl2")) |>
  group_by(role, app, truth, arm) |>
  summarise(rel = median(rel_rmse_median), .groups = "drop") |>
  pivot_wider(names_from = arm, values_from = rel) |>
  mutate(ratio = ncv_cl2 / reml_cl2)
ptr <- \(tr, r = c("counted", "stress test"))
  prel_truth$ratio[prel_truth$truth == tr & prel_truth$role %in% r]
fig4 <- prel |>
  filter(estimand == "beta", arm %in% c("reml_cl2", "ncv_cl2")) |>
  group_by(role, app, truth, arm) |>
  summarise(rel = median(rel_rmse_median), .groups = "drop") |>
  mutate(
    fit = factor(if_else(arm == "reml_cl2", "REML", "NCV"), c("REML", "NCV")),
    truth = factor(
      c(REML = "REML-fitted", MID = "intermediate", NCV = "NCV-fitted")[truth],
      c("REML-fitted", "intermediate", "NCV-fitted")
    ),
    app = factor(APP_LABELS[app], rev(APP_LABELS))
  )
ptab <- prel_app |>
  filter(estimand == "beta") |>
  select(role, app, arm, coverage, rel, hw, det) |>
  pivot_wider(names_from = arm, values_from = c(coverage, rel, hw, det)) |>
  left_join(
    POOL_MEANS |>
      filter(arm == "reml_cl2", estimand == "beta") |>
      select(app, q05 = q05_mean),
    by = "app"
  ) |>
  mutate(dataset = factor(APP_LABELS[app], APP_LABELS)) |>
  arrange(desc(role == "counted"), dataset) |>
  select(
    dataset,
    rel_reml_cl2,
    rel_ncv_cl2,
    coverage_reml_cl2,
    q05,
    hw_reml_cl2,
    det_reml_cl2,
    coverage_ncv_cl2,
    det_ncv_cl2
  )
# mean, alpha, gamma relative errors per dataset (REML vs NCV)
pr_other <- \(e, a) pv(e, a, "rel")

# ---- Section 4.4 / 5: other remedies (details in the supplement) -----------------------------
# combined bias correction + smoothing-parameter term around the NCV fit (ncv-rbc/)
CB_P <- read_pcells("ncv-rbc/summaries/combo-plasmode-cells.csv")
RB_P <- read_pcells("ncv-rbc/summaries/rbc-plasmode-cells.csv")
cb_app <- CB_P |>
  filter(estimand == "beta") |>
  group_by(role, app, arm) |>
  summarise(cov = mean(coverage), q05 = mean(coverage_q05), .groups = "drop")
cbv <- \(a, r = c("counted", "stress test"))
  cb_app$cov[cb_app$arm == a & cb_app$role %in% r]
rb_app <- RB_P |>
  filter(estimand == "beta") |>
  group_by(role, app, arm) |>
  summarise(cov = mean(coverage), .groups = "drop")
rbv <- \(a, r = c("counted", "stress test"))
  rb_app$cov[rb_app$arm == a & rb_app$role %in% r]
CB_S <- read_scells("ncv-rbc/summaries/combo-syn-cells.csv")
cbs <- \(a, blk, col = "coverage")
  mean(CB_S[[col]][CB_S$arm == a & CB_S$block == blk & CB_S$estimand == "beta"])
# NCV estimate with REML + CL2's half-width (summaries/hybrid/), plasmode, REML-fitted truth
hyb_p <- HSUM_P |>
  filter(
    flip == "curve",
    error == "attached",
    G > 40,
    residual %in% c("REML", "NCV"),
    estimand == "beta",
    arm %in% c("reml_cl2", "hybrid")
  )
hyb_v <- \(a, tr, col = "coverage_mean")
  mean(hyb_p[[col]][hyb_p$arm == a & hyb_p$truth == tr])
# bias allowance (NCV + CL2 + delta): plasmode beta under the REML-fitted truth
ba_reml <- bt("ncv_cl2_bias", "REML", role_sel = "counted")
# mgcv jackknife covariance for curve-blocked NCV
JACK <- readr::read_csv(
  "ncv-jack-recheck/jack-reps.csv",
  show_col_types = FALSE
)
jk <- \(a, col = "coverage") {
  JACK |>
    filter(arm == a) |>
    group_by(family, estimand) |>
    summarise(v = mean(.data[[col]]), .groups = "drop") |>
    pull(v)
}
# smoothing-parameter uncertainty term for REML (lambda-correction/)
LC <- readr::read_csv(
  "lambda-correction/summaries/lc-cells.csv",
  show_col_types = FALSE
)
lc <- LC |>
  filter(
    estimand == "beta",
    crit == "z",
    error != "iid",
    arm %in% c("cl2", "ij_lin")
  ) |>
  select(cell, G, arm, coverage, score) |>
  pivot_wider(names_from = arm, values_from = c(coverage, score)) |>
  mutate(gain = coverage_ij_lin - coverage_cl2, sr = score_ij_lin / score_cl2)
# comparators (comparators/): pcre, GLS, curve bootstrap
CMP <- readr::read_csv(
  "comparators/summaries/comparators-cells.csv",
  show_col_types = FALSE
) |>
  mutate(
    setting = case_when(
      source == "plasmode" ~ "plasmode",
      error == "iid" ~ "independent",
      TRUE ~ "dependent"
    )
  )
CMPF <- readr::read_csv(
  "comparators/summaries/comparators-fits.csv",
  show_col_types = FALSE
)
cmp_cov <- \(a, set, e = "beta")
  CMP$coverage[CMP$arm == a & CMP$setting == set & CMP$estimand == e]
cmp_w <- CMP |>
  filter(estimand == "beta") |>
  group_by(setting, cell) |>
  mutate(wr = width / width[arm == "reml_cl2"]) |>
  ungroup()
cmp_wr <- \(a, set) geo(cmp_w$wr[cmp_w$arm == a & cmp_w$setting == set])
cmp_time <- CMPF |>
  group_by(arm) |>
  summarise(med = median(time), .groups = "drop")
# MSE of pcre / GLS relative to the NCV fit (same data; NCV from the study records)
ncv_mse_syn <- S |>
  filter(cell %in% CMP$cell[CMP$source == "synthetic"], arm == "ncv_cl2") |>
  transmute(source = "synthetic", cell, estimand, mse_ncv = mse)
ncv_mse_pm <- P |>
  filter(
    truth == "MID",
    residual == "REML",
    flip == "curve",
    error == "attached",
    arm == "ncv_cl2"
  ) |>
  select(app, G, estimand, mse_ncv = mse) |>
  inner_join(
    CMP |> filter(source == "plasmode") |> distinct(cell, app, G),
    by = c("app", "G")
  ) |>
  transmute(source = "plasmode", cell, estimand, mse_ncv)
cmp_mse <- CMP |>
  filter(arm %in% c("reml_cl2", "pcre", "gls_fpc")) |>
  inner_join(
    bind_rows(ncv_mse_syn, ncv_mse_pm),
    by = c("source", "cell", "estimand")
  ) |>
  mutate(ratio = mse / mse_ncv)
cmr <- \(a, set, e = "beta")
  cmp_mse$ratio[
    cmp_mse$arm == a & cmp_mse$setting == set & cmp_mse$estimand == e
  ]
# AR(1) working model: plasmode (intermediate truth, REML residuals) and synthetic
ar1_syn <- S |> filter(arm == "ar1_model", estimand %in% c("mean", "beta"))
# functional GEE (competitors/fgee/)
FG_B <- readr::read_csv(
  "competitors/fgee/summaries/B-summary.csv",
  show_col_types = FALSE
)
FG_BP <- readr::read_csv(
  "competitors/fgee/summaries/B-paired-ncv-summary.csv",
  show_col_types = FALSE
)
fgb <- FG_B$coverage[
  FG_B$method %in%
    c("fgee_independent", "fgee_ar1", "fgee_fpca") &
    FG_B$estimand == "beta" &
    FG_B$error != "iid"
]
fg_vs_ncv <- FG_BP$est[
  FG_BP$contrast == "fgee_ind_vs_ncv" &
    FG_BP$estimand == "beta" &
    FG_BP$error != "iid"
]
# misspecified truths (misspec/)
MS <- readr::read_csv(
  "misspec/summaries/misspec-cells.csv",
  show_col_types = FALSE
)
msv <- \(tr, b, a, e = "beta", err = c("fpc", "iid"), col = "coverage") {
  MS[[col]][
    MS$truth_label == tr &
      MS$basis == b &
      MS$arm == a &
      MS$estimand == e &
      MS$error %in% err
  ]
}
# timing (Gaussian, G = 100, D = 61)
med_time <- \(m) median(TIMING$fit_time[TIMING$method == m])
cl2_time <- median(TIMING$cl2_time, na.rm = TRUE)
# combined interval vs REML + CL2 per dataset (independent curves): coverage gap and width ratio
cb_cmp <- CB_P |>
  filter(estimand == "beta", arm == "rbc_p_ijx") |>
  group_by(role, app) |>
  summarise(cov = mean(coverage), hw = mean(halfwidth), .groups = "drop") |>
  inner_join(
    RB_P |>
      filter(estimand == "beta", arm == "reml_cl2") |>
      group_by(app) |>
      summarise(
        cov_cl2 = mean(coverage),
        hw_cl2 = mean(halfwidth),
        .groups = "drop"
      ),
    by = "app"
  ) |>
  mutate(gap = cov_cl2 - cov, wr = hw / hw_cl2)
cbc <- cb_cmp |> filter(role == "counted")

# ---- revision additions ----------------------------------------------------------------------
# REML + CL2 coverage of beta below 0.95 in how many plasmode combinations
p_cl2_beta <- P |> filter(arm == "reml_cl2", estimand == "beta")
n_p_below <- sum(p_cl2_beta$coverage_mean < LEVEL)
n_p_all <- nrow(p_cl2_beta)
# absolute shortfall with 40 subjects, datasets with independent curves
g40_abs <- P_G40 |>
  filter(arm == "reml_cl2", estimand %in% c("mean", "beta")) |>
  group_by(role, app, estimand) |>
  summarise(cov = mean(coverage_mean), .groups = "drop")
g40a <- \(e, r = "counted")
  LEVEL - g40_abs$cov[g40_abs$estimand == e & g40_abs$role == r]
# estimation by pcre and GLS with FPCA residual covariance relative to NCV (geometric mean of
# per-setting MSE ratios; plasmode: midpoint truth, REML residuals)
cmg <- \(a, set) geo(cmr(a, set))
pcre_time <- CMPF |>
  filter(arm == "PCRE") |>
  summarise(med = median(time), q90 = quantile(time, 0.9), max = max(time))
boot_time <- median(CMPF$time[CMPF$arm == "BOOT"])
# pcre / GLS intervals in the four synthetic Gaussian settings with dependent errors
cmp_rel <- CMP |>
  filter(setting == "dependent") |>
  group_by(cell, estimand) |>
  mutate(
    wr = width / width[arm == "reml_cl2"],
    sr = score / score[arm == "reml_cl2"]
  ) |>
  ungroup()
cmq <- \(a, e = "beta", col = "coverage")
  cmp_rel[[col]][cmp_rel$arm == a & cmp_rel$estimand == e]
n_cmp_dep <- n_distinct(cmp_rel$cell)
# NCV + CL2 on real residual curves: SE/SD and coverage with smoothing parameters fixed
FX <- read_pcells("ncv-rbc/summaries/fixlam-plasmode-cells.csv") |>
  filter(estimand == "beta", arm == "ncvfix_cl2") |>
  group_by(role, app) |>
  summarise(cov = mean(coverage), se_sd = median(se_sd), .groups = "drop")
ncv_sesd <- RB_P |>
  filter(estimand == "beta", arm == "ncv_cl2") |>
  group_by(role, app) |>
  summarise(se_sd = median(se_sd), .groups = "drop")
# plasmode mean: NCV / REML relative error per dataset
pm_mean_ratio <- prel |>
  filter(estimand == "mean", arm %in% c("reml_cl2", "ncv_cl2")) |>
  group_by(app, arm) |>
  summarise(r = median(rel_rmse_median), .groups = "drop") |>
  pivot_wider(names_from = arm, values_from = r) |>
  mutate(q = ncv_cl2 / reml_cl2)
# misregistration: Gaussian REML + CL2 (mean and beta)
warp_gau <- S |>
  filter(
    block == "warp",
    family == "gaussian",
    estimand %in% c("mean", "beta"),
    arm == "reml_cl2"
  )
# Satterthwaite: binary alpha at G = 40, dependent
satt_bin_alpha <- sv("cov_satt", 40, e = "alpha", fam = "binomial")

# ---- 20 curves, zero surface, Satterthwaite beyond the main factorial --------------------------
# main factorial at G = 20: coverage per setting (all signal levels), z critical value
lad <- SGL |> mutate(dependence = dep_lab(error))
l20 <- \(
  a,
  dep = "dependent",
  e = c("mean", "beta", "gamma"),
  fam = c("gaussian", "poisson", "binomial"),
  col = "coverage_mean",
  g = 20
)
  lad[[col]][
    lad$arm == a &
      lad$dependence == dep &
      lad$estimand %in% e &
      lad$family %in% fam &
      lad$G %in% g
  ]
# REML + CL2 with z, t_{G-1} and Satterthwaite, G = 20 (default signal), per setting
lads <- SGLS |> mutate(dependence = dep_lab(error))
ls20 <- \(
  cr,
  dep = "dependent",
  e = c("mean", "beta", "gamma"),
  fam = c("gaussian", "poisson", "binomial"),
  col = "coverage",
  g = 20
)
  lads[[col]][
    lads$crit == cr &
      lads$dependence == dep &
      lads$estimand %in% e &
      lads$family %in% fam &
      lads$G %in% g
  ]
# Satterthwaite width relative to z at G = 20 (default signal), per setting and estimand
ls20_wr <- lads |>
  filter(G == 20, crit %in% c("z", "satt")) |>
  select(cell, estimand, dependence, crit, width) |>
  pivot_wider(names_from = crit, values_from = width) |>
  mutate(wr = satt / z)
# REML / NCV relative error of beta at G = 20, per setting (dependent errors)
rel20 <- SGREL |>
  inner_join(DESIGN_X, by = "cell") |>
  filter(
    block == "small_G",
    estimand == "beta",
    arm %in% c("reml_cl2", "ncv_cl2")
  ) |>
  select(cell, family, error, signal, arm, rel_rmse_median) |>
  pivot_wider(names_from = arm, values_from = rel_rmse_median) |>
  mutate(ratio = reml_cl2 / ncv_cl2)
rr20 <- \(dep) rel20$ratio[(rel20$error != "iid") == dep]

# beta = 0: pointwise false-positive rate and size of the estimate
nl <- \(
  a,
  err = "fpc",
  fam = c("gaussian", "binomial"),
  g = c(40, 100),
  col = "fpr"
)
  SGN[[col]][
    SGN$arm == a & SGN$error %in% err & SGN$family %in% fam & SGN$G %in% g
  ]
null_tab <- SGN |>
  filter(arm %in% c("reml_model", "reml_cl2", "reml_cl2_satt", "ncv_cl2")) |>
  select(family, error, G, arm, fpr) |>
  pivot_wider(names_from = arm, values_from = fpr) |>
  left_join(
    SGN |>
      filter(arm %in% c("reml_cl2", "ncv_cl2")) |>
      mutate(fit = if_else(arm == "reml_cl2", "rms_reml", "rms_ncv")) |>
      select(family, error, G, fit, rms_rel) |>
      pivot_wider(names_from = fit, values_from = rms_rel),
    by = c("family", "error", "G")
  ) |>
  mutate(
    family = factor(FAMILY_LABELS[family], FAMILY_LABELS),
    error = factor(ERR3[error], ERR3)
  ) |>
  arrange(family, desc(error), G) |>
  mutate(
    dup_e = duplicated(paste(family, error)),
    family = if_else(duplicated(family), "", as.character(family)),
    error = if_else(dup_e, "", as.character(error))
  ) |>
  select(
    family,
    error,
    G,
    reml_model,
    reml_cl2,
    reml_cl2_satt,
    ncv_cl2,
    rms_reml,
    rms_ncv
  )
# NCV / REML size of the null estimate
null_rms <- SGN |>
  filter(arm %in% c("reml_cl2", "ncv_cl2")) |>
  select(family, error, G, arm, rms) |>
  pivot_wider(names_from = arm, values_from = rms) |>
  mutate(ratio = reml_cl2 / ncv_cl2)
nrr <- \(err) null_rms$ratio[null_rms$error == err]

# Satterthwaite in the synthetic extensions (outside the main factorial, G = 20 and beta = 0)
satx_ext <- SATX |> filter(!block %in% c("small_G", "null_beta"))
satxc_ext <- SATXC |>
  filter(!block %in% c("small_G", "null_beta"), contrast == "Satt - z")
sxc <- \(g, col = "diff") mean(satxc_ext[[col]][satxc_ext$G == g])
satx_low <- satx_ext |>
  filter(crit %in% c("z", "satt")) |>
  select(cell, block, family, error, G, estimand, crit, coverage) |>
  pivot_wider(names_from = crit, values_from = coverage) |>
  filter(satt < LEVEL - ADEQUATE)
# every such combination is a count or binary response: misregistration, or alpha / gamma
stopifnot(
  all(satx_low$family %in% c("poisson", "binomial")),
  all(
    satx_low$block == "warp" |
      (satx_low$family == "binomial" &
        satx_low$estimand %in% c("alpha", "gamma"))
  ),
  satt_min_gain > 0
)
n_satx <- nrow(satx_ext |> filter(crit == "satt"))
# per-point Satterthwaite df relative to G - 1: beta, gamma, mean over all synthetic settings
satt_df_rel <- c(
  SATD$df_median_rel[SATD$estimand %in% c("beta", "gamma", "mean")],
  SATXD$df_median_rel[SATXD$estimand %in% c("beta", "f", "gamma", "mean")]
)

# REML + CL2 with Satterthwaite: all estimands except weather beta, per dataset
sp_ok <- satp_app |>
  filter(crit == "satt", !(app == "aemet" & estimand == "beta"))
