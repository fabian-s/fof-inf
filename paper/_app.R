# Objects for the main paper's new results (CL2 on GLS/pcre, interval scores) and Appendix C
# (computing time). Sourced after paper/_main.R; the repository root is the working directory.

# ---- interval score relative to REML + CL2 -------------------------------------------------
score_rel <- S |>
  filter(estimand %in% c("beta", "mean", "alpha", "gamma"), arm %in% MAIN3) |>
  group_by(cell, estimand) |>
  mutate(sr = score_mean / score_mean[arm == "reml_cl2"]) |>
  ungroup() |>
  mutate(dependence = dep_lab(error)) |>
  group_by(estimand, dependence, arm) |>
  summarise(score_ratio = geo(sr), .groups = "drop")
trade_tab <- trade_tab |>
  left_join(
    score_rel |>
      filter(dependence == "dependent") |>
      mutate(
        estimand = factor(
          ESTIMAND_LABELS[estimand],
          levels(trade_tab$estimand)
        ),
        arm = factor(LAB3[arm], LAB3)
      ) |>
      select(estimand, arm, score_ratio),
    by = c("estimand", "arm")
  )
scr <- \(e, a, dep = "dependent")
  score_rel$score_ratio[
    score_rel$estimand == e & score_rel$arm == a & score_rel$dependence == dep
  ]
# plasmode: interval score of NCV + CL2 and model-based relative to REML + CL2, beta, all subjects
p_score <- P_POOL |>
  filter(estimand == "beta", arm %in% MAIN3, residual %in% c("REML", "NCV")) |>
  group_by(app, truth, residual) |>
  mutate(sr = score_mean / score_mean[arm == "reml_cl2"]) |>
  group_by(app, arm) |>
  summarise(sr = geo(sr), .groups = "drop")
ptab <- ptab |>
  left_join(
    p_score |>
      filter(arm == "ncv_cl2") |>
      transmute(dataset = factor(APP_LABELS[app], APP_LABELS), score_ncv = sr),
    by = "dataset"
  )

# ---- CL2 on the GLS and pcre fits (replicates 1-100, where all constructions exist) ---------
CA <- readr::read_csv(
  "comparators/summaries/cl2-arms-cells.csv",
  show_col_types = FALSE
) |>
  filter(window == "1-100") |>
  mutate(
    setting = case_when(
      source == "plasmode" ~ "plasmode",
      error == "iid" ~ "independent",
      TRUE ~ "dependent"
    )
  )
ca <- CA |>
  group_by(setting, estimand, cell) |>
  mutate(
    wr = width / width[arm == "reml_cl2"],
    sr = score / score[arm == "reml_cl2"]
  ) |>
  group_by(setting, estimand, arm) |>
  summarise(
    cov = mean(coverage),
    wr = geo(wr),
    sr = geo(sr),
    se_sd = mean(se_sd_ratio),
    .groups = "drop"
  )
cav <- \(a, set, col = "cov", e = "beta")
  ca[[col]][ca$arm == a & ca$setting == set & ca$estimand == e]
ca_pm_app <- CA |>
  filter(source == "plasmode", estimand == "beta", arm == "gls_fpc_cl2") |>
  group_by(app) |>
  summarise(v = mean(coverage), .groups = "drop")

# ---- Appendix C: computing time -------------------------------------------------------------
tsum <- function(x) {
  x <- x[is.finite(x)]
  tibble::tibble(
    n = length(x),
    median = median(x),
    q25 = quantile(x, 0.25),
    q75 = quantile(x, 0.75),
    q90 = quantile(x, 0.9),
    max = max(x)
  )
}
get_kv <- function(info, key)
  suppressWarnings(as.numeric(sub(
    paste0(".*", key, "=([0-9.eE+-]+).*"),
    "\\1",
    info
  )))
CAF <- readr::read_csv(
  "comparators/summaries/cl2-arms-fits.csv",
  show_col_types = FALSE
) |>
  mutate(cl2_time = get_kv(info, "cl2_time"))
LCT <- readr::read_csv(
  "lambda-correction/summaries/lc-times.csv",
  show_col_types = FALSE
)
LCPT <- readr::read_csv(
  "lambda-correction/summaries-plasmode/lcp-times.csv.gz",
  show_col_types = FALSE
)
NCVT <- readr::read_csv(
  "lambda-correction/summaries-ncv/ncv-times.csv",
  show_col_types = FALSE
)
NCVPT <- readr::read_csv(
  "lambda-correction/summaries-ncv-plasmode/ncvp-times.csv.gz",
  show_col_types = FALSE
)
NCVGT <- readr::read_csv(
  "lambda-correction/summaries-ncvg/ncv-times.csv",
  show_col_types = FALSE
)
satt_line <- grep(
  "df part median",
  readLines("satterthwaite/summaries/SUMMARY.md"),
  value = TRUE
)
satt_df_time <- as.numeric(sub(
  ".*df part median ([0-9.]+) s.*",
  "\\1",
  satt_line
))
satt_task_time <- as.numeric(sub(
  ".*per task: median ([0-9.]+) s.*",
  "\\1",
  satt_line
))
TB <- TIMING
time_tab <- bind_rows(
  tsum(TB$fit_time[TB$method == "REML"]) |>
    mutate(method = "REML fit", data = "A"),
  tsum(TB$fit_time[TB$method == "NCV"]) |>
    mutate(method = "curve-blocked NCV fit", data = "A"),
  tsum(TB$fit_time[TB$method == "AR1"]) |>
    mutate(method = "AR(1) working model (profiled)", data = "A"),
  tsum(TB$cl2_time[TB$method == "REML"]) |>
    mutate(method = "CL2 covariance (REML fit)", data = "A"),
  tsum(CMPF$time[CMPF$arm == "REML"]) |>
    mutate(method = "REML fit", data = "B"),
  tsum(CMPF$time[CMPF$arm == "GLSfpc"]) |>
    mutate(method = "GLS, FPCA residual covariance", data = "B"),
  tsum(CMPF$time[CMPF$arm == "GLSraw"]) |>
    mutate(method = "GLS, raw residual covariance", data = "B"),
  tsum(CMPF$time[CMPF$arm == "PCRE"]) |>
    mutate(method = "pcre fit", data = "B"),
  tsum(CMPF$time[CMPF$arm == "BOOT"]) |>
    mutate(method = "curve bootstrap (199 refits)", data = "B"),
  tsum(CAF$cl2_time[CAF$arm %in% c("GLSfpc", "GLSraw")]) |>
    mutate(method = "CL2 on a GLS fit", data = "B"),
  tsum(CAF$cl2_time[CAF$arm == "PCRE"]) |>
    mutate(method = "CL2 on a pcre fit", data = "B"),
  tsum(LCT$pieces + LCT$ij) |>
    mutate(method = "smoothing-parameter term, REML", data = "C"),
  tsum(LCPT$pieces + LCPT$ij) |>
    mutate(method = "smoothing-parameter term, REML", data = "D"),
  tsum(NCVT$pieces + NCVT$hessian + NCVT$scores) |>
    mutate(method = "smoothing-parameter term, NCV", data = "C"),
  tsum(NCVPT$pieces + NCVPT$hessian + NCVPT$scores) |>
    mutate(method = "smoothing-parameter term, NCV", data = "D"),
  tsum(NCVGT$pieces + NCVGT$hessian_scores) |>
    mutate(
      method = "smoothing-parameter term, NCV (Poisson, binary)",
      data = "C"
    )
) |>
  select(method, data, n, median, q25, q75, q90, max)
tt <- \(m, d, col = "median")
  time_tab[[col]][time_tab$method == m & time_tab$data == d]

# ---- round-1 revision objects ----------------------------------------------------------------
# binary intercept: studentized error of the REML + CL2 interval (core, G = 40/100)
zb <- ZS |>
  filter(
    fit == "REML",
    block == "core",
    family == "binomial",
    estimand == "alpha"
  ) |>
  group_by(dep = error != "iid", G) |>
  summarise(
    m = mean(z_mean_cl2),
    s = mean(z_sd_cl2),
    sm = mean(z_sd_model),
    .groups = "drop"
  )
zbv <- \(dep, g, col) zb[[col]][zb$dep == dep & zb$G == g]
# CL2 cost under independence, binary mean at G = 40
bin40 <- iid40_cost |> filter(family == "binomial", estimand == "mean")
# NCV / REML MSE ratio for the mean (dependent, main factorial, G = 100) and for beta (iid, Gaussian)
mean_ratio_dep <- ratio_of(r_dep, "mean")
beta_ratio_iid_gau <- MSE_RATIO$ratio[
  MSE_RATIO$block == "core" &
    MSE_RATIO$error == "iid" &
    MSE_RATIO$estimand == "beta" &
    MSE_RATIO$family == "gaussian"
]
# plasmode: attached curve-flip combinations of the six datasets with independent source curves
p_ind <- P |>
  filter(
    arm == "reml_cl2",
    estimand == "beta",
    app %in% APP_COUNTED,
    error == "attached",
    flip == "curve"
  )
# REML + CL2 standard error vs replicate spread on the plasmode data (beta, all subjects)
cl2_sesd_p <- RB_P |>
  filter(estimand == "beta", arm == "reml_cl2") |>
  group_by(role, app) |>
  summarise(s = median(se_sd), .groups = "drop")
# scaled-t NCV convergence flags
scat_nc <- 1 -
  weighted.mean(
    FITS$converged[FITS$fit == "NCV" & FITS$fit_family == "scat"],
    FITS$n[FITS$fit == "NCV" & FITS$fit_family == "scat"]
  )
# subject-clustered CL2 in the function-on-scalar design with five curves per subject
FG_A <- readr::read_csv(
  "competitors/fgee/summaries/A-pooled.csv",
  show_col_types = FALSE
)
fga_cl2 <- \(n)
  FG_A$coverage[FG_A$method == "refund_cl2" & FG_A$n_i == 5 & FG_A$N %in% n]
# midpoint truth between the two fits
n_mid_interior <- sum(DATASETS$mid_interior)
mid_above <- DATASETS |>
  filter(edf_ff_MID > pmax(edf_ff_REML, edf_ff_NCV)) |>
  pull(app) |>
  as.character()
# the interval constructions evaluated
constructions <- tibble::tribble(
  ~interval,
  ~fit,
  ~covariance,
  ~critical,
  ~evaluated,
  "model-based",
  "REML",
  "mgcv's $\\mathbf V_c$",
  "z",
  "all settings, plasmode",
  "model-based",
  "NCV",
  "$\\mathbf V_p$",
  "z",
  "all settings (supplement)",
  "REML + CL2",
  "REML",
  "CL2, full block, + $(\\mathbf V_p - \\mathbf V_e)$",
  "z",
  "all settings, plasmode",
  "REML + CL2, Satterthwaite",
  "REML",
  "as REML + CL2",
  "$t_\\nu$ (Satterthwaite)",
  "main factorial (default signal), extensions, plasmode",
  "NCV + CL2",
  "NCV",
  "CL2, full block, + $(\\mathbf V_p - \\mathbf V_e)$",
  "z",
  "all settings, plasmode",
  "CR1, CL2 shortcut, frequentist CL2, $t_{G-1}$",
  "REML",
  "variants of CL2",
  "$z$ or $t_{G-1}$",
  "main factorial, default signal"
)

# ---- round-2 revision objects ----------------------------------------------------------------
pm9 <- \(e = c("mean", "beta"))
  c(pm("reml_cl2", e), pm("reml_cl2", e, "stress test"))
wthr <- P_POOL |>
  filter(
    app == "aemet",
    arm == "reml_cl2",
    estimand == "beta",
    residual %in% c("REML", "NCV")
  )
wthr40 <- P_G40 |>
  filter(
    app == "aemet",
    arm == "reml_cl2",
    estimand == "beta",
    residual %in% c("REML", "NCV")
  )
others_beta <- P_POOL |>
  filter(
    app != "aemet",
    arm == "reml_cl2",
    estimand == "beta",
    residual %in% c("REML", "NCV")
  ) |>
  group_by(app) |>
  summarise(c = mean(coverage_mean), .groups = "drop")
wds <- DATASETS |> filter(app == "aemet")
illus90 <- 2 * pnorm(1.96 * 0.9) - 1
# CL2 cost under independence, binary, G = 100, default signal, per estimand
bin100 <- S |>
  filter(
    block == "core",
    G == 100,
    error == "iid",
    family == "binomial",
    signal == "mid",
    arm %in% c("reml_model", "reml_cl2"),
    estimand %in% c("mean", "alpha", "gamma")
  ) |>
  select(estimand, arm, coverage_mean) |>
  pivot_wider(names_from = arm, values_from = coverage_mean) |>
  mutate(cost = reml_model - reml_cl2)
# binary intercept, independent errors, G = 40: studentized SD per signal level
zb_sig <- ZS |>
  filter(
    fit == "REML",
    block == "core",
    family == "binomial",
    estimand == "alpha",
    error == "iid",
    G == 40
  )
# fixed-lambda diagnostic (covariance components), binary, G = 40
ABL <- readr::read_csv(
  "ablation/summaries/ablation-cells.csv",
  show_col_types = FALSE
)
abl_bin <- \(e, err, a)
  mean(ABL$coverage[
    ABL$family == "binomial" &
      ABL$G == 40 &
      ABL$estimand == e &
      ABL$error %in% err &
      ABL$arm == a &
      ABL$crit == "z"
  ])
n_blocks_weather <- PD_BLOCKS$n_blocks[PD_BLOCKS$app == "aemet"]
blk_weather <- PB |>
  filter(
    app == "aemet",
    name == "coverage:reml_blockcl2",
    cells == "block_flips"
  )
# corrected neighbourhood (jackknife) covariance for curve-blocked NCV (ncv-jack-recheck/)
JF_S <- readr::read_csv(
  "ncv-jack-recheck/summaries/jackfix-syn-cells.csv",
  show_col_types = FALSE
)
JF_P <- readr::read_csv(
  "ncv-jack-recheck/summaries/jackfix-pm-cells.csv",
  show_col_types = FALSE
)
JF_A <- readr::read_csv(
  "ncv-jack-recheck/summaries/jackfix-pm-apps.csv",
  show_col_types = FALSE
)
JF_I <- readr::read_csv(
  "ncv-jack-recheck/summaries/jackfix-syn-info.csv",
  show_col_types = FALSE
)
jfs <- \(a, col = "coverage", e = "beta")
  JF_S[[col]][JF_S$arm == a & JF_S$estimand == e]
jfp <- \(a, col = "coverage", e = "beta")
  JF_P[[col]][JF_P$arm == a & JF_P$estimand == e]
jf_below <- JF_A |>
  filter(estimand == "beta", arm %in% c("jk_nei", "reml_cl2")) |>
  select(app, arm, coverage) |>
  tidyr::pivot_wider(names_from = arm, values_from = coverage)
jf_truth <- JF_P |>
  filter(estimand == "beta", arm == "jk_nei") |>
  group_by(truth) |>
  summarise(c = mean(coverage), .groups = "drop")
jf_app_mean <- JF_P |>
  filter(estimand == "beta", arm == "jk_nei") |>
  group_by(app) |>
  summarise(c = mean(coverage), .groups = "drop")
jf_bin <- \(e, a = "jk_nei")
  JF_S$coverage[JF_S$family == "binomial" & JF_S$estimand == e & JF_S$arm == a]
# scaled-t replicates whose NCV fit is flagged as not converged against the others, per interval
# (analysis/convergence-flag.R)
CONV_FLAG <- read_summary("final", "convergence-flag")
cflag <- \(arm, flagged, err = "fpc", e = "beta", col = "coverage") {
  CONV_FLAG[[col]][
    CONV_FLAG$arm == arm &
      CONV_FLAG$flagged == flagged &
      CONV_FLAG$error == err &
      CONV_FLAG$estimand == e
  ]
}
cflag_gap <- CONV_FLAG |>
  dplyr::filter(arm == "ncv_cl2") |>
  select(error, estimand, flagged, coverage) |>
  tidyr::pivot_wider(
    names_from = flagged,
    values_from = coverage,
    names_prefix = "f"
  ) |>
  mutate(gap = fFALSE - fTRUE)
cflag_arm_gap <- \(arm, err = "fpc", e = "beta") {
  cflag(arm, FALSE, err, e) - cflag(arm, TRUE, err, e)
}
cflag_rate <- \(err) {
  x <- CONV_FLAG[
    CONV_FLAG$arm == "ncv_cl2" &
      CONV_FLAG$estimand == "beta" &
      CONV_FLAG$error == err,
  ]
  sum(x$n_rep[x$flagged]) / sum(x$n_rep)
}
# CL2 against CR1 with the same allowance: coverage gain and width ratio per setting, dependent
# main-factorial settings at the default signal, normal critical value
CR1_CONTR <- readr::read_csv(
  "ablation/summaries/ablation-contrasts.csv",
  show_col_types = FALSE
) |>
  dplyr::filter(
    contrast == "CL2 - CR1",
    error != "iid",
    estimand %in% c("mean", "beta")
  )
CR1_WIDTH <- readr::read_csv(
  "ablation/summaries/ablation-cells.csv",
  show_col_types = FALSE
) |>
  dplyr::filter(
    crit == "z",
    error != "iid",
    arm %in% c("cr1_b", "cl2e_b"),
    estimand %in% c("mean", "beta")
  ) |>
  select(cell, G, estimand, arm, width) |>
  tidyr::pivot_wider(names_from = arm, values_from = width) |>
  mutate(ratio = cl2e_b / cr1_b)
cr1_gain <- \(g) CR1_CONTR$diff[CR1_CONTR$G == g]
cr1_wr <- \(g) CR1_WIDTH$ratio[CR1_WIDTH$G == g]
# REML / NCV ratio of the relative error of beta under dependence, per block of settings
rr_block <- function(blk) {
  REL |>
    dplyr::filter(
      block == blk,
      estimand == "beta",
      error != "iid",
      arm %in% c("reml_cl2", "ncv_cl2")
    ) |>
    select(cell, arm, rel_rmse_median) |>
    tidyr::pivot_wider(names_from = arm, values_from = rel_rmse_median) |>
    mutate(ratio = reml_cl2 / ncv_cl2) |>
    pull(ratio)
}

# restarts of flagged scaled-t NCV fits: first flagged fits (independent errors) and a random
# sample under smooth errors (analysis/scat-convergence.R)
SCAT_SMOOTH <- read_summary("final", "scat-convergence-smooth")
scat_spread <- SCAT_SMOOTH |>
  group_by(task) |>
  summarise(
    crit = diff(range(criterion)),
    # study fit above the best refit (> 0) or best-agreeing: which fit is the better one
    study_above = criterion[start == "production"] - min(criterion),
    worst_refit_above = max(criterion) - criterion[start == "production"],
    surf = max(beta_max_diff),
    se = max(se_cl2_rel_diff),
    .groups = "drop"
  ) |>
  mutate(other_fit = se > 0.5)
stopifnot(sum(scat_spread$other_fit) == 2, nrow(scat_spread) == 20)
scat_same <- dplyr::filter(scat_spread, !other_fit, crit < 0.01)
scat_near <- dplyr::filter(scat_spread, !other_fit, crit >= 0.01)
stopifnot(nrow(scat_same) == 17, nrow(scat_near) == 1)
scat_study_worse <- dplyr::filter(scat_spread, other_fit, study_above > 1)
scat_refit_worse <- dplyr::filter(scat_spread, other_fit, study_above <= 1)
stopifnot(nrow(scat_study_worse) == 1, nrow(scat_refit_worse) == 1)
scat_refit_flagged <- SCAT_SMOOTH |>
  dplyr::filter(start != "production") |>
  summarise(n = n(), f = sum(conv != "full convergence"))
scat_truth_rms <- sqrt(REL$denom[REL$cell == 88 & REL$estimand == "beta"][1])
# no REML fit flagged, in any family (settings with fit records)
stopifnot(all(FITS$converged[FITS$fit == "REML"] == 1))
pm_flag_max <- 1 -
  min(
    read_summary("plasmode", "fits") |>
      dplyr::filter(fit == "NCV") |>
      pull(converged)
  )
sci_tex <- function(x) {
  e <- floor(log10(x))
  sprintf("$%d \\times 10^{%d}$", round(x / 10^e), e)
}

# Low-smoothing NCV fits and the REML fit of the same replicates (analysis/ncv-low-smoothing.R)
LOWSM <- read_summary("final", "ncv-low-smoothing")
lowsm <- function(
  col,
  dep = "dependent",
  fam = "no warp",
  err = "all",
  thr = 3
) {
  x <- LOWSM[[col]][
    LOWSM$dependence == dep &
      LOWSM$family == fam &
      LOWSM$error == err &
      LOWSM$threshold == thr &
      is.na(LOWSM$G)
  ]
  stopifnot(length(x) == 1)
  x
}
LOWSM_WARP <- dplyr::filter(
  LOWSM,
  threshold == 3,
  dependence == "dependent",
  error == "warp"
)
LOWSM_DEP <- dplyr::filter(
  LOWSM,
  threshold == 3,
  dependence == "dependent",
  error %in% c("fpc", "ou", "ou0")
)
lowsm_flag <- LOWSM_DEP |>
  dplyr::filter(family != "scat") |>
  summarise(
    n = sum(round(share_low * n_rep)),
    f = sum(round(share_low * n_rep * flagged_low))
  )
stopifnot(
  lowsm("n_cells") +
    lowsm("n_cells", "independent") +
    sum(LOWSM_WARP$n_cells) ==
    86,
  sum(LOWSM_DEP$n_cells) == lowsm("n_cells"),
  lowsm("share_low_reml", "independent") == 0,
  lowsm("ncv_beats_reml_low", "independent") < 0.05,
  lowsm("smoother_than_reml_low") == 1,
  lowsm("ncv_beats_reml_low") > 0.95,
  # neither fit is uniformly more stable under dependence; REML is under independence, and under misregistration
  lowsm("sd_ncv") < lowsm("sd_reml"),
  lowsm("sd2_ncv") > lowsm("sd2_reml"),
  lowsm("sd_ncv", "independent") > lowsm("sd_reml", "independent"),
  lowsm("sd2_ncv", "independent") > lowsm("sd2_reml", "independent"),
  all(LOWSM_WARP$sd_ncv > LOWSM_WARP$sd_reml),
  lowsm("within_ncv_remltail") < 1.1,
  -lowsm("cov_diff_ncv_remltail") < 0.01,
  lowsm_flag$n == 542,
  lowsm_flag$f == 1
)

# ---- review revision ---------------------------------------------------------------------------
# REML / NCV ratio of the median relative errors of beta, range per block of settings (tbl-ratio)
ratio_tab <- tibble::tribble(
  ~settings,
  ~ratio,
  "main factorial, dependent errors",
  rng(rr(TRUE), 1),
  "main factorial, independent errors",
  rng(rr(FALSE), 1),
  "20 curves, dependent errors",
  rng(rr20(TRUE), 1),
  "20 curves, independent errors",
  rng(rr20(FALSE), 1),
  "four times denser grid",
  rng(rr_block("dense_grid"), 1),
  "larger bases",
  rng(rr_block("basis_size"), 1),
  "rough truth",
  rng(rr_block("rough_truth"), 1),
  "misregistration",
  rng(rr_block("warp"), 1),
  "smooth effect in place of the functional term",
  rng(rr_est("f"), 1)
)
# REML / NCV ratio of mean squared errors of beta under independent errors, main factorial (smooth
# truth, by family) and rough-truth block
mse_iid <- \(blk, fam = unique(MSE_RATIO$family))
  1 /
    MSE_RATIO$ratio[
      MSE_RATIO$block == blk &
        MSE_RATIO$error == "iid" &
        MSE_RATIO$estimand == "beta" &
        MSE_RATIO$family %in% fam
    ]
stopifnot(
  length(mse_iid("core", "binomial")) == 4,
  length(mse_iid("rough_truth")) == 7
)
