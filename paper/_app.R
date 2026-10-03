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
