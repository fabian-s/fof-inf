# Helpers for writeup/RESULTS.qmd: load the summary tables, label them, and
# provide formatting, figure and table writers. Everything here reads the committed
# summary CSVs (summaries/final, summaries/plasmode, summaries/hybrid, summaries/dti)
# plus writeup/tables/synthetic-design-effects.csv (analysis/describe-cells.R output).
# Sourced with the repository root as working directory.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

RESTART <- normalizePath(".")
FIG_DIR <- file.path(RESTART, "writeup", "figures")
TAB_DIR <- file.path(RESTART, "writeup", "tables")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TAB_DIR, showWarnings = FALSE, recursive = TRUE)
LEVEL <- 0.95
ADEQUATE <- 0.02
PRACTICAL_RATIO <- 0.95

# The largest tables are stored gzipped (<name>.csv.gz).
read_summary <- function(dir, name) {
  path <- file.path(RESTART, "summaries", dir, paste0(name, ".csv"))
  if (!file.exists(path)) path <- paste0(path, ".gz")
  readr::read_csv(path, show_col_types = FALSE, guess_max = 100000)
}

# Labels ---------------------------------------------------------------------------------

ARM_LABELS <- c(
  reml_model = "REML, model-based",
  reml_cl2 = "REML + CL2",
  ncv_model = "NCV, model-based",
  ncv_cl2 = "NCV + CL2",
  ncv_cl2_bias = "NCV + CL2, bias-aware",
  ncv_cl2f_bias = "NCV + CL2 (freq.), bias-aware",
  ar1_model = "AR(1) working model",
  hybrid = "NCV estimate + REML CL2 SE"
)
ARM_COLOURS <- setNames(
  c(
    "#999999",
    "#0072B2",
    "#BBBBBB",
    "#E69F00",
    "#D55E00",
    "#CC79A7",
    "#009E73",
    "#56B4E9"
  ),
  ARM_LABELS
)
MAIN_ARMS <- c("reml_model", "reml_cl2", "ncv_cl2", "ncv_cl2_bias")
ESTIMAND_LABELS <- c(
  mean = "E(Y | X)",
  beta = "beta(s,t)",
  f = "f(x,t)",
  alpha = "alpha(t)",
  gamma = "gamma(t)"
)
ERROR_LABELS <- c(
  iid = "independent",
  ou = "rough stationary (OU)",
  fpc = "smooth (DTI FPCs)",
  ou0 = "exact AR(1)",
  osc = "sign-changing",
  fpc_raw = "smooth, variance along t",
  fpc_het = "smooth, variance by covariate",
  fpc_raw_het = "smooth, both",
  warp = "misregistration"
)
FAMILY_LABELS <- c(
  gaussian = "Gaussian",
  poisson = "Poisson",
  binomial = "binary",
  scat = "scaled t",
  betar = "beta",
  nb = "negative binomial"
)
APP_LABELS <- c(
  ecg_gls4b = "ECG strain",
  af = "AF trial",
  running = "running",
  dti = "DTI",
  gait_vc = "gait",
  ecg_roahd = "ECG 8-lead",
  ocean = "ocean",
  aemet = "weather",
  adelaide = "electricity"
)
APP_COUNTED <- c("ecg_gls4b", "af", "running", "dti", "gait_vc", "ecg_roahd")
APP_STRESS <- c("ocean", "aemet", "adelaide")
TRUTH_LEVELS <- c("REML", "MID", "NCV")
RESIDUAL_LABELS <- c(
  REML = "REML residuals",
  NCV = "NCV residuals",
  beat_m = "beat differences (variance-matched)",
  beat_1 = "beat differences (own scale)"
)

label_arm <- function(x) factor(ARM_LABELS[x], levels = ARM_LABELS)
label_estimand <- function(x)
  factor(ESTIMAND_LABELS[x], levels = ESTIMAND_LABELS)
label_error <- function(x) factor(ERROR_LABELS[x], levels = ERROR_LABELS)
SHORT_ERROR <- c(
  iid = "iid",
  ou = "OU",
  fpc = "smooth",
  ou0 = "AR(1)",
  osc = "sign-chg.",
  fpc_raw = "var(t)",
  fpc_het = "var(z)",
  fpc_raw_het = "var(t,z)",
  warp = "misreg."
)
label_error_short <- function(x) factor(SHORT_ERROR[x], levels = SHORT_ERROR)
label_family <- function(x) factor(FAMILY_LABELS[x], levels = FAMILY_LABELS)
label_app <- function(x) factor(APP_LABELS[x], levels = APP_LABELS)
label_signal <- function(x) factor(x, levels = c("low", "mid", "high"))

# Formatting -----------------------------------------------------------------------------

f <- function(x, d = 3) formatC(x, digits = d, format = "f")
pp <- function(x, d = 1) sprintf("%s pp", f(100 * x, d))
rng <- function(x, d = 2, sep = "–") {
  x <- x[!is.na(x)]
  if (!length(x)) return("NA")
  lo <- f(min(x), d)
  hi <- f(max(x), d)
  if (lo == hi) return(lo)
  if (min(x) < 0) paste(lo, "to", hi) else paste0(lo, sep, hi)
}
rng_pp <- function(x, d = 1) {
  x <- x[!is.na(x)]
  if (!length(x)) return("NA")
  lo <- f(100 * min(x), d)
  hi <- f(100 * max(x), d)
  if (lo == hi) return(paste(lo, "pp"))
  if (min(x) < 0) paste(lo, "to", hi, "pp") else paste0(lo, "–", hi, " pp")
}
plural <- function(n, word) sprintf("%d %s%s", n, word, if (n == 1) "" else "s")
relabel_apps <- function(x) {
  for (a in names(APP_LABELS))
    x <- gsub(paste0("\\b", a, "\\b"), APP_LABELS[a], x)
  x
}
ci <- function(est, lo, hi, d = 3) {
  sprintf("%s [%s, %s]", f(est, d), f(lo, d), f(hi, d))
}
# Spearman rank correlation with a two-sided permutation p-value (N_PERM random
# permutations, fixed seed); at n = 9 only |rho| above about 0.6 reaches p < 0.05.
N_PERM <- 4000
rho_perm <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  n <- sum(ok)
  if (n < 3) return(list(rho = NA_real_, p = NA_real_, n = n))
  x <- x[ok]
  y <- y[ok]
  r <- cor(x, y, method = "spearman")
  perm <- withr::with_seed(
    1,
    replicate(N_PERM, cor(x, sample(y), method = "spearman"))
  )
  list(rho = r, p = mean(abs(perm) >= abs(r) - 1e-12), n = n)
}
rho_text <- function(x, y) {
  z <- rho_perm(x, y)
  if (is.na(z$rho)) return(sprintf("n = %d, no correlation computed", z$n))
  sprintf("rho = %s, p = %s, n = %d", f(z$rho, 2), f(z$p, 2), z$n)
}
rho_val <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3) return(NA_real_)
  cor(x[ok], y[ok], method = "spearman")
}

# Output ---------------------------------------------------------------------------------

theme_doc <- theme_bw(base_size = 9) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95", colour = NA),
    plot.title = element_text(size = 9, face = "bold")
  )

save_fig <- function(plot, name, width = 6.5, height = 4) {
  for (ext in c("pdf", "png")) {
    ggsave(
      file.path(FIG_DIR, paste0(name, ".", ext)),
      plot,
      width = width,
      height = height,
      dpi = 200,
      device = if (ext == "pdf") grDevices::cairo_pdf else NULL
    )
  }
  invisible(plot)
}

save_tab <- function(df, name) {
  readr::write_csv(df, file.path(TAB_DIR, paste0(name, ".csv")))
  invisible(df)
}

# A kable with numeric columns rounded; used for every table in the document.
show <- function(df, digits = 3, caption = NULL, ...) {
  df <- mutate(
    df,
    across(
      where(is.numeric),
      \(x) if (all(x == round(x), na.rm = TRUE)) x else round(x, digits)
    )
  )
  knitr::kable(df, caption = caption, ...)
}

# Coverage dot plot: point = grid-averaged coverage, bar = +/- 2 MC SE.
coverage_plot <- function(
  d,
  y,
  rows = NULL,
  cols = NULL,
  x = "coverage_mean",
  se = "coverage_se",
  xlab = "coverage (grid average, nominal 0.95)",
  dodge = 0.7
) {
  ggplot(d, aes(.data[[x]], .data[[y]], colour = arm)) +
    geom_vline(xintercept = LEVEL, linetype = 2, colour = "grey40") +
    geom_pointrange(
      aes(
        xmin = .data[[x]] - 2 * .data[[se]],
        xmax = .data[[x]] + 2 * .data[[se]]
      ),
      position = position_dodge(width = dodge),
      size = 0.2
    ) +
    facet_grid(rows = rows, cols = cols, scales = "free_y", space = "free_y") +
    scale_colour_manual(values = ARM_COLOURS, name = NULL, drop = TRUE) +
    guides(colour = guide_legend(nrow = 2)) +
    labs(x = xlab, y = NULL) +
    theme_doc
}

# Synthetic study ------------------------------------------------------------------------

S <- read_summary("final", "cells")
DESIGN <- S |>
  distinct(
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
    extra
  )
W <- read_summary("final", "within-cell-contrasts") |>
  left_join(DESIGN, by = "cell")
PROF <- read_summary("final", "profiles")
Q1 <- read_summary("final", "question1")
REC <- read_summary("final", "recommendation")
RECSIG <- read_summary("final", "recommendation-by-signal")
RECCRIT <- read_summary("final", "recommendation-criteria")
RV <- read_summary("final", "rule-variants")
MC <- read_summary("final", "matched-contrasts") |>
  left_join(DESIGN, by = "cell")
EDF <- read_summary("final", "term-edf") |>
  left_join(DESIGN, by = "cell")
FITS <- read_summary("final", "fit-diagnostics") |>
  left_join(DESIGN, by = "cell")
TIMING <- read_summary("final", "timing-benchmark")
CONC <- read_summary("final", "concurvity")
SCAT <- read_summary("final", "scat-convergence")
DE_FILE <- file.path(TAB_DIR, "synthetic-design-effects.csv")
SDE <- if (file.exists(DE_FILE))
  readr::read_csv(DE_FILE, show_col_types = FALSE) else NULL

# The recommendation pool: core cells at G = 100 with dependence, all
# signal levels, plus the Gaussian smooth-truth misregistration cell at G = 100.
pool_cells <- function(family) {
  DESIGN$cell[
    DESIGN$family == family &
      DESIGN$G == 100 &
      ((DESIGN$block == "core" & DESIGN$error %in% c("ou", "fpc")) |
        (DESIGN$block == "warp" &
          DESIGN$family == "gaussian" &
          DESIGN$truth == "smooth"))
  ]
}

# Wide coverage table for a set of cells and arms.
cov_wide <- function(
  cells_sel,
  estimands = c("mean", "beta"),
  arms = MAIN_ARMS,
  value = "coverage_mean"
) {
  S |>
    filter(cell %in% cells_sel, estimand %in% estimands, arm %in% arms) |>
    select(
      cell,
      family,
      error,
      G,
      signal,
      truth,
      D,
      basis,
      estimand,
      arm,
      all_of(value)
    ) |>
    pivot_wider(names_from = arm, values_from = all_of(value)) |>
    arrange(family, error, G, signal, truth, estimand)
}
# Readable labels for a table built from the design columns / arm codes.
pretty <- function(df) {
  if ("arm" %in% names(df) && !is.factor(df$arm)) df$arm <- ARM_LABELS[df$arm]
  # unknown codes (e.g. the plasmode error column: attached / detached) pass through unchanged
  if ("error" %in% names(df) && !is.factor(df$error))
    df$error <- ifelse(
      df$error %in% names(SHORT_ERROR),
      unname(SHORT_ERROR[df$error]),
      df$error
    )
  if ("family" %in% names(df) && !is.factor(df$family))
    df$family <- unname(FAMILY_LABELS[df$family])
  if ("fit_family" %in% names(df))
    df$fit_family <- unname(FAMILY_LABELS[df$fit_family])
  if ("estimand" %in% names(df) && !is.factor(df$estimand))
    df$estimand <- unname(ESTIMAND_LABELS[df$estimand])
  if ("model" %in% names(df))
    df$model <- c(M1 = "ff(X)", M2 = "f(x,t)")[df$model]
  if ("app" %in% names(df) && !is.factor(df$app))
    df$app <- unname(APP_LABELS[df$app])
  if ("residual" %in% names(df) && !is.factor(df$residual))
    df$residual <- unname(RESIDUAL_LABELS[df$residual])
  arms_in <- intersect(names(df), names(ARM_LABELS))
  if (length(arms_in))
    names(df)[match(arms_in, names(df))] <- ARM_LABELS[arms_in]
  df
}

# The recommendation rule on point estimates (no bootstrap flags) for a set of cells:
# under-only and symmetric calibration criteria.
symmetric_error <- function(coverage, coverage_se) {
  sqrt(max(0, mean((coverage - LEVEL)^2 - coverage_se^2)))
}
calibration_error <- function(coverage, coverage_se) {
  under <- coverage < LEVEL
  sqrt(max(0, mean(under * ((coverage - LEVEL)^2 - coverage_se^2))))
}
rule_point <- function(ce_p, ce_f, ratio) {
  adequate <- c(ce_p, ce_f) <= ADEQUATE
  if (adequate[1] && !adequate[2]) return("P")
  if (!adequate[1] && adequate[2]) return("F")
  if (!any(adequate)) return("N")
  if (ratio <= PRACTICAL_RATIO) return("P")
  if (ratio >= 1 / PRACTICAL_RATIO) return("F")
  "E"
}
rule_both <- function(d, proposal = "ncv_cl2_bias", fallback = "reml_cl2") {
  p <- d |> filter(arm == proposal) |> arrange(cell)
  q <- d |> filter(arm == fallback) |> arrange(cell)
  stopifnot(identical(p$cell, q$cell))
  ratio <- exp(mean(log(p$score_mean / q$score_mean)))
  tibble::tibble(
    n_cells = nrow(p),
    ce_proposal_under = calibration_error(p$coverage_mean, p$coverage_se),
    ce_fallback_under = calibration_error(q$coverage_mean, q$coverage_se),
    ce_proposal_sym = symmetric_error(p$coverage_mean, p$coverage_se),
    ce_fallback_sym = symmetric_error(q$coverage_mean, q$coverage_se),
    score_ratio = ratio,
    decision_under = rule_point(ce_proposal_under, ce_fallback_under, ratio),
    decision_sym = rule_point(ce_proposal_sym, ce_fallback_sym, ratio)
  )
}

# Per-cell NCV/REML MSE ratios with bootstrap intervals.
MSE_RATIO <- W |>
  filter(contrast == "ncv vs reml: squared error") |>
  select(
    cell,
    estimand,
    ratio,
    ratio_lower,
    ratio_upper,
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
    basis
  )

# Plasmode study -------------------------------------------------------------------------

P <- read_summary("plasmode", "coverage") |>
  mutate(
    cell = paste(app, truth, residual, flip, error, G),
    role = if_else(app %in% APP_COUNTED, "counted", "stress test"),
    G_label = if_else(G > 40, "all", "40")
  )
PM <- read_summary("plasmode", "mse-ratio") |>
  mutate(
    role = if_else(app %in% APP_COUNTED, "counted", "stress test"),
    G_label = if_else(G > 40, "all", "40")
  )
PR <- read_summary("plasmode", "rule")
PC <- read_summary("plasmode", "consistency")
PF <- read_summary("plasmode", "flips")
PA <- read_summary("plasmode", "attribution")
PB <- read_summary("plasmode", "block-arm")
PFIT <- read_summary("plasmode", "fits")
PD <- read_summary("plasmode", "diagnostics")

pd_wide <- function(src) {
  PD |>
    filter(source == src) |>
    select(-source) |>
    pivot_wider(names_from = metric, values_from = value)
}
PD_MODEL <- pd_wide("model (REML)")
PD_RESID <- PD |>
  filter(startsWith(source, "residual")) |>
  mutate(residual = sub("residual ", "", source)) |>
  select(-source) |>
  pivot_wider(names_from = metric, values_from = value)
PD_TRUTH <- PD |>
  filter(startsWith(source, "truth")) |>
  mutate(truth = sub("truth ", "", source)) |>
  select(-source) |>
  pivot_wider(names_from = metric, values_from = value)
PD_BLOCKS <- pd_wide("blocks (REML residuals)")

# Phase variability of each dataset (share of response variance removed by SRVF
# registration, from the dataset screening, table column
# "SRVF removed"). These are screening facts, not results.
PHASE_SHARE <- tibble::tribble(
  ~app,
  ~phase_share,
  "ecg_gls4b",
  0.63,
  "running",
  0.30,
  "dti",
  0.22,
  "gait_vc",
  0.03,
  "ecg_roahd",
  0.33,
  "af",
  0.43,
  "ocean",
  0,
  "aemet",
  0.06,
  "adelaide",
  0
)

# One row per dataset: design properties and the frozen fits' EDFs.
DATASETS <- PD_MODEL |>
  select(app, G, D, R2, rank_X) |>
  left_join(
    PD_RESID |>
      filter(residual == "REML") |>
      select(
        app,
        DE,
        neg_share,
        min_r,
        sd_ratio,
        lag1,
        lag25,
        fpc90,
        eigen_share1
      ),
    by = "app"
  ) |>
  left_join(
    PD_TRUTH |>
      select(
        app,
        truth,
        edf_ff = `edf:te(X.smat,X.tmat):L.X`,
        edf_total,
        mid_interior
      ) |>
      pivot_wider(
        names_from = truth,
        values_from = c(edf_ff, edf_total, mid_interior)
      ),
    by = "app"
  ) |>
  mutate(
    DE_per_D = DE / D,
    edf_gap = edf_ff_REML - edf_ff_NCV,
    edf_ratio = edf_ff_REML / edf_ff_NCV,
    mid_interior = mid_interior_MID == 1,
    role = if_else(app %in% APP_COUNTED, "counted", "stress test")
  ) |>
  select(-mid_interior_REML, -mid_interior_NCV, -mid_interior_MID) |>
  left_join(PHASE_SHARE, by = "app") |>
  mutate(app = factor(app, levels = names(APP_LABELS))) |>
  arrange(app)

# The per-dataset pool: attached, curve flips, all subjects.
P_POOL <- P |>
  filter(error == "attached", flip == "curve", G_label == "all")
P_G40 <- P |>
  filter(error == "attached", flip == "curve", G_label == "40")

# Pool means per dataset, estimand and arm (equal weight per cell), with the minimum
# 5% pointwise quantile over the pool's cells.
POOL_MEANS <- P_POOL |>
  group_by(app, role, estimand, arm) |>
  summarise(
    n_cells = n(),
    coverage = mean(coverage_mean),
    coverage_min = min(coverage_mean),
    coverage_max = max(coverage_mean),
    q05_min = min(pointwise_q05),
    q05_mean = mean(pointwise_q05),
    .groups = "drop"
  )

# NCV / REML MSE ratio per dataset, truth source and estimand (geometric mean over the
# pool's residual sources).
RATIO_BY_TRUTH <- PM |>
  filter(error == "attached", flip == "curve", G_label == "all") |>
  group_by(app, role, truth, estimand) |>
  summarise(
    n_cells = n(),
    ratio = exp(mean(log(mse_ratio_ncv_reml))),
    ratio_min = min(mse_ratio_ncv_reml),
    ratio_max = max(mse_ratio_ncv_reml),
    .groups = "drop"
  )

# Hybrid interval analysis ----------------------------------------------------------------

HS <- read_summary("hybrid", "synthetic-hybrid-digest")
HP <- read_summary("hybrid", "plasmode-hybrid-digest")
HMC_S <- read_summary("hybrid", "synthetic-hybrid-map-correlations")
HMC_P <- read_summary("hybrid", "plasmode-hybrid-map-correlations")
HSUM_S <- read_summary("hybrid", "synthetic-hybrid-summary")
HSUM_P <- read_summary("hybrid", "plasmode-hybrid-summary")
HPAIR_S <- read_summary("hybrid", "synthetic-hybrid-paired")
HPAIR_P <- read_summary("hybrid", "plasmode-hybrid-paired")
read_maps <- function(which)
  read_summary("hybrid", paste0(which, "-hybrid-maps"))

# Studentised-error diagnostics (writeup/lrz-zstats.R on LRZ; per cell x estimand x fit:
# mean and SD of (estimate - truth) / SE over grid points and replicates, and the share of
# |z| > 1.96, with the CL2 SE and the model-based SE).
zs_file <- \(name) file.path(TAB_DIR, paste0("zstats-", name, ".csv"))
ZS <- if (file.exists(zs_file("synthetic")))
  readr::read_csv(zs_file("synthetic"), show_col_types = FALSE) else NULL
ZP <- if (file.exists(zs_file("plasmode")))
  readr::read_csv(zs_file("plasmode"), show_col_types = FALSE) else NULL

# DTI application -------------------------------------------------------------------------

DTI_RECIPES <- read_summary("dti", "recipes")
