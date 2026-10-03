# Shared set-up for the supplement and the paper's appendices (sourced after paper/_common.R).
# supplement labels (set before the frozen write-up code builds its factors)
ARM_LABELS[["ncv_cl2_bias"]] <- "NCV + CL2, widened by NCV–REML difference"
ARM_LABELS[[
  "ncv_cl2f_bias"
]] <- "NCV + CL2 (freq.), widened by NCV–REML difference"
ARM_COLOURS <- setNames(unname(ARM_COLOURS), ARM_LABELS)
SHORT_ERROR[["ou"]] <- "rough"
SHORT_ERROR[["iid"]] <- "indep."
invisible(capture.output(source("paper/supp/_objects.R")))
source("paper/_supp.R")
options(knitr.kable.NA = "", scipen = 100)
theme_paper <- theme_bw(base_size = 9) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey95")
  )
# readable labels for tables built by the frozen write-up code
decode <- function(x) {
  codes <- c(
    ncv_cl2f_bias = "NCV + CL2 (freq.), NCV–REML widening",
    ncv_cl2_bias = "NCV + CL2, NCV–REML widening",
    ncv_cl2 = "NCV + CL2",
    reml_cl2 = "REML + CL2",
    reml_model = "REML, model-based",
    ncv_model = "NCV, model-based",
    ar1_model = "AR(1)",
    reml_curvecl2 = "REML + curve CL2",
    reml_blockcl2 = "REML + block CL2",
    reml_altcl2 = "REML + alt. block CL2",
    ncv_blockcl2 = "NCV + block CL2",
    ncv_curvecl2 = "NCV + curve CL2",
    mse_ratio_ncv_reml = "MSE NCV/REML",
    `IS_ncv_cl2_bias/reml_cl2` = "IS widened/REML+CL2",
    `IS_hybrid/reml_cl2` = "IS NCV est./REML+CL2",
    z_mean_cl2 = "mean z (CL2)",
    z_sd_cl2 = "SD z (CL2)",
    share_big_cl2 = "share |z|>1.96 (CL2)",
    z_mean_model = "mean z (model)",
    z_sd_model = "SD z (model)",
    share_big_model = "share |z|>1.96 (model)",
    log_mse_ratio = "log MSE NCV/REML",
    mse_ratio_block = "MSE block-NCV/NCV",
    mse_ratio_mismatch = "MSE mismatched blocks/NCV",
    coverage_loss_recipe = "coverage loss (method)",
    coverage_loss = "coverage loss",
    curve_flips = "curve flips",
    block_flips = "block flips",
    n_cells = "combinations",
    wiggly = "rough",
    basis_size = "basis size",
    rough_truth = "rough truth",
    dense_grid = "dense grid",
    term_type = "smooth effect",
    lowrank_covariate = "low-rank covariate",
    warp_ar1 = "misreg. (AR(1))",
    ar1_home = "exact AR(1)",
    heteroskedastic = "heteroskedastic",
    oscillating = "sign-changing",
    `coverage:` = "coverage ",
    `cov_` = "cov ",
    `q05_` = "q05 ",
    MID = "midpoint",
    xlarge = "4× finer knots",
    `ff EDF` = "surface EDF",
    `ff basis` = "surface basis",
    `ff REML` = "surface REML",
    hybrid = "NCV est., REML+CL2 width",
    standard = "NCV + CL2, NCV–REML widening",
    `bias-aware` = "NCV–REML widening",
    `bias allowance` = "NCV–REML widening",
    coverage_min = "min coverage",
    rel_error = "rel. error",
    mse_vs_REML = "MSE/REML",
    mse_vs_NCV = "MSE/NCV",
    vs_NCV_min = "MSE/NCV min",
    vs_NCV_max = "MSE/NCV max",
    mse_ratio_beta = "MSE ratio NCV/REML (β)",
    fit_family = "fitted family",
    beat_m = "beat differences (variance-matched)",
    beat_1 = "beat differences (own scale)"
  )
  for (k in names(codes)) x <- gsub(k, codes[[k]], x, fixed = TRUE)
  x
}
relabel <- function(df) {
  df <- as.data.frame(df)
  names(df) <- decode(names(df))
  names(df)[names(df) == "arm"] <- "interval"
  if (
    "cells" %in%
      names(df) &&
      all(
        as.character(df$cells) %in%
          c("curve_flips", "block_flips", "curve flips", "block flips")
      )
  )
    names(df)[names(df) == "cells"] <- "flips"
  nm <- c(iid = "independent", ou = "rough stationary", fpc = "smooth")
  names(df)[names(df) %in% names(nm)] <- nm[names(df)[names(df) %in% names(nm)]]
  for (j in names(df))
    if (is.character(df[[j]]) || is.factor(df[[j]]))
      df[[j]] <- decode(as.character(df[[j]]))
  if ("estimand" %in% names(df))
    df$estimand <- dplyr::coalesce(
      unname(ESTIMAND_LABELS[as.character(df$estimand)]),
      as.character(df$estimand)
    )
  if ("error" %in% names(df))
    df$error <- dplyr::coalesce(
      unname(SHORT_ERROR[as.character(df$error)]),
      as.character(df$error)
    )
  if ("family" %in% names(df))
    df$family <- dplyr::coalesce(
      unname(FAMILY_LABELS[as.character(df$family)]),
      as.character(df$family)
    )
  if ("role" %in% names(df))
    df$role <- dplyr::recode(
      as.character(df$role),
      counted = "indep. curves",
      `stress test` = "dep. source"
    )
  if ("app" %in% names(df))
    df$app <- dplyr::coalesce(
      unname(APP_LABELS[as.character(df$app)]),
      as.character(df$app)
    )
  if ("block" %in% names(df))
    df$block <- dplyr::recode(
      as.character(df$block),
      core = "main factorial",
      warp = "misreg.",
      ar1_home = "exact AR(1)"
    )
  for (j in intersect(c("signal", "basis"), names(df))) {
    v <- as.character(df[[j]])
    v <- gsub("\\bxlarge\\b", "4× finer knots", v)
    v <- gsub("\\blarge\\b", "2× finer knots", v)
    v <- gsub("\\bmid\\b", "default", v)
    df[[j]] <- v
  }
  names(df)[names(df) == "app"] <- "dataset"
  df[, setdiff(names(df), c("cell", "cell_role", "incl_note")), drop = FALSE]
}
stab <- function(
  df,
  caption,
  long = FALSE,
  digits = 3,
  font = 7,
  landscape = FALSE
) {
  df <- relabel(df)
  caption <- gsub("(?<!\\\\)%", "\\\\%", caption, perl = TRUE)
  k <- knitr::kable(
    df,
    format = "latex",
    booktabs = TRUE,
    longtable = long,
    digits = digits,
    caption = caption,
    linesep = ""
  )
  opts <- if (long) c("repeat_header") else c("HOLD_position", "scale_down")
  k <- kableExtra::kable_styling(k, font_size = font, latex_options = opts)
  if (landscape) k <- kableExtra::landscape(k)
  k
}
truth_lab <- c(REML = "REML-fitted", MID = "midpoint", NCV = "NCV-fitted")
relabel_truth <- function(p) {
  if ("truth" %in% names(p$data))
    p$data$truth <- factor(truth_lab[as.character(p$data$truth)], truth_lab)
  p
}

# ---- shared by the paper and the supplement (round 3) ------------------------------------------
# binary intercept, dependent errors: mean studentized CL2 error per error process and G,
# averaged over the signal levels
zbe <- ZS |>
  filter(
    fit == "REML",
    block == "core",
    family == "binomial",
    estimand == "alpha",
    error != "iid"
  ) |>
  group_by(error, G) |>
  summarise(m = mean(z_mean_cl2), .groups = "drop")
zbe_rng <- \(g) rng(zbe$m[zbe$G == g])
# smoothing-parameter term and t_{G-1} reference for REML + CL2 on the plasmode data
# (lambda-correction/summaries-plasmode/lcp-cells.csv; attached errors, curve flips)
LCP <- readr::read_csv(
  "lambda-correction/summaries-plasmode/lcp-cells.csv",
  show_col_types = FALSE
) |>
  select(-G) |>
  left_join(PCELL, by = "cell") |>
  filter(flip == "curve", error.y == "attached", estimand == "beta")
lcp_app <- LCP |>
  group_by(app, G_label, arm, crit) |>
  summarise(c = mean(coverage), .groups = "drop")
lcpv <- \(a, cr, g = "all", weather = FALSE) {
  x <- lcp_app[lcp_app$arm == a & lcp_app$crit == cr & lcp_app$G_label == g, ]
  if (weather) x$c[x$app == "aemet"] else x$c[x$app != "aemet"]
}
# corrected neighbourhood covariance on the plasmode data: SE/SD per dataset (beta)
# (median over the dataset's truth x residual combinations, as in the plasmode table of Supplement S4.2)
jf_sesd <- readr::read_csv(
  "ncv-jack-recheck/summaries/jackfix-pm-cells.csv",
  show_col_types = FALSE
) |>
  filter(estimand == "beta", arm == "jk_nei") |>
  group_by(app) |>
  summarise(se_sd = median(se_sd), .groups = "drop")
jfsd <- \(a) jf_sesd$se_sd[jf_sesd$app == a]
# ffpc() term (regression on the covariate's leading principal component scores) next to ff,
# REML fits, eight synthetic settings with smooth errors, 200 replicates on the study's simulated
# data (ffpc-pilot/ffpc-lrz-summary.csv): relative RMSE of beta pooled over replicates and grid,
# grid-averaged coverage of the true beta by 95% intervals with the exact Bayesian CL2 covariance
# and the normal critical value. Rules: k<K> fixed number of components, p<100 pve> share of
# variance explained. The 0.99 rule only in the 20-replicate pilot (ffpc-pilot-summary.csv).
FFPC <- readr::read_csv(
  "ffpc-pilot/ffpc-lrz-summary.csv",
  show_col_types = FALSE
)
FFPC_PILOT <- readr::read_csv(
  "ffpc-pilot/ffpc-pilot-summary.csv",
  show_col_types = FALSE
)
stopifnot(all(FFPC$n_rep == 200), all(FFPC$method == "REML"))
FFPC_DESIGN <- bind_rows(DESIGN, DESIGN_X) |>
  distinct(cell, .keep_all = TRUE) |>
  filter(cell %in% FFPC$cell) |>
  mutate(
    setting = paste0(
      FAMILY_LABELS[family],
      if_else(signal == "low", ", low signal", ""),
      if_else(truth == "wiggly", ", rough truth", ""),
      if_else(covariate == "lowrank", ", low-rank covariate", ""),
      ", G = ",
      G
    )
  ) |>
  arrange(
    covariate != "lowrank",
    truth != "smooth",
    family != "gaussian",
    desc(G),
    signal
  )
stopifnot(
  nrow(FFPC_DESIGN) == n_distinct(FFPC$cell),
  all(FFPC_DESIGN$error == "fpc"),
  sum(FFPC_DESIGN$covariate == "lowrank") == 1
)
FFPC_LR <- FFPC_DESIGN$cell[FFPC_DESIGN$covariate == "lowrank"]
FFPC_RICH <- setdiff(FFPC_DESIGN$cell, FFPC_LR)
# value for cells x rule x arm; the ff arms are the same under every rule
ffv <- \(cell, rule, arm = "ffpc_cl2", col = "rel_rmse", d = FFPC)
  d[[col]][d$cell %in% cell & d$pve == rule & d$arm == arm]
ffpc_tab <- FFPC |>
  filter(
    arm == "ffpc_cl2" | (pve == "k3" & arm %in% c("ff_reml_cl2", "ff_ncv_cl2"))
  ) |>
  mutate(key = if_else(arm == "ffpc_cl2", pve, arm)) |>
  select(cell, key, K = K_median, rel = rel_rmse, cov = coverage) |>
  pivot_wider(names_from = key, values_from = c(K, rel, cov)) |>
  right_join(select(FFPC_DESIGN, cell, setting), by = "cell") |>
  mutate(cell = factor(cell, FFPC_DESIGN$cell)) |>
  arrange(cell) |>
  select(
    setting,
    rel_ff_reml_cl2,
    cov_ff_reml_cl2,
    rel_ff_ncv_cl2,
    cov_ff_ncv_cl2,
    rel_k3,
    cov_k3,
    rel_k5,
    cov_k5,
    K_p095,
    rel_p095,
    cov_p095,
    K_p090,
    rel_p090,
    cov_p090
  )
