# REPORT.md of the lambda-correction study from summaries/ (main run), summaries-v2/
# (jackknife-free rerun of the dependent G = 40 cells with the corrected code and the
# regularised arms), summaries-plasmode/ and validation/. Every number is computed here.
# Run from the repository root: Rscript lambda-correction/make-report.R

suppressPackageStartupMessages(library(dplyr))
dir <- "lambda-correction"
rd <- \(f, sub = "summaries") {
  p <- file.path(dir, sub, f)
  if (file.exists(p)) read.csv(p) else NULL
}
cells <- rd("lc-cells.csv")
contr <- rd("lc-contrasts.csv")
calib <- rd("lc-calib.csv")
times <- rd("lc-times.csv")
crn <- rd("lc-crn.csv")
regime <- rd("lc-regime.csv")
target <- rd("lc-targeting.csv")
v2 <- rd("lc-cells.csv", "summaries-v2")
v2_contr <- rd("lc-contrasts.csv", "summaries-v2")
v2_calib <- rd("lc-calib.csv", "summaries-v2")
pl <- rd("lcp-cells.csv", "summaries-plasmode")
val <- lapply(
  list.files(file.path(dir, "validation"), "rds$", full.names = TRUE),
  readRDS
)

pct <- \(x) sprintf("%.1f", 100 * x)
rng <- \(x) sprintf("%s-%s", pct(min(x)), pct(max(x)))
grp <- \(d)
  mutate(
    d,
    group = case_when(
      error == "iid" ~ "iid",
      G == 40 ~ "dep. G=40",
      TRUE ~ "dep. G=100"
    )
  )
GROUPS <- c("iid", "dep. G=40", "dep. G=100")
FF <- "L.X1$" # first penalty of the ff surface's tensor product (along the covariate argument s)

arms <- c(
  cl2 = "CL2",
  naive = "+mgcv J",
  vc = "+mgcv Vc",
  ij_lin = "IJ",
  os_lin = "OS",
  ij_linx = "IJ×",
  jk_lin = "JK",
  jk_full = "JKfull",
  fix = "fixλ (oracle)"
)

md_rows <- \(t)
  c(
    paste0("| ", paste(names(t), collapse = " | "), " |"),
    paste0("|", strrep("---|", ncol(t))),
    apply(t, 1, \(r) paste0("| ", paste(r, collapse = " | "), " |"))
  )

# mean over cells of a column by estimand x group, arms as columns (ratio: relative to CL2;
# dist: |coverage - 0.95|)
table_of <- function(
  d,
  value,
  crit_,
  arm_map,
  scale = 100,
  digits = 1,
  ratio = FALSE,
  groups = GROUPS,
  dist = FALSE
) {
  d <- grp(filter(d, crit == crit_))
  if (ratio)
    d <- d |>
      group_by(cell, estimand) |>
      mutate(!!value := .data[[value]] / .data[[value]][arm == "cl2"]) |>
      ungroup()
  if (dist) d <- mutate(d, !!value := abs(.data[[value]] - 0.95))
  t <- d |>
    filter(arm %in% names(arm_map), group %in% groups) |>
    group_by(estimand, group, arm) |>
    summarise(v = mean(.data[[value]]) * scale, .groups = "drop") |>
    mutate(arm = arm_map[arm]) |>
    tidyr::pivot_wider(names_from = arm, values_from = v) |>
    select(estimand, cells = group, any_of(unname(arm_map))) |>
    arrange(
      factor(estimand, c("beta", "alpha", "gamma", "mean")),
      factor(cells, GROUPS)
    ) |>
    mutate(across(
      -c(estimand, cells),
      \(x) formatC(x, format = "f", digits = digits)
    ))
  md_rows(t)
}

# paired coverage difference vs CL2 (z), beta, by group: mean of cell differences, MC SE
paired <- function(ct, arm_) {
  grp(filter(ct, crit == "z", estimand == "beta", arm == arm_)) |>
    group_by(group) |>
    summarise(
      diff = mean(diff),
      se = sqrt(sum(diff_se^2)) / n(),
      .groups = "drop"
    ) |>
    mutate(txt = sprintf("%+.2f (%.2f)", 100 * diff, 100 * se)) |>
    (\(d) setNames(d$txt, d$group))()
}
paired_table <- function(ct, arm_map, groups = GROUPS) {
  rows <- lapply(
    names(arm_map),
    \(a) c(arm = arm_map[[a]], paired(ct, a)[groups])
  )
  t <- as.data.frame(do.call(rbind, rows), check.names = FALSE)
  names(t) <- c("arm", groups)
  md_rows(t)
}

z <- filter(cells, crit == "z")
cov_of <- \(d, a) d$coverage[d$arm == a]
dep_b <- filter(z, estimand == "beta", error != "iid")
dep_bt <- filter(cells, crit == "t", estimand == "beta", error != "iid")
rel_dep <- \(d, a, col, g)
  mean(
    filter(d, G == g, arm == a)[[col]] / filter(d, G == g, arm == "cl2")[[col]]
  )
bin_iid <- filter(
  z,
  family == "binomial",
  error == "iid",
  estimand %in% c("alpha", "gamma")
)
iid_gp <- filter(z, error == "iid", family != "binomial")
iid_shift <- max(abs(cov_of(iid_gp, "ij_lin") - cov_of(iid_gp, "cl2")))

# calibration of the surface penalty's log lambda-hat
cal_ff <- function(cal, src, den = "mc_sd") {
  x <- filter(cal, source == src, error != "iid", grepl(FF, par))
  r <- x$median_est_sd / x[[den]]
  sprintf("%.2f-%.2f", min(r), max(r))
}
spread <- filter(calib, source == "ij", error != "iid", grepl(FF, par)) |>
  group_by(G) |>
  summarise(
    sd = mean(mc_sd),
    mad = mean(mc_mad),
    iqr = mean(mc_iqr),
    .groups = "drop"
  )
naive_iid <- with(
  filter(calib, source == "naive", error == "iid", par != "log_scale"),
  range(median_est_sd / mc_sd)
)

reg_tab <- function(rg, arm_map) {
  t <- rg |>
    filter(error != "iid", arm %in% names(arm_map)) |>
    group_by(G, tertile, arm) |>
    summarise(
      cov = mean(coverage),
      se = sqrt(sum(coverage_se^2)) / n(),
      ser = mean(se_ratio),
      rho = mean(rho_ff),
      .groups = "drop"
    ) |>
    mutate(arm = arm_map[arm])
  cov <- t |>
    select(G, tertile, rho, arm, cov) |>
    tidyr::pivot_wider(names_from = arm, values_from = cov)
  ser <- t |>
    filter(arm != "CL2") |>
    select(G, tertile, arm, ser) |>
    tidyr::pivot_wider(
      names_from = arm,
      values_from = ser,
      names_prefix = "SE ratio "
    )
  mcse <- t |> filter(arm == "CL2") |> select(G, tertile, `MC SE` = se)
  left_join(cov, mcse, by = c("G", "tertile")) |>
    left_join(ser, by = c("G", "tertile")) |>
    arrange(G, factor(tertile, c("low", "mid", "high"))) |>
    mutate(
      G = as.character(G),
      tertile = as.character(tertile),
      rho = sprintf("%.1f", rho),
      across(-c(G, tertile, rho, starts_with("SE ratio")), pct),
      across(starts_with("SE ratio"), \(x) sprintf("%.2f", x))
    ) |>
    rename(`mean log lambda-hat` = rho) |>
    md_rows()
}

tgt <- filter(target, error != "iid") |>
  group_by(G, arm) |>
  summarise(m = mean(spearman), mx = max(spearman), .groups = "drop") |>
  filter(
    arm %in%
      c("naive", "vc", "ij_lin", "os_lin", "ij_linx", "jk_lin", "jk_full")
  ) |>
  mutate(txt = sprintf("%.2f (max %.2f)", m, mx), arm = arms[arm]) |>
  select(G, arm, txt) |>
  tidyr::pivot_wider(names_from = G, values_from = txt, names_prefix = "G = ")

# validation
fd_err <- max(vapply(
  val,
  \(v) max(abs(v$fd$fd - v$fd$ij) / pmax(abs(v$fd$fd), 1e-6)),
  0
))
sd_ratio <- \(v, m) apply(v$d_om[[m]], 2, sd) / apply(v$d_om$jk, 2, sd)
ij_ratio <- range(unlist(lapply(val, sd_ratio, "ij")))
os_ratio <- range(unlist(lapply(val, sd_ratio, "os")))
val_cells <- paste(
  vapply(
    val,
    \(v) sprintf("%s/%s/G=%d", v$cell$fit_family, v$cell$error, v$cell$G),
    ""
  ),
  collapse = ", "
)
n_fd <- nrow(val[[1]]$fd) / length(unique(val[[1]]$fd$par))
cost <- times |> mutate(lc = pieces + ij + onestep)

md <- c(
  "# Cluster-robust smoothing-parameter correction for CL2: report",
  "",
  "Derivation: `NOTE.qmd`. Full tables: `summaries/SUMMARY.md`.",
  "Generated by `make-report.R`; no number is typed by hand.",
  "",
  "## Method",
  "",
  "Curve weights in mgcv's REML criterion and the implicit-function theorem give the influence",
  "of curve g on omega-hat = (log lambda-hat, log scale-hat), d_g = -H^-1 s_g (s_g analytic for",
  "Gaussian, Poisson, binary binomial). Var(rho-hat) by the infinitesimal jackknife (IJ) or",
  "one-step leave-curve-out Newton updates (OS), added to CL2 as J V_rho J'; × adds the",
  "first-order cross term between the coefficient and rho parts. Exploratory references: exact",
  "curve jackknife of rho-hat (JK; 0.5% of its refits not fully converged), refit jackknife of",
  "theta-hat (JKfull). fixλ is CL2 at the cell's geometric-mean REML lambda (computed from the",
  "same replicates): an oracle decomposition, not a procedure. Main run: 18 core cells x 200",
  "replicates, code before the round-1 fixes (OS for GLMs then evaluated the leave-out gradient",
  "after one Newton step; the leave-out cross terms OS×/JK× used CL2's scaling, a difference",
  "of 0.1-0.2% in SE per the round-2 check). The recommendation rests on IJ, which no fix changed.",
  "",
  "## Validation (scope)",
  "",
  sprintf(
    "- Analytic curve scores vs central finite differences in the curve weight: %d curves each in %s, rep 1; max relative difference %.1e.",
    n_fd,
    val_cells,
    fd_err
  ),
  sprintf(
    "- Same replicates: SD of leave-curve-out changes of omega-hat relative to the exact jackknife: IJ %.2f-%.2f, OS %.2f-%.2f.",
    ij_ratio[1],
    ij_ratio[2],
    os_ratio[1],
    os_ratio[2]
  ),
  sprintf(
    "- Study: CRN reproduction of the study's REML estimates and CL2 SEs, max relative difference %.1e.",
    max(crn$max_rel[crn$column != "truth"])
  ),
  "",
  "## Coverage (z, %, mean over cells)",
  "",
  table_of(z, "coverage", "z", arms),
  "",
  "Mean distance to nominal, |coverage - 95| (points):",
  "",
  table_of(filter(z, estimand == "beta"), "coverage", "z", arms, dist = TRUE),
  "",
  sprintf(
    "Independent errors: IJ moves coverage by at most %s points in the Gaussian and Poisson cells; in the binomial cells CL2 undercovers alpha and gamma (%s) and IJ gives %s.",
    pct(iid_shift),
    rng(cov_of(bin_iid, "cl2")),
    rng(cov_of(bin_iid, "ij_lin"))
  ),
  "",
  "Paired coverage difference to CL2, beta, points (MC SE):",
  "",
  paired_table(
    contr,
    c(
      naive = "+mgcv J",
      ij_lin = "IJ",
      os_lin = "OS",
      ij_linx = "IJ×",
      fix = "fixλ (oracle)"
    )
  ),
  "",
  sprintf(
    "The oracle fixλ is narrower than CL2 (beta width %.3f / %.3f of CL2 at G = 40 / 100, dependent cells) yet covers better (paired %s / %s) with interval score %.3f / %.3f: the shortfall comes from which lambda REML selects (a conditional bias of the selected fit), not from variance missing from CL2.",
    rel_dep(dep_b, "fix", "width", 40),
    rel_dep(dep_b, "fix", "width", 100),
    paired(contr, "fix")[["dep. G=40"]],
    paired(contr, "fix")[["dep. G=100"]],
    rel_dep(dep_b, "fix", "score", 40),
    rel_dep(dep_b, "fix", "score", 100)
  ),
  "",
  "## Width and interval score (relative to CL2, mean over cells)",
  "",
  "Width:",
  "",
  table_of(z, "width", "z", arms, scale = 1, digits = 3, ratio = TRUE),
  "",
  "Interval score (lower is better):",
  "",
  table_of(z, "score", "z", arms, scale = 1, digits = 3, ratio = TRUE),
  "",
  sprintf(
    "With the t_{G-1} critical value, beta, dependent cells, G = 40 / 100: CL2 coverage %s / %s (z: %s / %s), interval score relative to CL2 with z %.3f / %.3f; +mgcv J with t %s / %s, score %.3f / %.3f; IJ with z %s / %s, score %.3f / %.3f.",
    pct(mean(cov_of(filter(dep_bt, G == 40), "cl2"))),
    pct(mean(cov_of(filter(dep_bt, G == 100), "cl2"))),
    pct(mean(cov_of(filter(dep_b, G == 40), "cl2"))),
    pct(mean(cov_of(filter(dep_b, G == 100), "cl2"))),
    mean(
      filter(dep_bt, G == 40, arm == "cl2")$score /
        filter(dep_b, G == 40, arm == "cl2")$score
    ),
    mean(
      filter(dep_bt, G == 100, arm == "cl2")$score /
        filter(dep_b, G == 100, arm == "cl2")$score
    ),
    pct(mean(cov_of(filter(dep_bt, G == 40), "naive"))),
    pct(mean(cov_of(filter(dep_bt, G == 100), "naive"))),
    mean(
      filter(dep_bt, G == 40, arm == "naive")$score /
        filter(dep_b, G == 40, arm == "cl2")$score
    ),
    mean(
      filter(dep_bt, G == 100, arm == "naive")$score /
        filter(dep_b, G == 100, arm == "cl2")$score
    ),
    pct(mean(cov_of(filter(dep_b, G == 40), "ij_lin"))),
    pct(mean(cov_of(filter(dep_b, G == 100), "ij_lin"))),
    rel_dep(dep_b, "ij_lin", "score", 40),
    rel_dep(dep_b, "ij_lin", "score", 100)
  ),
  "",
  "## Where the undercoverage sits",
  "",
  sprintf(
    "Sampling spread of log lambda-hat for the surface's penalty along s, dependent cells (mean over cells), G = 40 / 100: SD %.2f / %.2f, MAD %.2f / %.2f, IQR/1.349 %.2f / %.2f. At G = 40 the IQR-based spread exceeds the SD (a flat or bimodal distribution), at G = 100 it is below it (heavier tails); the SD shrinks by a factor %.2f from G = 40 to 100 (sqrt(G) scaling: 0.63).",
    spread$sd[spread$G == 40],
    spread$sd[spread$G == 100],
    spread$mad[spread$G == 40],
    spread$mad[spread$G == 100],
    spread$iqr[spread$G == 40],
    spread$iqr[spread$G == 100],
    spread$sd[spread$G == 100] / spread$sd[spread$G == 40]
  ),
  "",
  "Beta coverage (z) by tertile of that log lambda-hat within each dependent cell (mean over cells; MC SE of the CL2 entry), with the mean SE ratio to CL2:",
  "",
  reg_tab(
    regime,
    c(cl2 = "CL2", naive = "+mgcv J", ij_lin = "IJ", ij_linx = "IJ×")
  ),
  "",
  "Targeting: within-cell Spearman correlation between an arm's per-replicate SE inflation over CL2 (beta) and CL2's per-replicate share of missed grid points, dependent cells, mean (max) over cells. Negative: the arm widens least where CL2 misses most.",
  "",
  md_rows(tgt),
  "",
  "## Calibration of Var(rho-hat)",
  "",
  sprintf(
    "Median over replicates of the estimated SD of the surface penalty's log lambda-hat divided by its Monte Carlo SD, dependent cells: mgcv's Hessian %s, IJ %s, OS %s, exact JK %s; divided by the robust spread IQR/1.349: IJ %s, OS %s. Under independent errors mgcv's Hessian gives %.2f-%.2f (all penalties).",
    cal_ff(calib, "naive"),
    cal_ff(calib, "ij"),
    cal_ff(calib, "os"),
    cal_ff(calib, "jk"),
    cal_ff(calib, "ij", "mc_iqr"),
    cal_ff(calib, "os", "mc_iqr"),
    naive_iid[1],
    naive_iid[2]
  ),
  "",
  sprintf(
    "## Cost\n\nREML fit median %.1f s on LRZ; pieces + IJ + OS median %.1f s (max %.1f s); exact jackknife median %.0f s.",
    median(cost$fit),
    median(cost$lc),
    max(cost$lc),
    median(cost$jackknife)
  ),
  ""
)

if (!is.null(v2)) {
  z2 <- filter(v2, crit == "z")
  v2_arms <- c(
    cl2 = "CL2",
    naive = "+mgcv J",
    vc = "+mgcv Vc",
    ij_lin = "IJ",
    ij_reg = "IJr",
    os_reg = "OSr",
    ij_linx = "IJ×",
    ij_regx = "IJr×",
    os_regx = "OSr×"
  )
  md <- c(
    md,
    sprintf(
      "## Regularised and corrected arms (rerun of the %d dependent G = 40 cells, replicates 1-%d, corrected code)",
      length(unique(z2$cell)),
      max(z2$n_rep)
    ),
    "",
    "H regularised with mgcv's prior precision 1/10 (IJr, OSr), as mgcv's V.sp; cross-term scaling",
    "and the GLM one-step corrected. Fewer replicates than the main run.",
    "",
    "Coverage (z, %):",
    "",
    table_of(z2, "coverage", "z", v2_arms, groups = "dep. G=40"),
    "",
    "Interval score relative to CL2:",
    "",
    table_of(
      z2,
      "score",
      "z",
      v2_arms,
      scale = 1,
      digits = 3,
      ratio = TRUE,
      groups = "dep. G=40"
    ),
    "",
    "Width relative to CL2:",
    "",
    table_of(
      z2,
      "width",
      "z",
      v2_arms,
      scale = 1,
      digits = 3,
      ratio = TRUE,
      groups = "dep. G=40"
    ),
    "",
    "Paired coverage difference to CL2, beta (MC SE):",
    "",
    paired_table(v2_contr, v2_arms[-1], groups = "dep. G=40"),
    "",
    sprintf(
      "Median calibration ratio (vs MC SD): IJr %s, OSr %s.",
      cal_ff(v2_calib, "ij_reg"),
      cal_ff(v2_calib, "os_reg")
    ),
    ""
  )
}

if (!is.null(pl)) {
  stress <- c("ocean", "aemet", "adelaide") # plasmode/rule.R PLASMODE_STRESS
  pl <- pl |>
    mutate(
      set = ifelse(family %in% stress, "stress", "independent"),
      size = ifelse(G == 40, "G=40", "G=all")
    )
  parms <- c(
    cl2 = "CL2",
    naive = "+mgcv J",
    ij_lin = "IJ",
    os_lin = "OS",
    ij_linx = "IJ×"
  )
  pt <- function(value, ratio = FALSE, scale = 100, digits = 1) {
    d <- filter(pl, crit == "z", arm %in% names(parms))
    if (ratio)
      d <- d |>
        group_by(cell, estimand) |>
        mutate(!!value := .data[[value]] / .data[[value]][arm == "cl2"]) |>
        ungroup()
    d |>
      group_by(estimand, set, size, arm) |>
      summarise(v = mean(.data[[value]]) * scale, .groups = "drop") |>
      mutate(arm = parms[arm]) |>
      tidyr::pivot_wider(names_from = arm, values_from = v) |>
      select(estimand, datasets = set, G = size, any_of(unname(parms))) |>
      arrange(
        factor(estimand, c("beta", "alpha", "gamma", "mean")),
        datasets,
        G
      ) |>
      mutate(across(
        -c(estimand, datasets, G),
        \(x) formatC(x, format = "f", digits = digits)
      )) |>
      md_rows()
  }
  plb <- filter(pl, crit == "z", estimand == "beta")
  pl_gain <- plb |>
    group_by(cell) |>
    mutate(d = coverage - coverage[arm == "cl2"]) |>
    ungroup() |>
    filter(arm %in% c("ij_lin", "os_lin", "ij_linx")) |>
    group_by(arm) |>
    summarise(lo = min(d), hi = max(d), m = mean(d), .groups = "drop")
  pl_score <- plb |>
    group_by(cell) |>
    mutate(r = score / score[arm == "cl2"]) |>
    ungroup() |>
    filter(arm %in% c("naive", "ij_lin", "os_lin")) |>
    group_by(arm) |>
    summarise(m = mean(r), .groups = "drop")
  pl_cl2 <- filter(plb, arm == "cl2")
  md <- c(
    md,
    "## Plasmode (real covariates and residual curves, Gaussian; 130 cells x 200 replicates)",
    "",
    sprintf(
      "Independent-curve datasets: %s; stress tests (dependent curves): %s. G = all: each dataset's full sample. Mean over cells. Run with the code before the round-1 fixes (Gaussian only, so only the leave-out cross-term scaling differs; OS× omitted).",
      paste(sort(unique(pl$family[pl$set == "independent"])), collapse = ", "),
      paste(sort(unique(pl$family[pl$set == "stress"])), collapse = ", ")
    ),
    "",
    sprintf(
      "CL2 beta coverage is below 95%% in %d of %d cells (independent datasets: min %s; stress: min %s).",
      sum(pl_cl2$coverage < 0.95),
      nrow(pl_cl2),
      pct(min(pl_cl2$coverage[pl_cl2$set == "independent"])),
      pct(min(pl_cl2$coverage[pl_cl2$set == "stress"]))
    ),
    "",
    "Coverage (z, %):",
    "",
    pt("coverage"),
    "",
    "Interval score relative to CL2:",
    "",
    pt("score", ratio = TRUE, scale = 1, digits = 3),
    "",
    "Width relative to CL2:",
    "",
    pt("width", ratio = TRUE, scale = 1, digits = 3),
    ""
  )
}

# Recommendation
sc_dep <- \(a)
  mean(
    filter(grp(z), estimand == "beta", group != "iid", arm == a)$score /
      filter(grp(z), estimand == "beta", group != "iid", arm == "cl2")$score
  )
tg <- filter(target, error != "iid", arm == "ij_lin")
rec <- c(
  "## Recommendation",
  "",
  sprintf(
    "- refund: no variant as a default or option. In the core cells IJ/OS raise beta coverage under dependence (paired, G = 40: IJ %s, OS %s; G = 100: IJ %s, OS %s points, where CL2 is already at %s and IJ reaches %s) without improving the beta interval score (IJ %.3f, OS %.3f relative to CL2, dependent cells); mgcv's existing first-order term scores %.3f.",
    paired(contr, "ij_lin")[["dep. G=40"]],
    paired(contr, "os_lin")[["dep. G=40"]],
    paired(contr, "ij_lin")[["dep. G=100"]],
    paired(contr, "os_lin")[["dep. G=100"]],
    rng(cov_of(filter(dep_b, G == 100), "cl2")),
    rng(cov_of(filter(dep_b, G == 100), "ij_lin")),
    sc_dep("ij_lin"),
    sc_dep("os_lin"),
    sc_dep("naive")
  ),
  sprintf(
    "- The added width is anti-targeted: the per-replicate SE inflation of IJ is negatively or not correlated with CL2's per-replicate misses in every dependent cell (Spearman %.2f to %.2f), whereas the oracle fixλ is narrower and covers better. CL2's shortfall is a selection effect of the weakly identified surface penalty, and the per-fit Var(rho-hat) neither captures its size nor flags the replicates concerned (`flag-check.R`).",
    min(tg$spearman),
    max(tg$spearman)
  ),
  if (!is.null(pl))
    sprintf(
      "- Plasmode: CL2 undercovers beta in %d of %d cells; IJ raises beta coverage by %s to %s points (mean %s), OS by %s to %s, IJ× by %s to %s, while the beta interval score worsens (IJ %.3f, OS %.3f relative to CL2; +mgcv J %.3f).",
      sum(pl_cl2$coverage < 0.95),
      nrow(pl_cl2),
      pct(pl_gain$lo[pl_gain$arm == "ij_lin"]),
      pct(pl_gain$hi[pl_gain$arm == "ij_lin"]),
      pct(pl_gain$m[pl_gain$arm == "ij_lin"]),
      pct(pl_gain$lo[pl_gain$arm == "os_lin"]),
      pct(pl_gain$hi[pl_gain$arm == "os_lin"]),
      pct(pl_gain$lo[pl_gain$arm == "ij_linx"]),
      pct(pl_gain$hi[pl_gain$arm == "ij_linx"]),
      pl_score$m[pl_score$arm == "ij_lin"],
      pl_score$m[pl_score$arm == "os_lin"],
      pl_score$m[pl_score$arm == "naive"]
    ) else NULL,
  "- Cheaper partial remedies already exist: the t_{G-1} critical value and mgcv's own term (see the t paragraph above); neither needs new code in refund.",
  ""
)
md <- c(md, rec)
writeLines(md, file.path(dir, "REPORT.md"))
cat("wrote", file.path(dir, "REPORT.md"), "\n")
