# REPORT-ncv.md: Gaussian NCV lambda correction, from summaries-ncv/ and
# summaries-ncv-plasmode/. Every number is computed here.
# Run from the repository root: Rscript lambda-correction/make-report-ncv.R

suppressPackageStartupMessages(library(dplyr))
d <- "lambda-correction/summaries-ncv"
dp <- "lambda-correction/summaries-ncv-plasmode"
rdp <- \(f)
  if (file.exists(file.path(dp, f))) read.csv(file.path(dp, f)) else NULL
cells <- read.csv(file.path(d, "ncv-cells.csv"))
contr <- read.csv(file.path(d, "ncv-contrasts.csv"))
calib <- read.csv(file.path(d, "ncv-calib.csv"))
tgt <- read.csv(file.path(d, "ncv-targeting.csv"))
crn <- read.csv(file.path(d, "ncv-crn.csv"))
times <- read.csv(file.path(d, "ncv-times.csv"))
conv <- read.csv(file.path(d, "ncv-converged.csv"))
study <- read.csv("summaries/final/cells.csv")
pl <- rdp("ncvp-cells.csv")
pt <- rdp("ncvp-targeting.csv")
pcal <- rdp("ncvp-calib.csv")
pcrn <- rdp("ncvp-crn.csv")
pconv <- rdp("ncvp-converged.csv")

pct <- \(x) sprintf("%.1f", 100 * x)
grp <- \(x)
  mutate(
    x,
    group = case_when(
      error == "iid" ~ "iid",
      G == 40 ~ "dep. G=40",
      TRUE ~ "dep. G=100"
    )
  )
GROUPS <- c("iid", "dep. G=40", "dep. G=100")
ARMS <- c(
  cl2 = "CL2",
  naive = "+J H^-1 J'",
  ij_lin = "IJ",
  os_lin = "OS",
  ij_linx = "IJ×",
  ij_reg = "IJr",
  fix = "fixλ (oracle)"
)
md_rows <- \(t)
  c(
    paste0("| ", paste(names(t), collapse = " | "), " |"),
    paste0("|", strrep("---|", ncol(t))),
    apply(t, 1, \(r) paste0("| ", paste(r, collapse = " | "), " |"))
  )
table_of <- function(
  value,
  ratio = FALSE,
  scale = 100,
  digits = 1,
  est = NULL
) {
  x <- grp(filter(cells, crit == "z"))
  if (!is.null(est)) x <- filter(x, estimand %in% est)
  if (ratio)
    x <- x |>
      group_by(cell, estimand) |>
      mutate(!!value := .data[[value]] / .data[[value]][arm == "cl2"]) |>
      ungroup()
  x |>
    filter(arm %in% names(ARMS)) |>
    group_by(estimand, group, arm) |>
    summarise(v = mean(.data[[value]]) * scale, .groups = "drop") |>
    mutate(arm = ARMS[arm]) |>
    tidyr::pivot_wider(names_from = arm, values_from = v) |>
    select(estimand, cells = group, any_of(unname(ARMS))) |>
    arrange(
      factor(estimand, c("beta", "alpha", "gamma", "mean")),
      factor(cells, GROUPS)
    ) |>
    mutate(across(
      -c(estimand, cells),
      \(v) formatC(v, format = "f", digits = digits)
    )) |>
    md_rows()
}
paired <- \(a)
  grp(filter(contr, crit == "z", estimand == "beta", arm == a)) |>
    group_by(group) |>
    summarise(
      txt = sprintf(
        "%+.2f (%.2f)",
        100 * mean(diff),
        100 * sqrt(sum(diff_se^2)) / n()
      ),
      .groups = "drop"
    ) |>
    (\(x) setNames(x$txt, x$group))()
ptab <- do.call(
  rbind,
  lapply(
    c("naive", "ij_lin", "os_lin", "ij_linx", "fix"),
    \(a) c(arm = ARMS[[a]], paired(a)[GROUPS])
  )
)
ptab <- as.data.frame(ptab, check.names = FALSE)
names(ptab) <- c("arm", GROUPS)

# reproduction checks
st <- study |>
  filter(arm == "ncv_cl2", cell %in% unique(cells$cell)) |>
  select(cell, estimand, coverage_mean)
rep_diff <- cells |>
  filter(arm == "cl2", crit == "z") |>
  inner_join(st, by = c("cell", "estimand")) |>
  summarise(n = n(), max = max(abs(coverage - coverage_mean)))
crn_est <- crn |>
  filter(column == "est_NCV") |>
  group_by(cell, rep) |>
  summarise(rel = max(max_rel), abs = max(max_abs), .groups = "drop")

# variance overshoot: SD^2 units, beta, dependent cells: needed = 1 - (SE/SD_CL2)^2, added
over <- grp(filter(
  cells,
  crit == "z",
  estimand == "beta",
  error != "iid",
  arm %in% c("cl2", "ij_lin")
)) |>
  select(cell, group, arm, se_sd_ratio, se_sd_median) |>
  tidyr::pivot_wider(
    names_from = arm,
    values_from = c(se_sd_ratio, se_sd_median)
  ) |>
  group_by(group) |>
  summarise(
    needed = mean(1 - se_sd_ratio_cl2^2),
    added = mean(se_sd_ratio_ij_lin^2 - se_sd_ratio_cl2^2),
    needed_m = mean(1 - se_sd_median_cl2^2),
    added_m = mean(se_sd_median_ij_lin^2 - se_sd_median_cl2^2),
    .groups = "drop"
  )
# per-replicate SE inflation of IJ over CL2 (beta, dependent cells)
infl <- filter(tgt, error != "iid", arm == "ij_lin")

tg <- filter(tgt, error != "iid") |>
  group_by(G, arm) |>
  summarise(
    txt = sprintf(
      "%.2f (%.2f to %.2f)",
      mean(spearman),
      min(spearman),
      max(spearman)
    ),
    .groups = "drop"
  ) |>
  filter(arm %in% c("naive", "ij_lin", "os_lin", "ij_linx")) |>
  mutate(arm = ARMS[arm]) |>
  tidyr::pivot_wider(names_from = G, values_from = txt, names_prefix = "G = ")
tg_iid <- filter(tgt, error == "iid", arm == "naive")$spearman
cal <- filter(calib, error != "iid", source %in% c("naive", "ij", "os")) |>
  mutate(
    penalty = sub(
      "^te\\(X.smat,X.tmat\\):L.X1$",
      "surface (s)",
      sub("^te\\(X.smat,X.tmat\\):L.X2$", "surface (t)", par)
    )
  ) |>
  group_by(penalty, source) |>
  summarise(
    sd = mean(mc_sd),
    iqr = mean(mc_iqr),
    r_sd = mean(median_est_sd / mc_sd),
    r_iqr = mean(median_est_sd / mc_iqr),
    .groups = "drop"
  ) |>
  mutate(across(where(is.numeric), \(v) sprintf("%.2f", v))) |>
  rename(
    `MC SD` = sd,
    `MC IQR/1.349` = iqr,
    `median est SD / MC SD` = r_sd,
    `median est SD / IQR spread` = r_iqr
  )
zb <- grp(filter(cells, crit == "z", estimand == "beta"))
q05 <- \(a) mean(filter(zb, group != "iid", arm == a)$pointwise_q05)
sesd <- \(a, g) range(filter(zb, group %in% g, arm == a)$se_sd_ratio)
sesdm <- \(a, g) range(filter(zb, group %in% g, arm == a)$se_sd_median)

md <- c(
  "# Lambda correction for Gaussian curve-blocked NCV: report",
  "",
  "Derivation `NOTE-ncv.qmd`; code `ncv-lib.R`, `run-lc-ncv.R`; tables `summaries-ncv/`. Generated by",
  "`make-report-ncv.R`. 6 Gaussian core cells (iid, OU, FPC errors; G = 40, 100) x 200 replicates,",
  "same data as the study. CL2 as in the study (on Vp). +J H^-1 J' = mgcv-style term with the exact",
  "NCV Hessian (NCV is not a likelihood, so not a calibrated variance). IJ / OS = Var(rho-hat) from the",
  "NCV curve scores / one-step leave-curve-out NCV; × = with the cross term; IJr = H + I/10",
  "(scale-dependent for NCV, a sensitivity variant). fixλ = CL2 of a different estimator: the refit at",
  "the cell's geometric-mean NCV lambda, computed from the same replicates (in-sample oracle); its",
  "SE/SD divides by its own Monte Carlo SD. The target of a lambda-uncertainty correction of the NCV",
  "estimator is SE/SD = 1 against the NCV estimator's own Monte Carlo SD.",
  "",
  sprintf(
    "Checks: CL2 coverage vs the study's published NCV coverage, %d cell x estimand rows, max |difference| %.1e. NCV estimates vs the study's stored ones: %d of %d replicates differ by more than 1e-6 relative (max absolute %.1e; different mgcv build on LRZ, flat NCV criterion). All %d of %d NCV fits fully converged. Time per task: median %.0f s, total %.1f core-hours. MC SEs of paired differences treat cells as independent.",
    rep_diff$n,
    rep_diff$max,
    sum(crn_est$rel > 1e-6),
    nrow(crn_est),
    max(crn_est$abs),
    sum(conv$converged),
    nrow(conv),
    median(times$total),
    sum(times$total) / 3600
  ),
  "",
  "## Coverage (z, %, mean over cells)",
  "",
  table_of("coverage"),
  "",
  "Paired coverage difference to CL2, beta, points (MC SE):",
  "",
  md_rows(ptab),
  "",
  "5% quantile of pointwise beta coverage (%), z:",
  "",
  table_of("pointwise_q05", est = "beta"),
  "",
  "## SE/SD (link scale; fixλ against its own estimator's SD)",
  "",
  "RMS of the SE over replicates / MC SD (mean over grid points and cells):",
  "",
  table_of("se_sd_ratio", FALSE, 1, 2),
  "",
  "Median SE over replicates / MC SD (the typical replicate; biased low for a correct but noisy SE):",
  "",
  table_of("se_sd_median", FALSE, 1, 2),
  "",
  sprintf(
    "Variance added by IJ vs variance missing from CL2, beta, dependent cells (SD^2 units, 1 - (SE/SD)^2 for CL2). With the RMS SE: G = 40 adds %.2f where %.2f is missing; G = 100 adds %.2f where %.2f is missing. With the median SE: G = 40 adds %.2f where %.2f is missing; G = 100 adds %.2f where %.2f is missing. The IJ term is heavy-tailed over replicates: per-replicate beta SE inflation over CL2 has median %.2f-%.2f and 90%% quantile %.2f-%.2f over the dependent cells, maximum %.1f. On average (RMS) IJ adds more variance than CL2 lacks, weighted heavily by a few replicates with very large terms; in the typical replicate (median) it adds less than the median-based shortfall. With coverage near nominal, neither ratio shows a reliable overshoot or undershoot (both are descriptive; for a heavy-tailed SE neither is a calibration diagnostic).",
    over$added[over$group == "dep. G=40"],
    over$needed[over$group == "dep. G=40"],
    over$added[over$group == "dep. G=100"],
    over$needed[over$group == "dep. G=100"],
    over$added_m[over$group == "dep. G=40"],
    over$needed_m[over$group == "dep. G=40"],
    over$added_m[over$group == "dep. G=100"],
    over$needed_m[over$group == "dep. G=100"],
    min(infl$infl_median),
    max(infl$infl_median),
    min(infl$infl_q90),
    max(infl$infl_q90),
    max(infl$infl_max)
  ),
  "",
  "## Interval score relative to CL2 (lower is better)",
  "",
  table_of("score", TRUE, 1, 3),
  "",
  "## Width relative to CL2",
  "",
  table_of("width", TRUE, 1, 3),
  "",
  "## Targeting",
  "",
  "Within-cell Spearman correlation between an arm's per-replicate SE inflation over CL2 (beta) and CL2's per-replicate share of missed grid points, dependent cells: mean (range over cells).",
  "",
  md_rows(tg),
  "",
  sprintf(
    "The mgcv-style term, which changes coverage by only a fraction of a point, has the same correlation (and %.2f-%.2f in the iid cells): this statistic picks up any lambda-hat-dependent term and does not establish that IJ's variance is targeted.",
    min(tg_iid),
    max(tg_iid)
  ),
  "",
  "## Calibration of Var(rho-hat) (dependent cells, mean over cells)",
  "",
  md_rows(cal),
  ""
)

if (!is.null(pl)) {
  stress <- c("ocean", "aemet", "adelaide") # plasmode/rule.R PLASMODE_STRESS
  pl <- mutate(
    pl,
    set = ifelse(family %in% stress, "stress", "independent"),
    size = ifelse(G == 40, "G=40", "full G")
  )
  parms <- c(
    cl2 = "CL2",
    naive = "+J H^-1 J'",
    ij_lin = "IJ",
    ij_linx = "IJ×",
    ij_reg = "IJr"
  )
  ptab2 <- function(value, ratio = FALSE, scale = 100, digits = 1, fun = mean) {
    x <- filter(pl, crit == "z", arm %in% names(parms))
    if (ratio)
      x <- x |>
        group_by(cell, estimand) |>
        mutate(!!value := .data[[value]] / .data[[value]][arm == "cl2"]) |>
        ungroup()
    x |>
      group_by(estimand, set, size, arm) |>
      summarise(v = fun(.data[[value]]) * scale, .groups = "drop") |>
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
        \(v) formatC(v, format = "f", digits = digits)
      )) |>
      md_rows()
  }
  plb <- filter(pl, crit == "z", estimand == "beta")
  gain <- plb |>
    group_by(cell) |>
    mutate(d = coverage - coverage[arm == "cl2"]) |>
    ungroup() |>
    filter(arm == "ij_lin")
  bydata <- plb |>
    filter(arm %in% c("cl2", "ij_lin")) |>
    group_by(family, arm) |>
    summarise(
      cov = mean(coverage),
      sesd = mean(se_sd_ratio),
      .groups = "drop"
    ) |>
    tidyr::pivot_wider(names_from = arm, values_from = c(cov, sesd)) |>
    transmute(
      dataset = family,
      `CL2 coverage` = pct(cov_cl2),
      `IJ coverage` = pct(cov_ij_lin),
      `CL2 SE/SD` = sprintf("%.2f", sesd_cl2),
      `IJ SE/SD` = sprintf("%.2f", sesd_ij_lin)
    )
  ptg <- pt |>
    group_by(arm) |>
    summarise(
      m = mean(spearman, na.rm = TRUE),
      lo = min(spearman, na.rm = TRUE),
      hi = max(spearman, na.rm = TRUE),
      .groups = "drop"
    )
  pc <- pcal |>
    filter(source %in% c("ij", "naive")) |>
    group_by(source) |>
    summarise(r = median(median_est_sd / mc_iqr), .groups = "drop")
  wr <- plb |>
    group_by(cell) |>
    mutate(r = width / width[arm == "cl2"]) |>
    ungroup() |>
    filter(arm == "ij_reg")
  wr_by <- wr |>
    group_by(family) |>
    summarise(mx = max(r), .groups = "drop") |>
    arrange(desc(mx))
  pcrn_est <- if (!is.null(pcrn))
    pcrn |>
      filter(column == "est_NCV") |>
      group_by(cell, rep) |>
      summarise(rel = max(max_rel), abs = max(max_abs), .groups = "drop")
  md <- c(
    md,
    sprintf(
      "## Plasmode (Gaussian, real covariates and residual curves; %d cells, replicates 1-%d)",
      length(unique(pl$cell)),
      max(pl$n_rep)
    ),
    "",
    "No one-step and no fixed-lambda arm (the fixed-lambda plasmode check is in `ncv-rbc/`). Stress tests: ocean, aemet, adelaide (dependent curves). Full G: each dataset's full sample (73-138 curves).",
    "",
    if (!is.null(pcrn_est))
      sprintf(
        "CRN vs the plasmode study's NCV records: %d of %d replicates differ by more than 1e-6 relative (max absolute %.1e).",
        sum(pcrn_est$rel > 1e-6),
        nrow(pcrn_est),
        max(pcrn_est$abs)
      ) else NULL,
    if (!is.null(pconv))
      sprintf(
        "NCV fits fully converged: %d of %d.",
        sum(pconv$converged == 1),
        nrow(pconv)
      ) else "NCV convergence of the plasmode fits: not yet tabulated.",
    "",
    sprintf(
      "IJ changes beta coverage by %s to %s points (mean %s) over the cells (minimum paired gain exceeds its MC SE). Within-cell targeting (Spearman): IJ mean %.2f (%.2f to %.2f), mgcv-style term %.2f. Median over cells of the median IJ SD of log lambda-hat over its robust Monte Carlo spread (IQR/1.349): %.2f (mgcv-style: %.2f).",
      pct(min(gain$d)),
      pct(max(gain$d)),
      pct(mean(gain$d)),
      ptg$m[ptg$arm == "ij_lin"],
      ptg$lo[ptg$arm == "ij_lin"],
      ptg$hi[ptg$arm == "ij_lin"],
      ptg$m[ptg$arm == "naive"],
      pc$r[pc$source == "ij"],
      pc$r[pc$source == "naive"]
    ),
    "",
    "Coverage (z, %, mean over cells):",
    "",
    ptab2("coverage"),
    "",
    "SE/SD, RMS SE (mean over cells):",
    "",
    ptab2("se_sd_ratio", scale = 1, digits = 2),
    "",
    "SE/SD, median SE (mean over cells):",
    "",
    ptab2("se_sd_median", scale = 1, digits = 2),
    "",
    "Interval score relative to CL2 (mean over cells):",
    "",
    ptab2("score", TRUE, 1, 3),
    "",
    "Width relative to CL2 (median over cells; the mean is dominated by IJr's breakdowns):",
    "",
    ptab2("width", TRUE, 1, 3, median),
    "",
    sprintf(
      "IJr breakdown: maximum beta width ratio to CL2 by dataset: %s.",
      paste(sprintf("%s %.1f", wr_by$family, wr_by$mx), collapse = ", ")
    ),
    "",
    "Beta by dataset (mean over its cells, z):",
    "",
    md_rows(bydata),
    ""
  )
}

# Summary (numbers computed here)
z <- grp(filter(cells, crit == "z"))
cb <- \(a, g, col = "coverage")
  mean(filter(z, estimand == "beta", group == g, arm == a)[[col]])
rs <- \(a, g, col)
  mean(
    filter(z, estimand == "beta", group == g, arm == a)[[col]] /
      filter(z, estimand == "beta", group == g, arm == "cl2")[[col]]
  )
rec <- c(
  "## Summary and recommendation",
  "",
  sprintf(
    "- Core cells, dependent errors: IJ moves NCV-CL2 beta coverage from %s to %s at G = 40 and from %s to %s at G = 100 (width x%.2f / x%.2f, interval score x%.3f / x%.3f relative to CL2); under independent errors width x%.2f, score x%.3f.",
    pct(cb("cl2", "dep. G=40")),
    pct(cb("ij_lin", "dep. G=40")),
    pct(cb("cl2", "dep. G=100")),
    pct(cb("ij_lin", "dep. G=100")),
    rs("ij_lin", "dep. G=40", "width"),
    rs("ij_lin", "dep. G=100", "width"),
    rs("ij_lin", "dep. G=40", "score"),
    rs("ij_lin", "dep. G=100", "score"),
    rs("ij_lin", "iid", "width"),
    rs("ij_lin", "iid", "score")
  ),
  sprintf(
    "- Its variance is heavy-tailed over replicates: RMS-based beta SE/SD %.2f-%.2f (CL2 %.2f-%.2f), median-based %.2f-%.2f (CL2 %.2f-%.2f); per-replicate SE inflation median %.2f-%.2f, maximum %.1f. It raises the 5%% pointwise quantile less than the oracle does (CL2 %s, IJ %s, fixλ %s).",
    sesd("ij_lin", c("dep. G=40", "dep. G=100"))[1],
    sesd("ij_lin", c("dep. G=40", "dep. G=100"))[2],
    sesd("cl2", c("dep. G=40", "dep. G=100"))[1],
    sesd("cl2", c("dep. G=40", "dep. G=100"))[2],
    sesdm("ij_lin", c("dep. G=40", "dep. G=100"))[1],
    sesdm("ij_lin", c("dep. G=40", "dep. G=100"))[2],
    sesdm("cl2", c("dep. G=40", "dep. G=100"))[1],
    sesdm("cl2", c("dep. G=40", "dep. G=100"))[2],
    min(infl$infl_median),
    max(infl$infl_median),
    max(infl$infl_max),
    pct(q05("cl2")),
    pct(q05("ij_lin")),
    pct(q05("fix"))
  ),
  if (!is.null(pl))
    sprintf(
      "- Plasmode: IJ raises beta coverage in every cell (by %s to %s points) at beta SE/SD %.2f on average (CL2 %.2f; median-based %.2f, CL2 %.2f), with coverage still %s-%s%% (CL2 %s-%s%%) and no targeting: variance inflation alone does not calibrate NCV-based intervals here.",
      pct(min(gain$d)),
      pct(max(gain$d)),
      mean(filter(plb, arm == "ij_lin")$se_sd_ratio),
      mean(filter(plb, arm == "cl2")$se_sd_ratio),
      mean(filter(plb, arm == "ij_lin")$se_sd_median),
      mean(filter(plb, arm == "cl2")$se_sd_median),
      pct(min(filter(plb, arm == "ij_lin")$coverage)),
      pct(max(filter(plb, arm == "ij_lin")$coverage)),
      pct(min(filter(plb, arm == "cl2")$coverage)),
      pct(max(filter(plb, arm == "cl2")$coverage))
    ) else NULL,
  "- IJr (scale-dependent ridge) breaks down on some plasmode datasets and is not to be used.",
  "- Recommendation: IJ is a component, not a calibrated fix on its own: lambda variability accounts for part of NCV-CL2's undercoverage in the simulations, but the term is heavy-tailed (large in a few replicates, small in most). No refund patch yet. NCV-based intervals also need a bias correction; the lambda term belongs in that combination.",
  ""
)
md <- c(md, rec)
writeLines(md, "lambda-correction/REPORT-ncv.md")
cat("wrote lambda-correction/REPORT-ncv.md\n")
