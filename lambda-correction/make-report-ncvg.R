# REPORT-ncvg.md: NCV lambda correction for Poisson and binomial core cells, from
# summaries-ncvg/. Every number is computed here.
# Run from the repository root: Rscript lambda-correction/make-report-ncvg.R

suppressPackageStartupMessages(library(dplyr))
d <- "lambda-correction/summaries-ncvg"
cells <- read.csv(file.path(d, "ncv-cells.csv"))
contr <- read.csv(file.path(d, "ncv-contrasts.csv"))
calib <- read.csv(file.path(d, "ncv-calib.csv"))
tgt <- read.csv(file.path(d, "ncv-targeting.csv"))
crn <- read.csv(file.path(d, "ncv-crn.csv"))
times <- read.csv(file.path(d, "ncv-times.csv"))
fails <- if (file.size(file.path(d, "ncv-failures.csv")) > 3)
  read.csv(file.path(d, "ncv-failures.csv")) else data.frame()
study <- read.csv("summaries/final/cells.csv")
bnd <- read.csv(file.path(d, "ncv-boundary.csv"))

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
FAMS <- c("poisson", "binomial")
ARMS <- c(
  cl2 = "CL2",
  naive = "+J H^-1 J'",
  ij_lin = "IJ",
  ij_linx = "IJ×",
  fix = "fixλ (geo-mean λ̂)"
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
    group_by(estimand, family, group, arm) |>
    summarise(v = mean(.data[[value]]) * scale, .groups = "drop") |>
    mutate(arm = ARMS[arm]) |>
    tidyr::pivot_wider(names_from = arm, values_from = v) |>
    select(estimand, family, cells = group, any_of(unname(ARMS))) |>
    arrange(
      factor(estimand, c("beta", "alpha", "gamma", "mean")),
      factor(family, FAMS),
      factor(cells, GROUPS)
    ) |>
    mutate(across(
      -c(estimand, family, cells),
      \(v) formatC(v, format = "f", digits = digits)
    )) |>
    md_rows()
}
paired <- \(a)
  grp(filter(contr, crit == "z", estimand == "beta", arm == a)) |>
    group_by(family, group) |>
    summarise(
      txt = sprintf(
        "%+.2f (%.2f)",
        100 * mean(diff),
        100 * sqrt(sum(diff_se^2)) / n()
      ),
      .groups = "drop"
    )
ptab <- bind_rows(lapply(
  c("naive", "ij_lin", "ij_linx", "fix"),
  \(a) mutate(paired(a), arm = ARMS[[a]])
)) |>
  tidyr::pivot_wider(names_from = group, values_from = txt) |>
  arrange(factor(family, FAMS)) |>
  select(family, arm, any_of(GROUPS))

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

# variance overshoot: SD^2 units, beta, by family x group: needed = 1 - (SE/SD_CL2)^2, added
over <- grp(filter(
  cells,
  crit == "z",
  estimand == "beta",
  arm %in% c("cl2", "ij_lin")
)) |>
  select(cell, family, group, arm, se_sd_ratio, se_sd_median) |>
  tidyr::pivot_wider(
    names_from = arm,
    values_from = c(se_sd_ratio, se_sd_median)
  ) |>
  group_by(family, group) |>
  summarise(
    needed = mean(1 - se_sd_ratio_cl2^2),
    added = mean(se_sd_ratio_ij_lin^2 - se_sd_ratio_cl2^2),
    needed_m = mean(1 - se_sd_median_cl2^2),
    added_m = mean(se_sd_median_ij_lin^2 - se_sd_median_cl2^2),
    .groups = "drop"
  ) |>
  arrange(factor(family, FAMS), factor(group, GROUPS)) |>
  transmute(
    family,
    cells = group,
    `missing, RMS SE` = sprintf("%.2f", needed),
    `added by IJ, RMS SE` = sprintf("%.2f", added),
    `missing, median SE` = sprintf("%.2f", needed_m),
    `added by IJ, median SE` = sprintf("%.2f", added_m)
  )

tg <- grp(tgt) |>
  group_by(family, group, arm) |>
  summarise(
    txt = sprintf(
      "%.2f (%.2f to %.2f)",
      mean(spearman),
      min(spearman),
      max(spearman)
    ),
    .groups = "drop"
  ) |>
  filter(arm %in% c("naive", "ij_lin", "ij_linx")) |>
  mutate(arm = ARMS[arm]) |>
  tidyr::pivot_wider(names_from = group, values_from = txt) |>
  arrange(factor(family, FAMS)) |>
  select(family, arm, any_of(GROUPS))
cal <- grp(filter(calib, source %in% c("naive_mgcv", "naive", "ij"))) |>
  mutate(
    penalty = sub(
      "^te\\(X.smat,X.tmat\\):L.X1$",
      "surface (s)",
      sub("^te\\(X.smat,X.tmat\\):L.X2$", "surface (t)", par)
    )
  ) |>
  group_by(family, penalty, source) |>
  summarise(
    sd = mean(mc_sd),
    iqr = mean(mc_iqr),
    r_sd = mean(median_est_sd / mc_sd),
    r_iqr = mean(median_est_sd / mc_iqr),
    .groups = "drop"
  ) |>
  arrange(factor(family, FAMS), penalty, source) |>
  mutate(across(where(is.numeric), \(v) sprintf("%.2f", v))) |>
  rename(
    `MC SD` = sd,
    `MC IQR/1.349` = iqr,
    `median est SD / MC SD` = r_sd,
    `median est SD / IQR spread` = r_iqr
  )

infl <- grp(filter(tgt, arm %in% c("ij_lin", "ij_linx"))) |>
  group_by(family, group, arm) |>
  summarise(
    med = mean(infl_median),
    q90 = mean(infl_q90),
    mx = max(infl_max),
    .groups = "drop"
  ) |>
  arrange(factor(family, FAMS), factor(group, GROUPS), arm) |>
  transmute(
    family,
    cells = group,
    arm = ARMS[arm],
    `median` = sprintf("%.2f", med),
    `90% quantile` = sprintf("%.2f", q90),
    `max` = sprintf("%.1f", mx)
  )

# summary numbers
z <- grp(filter(cells, crit == "z", estimand == "beta"))
cv <- \(a, f, g, col = "coverage")
  mean(filter(z, family == f, group == g, arm == a)[[col]])
rel <- \(a, f, g, col)
  mean(
    filter(z, family == f, group == g, arm == a)[[col]] /
      filter(z, family == f, group == g, arm == "cl2")[[col]]
  )
sumline <- \(f, g)
  sprintf(
    "%s, %s: CL2 %s, IJ %s, IJ× %s, fixλ %s; width IJ x%.2f, IJ× x%.2f; interval score IJ x%.3f, IJ× x%.3f; median SE/SD CL2 %.2f, IJ %.2f, IJ× %.2f.",
    f,
    g,
    pct(cv("cl2", f, g)),
    pct(cv("ij_lin", f, g)),
    pct(cv("ij_linx", f, g)),
    pct(cv("fix", f, g)),
    rel("ij_lin", f, g, "width"),
    rel("ij_linx", f, g, "width"),
    rel("ij_lin", f, g, "score"),
    rel("ij_linx", f, g, "score"),
    cv("cl2", f, g, "se_sd_median"),
    cv("ij_lin", f, g, "se_sd_median"),
    cv("ij_linx", f, g, "se_sd_median")
  )

btab <- grp(bnd) |>
  mutate(
    penalty = sub(
      "^te\\(X.smat,X.tmat\\):L.X1$",
      "surface (s)",
      sub("^te\\(X.smat,X.tmat\\):L.X2$", "surface (t)", par)
    )
  ) |>
  group_by(family, group, penalty) |>
  summarise(rate = pct(mean(rate)), .groups = "drop") |>
  tidyr::pivot_wider(names_from = penalty, values_from = rate) |>
  arrange(factor(family, FAMS), factor(group, GROUPS)) |>
  rename(cells = group)

# findings (computed)
b40 <- filter(grp(bnd), family == "binomial", group == "dep. G=40")
cal_r <- grp(filter(calib, source == "ij")) |>
  group_by(family) |>
  summarise(
    lo = min(median_est_sd / mc_sd),
    hi = max(median_est_sd / mc_sd),
    sd_lo = min(mc_sd),
    sd_hi = max(mc_sd),
    .groups = "drop"
  )
tg_r <- filter(tgt, error != "iid") |>
  group_by(arm) |>
  summarise(lo = min(spearman), hi = max(spearman), .groups = "drop")
findings <- c(
  sprintf(
    "- fixλ in the binomial G = 40 dependent cells: coverage %s at median SE/SD %.2f against its own SD and width x%.2f: the geometric mean of a lambda-hat distribution with frequent boundary solutions (boundary rate per penalty %s%% on average, up to %s%% in a cell) oversmooths, so its undercoverage is bias. It is no benchmark there.",
    pct(cv("fix", "binomial", "dep. G=40")),
    cv("fix", "binomial", "dep. G=40", "se_sd_median"),
    rel("fix", "binomial", "dep. G=40", "width"),
    pct(mean(b40$rate)),
    pct(max(b40$rate))
  ),
  sprintf(
    "- Var(rho-hat): median IJ SD / MC SD of log lambda-hat %.2f-%.2f (Poisson) and %.2f-%.2f (binomial; MC SD %.1f-%.1f on the log scale): for binomial, lambda-hat jumps between interior and boundary solutions, which a first-order term cannot represent; boundary components get zero variance by construction.",
    cal_r$lo[cal_r$family == "poisson"],
    cal_r$hi[cal_r$family == "poisson"],
    cal_r$lo[cal_r$family == "binomial"],
    cal_r$hi[cal_r$family == "binomial"],
    min(filter(grp(calib), family == "binomial", source == "ij")$mc_sd),
    max(filter(grp(calib), family == "binomial", source == "ij")$mc_sd)
  ),
  sprintf(
    "- Targeting (dependent cells, Spearman): IJ %.2f to %.2f, mgcv-style term %.2f to %.2f: the gain comes from the size of the term, not from better targeting.",
    tg_r$lo[tg_r$arm == "ij_lin"],
    tg_r$hi[tg_r$arm == "ij_lin"],
    tg_r$lo[tg_r$arm == "naive"],
    tg_r$hi[tg_r$arm == "naive"]
  ),
  sprintf(
    "- Independent errors: CL2 already overcovers beta (%s / %s); IJ costs interval score x%.3f / x%.3f (Poisson / binomial).",
    pct(cv("cl2", "poisson", "iid")),
    pct(cv("cl2", "binomial", "iid")),
    rel("ij_lin", "poisson", "iid", "score"),
    rel("ij_lin", "binomial", "iid", "score")
  )
)

md <- c(
  "# Lambda correction for curve-blocked NCV, Poisson and binomial: report",
  "",
  "Derivation `NOTE-ncvg.qmd`; code `ncvg-lib.R`, `run-lc-ncvg.R`; tables `summaries-ncvg/`. Generated by",
  "`make-report-ncvg.R`. 12 core cells (Poisson, binomial; iid, OU, FPC errors; G = 40, 100; mid signal) x",
  "200 replicates, same data as the study. CL2 as in the study (on Vp). +J H^-1 J' = mgcv-style term with",
  "the exact NCV Hessian (not a calibrated variance). IJ = Var(rho-hat) from the NCV curve scores; × = with",
  "the cross term. fixλ = CL2 of a different estimator: the refit at the cell's geometric-mean NCV lambda,",
  "computed from the same replicates (a hindsight reference, not an oracle; it breaks down where lambda-hat",
  "is often on a boundary, see below); its SE/SD divides by its own Monte Carlo SD. Nominal level 95%;",
  "z = normal critical value (t_{G-1} rows are in the CSV). Coverage and SE/SD on the link scale. Primary",
  "criteria: coverage and interval score; the SE/SD ratios are descriptive (neither is a calibration",
  "diagnostic for an SE whose inflation is heavy-tailed over replicates).",
  "",
  sprintf(
    "Checks: CL2 coverage vs the study's published NCV coverage, %d cell x estimand rows, max |difference| %.1e. NCV estimates vs the study's stored ones: %d of %d replicates differ by more than 1e-6 relative (max absolute %.1e). Failed tasks: %d. NCV fits fully converged: %d of %d. Time per task: median %.0f s, total %.1f core-hours. MC SEs of paired differences treat cells as independent.",
    rep_diff$n,
    rep_diff$max,
    sum(crn_est$rel > 1e-6),
    nrow(crn_est),
    max(crn_est$abs),
    nrow(fails),
    sum(times$converged),
    nrow(times),
    median(times$total),
    sum(times$total) / 3600
  ),
  "",
  "## Summary (beta, z, mean over cells)",
  "",
  paste(
    "-",
    c(
      sumline("poisson", "dep. G=40"),
      sumline("poisson", "dep. G=100"),
      sumline("binomial", "dep. G=40"),
      sumline("binomial", "dep. G=100"),
      sumline("poisson", "iid"),
      sumline("binomial", "iid")
    )
  ),
  "",
  findings,
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
  "## SE/SD (fixλ against its own estimator's SD)",
  "",
  "RMS of the SE over replicates / MC SD of the estimate (mean over grid points and cells):",
  "",
  table_of("se_sd_ratio", FALSE, 1, 2),
  "",
  "Median SE over replicates / MC SD (the typical replicate; biased low for a correct but noisy SE):",
  "",
  table_of("se_sd_median", FALSE, 1, 2),
  "",
  "Variance added by IJ vs variance missing from CL2, beta (SD^2 units, 1 - (SE/SD)^2 for CL2, mean over cells). The two metrics disagree in direction; with coverage near nominal, neither shows a reliable overshoot or undershoot:",
  "",
  md_rows(over),
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
  "Within-cell Spearman correlation between an arm's per-replicate SE inflation over CL2 (beta) and CL2's per-replicate share of missed grid points: mean (range over cells). Any lambda-hat-dependent term can correlate; compare with the mgcv-style row.",
  "",
  md_rows(tg),
  "",
  "Per-replicate SE inflation over CL2 (beta; mean over grid points of SE / SE_CL2), distribution over replicates: median and 90% quantile averaged over cells, maximum over cells. The IJ term is small in most replicates and very large in a few, which the RMS-based SE/SD weights heavily.",
  "",
  md_rows(infl),
  "",
  "## Calibration of Var(rho-hat) (mean over cells; fails for binomial, see Summary)",
  "",
  md_rows(cal),
  "",
  "## Boundary smoothing parameters",
  "",
  "Share of fits (%) in which a smoothing parameter is on a flat boundary of the NCV criterion (curvature < 1e-3 x the largest; treated as fixed in IJ, `NOTE-ncvg.qmd`):",
  "",
  md_rows(btab),
  ""
)
writeLines(md, "lambda-correction/REPORT-ncvg.md")
cat("wrote lambda-correction/REPORT-ncvg.md\n")
