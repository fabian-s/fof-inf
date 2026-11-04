# Plotting functions for paper.qmd, app/ and supplement.qmd (sourced in the setup chunk
# of those documents). Plots existing results objects; computes no new results.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

# ---- common encoding -------------------------------------------------------------------------
# One colour and one shape per interval in both documents, so that intervals stay distinct in
# grayscale and under colour-vision deficiency. Filled = default; open = second level of a factor
# named in the caption.
ARM_STYLE <- tibble::tribble(
  ~key,
  ~col,
  ~shape,
  "model",
  "#666666",
  21L,
  "cl2",
  "#0072B2",
  21L,
  "ncv",
  "#D55E00",
  24L,
  "wide",
  "#CC79A7",
  23L,
  "widef",
  "#882255",
  25L,
  "ar1",
  "#009E73",
  22L,
  "hyb",
  "#56B4E9",
  23L
)
arm_key <- function(x) {
  x <- as.character(x)
  dplyr::case_when(
    grepl("freq", x) & grepl("bias|widen|NCV–REML", x) ~ "widef",
    grepl("bias|widen|NCV–REML", x) ~ "wide",
    grepl("AR\\(1\\)", x) ~ "ar1",
    grepl("REML \\+ CL2 width|half-width", x) ~ "hyb",
    grepl("model", x) ~ "model",
    grepl("^NCV", x) ~ "ncv",
    TRUE ~ "cl2"
  )
}
arm_col <- function(x) ARM_STYLE$col[match(arm_key(x), ARM_STYLE$key)]
# model-based intervals: a solid grey square, every other interval its own filled symbol
arm_shape <- function(x) {
  s <- ARM_STYLE$shape[match(arm_key(x), ARM_STYLE$key)]
  s[arm_key(x) == "model"] <- 22L
  s[arm_key(x) == "ar1"] <- 25L
  s
}
scale_arm <- function(levels, name = NULL, nrow = NULL) {
  v_col <- setNames(arm_col(levels), levels)
  v_shape <- setNames(arm_shape(levels), levels)
  list(
    scale_colour_manual(
      values = v_col,
      name = name,
      breaks = levels,
      drop = FALSE
    ),
    scale_shape_manual(
      values = v_shape,
      name = name,
      breaks = levels,
      drop = FALSE
    ),
    scale_fill_manual(
      values = v_col,
      name = name,
      breaks = levels,
      drop = FALSE
    ),
    if (!is.null(nrow))
      guides(
        colour = guide_legend(nrow = nrow),
        shape = guide_legend(nrow = nrow),
        fill = guide_legend(nrow = nrow)
      )
  )
}
FIT_COL <- c(REML = "#0072B2", NCV = "#D55E00")
FIT_SHAPE <- c(REML = 21L, NCV = 24L)
FAM_COL <- c(
  Gaussian = "#000000",
  Poisson = "#56B4E9",
  binary = "#CC79A7",
  `scaled t` = "#009E73",
  beta = "#E69F00",
  `negative binomial` = "#999933"
)
FAM_SHAPE <- c(
  Gaussian = 21L,
  Poisson = 24L,
  binary = 22L,
  `scaled t` = 23L,
  beta = 25L,
  `negative binomial` = 3L
)
REF <- "grey35"
theme_fig <- function(base = 9) {
  theme_bw(base_size = base) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey90", linewidth = 0.3),
      panel.spacing = unit(0.35, "lines"),
      legend.position = "bottom",
      legend.margin = margin(0, 0, 0, 0),
      legend.box.margin = margin(-3, 0, 0, 0),
      legend.key.height = unit(0.32, "cm"),
      legend.key.width = unit(0.4, "cm"),
      legend.spacing.y = unit(0, "cm"),
      strip.background = element_rect(fill = "grey93", colour = NA),
      strip.text = element_text(size = base - 1, margin = margin(2, 2, 2, 2)),
      axis.text = element_text(size = base - 1.5, colour = "grey20"),
      axis.title = element_text(size = base - 0.5)
    )
}
num_lab <- function(x) sub("\\.?0+$", "", formatC(x, format = "f", digits = 3))
DODGE <- function(w = 0.72) position_dodge(width = w, orientation = "y")
# horizontal point + interval bar, dodged along y
hpoint <- function(
  xmin = NULL,
  xmax = NULL,
  w = 0.72,
  size = 1.5,
  stroke = 0.35,
  lw = 0.4,
  open = NULL
) {
  c(
    if (!rlang::quo_is_null(rlang::enquo(xmin)))
      list(geom_linerange(
        aes(xmin = {{ xmin }}, xmax = {{ xmax }}),
        position = DODGE(w),
        linewidth = lw,
        orientation = "y"
      )),
    list(geom_point(position = DODGE(w), size = size, stroke = stroke))
  )
}
vref <- function(x, lt = 2)
  geom_vline(xintercept = x, linetype = lt, colour = REF, linewidth = 0.35)
COV_BREAKS <- function(by = 0.1) seq(0, 1, by)

# ---- main text: estimation ---------------------------------------------------------------------
# scatter of per-setting relative errors with equal axes; open symbols for one level of `open_if`
LOG_REL <- c(0.05, 0.1, 0.2, 0.5, 1, 2, 5)
# solid_alpha: transparency of the solid symbols' fill only (outlines and open symbols stay opaque);
# NULL keeps one alpha for the whole symbol.
rel_scatter <- function(
  d,
  x,
  y,
  open_if = FALSE,
  xlab = "relative error, NCV",
  ylab = "relative error, REML",
  size = 1.5,
  solid_alpha = NULL
) {
  solid <- \(col)
    if (is.null(solid_alpha)) col else scales::alpha(col, solid_alpha)
  d <- d |>
    mutate(
      family = factor(as.character(family), names(FAM_COL)),
      .fill = if_else(
        rep_len({{ open_if }}, n()),
        "white",
        solid(unname(FAM_COL[as.character(family)]))
      )
    )
  lim <- range(c(pull(d, {{ x }}), pull(d, {{ y }})))
  ggplot(d, aes({{ x }}, {{ y }}, colour = family, shape = family)) +
    geom_abline(linetype = 2, colour = REF, linewidth = 0.35) +
    geom_hline(yintercept = 1, linetype = 3, colour = REF, linewidth = 0.4) +
    geom_point(
      aes(fill = I(.fill)),
      size = size,
      stroke = 0.4,
      alpha = if (is.null(solid_alpha)) 0.9 else NA
    ) +
    scale_x_log10(
      breaks = LOG_REL,
      labels = num_lab,
      limits = lim,
      minor_breaks = NULL
    ) +
    scale_y_log10(
      breaks = LOG_REL,
      labels = num_lab,
      limits = lim,
      minor_breaks = NULL
    ) +
    scale_colour_manual(values = FAM_COL, name = NULL, drop = TRUE) +
    scale_shape_manual(values = FAM_SHAPE, name = NULL, drop = TRUE) +
    coord_fixed() +
    labs(x = xlab, y = ylab) +
    theme_fig()
}
FAM_PANEL <- c(
  Gaussian = "Gaussian",
  Poisson = "Poisson",
  binary = "binary",
  `scaled t` = "scaled t, beta, negative binomial",
  beta = "scaled t, beta, negative binomial",
  `negative binomial` = "scaled t, beta, negative binomial"
)
fig_estimation <- function(fig2a) {
  d <- fig2a |>
    mutate(panel = factor(FAM_PANEL[as.character(family)], unique(FAM_PANEL)))
  rel_scatter(
    d,
    ncv_cl2,
    reml_cl2,
    truth == "rough truth",
    size = 1.1,
    solid_alpha = 0.65
  ) +
    facet_wrap(~panel, nrow = 1, labeller = label_wrap_gen(20)) +
    guides(
      colour = guide_legend(
        nrow = 1,
        override.aes = list(
          size = 1.8,
          fill = scales::alpha(unname(FAM_COL), 0.65)
        )
      ),
      shape = guide_legend(nrow = 1)
    )
}

# ---- main text: coverage in the main factorial ---------------------------------------------------
fig_coverage <- function(fig1) {
  d <- fig1 |>
    mutate(
      G0 = as.character(G),
      panel = factor(paste0(as.character(estimand), ", ", G0)),
      panel = factor(
        panel,
        paste0(
          rep(levels(droplevels(estimand)), each = 2),
          ", ",
          c("G = 40", "G = 100")
        )
      )
    )
  lv <- levels(d$arm)
  ggplot(d, aes(coverage_mean, error, colour = arm, shape = arm, fill = arm)) +
    vref(LEVEL) +
    hpoint(
      coverage_mean - 2 * coverage_se,
      coverage_mean + 2 * coverage_se,
      w = 0.92
    ) +
    facet_grid(panel ~ family) +
    scale_x_continuous(breaks = COV_BREAKS(0.1), labels = num_lab) +
    scale_y_discrete(expand = expansion(add = 0.55)) +
    scale_arm(lv) +
    labs(x = "coverage (grid average)", y = NULL) +
    theme_fig()
}

# ---- main text: Satterthwaite against the normal critical value ----------------------------------
# dumbbells within family x number of curves; the width ratio is printed at the right edge
fig_satt <- function(m_sat_tab) {
  d <- m_sat_tab |>
    pivot_longer(
      -c(family, estimand),
      names_to = c(".value", "G"),
      names_pattern = "(cov_z|cov_satt|wr)_(\\d+)"
    ) |>
    mutate(
      G = factor(paste0("G = ", G), paste0("G = ", c(20, 40, 100))),
      estimand = factor(
        as.character(estimand),
        rev(c("E(Y | X)", "β(s,t)", "α(t)", "γ(t)"))
      )
    )
  ggplot(d, aes(y = estimand)) +
    vref(LEVEL) +
    geom_segment(
      aes(x = cov_z, xend = cov_satt, yend = estimand),
      colour = FIT_COL[["REML"]],
      linewidth = 0.5
    ) +
    geom_point(
      aes(cov_z, shape = "normal"),
      colour = FIT_COL[["REML"]],
      fill = "white",
      size = 1.7,
      stroke = 0.5
    ) +
    geom_point(
      aes(cov_satt, shape = "Satterthwaite"),
      colour = FIT_COL[["REML"]],
      fill = FIT_COL[["REML"]],
      size = 1.7,
      stroke = 0.5
    ) +
    geom_text(
      aes(x = 0.998, label = sprintf("%.2f", wr)),
      size = 2.3,
      hjust = 1,
      colour = "grey25"
    ) +
    facet_grid(G ~ family) +
    scale_x_continuous(
      limits = c(0.85, 1),
      breaks = seq(0.85, 0.95, 0.05),
      labels = num_lab
    ) +
    scale_shape_manual(
      values = c(normal = 21, Satterthwaite = 21),
      name = "critical value"
    ) +
    guides(
      shape = guide_legend(
        override.aes = list(
          fill = c("white", FIT_COL[["REML"]]),
          colour = FIT_COL[["REML"]],
          size = 1.8
        )
      )
    ) +
    labs(x = "coverage (grid average)", y = NULL) +
    theme_fig()
}

# ---- dot-plot grid: one panel per performance measure --------------------------------------------
MEASURES <- c(
  rel = "relative error",
  coverage = "coverage",
  q05 = "lower tail",
  hw = "half-width",
  det = "detection",
  score = "interval-score ratio"
)
measure_refs <- function(levels_m) {
  tibble::tibble(
    measure = factor(c("coverage", "interval-score ratio"), levels_m),
    x = c(LEVEL, 1)
  )
}
fig_tradeoff_summary <- function(trade_tab) {
  d <- trade_tab |>
    rename(score = score_ratio) |>
    mutate(
      strip = factor(
        sprintf("%s\n%d settings", as.character(estimand), settings),
        sprintf(
          "%s\n%d settings",
          levels(estimand),
          tapply(settings, estimand, first)[levels(estimand)]
        )
      ),
      arm = factor(as.character(arm), rev(levels(arm)))
    ) |>
    select(strip, arm, all_of(names(MEASURES))) |>
    pivot_longer(all_of(names(MEASURES)), names_to = "m", values_to = "x") |>
    filter(!is.na(x)) |>
    mutate(measure = factor(MEASURES[m], MEASURES))
  lv <- rev(levels(d$arm))
  ggplot(d, aes(x, arm, colour = arm, shape = arm, fill = arm)) +
    geom_vline(
      data = measure_refs(MEASURES),
      aes(xintercept = x),
      linetype = 2,
      colour = REF,
      linewidth = 0.35
    ) +
    geom_point(size = 1.8, stroke = 0.35) +
    facet_grid(
      strip ~ measure,
      scales = "free_x",
      labeller = labeller(measure = label_wrap_gen(10))
    ) +
    scale_x_continuous(breaks = \(l) base::pretty(l, 3), labels = num_lab) +
    scale_arm(lv) +
    guides(colour = guide_legend(reverse = FALSE)) +
    labs(x = NULL, y = NULL) +
    theme_fig() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank()
    )
}
fig_tradeoff <- function(fig2b) {
  d <- fig2b |>
    mutate(
      .fill = if_else(truth == "rough truth", "white", unname(arm_col(arm))),
      arm = factor(as.character(arm), c("REML + CL2", "NCV + CL2")),
      family = FAMILY_LABELS[as.character(family)],
      panel = factor(FAM_PANEL[family], unique(FAM_PANEL))
    )
  ggplot(d, aes(coverage, detect, colour = arm, shape = arm)) +
    geom_vline(
      xintercept = LEVEL,
      linetype = 2,
      colour = REF,
      linewidth = 0.35
    ) +
    geom_point(aes(fill = I(.fill)), size = 1.5, stroke = 0.45, alpha = 0.85) +
    facet_wrap(~panel, ncol = 2) +
    scale_colour_manual(
      values = setNames(arm_col(levels(d$arm)), levels(d$arm)),
      name = NULL
    ) +
    scale_shape_manual(
      values = setNames(arm_shape(levels(d$arm)), levels(d$arm)),
      name = NULL
    ) +
    scale_x_continuous(breaks = seq(0.75, 1, 0.05), labels = num_lab) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.2),
      labels = num_lab
    ) +
    guides(
      colour = guide_legend(
        override.aes = list(size = 1.8, fill = unname(arm_col(levels(d$arm))))
      )
    ) +
    labs(x = "coverage (grid average)", y = "detection rate") +
    theme_fig()
}

# ---- main text: zero surface --------------------------------------------------------------------
fig_null <- function(SGN, ERR3) {
  arms <- c(
    reml_model = "model-based (REML)",
    reml_cl2 = "REML + CL2, normal",
    reml_cl2_satt = "REML + CL2, Satterthwaite",
    ncv_cl2 = "NCV + CL2"
  )
  col <- c("#666666", "#0072B2", "#56B4E9", "#D55E00")
  shp <- c(22L, 21L, 23L, 24L)
  names(col) <- names(shp) <- arms
  d <- SGN |>
    filter(arm %in% names(arms)) |>
    mutate(
      arm = factor(arms[arm], arms),
      family = factor(FAMILY_LABELS[family], c("Gaussian", "binary")),
      error = factor(ERR3[error], c("smooth", "independent")),
      G = factor(paste0("G = ", G), paste0("G = ", c(100, 40)))
    )
  p1 <- ggplot(d, aes(fpr, G, colour = arm, shape = arm, fill = arm)) +
    geom_vline(
      xintercept = 0.05,
      linetype = 2,
      colour = REF,
      linewidth = 0.35
    ) +
    geom_point(position = DODGE(0.6), size = 1.6, stroke = 0.35) +
    facet_grid(family + error ~ .) +
    scale_x_log10(
      breaks = c(0.01, 0.02, 0.05, 0.1, 0.2, 0.4),
      labels = num_lab
    ) +
    scale_colour_manual(values = col, name = NULL) +
    scale_shape_manual(values = shp, name = NULL) +
    scale_fill_manual(values = col, name = NULL) +
    guides(
      colour = guide_legend(ncol = 2),
      shape = guide_legend(ncol = 2),
      fill = guide_legend(ncol = 2)
    ) +
    labs(x = "pointwise false-positive rate (log scale)", y = NULL) +
    theme_fig() +
    theme(strip.text.y = element_blank())
  ds <- d |>
    filter(arm %in% arms[c("reml_cl2", "ncv_cl2")]) |>
    mutate(
      fit = factor(
        if_else(grepl("^REML", arm), "REML", "NCV"),
        c("REML", "NCV")
      )
    )
  p2 <- ggplot(ds, aes(rms_rel, G, colour = fit, shape = fit, fill = fit)) +
    geom_vline(xintercept = 1, linetype = 3, colour = REF, linewidth = 0.4) +
    geom_point(position = DODGE(0.5), size = 1.6, stroke = 0.35) +
    facet_grid(family + error ~ .) +
    scale_x_log10(breaks = c(0.01, 0.03, 0.1, 0.3, 1, 3), labels = num_lab) +
    scale_colour_manual(values = FIT_COL, name = "estimate") +
    scale_shape_manual(values = FIT_SHAPE, name = "estimate") +
    scale_fill_manual(values = FIT_COL, name = "estimate") +
    labs(x = "root-mean-square size (log scale)", y = NULL) +
    theme_fig() +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  p1 + p2 + plot_layout(widths = c(1.25, 1))
}

# ---- main text: plasmode ---------------------------------------------------------------------------
fig_plasmode_est <- function(fig4) {
  d <- fig4 |>
    mutate(
      group = factor(
        if_else(role == "counted", "indep. source", "dep. source"),
        c("indep. source", "dep. source")
      ),
      truth = factor(
        c(
          `REML-fitted` = "REML-fitted truth",
          intermediate = "midpoint truth",
          `NCV-fitted` = "NCV-fitted truth"
        )[as.character(truth)],
        c("REML-fitted truth", "midpoint truth", "NCV-fitted truth")
      )
    )
  w <- d |>
    select(group, app, truth, fit, rel) |>
    pivot_wider(names_from = fit, values_from = rel)
  ggplot(d, aes(rel, app)) +
    geom_vline(xintercept = 1, linetype = 3, colour = REF, linewidth = 0.4) +
    geom_segment(
      data = w,
      aes(x = REML, xend = NCV, y = app, yend = app),
      colour = "grey60",
      linewidth = 0.5
    ) +
    geom_point(
      aes(colour = fit, shape = fit, fill = fit),
      size = 1.8,
      stroke = 0.35
    ) +
    facet_grid(group ~ truth, scales = "free_y", space = "free_y") +
    scale_x_log10(breaks = c(0.5, 0.7, 1, 2, 3, 5), labels = num_lab) +
    scale_colour_manual(values = FIT_COL, name = "fit") +
    scale_shape_manual(values = FIT_SHAPE, name = "fit") +
    scale_fill_manual(values = FIT_COL, name = "fit") +
    labs(x = "relative error of β(s,t) (log scale)", y = NULL) +
    theme_fig()
}
fig_plasmode <- function(fig3) {
  lv <- levels(fig3$arm)
  ggplot(fig3, aes(coverage, app, colour = arm, shape = arm, fill = arm)) +
    vref(LEVEL) +
    hpoint(coverage_min, coverage_max, w = 0.7, size = 1.5) +
    facet_grid(role ~ estimand, scales = "free_y", space = "free_y") +
    scale_x_continuous(breaks = COV_BREAKS(0.1), labels = num_lab) +
    scale_arm(lv) +
    labs(x = "coverage (grid average)", y = NULL) +
    theme_fig()
}
fig_plasmode_info <- function(ptab) {
  d <- ptab |>
    mutate(
      group = factor(
        c(rep("indep. source", 6), rep("dep. source", 3)),
        c("indep. source", "dep. source")
      ),
      dataset = factor(as.character(dataset), rev(as.character(dataset)))
    )
  long <- bind_rows(
    d |>
      transmute(
        group,
        dataset,
        arm = "REML + CL2",
        rel = rel_reml_cl2,
        coverage = coverage_reml_cl2,
        q05 = q05,
        hw = hw_reml_cl2,
        det = det_reml_cl2
      ),
    d |>
      transmute(
        group,
        dataset,
        arm = "NCV + CL2",
        rel = rel_ncv_cl2,
        coverage = coverage_ncv_cl2,
        det = det_ncv_cl2,
        score = score_ncv
      )
  ) |>
    pivot_longer(
      c(rel, coverage, q05, hw, det, score),
      names_to = "m",
      values_to = "x"
    ) |>
    filter(!is.na(x)) |>
    mutate(arm = factor(arm, c("REML + CL2", "NCV + CL2")))
  lv <- levels(long$arm)
  panel <- function(
    m,
    title,
    breaks,
    log = FALSE,
    ref = NULL,
    first = FALSE,
    last = FALSE
  ) {
    dm <- long |> filter(m == !!m) |> mutate(title = title)
    p <- ggplot(dm, aes(x, dataset)) +
      (if (!is.null(ref))
        geom_vline(
          xintercept = ref,
          linetype = 2,
          colour = REF,
          linewidth = 0.35
        )) +
      geom_line(
        aes(group = dataset),
        colour = "grey65",
        linewidth = 0.5,
        orientation = "y"
      ) +
      geom_point(
        aes(colour = arm, shape = arm, fill = arm),
        size = 1.7,
        stroke = 0.3
      ) +
      facet_grid(group ~ title, scales = "free_y", space = "free_y") +
      (if (log) scale_x_log10(breaks = breaks, labels = num_lab) else
        scale_x_continuous(breaks = breaks, labels = num_lab)) +
      scale_arm(lv) +
      labs(x = NULL, y = NULL) +
      theme_fig() +
      theme(
        plot.margin = margin(2, 1.5, 2, 1.5),
        strip.text.x = element_text(size = 7)
      )
    if (!first)
      p <- p +
        theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
    if (!last)
      p <- p +
        theme(
          strip.text.y = element_blank(),
          strip.background.y = element_blank()
        )
    if (!first) p <- p + guides(colour = "none", shape = "none", fill = "none")
    p
  }
  wrap_plots(
    panel(
      "rel",
      "relative error\n(log scale)",
      c(0.5, 1, 2, 3),
      log = TRUE,
      first = TRUE
    ),
    panel("coverage", "coverage", seq(0.6, 0.9, 0.1), ref = LEVEL),
    panel("q05", "lower\ntail", seq(0.85, 0.9, 0.05)),
    panel("hw", "half-width\n(log scale)", c(1, 2, 5), log = TRUE),
    panel("det", "detection", seq(0.2, 0.8, 0.2)),
    panel("score", "interval-score\nratio", c(1, 1.5), ref = 1, last = TRUE),
    nrow = 1,
    widths = c(1.45, 1.25, 0.8, 0.9, 1.1, 1)
  ) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
}

# ---- Appendix A.6: covariance components -----------------------------------------------------------
S1_LAB <- c(
  "model" = "model-based",
  "CR1" = "CR1",
  "CL2 $I - H$" = "CL2, shortcut",
  "CL2" = "CL2",
  "CL2 freq." = "CL2, frequentist",
  "CL2 $t_{G-1}$" = "CL2, t(G-1)",
  "CL2 fixed $\\lambda$" = "CL2, fixed λ"
)
S1_ERR <- c(iid = "independent", ou = "rough stationary", fpc = "smooth")
s1_long <- function(s1_cells, s1_arms, e) {
  s1_cells |>
    filter(estimand %in% e) |>
    inner_join(s1_arms, by = c("arm", "crit")) |>
    mutate(
      label = factor(S1_LAB[label], rev(unname(S1_LAB))),
      family = label_family(family),
      error = factor(S1_ERR[error], S1_ERR),
      est = label_estimand(estimand),
      G = factor(paste0("G = ", G), paste0("G = ", c(40, 100)))
    ) |>
    droplevels()
}
G_COL <- c(`G = 40` = "#D55E00", `G = 100` = "#0072B2")
G_SHAPE <- c(`G = 40` = 24L, `G = 100` = 21L)
# coverage by interval with G as second level; values off the axis are printed at its left edge
fig_s1_cov <- function(d, rows, lo = 0.84) {
  off <- d |>
    filter(coverage < lo) |>
    arrange(G) |>
    group_by(pick(all_of(c(rows, "family", "label")))) |>
    summarise(
      txt = paste0(G, ": ", sprintf("%.2f", coverage), collapse = "; "),
      .groups = "drop"
    )
  on <- d |> filter(coverage >= lo)
  ggplot(on, aes(coverage, label, colour = G, shape = G, fill = G)) +
    vref(LEVEL) +
    geom_line(
      aes(group = label),
      colour = "grey70",
      linewidth = 0.4,
      orientation = "y"
    ) +
    geom_point(size = 1.5, stroke = 0.3) +
    geom_segment(
      data = off,
      aes(x = lo + 0.012, xend = lo + 0.001, y = label, yend = label),
      inherit.aes = FALSE,
      colour = "grey25",
      linewidth = 0.35,
      arrow = arrow(length = unit(0.05, "in"), type = "closed")
    ) +
    geom_text(
      data = off,
      aes(x = lo + 0.015, y = label, label = txt),
      inherit.aes = FALSE,
      hjust = 0,
      size = 2.3,
      colour = "grey15"
    ) +
    facet_grid(rows = vars(!!!rlang::syms(rows)), cols = vars(family)) +
    scale_x_continuous(
      limits = c(lo, 1),
      expand = expansion(mult = c(0, 0.04)),
      breaks = seq(0.85, 1, 0.05),
      labels = num_lab
    ) +
    scale_y_discrete(drop = FALSE) +
    scale_colour_manual(values = G_COL, name = NULL) +
    scale_shape_manual(values = G_SHAPE, name = NULL) +
    scale_fill_manual(values = G_COL, name = NULL) +
    labs(x = "coverage (grid average)", y = NULL) +
    theme_fig()
}
S1_CT <- c(
  "CL2 - CR1" = "CL2 minus CR1",
  "exact - shortcut" = "CL2 minus CL2, shortcut",
  "Bayes - freq (CL2)" = "Bayesian minus frequentist, CL2",
  "Bayes - freq (CR1)" = "Bayesian minus frequentist, CR1",
  "t - z (CL2)" = "t(G-1) minus normal",
  "fixed - reselected (CL2)" = "fixed minus reselected λ, CL2",
  "fixed - reselected (model)" = "fixed minus reselected λ, model-based"
)
fig_s1_contrasts <- function(s1_ct) {
  d <- tidyr::expand_grid(
    contrast = names(S1_CT),
    dep = c(FALSE, TRUE),
    G = c(40, 100)
  ) |>
    mutate(
      diff = 100 *
        purrr::pmap_dbl(
          list(contrast, dep, G),
          \(c, de, g) s1_ct(c, de, g)[["diff"]]
        )
    ) |>
    mutate(
      contrast = factor(S1_CT[contrast], rev(S1_CT)),
      dep = factor(
        if_else(dep, "dependent errors", "independent errors"),
        c("independent errors", "dependent errors")
      ),
      G = factor(paste0("G = ", G), paste0("G = ", c(40, 100)))
    )
  ggplot(d, aes(diff, contrast, colour = G, shape = G, fill = G)) +
    geom_vline(xintercept = 0, linetype = 2, colour = REF, linewidth = 0.35) +
    geom_line(
      aes(group = contrast),
      colour = "grey70",
      linewidth = 0.4,
      orientation = "y"
    ) +
    geom_point(size = 1.6, stroke = 0.3) +
    facet_grid(. ~ dep) +
    scale_x_continuous(breaks = seq(-2, 10, 2)) +
    scale_colour_manual(values = G_COL, name = NULL) +
    scale_shape_manual(values = G_SHAPE, name = NULL) +
    scale_fill_manual(values = G_COL, name = NULL) +
    labs(x = "difference in coverage (pp)", y = NULL) +
    theme_fig()
}
# critical values: coverage, one panel per family x dependence and number of curves
CRIT_STYLE <- list(
  col = c("#D55E00", "#999933", "#0072B2", "#CC79A7"),
  shape = c(24L, 22L, 21L, 23L)
)
fig_s1_satt <- function(sat_w, crit_lab) {
  crit_lab[["z"]] <- "normal"
  lv <- unname(crit_lab)
  col <- setNames(CRIT_STYLE$col, lv)
  shp <- setNames(CRIT_STYLE$shape, lv)
  d <- sat_w |>
    mutate(
      crit = factor(crit_lab[crit], lv),
      family = label_family(family) |> droplevels(),
      est = factor(
        ESTIMAND_LABELS[estimand],
        rev(ESTIMAND_LABELS[c("mean", "beta", "alpha", "gamma")])
      ),
      dep = factor(dependence, c("dependent", "independent")),
      G = factor(paste0("G = ", G), paste0("G = ", c(40, 100)))
    )
  ggplot(d, aes(cov, est, colour = crit, shape = crit, fill = crit)) +
    vref(LEVEL) +
    geom_point(position = DODGE(0.55), size = 1.35, stroke = 0.3) +
    facet_grid(dep + G ~ family) +
    scale_x_continuous(breaks = seq(0.85, 1, 0.05), labels = num_lab) +
    scale_colour_manual(values = col, name = NULL) +
    scale_shape_manual(values = shp, name = NULL) +
    scale_fill_manual(values = col, name = NULL) +
    labs(x = "coverage (grid average)", y = NULL) +
    theme_fig()
}

# ---- Appendix C: computing time ------------------------------------------------------------------------
TIME_SRC <- c(
  A = "A: local benchmark",
  B = "B: comparator runs",
  C = "C: main factorial",
  D = "D: plasmode"
)
fig_time <- function(time_tab) {
  d <- time_tab |>
    mutate(
      src = factor(TIME_SRC[data], TIME_SRC),
      lab = sprintf(
        "%s (n = %s)",
        method,
        format(n, big.mark = ",", trim = TRUE)
      ),
      lab = reorder(lab, median)
    )
  ggplot(d, aes(y = lab)) +
    geom_linerange(
      aes(xmin = q25, xmax = q90),
      linewidth = 0.3,
      colour = "grey45"
    ) +
    geom_linerange(
      aes(xmin = q25, xmax = q75),
      linewidth = 1.4,
      colour = "#7FB2D5"
    ) +
    geom_point(aes(x = median), size = 1.5, colour = "#0072B2") +
    geom_point(
      aes(x = max),
      shape = 4,
      size = 1.2,
      stroke = 0.4,
      colour = "grey25"
    ) +
    facet_grid(
      src ~ .,
      scales = "free_y",
      space = "free_y",
      switch = "y",
      labeller = label_wrap_gen(11)
    ) +
    scale_x_log10(
      breaks = 10^(-1:4),
      labels = c("0.1", "1", "10", "100", "1,000", "10,000"),
      minor_breaks = rep(1:9, 6) * rep(10^(-1:4), each = 9),
      guide = guide_axis(minor.ticks = TRUE),
      expand = expansion(mult = c(0.03, 0.08))
    ) +
    labs(x = "seconds (log scale)", y = NULL) +
    theme_fig() +
    theme(
      strip.placement = "outside",
      strip.text.y.left = element_text(angle = 0, hjust = 0),
      panel.grid.major.y = element_line(colour = "grey94")
    )
}

# ---- supplement: restyled figures ----------------------------------------------------------------
ERR_LONG <- c(
  iid = "independent",
  indep. = "independent",
  OU = "rough stationary",
  rough = "rough stationary",
  smooth = "smooth",
  ou = "rough stationary",
  fpc = "smooth"
)
ARM_MAIN <- c(
  "REML, model-based",
  "REML + CL2",
  "NCV + CL2",
  "NCV + CL2, widened by NCV–REML difference"
)
# the supplement's interval labels, in the fixed order and with the common encoding
arm_levels <- function(x) {
  x <- as.character(x)
  c(
    intersect(
      c(
        ARM_MAIN,
        "NCV + CL2 (freq.), widened by NCV–REML difference",
        "AR(1) working model"
      ),
      x
    ),
    setdiff(
      unique(x),
      c(
        ARM_MAIN,
        "NCV + CL2 (freq.), widened by NCV–REML difference",
        "AR(1) working model"
      )
    )
  )
}
as_arm <- function(x) factor(as.character(x), arm_levels(x))
sfig_signal <- function(p) {
  d <- p$data |> mutate(arm = as_arm(arm))
  ggplot(d, aes(coverage_mean, signal, colour = arm, shape = arm, fill = arm)) +
    vref(LEVEL) +
    hpoint(
      coverage_mean - 2 * coverage_se,
      coverage_mean + 2 * coverage_se,
      w = 0.75,
      size = 1.4
    ) +
    facet_grid(truth + error ~ family) +
    scale_x_continuous(breaks = seq(0.8, 1, 0.05), labels = num_lab) +
    scale_arm(levels(d$arm), nrow = 2) +
    labs(x = "coverage (grid average)", y = "signal level") +
    theme_fig()
}
sfig_grid <- function(p) {
  d <- p$data |>
    mutate(
      arm = as_arm(arm),
      error = factor(ERR_LONG[as.character(error)], unique(ERR_LONG))
    )
  pd <- position_dodge(width = 0.18)
  ggplot(
    d,
    aes(D, coverage_mean, colour = arm, shape = arm, fill = arm, group = arm)
  ) +
    geom_hline(
      yintercept = LEVEL,
      linetype = 2,
      colour = REF,
      linewidth = 0.35
    ) +
    geom_line(linewidth = 0.4, position = pd) +
    geom_linerange(
      aes(
        ymin = coverage_mean - 2 * coverage_se,
        ymax = coverage_mean + 2 * coverage_se
      ),
      linewidth = 0.4,
      position = pd
    ) +
    geom_point(size = 1.4, stroke = 0.3, position = pd) +
    facet_grid(family ~ estimand + error) +
    scale_y_continuous(breaks = seq(0.4, 1, 0.2), labels = num_lab) +
    scale_arm(levels(d$arm)) +
    labs(x = "grid points per curve D", y = "coverage (grid average)") +
    theme_fig() +
    theme(strip.text.x = element_text(size = 7))
}
sfig_mse <- function(p) {
  d <- p$data |>
    mutate(
      estimand = factor(
        recode(
          as.character(estimand),
          `E(Y | X)` = "E(Y | X)",
          `gamma(t)` = "γ(t)",
          `alpha(t)` = "α(t)",
          `surface: beta(s,t) or f(x,t)` = "β(s,t) or f(x,t)"
        ),
        c("E(Y | X)", "β(s,t) or f(x,t)", "α(t)", "γ(t)")
      ) |>
        droplevels()
    )
  fam <- c(Gaussian = "#000000", Poisson = "#56B4E9", binary = "#CC79A7")
  ggplot(
    d,
    aes(
      ratio,
      error,
      colour = family,
      shape = signal,
      group = interaction(family, signal)
    )
  ) +
    vref(1) +
    geom_linerange(
      aes(xmin = ratio_lower, xmax = ratio_upper),
      position = DODGE(0.82),
      linewidth = 0.35
    ) +
    geom_point(position = DODGE(0.82), size = 1.2, stroke = 0.4) +
    facet_grid(setting ~ estimand, scales = "free_y", space = "free_y") +
    scale_shape_manual(
      values = c(low = 1, mid = 16, high = 17),
      name = "signal"
    ) +
    scale_x_log10(breaks = c(0.05, 0.1, 0.2, 0.5, 1, 2), labels = num_lab) +
    scale_colour_manual(values = fam, name = NULL) +
    labs(x = "MSE of NCV / MSE of REML (log scale; G = 100)", y = NULL) +
    theme_fig()
}
sfig_bias <- function(p) {
  d <- p$data |>
    select(family, error, estimand, fit, truth, component, value) |>
    pivot_wider(names_from = component, values_from = value) |>
    mutate(
      share = `squared bias` / (`squared bias` + variance),
      fit = factor(fit, c("REML", "NCV")),
      truth = factor(paste(truth, "truth"), c("smooth truth", "rough truth")),
      error = factor(ERR_LONG[as.character(error)], rev(unique(ERR_LONG)))
    )
  ggplot(d, aes(share, error, colour = fit, shape = fit, fill = fit)) +
    geom_point(position = DODGE(0.5), size = 1.6, stroke = 0.3) +
    facet_grid(estimand + truth ~ family, scales = "free_y", space = "free_y") +
    scale_x_continuous(
      limits = c(0, NA),
      breaks = seq(0, 1, 0.1),
      labels = num_lab
    ) +
    scale_colour_manual(values = FIT_COL, name = "fit") +
    scale_shape_manual(values = FIT_SHAPE, name = "fit") +
    scale_fill_manual(values = FIT_COL, name = "fit") +
    labs(
      x = "share of squared bias in the integrated squared error",
      y = NULL
    ) +
    theme_fig()
}
sfig_rel <- function(p) {
  rel_scatter(
    p$data,
    NCV,
    REML,
    FALSE,
    xlab = "relative error of β(s,t), NCV",
    ylab = "relative error of β(s,t), REML"
  ) +
    geom_vline(xintercept = 1, linetype = 3, colour = REF, linewidth = 0.4) +
    facet_wrap(~dependence) +
    guides(
      colour = guide_legend(
        nrow = 1,
        override.aes = list(size = 1.8, fill = unname(FAM_COL))
      ),
      shape = guide_legend(nrow = 1)
    )
}
TRUTH_LAB <- c(REML = "REML-fitted", MID = "midpoint", NCV = "NCV-fitted")
sfig_plasmode_cov <- function(p) {
  d <- p$data |>
    mutate(
      arm = as_arm(arm),
      truth = factor(TRUTH_LAB[as.character(truth)], rev(TRUTH_LAB)),
      residual = factor(
        as.character(residual),
        c("REML residuals", "NCV residuals")
      )
    )
  lv <- levels(d$arm)
  ggplot(d, aes(coverage_mean, truth, colour = arm, shape = arm, fill = arm)) +
    vref(LEVEL) +
    hpoint(
      coverage_mean - 2 * coverage_se,
      coverage_mean + 2 * coverage_se,
      w = 0.8,
      size = 1.45,
      lw = 0.4
    ) +
    facet_grid(app ~ estimand + residual) +
    scale_x_continuous(breaks = seq(0.2, 1, 0.2), labels = num_lab) +
    scale_arm(lv, nrow = 2) +
    labs(x = "coverage (grid average)", y = "truth source") +
    theme_fig() +
    theme(
      strip.text.y = element_text(size = 7),
      strip.text.x = element_text(size = 7)
    )
}
TRUTH_COL <- c(
  `REML-fitted` = "#0072B2",
  midpoint = "#009E73",
  `NCV-fitted` = "#D55E00"
)
TRUTH_SHAPE <- c(`REML-fitted` = 21L, midpoint = 22L, `NCV-fitted` = 24L)
sfig_plasmode_mse <- function(p) {
  d <- p$data |>
    mutate(
      truth = factor(TRUTH_LAB[as.character(truth)], TRUTH_LAB),
      app = factor(as.character(app), rev(levels(app))),
      .fill = if_else(
        grepl("^NCV", as.character(residual)),
        "white",
        unname(TRUTH_COL[as.character(truth)])
      )
    )
  pd <- position_dodge(width = 0.8, orientation = "y")
  ggplot(
    d,
    aes(
      mse_ratio_ncv_reml,
      app,
      colour = truth,
      shape = truth,
      group = interaction(truth, residual)
    )
  ) +
    vref(1) +
    geom_linerange(
      aes(xmin = ratio_lower, xmax = ratio_upper),
      position = pd,
      linewidth = 0.35
    ) +
    geom_point(
      aes(fill = I(.fill)),
      position = pd,
      size = 1.15,
      stroke = 0.35
    ) +
    facet_wrap(~estimand, nrow = 1) +
    scale_x_log10(breaks = c(0.03, 0.1, 0.3, 1, 3), labels = num_lab) +
    scale_colour_manual(values = TRUTH_COL, name = "truth") +
    scale_shape_manual(values = TRUTH_SHAPE, name = "truth") +
    guides(
      colour = guide_legend(
        override.aes = list(fill = unname(TRUTH_COL), size = 1.6)
      )
    ) +
    labs(x = "MSE of NCV / MSE of REML (log scale)", y = NULL) +
    theme_fig()
}
# strip plots: every setting one open symbol, jittered vertically; bar at the median
sfig_strip <- function(p, x, xlab, ref, log = FALSE, breaks, y = "arm") {
  med <- p$data |>
    group_by(.y = .data[[y]], estimand, setting) |>
    summarise(.med = median({{ x }}, na.rm = TRUE), .groups = "drop")
  ggplot(p$data, aes({{ x }}, .data[[y]])) +
    vref(ref) +
    geom_point(
      position = position_jitter(height = 0.22, width = 0, seed = 1),
      shape = 1,
      size = 0.9,
      stroke = 0.3,
      colour = "grey20"
    ) +
    geom_point(
      data = med,
      aes(.med, .y),
      inherit.aes = FALSE,
      shape = "|",
      size = 3.2,
      colour = "#D55E00"
    ) +
    facet_grid(estimand ~ setting) +
    (if (log) scale_x_log10(breaks = breaks, labels = num_lab) else
      scale_x_continuous(breaks = breaks, labels = num_lab)) +
    labs(x = xlab, y = NULL) +
    theme_fig()
}
sfig_misspec <- function(p) {
  d <- p$data |>
    mutate(
      setting = sub("^iid", "independent", as.character(setting)),
      setting = factor(setting, unique(setting)),
      basis = factor(
        c(default = "default", xlarge = "4× finer knots")[as.character(basis)],
        c("default", "4× finer knots")
      )
    )
  col <- c(default = "#D55E00", `4× finer knots` = "#0072B2")
  shp <- c(default = 24L, `4× finer knots` = 21L)
  db <- d |> filter(truth_label == "T1", estimand == ESTIMAND_LABELS[["beta"]])
  pd <- position_dodge(width = 0.6, orientation = "y")
  ggplot(
    d,
    aes(coverage, arm_l, colour = basis, shape = basis, group = basis)
  ) +
    vref(LEVEL) +
    geom_linerange(
      aes(xmin = coverage - 2 * coverage_se, xmax = coverage + 2 * coverage_se),
      position = pd,
      linewidth = 0.4
    ) +
    geom_point(aes(fill = basis), position = pd, size = 1.5, stroke = 0.3) +
    geom_point(
      data = db,
      aes(x = coverage_bump),
      fill = "white",
      position = pd,
      size = 1.5,
      stroke = 0.4,
      show.legend = FALSE
    ) +
    facet_grid(truth + estimand ~ setting) +
    scale_x_continuous(breaks = seq(0.2, 1, 0.2), labels = num_lab) +
    scale_colour_manual(values = col, name = "basis") +
    scale_shape_manual(values = shp, name = "basis") +
    scale_fill_manual(values = col, name = "basis") +
    labs(x = "coverage of the exact truth (grid average)", y = NULL) +
    theme_fig() +
    theme(strip.text.y = element_text(size = 7))
}
sfig_maps <- function(fig) {
  d <- fig[[1]]$data
  p_cov <- ggplot(d, aes(s, t)) +
    geom_raster(aes(fill = coverage)) +
    geom_contour(
      aes(z = coverage),
      breaks = 0.9,
      colour = "black",
      linewidth = 0.45
    ) +
    facet_grid(truth ~ method) +
    scale_fill_viridis_c(
      limits = c(0.7, 1),
      oob = scales::oob_squish,
      name = "coverage",
      breaks = c(0.7, 0.8, 0.9, 0.95, 1),
      labels = c("≤ 0.7", "0.8", "0.9", "0.95", "1")
    ) +
    coord_equal() +
    theme_doc
  (p_cov + fig[[2]] + plot_layout(widths = c(4, 1.15))) &
    scale_x_continuous(
      breaks = c(0, 0.5, 1),
      labels = num_lab,
      expand = c(0, 0)
    ) &
    scale_y_continuous(
      breaks = c(0, 0.5, 1),
      labels = num_lab,
      expand = c(0, 0)
    ) &
    theme(
      strip.text = element_text(size = 6.5),
      strip.background = element_rect(fill = "grey93", colour = NA),
      legend.position = "bottom",
      legend.title.position = "top",
      legend.key.height = unit(0.25, "cm"),
      legend.key.width = unit(0.75, "cm"),
      legend.margin = margin(0, 0, 0, 0),
      axis.text = element_text(size = 6.5),
      panel.spacing = unit(0.3, "lines")
    )
}
BASIS_LAB <- c(default = "default", large = "2× finer", xlarge = "4× finer")
sfig_basis <- function(fig) {
  style <- list(
    scale_x_discrete(labels = BASIS_LAB),
    theme_fig(),
    theme(strip.text = element_text(size = 7))
  )
  p1 <- fig[[1]] +
    aes(shape = fit) +
    style +
    scale_shape_manual(values = c(NCV = 17, REML = 16), name = NULL) +
    labs(y = "EDF of the surface term")
  p2 <- fig[[2]] +
    aes(shape = estimand) +
    style +
    scale_shape_manual(values = c(17, 16), name = NULL) +
    scale_y_log10(breaks = c(0.01, 0.03, 0.1, 0.3, 1), labels = num_lab) +
    labs(y = "MSE of NCV / MSE of REML")
  lv <- levels(droplevels(as_arm(fig[[3]]$data$arm)))
  p3 <- fig[[3]]
  p3$data <- p3$data |> mutate(arm = as_arm(arm))
  p3 <- suppressMessages(
    p3 +
      aes(shape = arm, fill = arm) +
      scale_arm(lv, nrow = 2) +
      style +
      labs(y = "coverage of β(s,t)")
  )
  p1 / p2 / p3
}

# ---- supplement: figures replacing tables ---------------------------------------------------------
sfig_tail <- function(d) {
  d <- d |>
    mutate(
      arm = factor(as.character(arm), rev(arm_levels(arm))),
      dependence = factor(dependence, c("dependent", "independent")),
      Gl = factor(paste0("G = ", G), paste0("G = ", c(40, 100)))
    )
  lv <- rev(levels(d$arm))
  ggplot(d, aes(y = arm, colour = arm, fill = arm)) +
    vref(LEVEL) +
    geom_segment(
      aes(x = q05_min, xend = q05_median, yend = arm),
      linewidth = 0.5
    ) +
    geom_point(aes(x = q05_median, shape = arm), size = 1.7, stroke = 0.3) +
    geom_point(aes(x = q05_min), shape = "|", size = 2.2) +
    geom_point(aes(x = coverage_min), shape = 4, size = 1.3, stroke = 0.45) +
    facet_grid(dependence + Gl ~ estimand) +
    scale_x_continuous(breaks = COV_BREAKS(0.1), labels = num_lab) +
    scale_y_discrete(labels = NULL) +
    scale_arm(lv, nrow = 2) +
    labs(
      x = "lower tail (segment and symbol); minimum grid-averaged coverage (cross)",
      y = NULL
    ) +
    theme_fig() +
    theme(axis.ticks.y = element_blank(), strip.text.y = element_text(size = 7))
}
sfig_ladder <- function(d) {
  d <- d |>
    select(family, errors, estimand, interval, cov_20, cov_40, cov_100) |>
    pivot_longer(
      starts_with("cov_"),
      names_to = "G",
      names_prefix = "cov_",
      values_to = "cov"
    ) |>
    mutate(
      G = factor(G, c("20", "40", "100")),
      interval = as_arm(interval),
      estimand = factor(
        as.character(estimand),
        c("E(Y | X)", "β(s,t)", "α(t)", "γ(t)")
      ),
      family = factor(as.character(family), c("Gaussian", "Poisson", "binary")),
      errors = factor(errors, c("dependent", "independent"))
    )
  ggplot(
    d,
    aes(
      G,
      cov,
      colour = interval,
      shape = interval,
      fill = interval,
      group = interval
    )
  ) +
    geom_hline(
      yintercept = LEVEL,
      linetype = 2,
      colour = REF,
      linewidth = 0.35
    ) +
    geom_line(linewidth = 0.4) +
    geom_point(size = 1.3, stroke = 0.3) +
    facet_grid(estimand ~ family + errors) +
    scale_y_continuous(breaks = seq(0.5, 1, 0.1), labels = num_lab) +
    scale_arm(levels(d$interval), nrow = 2) +
    labs(
      x = "number of curves G",
      y = "coverage (grid average; mean over settings)"
    ) +
    theme_fig() +
    theme(strip.text.x = element_text(size = 6.5))
}
HET_LAB <- c(
  smooth = "smooth",
  `var(z)` = "variance by covariate",
  `var(t)` = "variance along t",
  `var(t,z)` = "both"
)
sfig_het <- function(d) {
  long <- d |>
    select(-cell) |>
    pivot_longer(-c(error, G, estimand), names_to = "arm", values_to = "cov") |>
    filter(!is.na(cov)) |>
    mutate(
      arm = as_arm(arm),
      error = factor(HET_LAB[as.character(error)], rev(HET_LAB)),
      G = factor(paste0("G = ", G), paste0("G = ", c(40, 100))),
      estimand = factor(as.character(estimand), c("E(Y | X)", "β(s,t)", "γ(t)"))
    )
  ggplot(long, aes(cov, error, colour = arm, shape = arm, fill = arm)) +
    vref(LEVEL) +
    geom_point(position = DODGE(0.72), size = 1.4, stroke = 0.3) +
    facet_grid(G ~ estimand) +
    scale_x_continuous(breaks = COV_BREAKS(0.1), labels = num_lab) +
    scale_arm(levels(long$arm), nrow = 2) +
    labs(x = "coverage (grid average)", y = NULL) +
    theme_fig()
}
# "0.97 [0.96, 0.98]" -> estimate, lower, upper
parse_ci <- function(x) {
  m <- regmatches(x, regexec("^(-?[0-9.]+) \\[(-?[0-9.]+), (-?[0-9.]+)\\]", x))
  tibble::tibble(
    est = as.numeric(vapply(m, `[`, "", 2)),
    lo = as.numeric(vapply(m, `[`, "", 3)),
    hi = as.numeric(vapply(m, `[`, "", 4))
  )
}
sfig_prof <- function(d) {
  long <- d |>
    pivot_longer(c(iid, ou, fpc), names_to = "error", values_to = "txt") |>
    bind_cols(parse_ci(
      d |> pivot_longer(c(iid, ou, fpc), values_to = "txt") |> pull(txt)
    )) |>
    mutate(
      error = factor(S1_ERR[error], rev(S1_ERR)),
      family = factor(as.character(family), c("Gaussian", "Poisson", "binary")),
      estimand = factor(as.character(estimand), c("E(Y | X)", "β(s,t)"))
    )
  ggplot(long, aes(est, error)) +
    vref(1) +
    geom_linerange(
      aes(xmin = lo, xmax = hi),
      linewidth = 0.5,
      colour = "#0072B2"
    ) +
    geom_point(size = 1.5, colour = "#0072B2") +
    facet_grid(family ~ estimand) +
    scale_x_log10(breaks = c(0.1, 0.2, 0.3, 0.5, 0.7, 1), labels = num_lab) +
    labs(x = "NCV/REML ratio of mean squared errors (log scale)", y = NULL) +
    theme_fig()
}
# REML against NCV: dumbbells of one measure
sfig_edf <- function(d) {
  w <- d |>
    mutate(
      family = factor(
        FAMILY_LABELS[as.character(family)],
        c("Gaussian", "Poisson", "binary")
      ),
      error = factor(S1_ERR[as.character(error)], rev(S1_ERR)),
      term = factor(
        recode(
          as.character(term),
          `alpha(t)` = "α(t)",
          `beta(s,t)` = "β(s,t)",
          `gamma(t)` = "γ(t)"
        ),
        c("β(s,t)", "α(t)", "γ(t)")
      )
    )
  long <- w |>
    pivot_longer(c(REML, NCV), names_to = "fit", values_to = "x") |>
    mutate(fit = factor(fit, c("REML", "NCV")))
  ggplot(long, aes(x, error)) +
    geom_segment(
      data = w,
      aes(x = REML, xend = NCV, y = error, yend = error),
      colour = "grey60",
      linewidth = 0.5
    ) +
    geom_point(
      aes(colour = fit, shape = fit, fill = fit),
      size = 1.7,
      stroke = 0.3
    ) +
    facet_grid(family ~ term, scales = "free_x") +
    scale_x_continuous(breaks = \(l) base::pretty(l, 4)) +
    scale_colour_manual(values = FIT_COL, name = "fit") +
    scale_shape_manual(values = FIT_SHAPE, name = "fit") +
    scale_fill_manual(values = FIT_COL, name = "fit") +
    labs(x = "mean EDF per term", y = NULL) +
    theme_fig()
}
sfig_lsp <- function(d) {
  w <- d |>
    mutate(
      family = factor(
        FAMILY_LABELS[as.character(family)],
        c("Gaussian", "Poisson", "binary")
      ),
      error = factor(S1_ERR[as.character(error)], rev(S1_ERR)),
      G = factor(paste0("G = ", G), paste0("G = ", c(40, 100)))
    )
  long <- w |>
    pivot_longer(c(REML, NCV), names_to = "fit", values_to = "x") |>
    mutate(fit = factor(fit, c("REML", "NCV")))
  ggplot(long, aes(x, error)) +
    geom_segment(
      data = w,
      aes(x = REML, xend = NCV, y = error, yend = error),
      colour = "grey60",
      linewidth = 0.5
    ) +
    geom_point(
      aes(colour = fit, shape = fit, fill = fit),
      size = 1.7,
      stroke = 0.3
    ) +
    facet_grid(family ~ G) +
    scale_x_continuous(breaks = 0:6) +
    scale_colour_manual(values = FIT_COL, name = "fit") +
    scale_shape_manual(values = FIT_SHAPE, name = "fit") +
    scale_fill_manual(values = FIT_COL, name = "fit") +
    labs(x = "SD over replicates of the log smoothing parameters", y = NULL) +
    theme_fig()
}
# normal -> Satterthwaite dumbbells per row label, one panel per estimand
sfig_satt_rows <- function(d, group = NULL) {
  ests <- c(mean = "E(Y | X)", beta = "β(s,t)", alpha = "α(t)", gamma = "γ(t)")
  long <- d |>
    pivot_longer(
      matches("^(mean|beta|alpha|gamma) (z|satt)$"),
      names_to = c("est", ".value"),
      names_sep = " "
    ) |>
    filter(!is.na(z) | !is.na(satt)) |>
    mutate(est = factor(ests[est], ests))
  col <- FIT_COL[["REML"]]
  p <- ggplot(long, aes(y = row)) +
    vref(LEVEL) +
    geom_segment(
      aes(x = z, xend = satt, yend = row),
      colour = col,
      linewidth = 0.5
    ) +
    geom_point(
      aes(z, shape = "normal"),
      colour = col,
      fill = "white",
      size = 1.5,
      stroke = 0.45
    ) +
    geom_point(
      aes(satt, shape = "Satterthwaite"),
      colour = col,
      fill = col,
      size = 1.5,
      stroke = 0.45
    ) +
    scale_x_continuous(breaks = seq(0.85, 1, 0.05), labels = num_lab) +
    scale_shape_manual(
      values = c(normal = 21, Satterthwaite = 21),
      name = "critical value"
    ) +
    guides(
      shape = guide_legend(
        override.aes = list(fill = c("white", col), colour = col, size = 1.7)
      )
    ) +
    labs(x = "coverage (grid average)", y = NULL) +
    theme_fig()
  if (is.null(group)) p + facet_grid(. ~ est) else
    p +
      facet_grid(
        rlang::new_formula(rlang::sym(group), quote(est)),
        scales = "free_y",
        space = "free_y"
      )
}
sfig_satx <- function(d) {
  d <- d |>
    mutate(
      row = paste0(block, ", G = ", G),
      row = factor(row, rev(unique(row)))
    )
  sfig_satt_rows(d)
}
sfig_satp <- function(d) {
  d <- d |>
    mutate(
      dataset = factor(as.character(dataset), unique(as.character(dataset))),
      group = factor(
        if_else(
          dataset %in% APP_LABELS[APP_STRESS],
          "dep. source",
          "indep. source"
        ),
        c("indep. source", "dep. source")
      ),
      row = factor(paste0(dataset, ", G = ", G))
    ) |>
    arrange(dataset, desc(G)) |>
    mutate(row = factor(row, rev(unique(as.character(row)))))
  sfig_satt_rows(d, "group")
}
ROLE_LAB <- c(
  counted = "indep. source",
  `stress test` = "dep. source",
  `indep. curves` = "indep. source",
  `dep. source` = "dep. source",
  `indep. source` = "indep. source"
)
sfig_tails_p <- function(d) {
  long <- d |>
    pivot_longer(
      starts_with("q05 "),
      names_to = c("stat", "arm"),
      names_sep = ": ",
      values_to = "v"
    ) |>
    pivot_wider(names_from = stat, values_from = v) |>
    mutate(
      arm = as_arm(ARM_LABELS[arm]),
      estimand = factor(
        ESTIMAND_LABELS[as.character(estimand)],
        ESTIMAND_LABELS[c("mean", "beta")]
      ),
      app = factor(as.character(app), rev(APP_LABELS)),
      role = factor(
        ROLE_LAB[as.character(role)],
        c("indep. source", "dep. source")
      )
    )
  ggplot(long, aes(`q05 mean`, app, colour = arm, shape = arm, fill = arm)) +
    geom_linerange(
      aes(xmin = `q05 min`, xmax = `q05 mean`),
      position = DODGE(0.75),
      linewidth = 0.4
    ) +
    geom_point(position = DODGE(0.75), size = 1.4, stroke = 0.3) +
    facet_grid(role ~ estimand, scales = "free_y", space = "free_y") +
    scale_x_continuous(breaks = COV_BREAKS(0.1), labels = num_lab) +
    scale_arm(levels(long$arm), nrow = 2) +
    labs(
      x = "lower tail: mean over combinations (symbol), down to the minimum (bar)",
      y = NULL
    ) +
    theme_fig()
}
sfig_g40 <- function(d) {
  long <- d |>
    pivot_longer(
      -c(app, role, estimand),
      names_to = "arm",
      values_to = "loss"
    ) |>
    mutate(
      arm = as_arm(arm),
      estimand = factor(
        ESTIMAND_LABELS[as.character(estimand)],
        ESTIMAND_LABELS[c("mean", "beta")]
      ),
      app = factor(as.character(app), rev(APP_LABELS)),
      role = factor(
        ROLE_LAB[as.character(role)],
        c("indep. source", "dep. source")
      )
    )
  ggplot(long, aes(100 * loss, app, colour = arm, shape = arm, fill = arm)) +
    geom_vline(xintercept = 0, linetype = 2, colour = REF, linewidth = 0.35) +
    geom_point(position = DODGE(0.75), size = 1.4, stroke = 0.3) +
    facet_grid(role ~ estimand, scales = "free_y", space = "free_y") +
    scale_x_continuous(breaks = seq(-5, 40, 5)) +
    scale_arm(levels(long$arm), nrow = 2) +
    labs(
      x = "coverage with all subjects minus coverage with 40 subjects (pp)",
      y = NULL
    ) +
    theme_fig()
}
sfig_rbc <- function(d) {
  ms <- c(
    cov = "coverage",
    q05 = "lower tail",
    wr = "width relative to REML + CL2",
    det = "detection",
    rel = "relative error"
  )
  arms <- c(
    ncv_cl2 = "NCV + CL2",
    rbc_p = "bias correction",
    reml_cl2 = "REML + CL2"
  )
  col <- c(
    `NCV + CL2` = "#D55E00",
    `bias correction` = "#CC79A7",
    `REML + CL2` = "#0072B2"
  )
  shp <- c(`NCV + CL2` = 24L, `bias correction` = 23L, `REML + CL2` = 21L)
  long <- d |>
    mutate(
      arm = factor(arms[as.character(arm)], arms),
      group = factor(
        paste0(group, " (", n, ")"),
        rev(unique(paste0(group, " (", n, ")")))
      )
    ) |>
    pivot_longer(all_of(names(ms)), names_to = "m", values_to = "x") |>
    mutate(measure = factor(ms[m], ms))
  refs <- tibble::tibble(measure = factor("coverage", ms), x = LEVEL)
  ggplot(long, aes(x, group, colour = arm, shape = arm, fill = arm)) +
    geom_vline(
      data = refs,
      aes(xintercept = x),
      linetype = 2,
      colour = REF,
      linewidth = 0.35
    ) +
    geom_point(position = DODGE(0.6), size = 1.5, stroke = 0.3) +
    facet_grid(. ~ measure, scales = "free_x", labeller = label_wrap_gen(14)) +
    scale_x_continuous(breaks = \(l) base::pretty(l, 3), labels = num_lab) +
    scale_colour_manual(values = col, name = NULL) +
    scale_shape_manual(values = shp, name = NULL) +
    scale_fill_manual(values = col, name = NULL) +
    labs(x = NULL, y = NULL) +
    theme_fig()
}
# tables whose cells read "coverage (lower tail) third"
parse_cell <- function(x) {
  m <- regmatches(
    x,
    regexec("^([0-9.]+) \\(([0-9.]+)\\)(?: ([0-9.]+))?", x, perl = TRUE)
  )
  tibble::tibble(
    coverage = as.numeric(vapply(
      m,
      \(v) if (length(v)) v[2] else NA_character_,
      ""
    )),
    tail = as.numeric(vapply(
      m,
      \(v) if (length(v)) v[3] else NA_character_,
      ""
    )),
    third = as.numeric(vapply(
      m,
      \(v) if (length(v)) v[4] else NA_character_,
      ""
    ))
  )
}
sfig_cells <- function(d, third, id = "dataset", refline = NULL) {
  cols <- setdiff(names(d), id)
  stopifnot(length(cols) <= 5L)
  pal <- c("#D55E00", "#E69F00", "#CC79A7", "#009E73", "#0072B2")[seq_along(
    cols
  )]
  shp <- c(24L, 25L, 23L, 22L, 21L)[seq_along(cols)]
  names(pal) <- names(shp) <- cols
  long0 <- d |> pivot_longer(all_of(cols), names_to = "int", values_to = "txt")
  long <- bind_cols(long0, parse_cell(long0$txt)) |>
    mutate(
      int = factor(int, cols),
      app = factor(.data[[id]], rev(unique(.data[[id]]))),
      group = factor(
        if_else(
          .data[[id]] %in% APP_LABELS[APP_STRESS],
          "dep. source",
          "indep. source"
        ),
        c("indep. source", "dep. source")
      )
    ) |>
    pivot_longer(c(coverage, tail, third), names_to = "m", values_to = "x") |>
    filter(!is.na(x)) |>
    mutate(
      measure = factor(
        c(coverage = "coverage", tail = "lower tail", third = third)[m],
        c("coverage", "lower tail", third)
      )
    )
  refs <- tibble::tibble(
    measure = factor(c("coverage", third), levels(long$measure)),
    x = c(LEVEL, 1)
  )
  ggplot(long, aes(x, app, colour = int, shape = int, fill = int)) +
    geom_vline(
      data = refs,
      aes(xintercept = x),
      linetype = 2,
      colour = REF,
      linewidth = 0.35
    ) +
    geom_point(position = DODGE(0.72), size = 1.35, stroke = 0.3) +
    facet_grid(group ~ measure, scales = "free", space = "free_y") +
    scale_x_continuous(breaks = \(l) base::pretty(l, 4), labels = num_lab) +
    scale_colour_manual(values = pal, name = NULL) +
    scale_shape_manual(values = shp, name = NULL) +
    scale_fill_manual(values = pal, name = NULL) +
    guides(
      colour = guide_legend(nrow = 2),
      shape = guide_legend(nrow = 2),
      fill = guide_legend(nrow = 2)
    ) +
    labs(x = NULL, y = NULL) +
    theme_fig()
}
sfig_fgee_p <- function(s5_P, s5_PP, s5_apps) {
  ms <- c(
    reml_model = "REML, model-based",
    reml_cl2 = "REML + CL2",
    fgee_independent = "fGEE, independence",
    fgee_ar1 = "fGEE, AR(1)",
    fgee_fpca = "fGEE, FPCA"
  )
  col <- c("#666666", "#0072B2", "#E69F00", "#009E73", "#CC79A7")
  shp <- c(22L, 21L, 24L, 25L, 23L)
  names(col) <- names(shp) <- ms
  key <- s5_P |>
    distinct(app, G) |>
    mutate(app = factor(app, names(s5_apps))) |>
    arrange(app, desc(G)) |>
    mutate(
      row = paste0(s5_apps[as.character(app)], ", G = ", G),
      row = factor(row, rev(row)),
      app = as.character(app)
    )
  d1 <- s5_P |>
    filter(estimand == "beta", method %in% names(ms)) |>
    left_join(key, by = c("app", "G")) |>
    mutate(method = factor(ms[method], ms))
  p1 <- ggplot(
    d1,
    aes(coverage, row, colour = method, shape = method, fill = method)
  ) +
    vref(LEVEL) +
    geom_point(position = DODGE(0.7), size = 1.35, stroke = 0.3) +
    scale_x_continuous(breaks = COV_BREAKS(0.1), labels = num_lab) +
    scale_colour_manual(values = col, name = NULL) +
    scale_shape_manual(values = shp, name = NULL) +
    scale_fill_manual(values = col, name = NULL) +
    guides(
      colour = guide_legend(ncol = 2),
      shape = guide_legend(ncol = 2),
      fill = guide_legend(ncol = 2)
    ) +
    labs(x = "coverage (grid average)", y = NULL) +
    theme_fig()
  ct <- c(fgee_ind_vs_reml = "fGEE / REML", fgee_ind_vs_ncv = "fGEE / NCV")
  d2 <- s5_PP |>
    filter(estimand == "beta", contrast %in% names(ct)) |>
    left_join(key, by = c("app", "G")) |>
    mutate(contrast = factor(ct[contrast], ct))
  p2 <- ggplot(d2, aes(est, row)) +
    vref(1) +
    geom_linerange(
      aes(xmin = lower, xmax = upper),
      linewidth = 0.45,
      colour = "#E69F00"
    ) +
    geom_point(size = 1.3, colour = "#E69F00", shape = 24, fill = "#E69F00") +
    facet_grid(. ~ contrast) +
    scale_x_log10(breaks = c(0.2, 0.5, 1, 2, 5), labels = num_lab) +
    labs(x = "ratio of mean squared errors (log scale)", y = NULL) +
    theme_fig() +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  p1 + p2 + plot_layout(widths = c(1.25, 1))
}

# ---- S7: values of the paper's Appendix A.6 figures ------------------------------------------------
values_s1 <- function() {
  cells <- readr::read_csv(
    "ablation/summaries/ablation-cells.csv",
    show_col_types = FALSE
  )
  contr <- readr::read_csv(
    "ablation/summaries/ablation-contrasts.csv",
    show_col_types = FALSE
  )
  arms <- tibble::tribble(
    ~arm,
    ~crit,
    ~label,
    "model",
    "z",
    "model",
    "cr1_b",
    "z",
    "CR1",
    "cl2s_b",
    "z",
    "CL2, shortcut",
    "cl2e_b",
    "z",
    "CL2",
    "cl2e_f",
    "z",
    "CL2, freq.",
    "cl2e_b",
    "t",
    "CL2, t(G-1)",
    "fix_cl2e_b",
    "z",
    "CL2, fixed λ"
  )
  wide <- function(e) {
    cells |>
      filter(estimand == e) |>
      inner_join(arms, by = c("arm", "crit")) |>
      mutate(
        label = factor(label, arms$label),
        family = label_family(family),
        error = S1_ERR[error]
      ) |>
      select(family, error, G, label, coverage) |>
      pivot_wider(
        names_from = label,
        values_from = coverage,
        names_sort = TRUE
      ) |>
      arrange(family, factor(error, S1_ERR), G) |>
      mutate(family = as.character(family))
  }
  ag <- arms |>
    filter(arm %in% c("model", "cr1_b", "cl2e_b", "fix_cl2e_b"), crit == "z")
  agw <- cells |>
    filter(estimand %in% c("alpha", "gamma")) |>
    inner_join(ag, by = c("arm", "crit")) |>
    mutate(
      key = paste(c(alpha = "α(t)", gamma = "γ(t)")[estimand], label),
      key = factor(key, paste(rep(c("α(t)", "γ(t)"), each = 4), ag$label)),
      family = label_family(family),
      error = S1_ERR[error]
    ) |>
    select(family, error, G, key, coverage) |>
    pivot_wider(names_from = key, values_from = coverage, names_sort = TRUE) |>
    arrange(family, factor(error, S1_ERR), G) |>
    mutate(family = as.character(family))
  ct1 <- function(ct, dep, g) {
    d <- contr |>
      filter(
        contrast == ct,
        estimand %in% c("mean", "beta"),
        (error != "iid") == dep,
        G == g
      )
    sprintf(
      "%+.1f (%.2f–%.2f)",
      100 * mean(d$diff),
      100 * min(d$diff_se),
      100 * max(d$diff_se)
    )
  }
  cts <- tidyr::expand_grid(
    contrast = names(S1_CT),
    dep = c(FALSE, TRUE),
    G = c(40, 100)
  ) |>
    mutate(
      cell = purrr::pmap_chr(list(contrast, dep, G), ct1),
      col = paste0(if_else(dep, "dep., G = ", "ind., G = "), G)
    ) |>
    select(contrast, col, cell) |>
    pivot_wider(names_from = col, values_from = cell) |>
    mutate(contrast = unname(S1_CT[contrast]))
  list(
    beta = wide("beta"),
    mean = wide("mean"),
    alpha_gamma = agw,
    contrasts = cts
  )
}

# ---- supplement: misregistration ---------------------------------------------------------------------
sfig_warp <- function(d) {
  long <- d |>
    pivot_longer(
      -c(family, G, signal, truth, estimand),
      names_to = "arm",
      values_to = "cov"
    ) |>
    filter(!is.na(cov)) |>
    mutate(
      arm = as_arm(arm),
      truth = factor(
        c(smooth = "smooth truth", wiggly = "rough truth")[truth],
        c("rough truth", "smooth truth")
      ),
      G = factor(paste0("G = ", G), paste0("G = ", c(40, 100))),
      family = factor(as.character(family), c("Gaussian", "Poisson", "binary")),
      estimand = factor(as.character(estimand), c("E(Y | X)", "β(s,t)"))
    )
  ggplot(long, aes(cov, truth, colour = arm, shape = arm, fill = arm)) +
    vref(LEVEL) +
    geom_point(position = DODGE(0.78), size = 1.5, stroke = 0.3) +
    facet_grid(family + estimand ~ G) +
    scale_x_continuous(breaks = seq(0.5, 1, 0.1), labels = num_lab) +
    scale_arm(levels(long$arm), nrow = 3) +
    labs(x = "coverage (grid average)", y = NULL) +
    theme_fig()
}

# ---- S7: tables at their stated size, without scaling to the text width -------------------------------
ARM_SHORT <- c(
  "REML, model-based" = "REML model",
  "model-based (REML)" = "REML model",
  "NCV + CL2, widened by NCV–REML difference" = "NCV widened",
  "NCV + CL2, NCV–REML widening" = "NCV widened",
  "NCV + CL2 (freq.), widened by NCV–REML difference" = "NCV widened (freq.)",
  "AR(1) working model" = "AR(1)"
)
short_arm <- function(x) {
  x <- as.character(x)
  i <- x %in% names(ARM_SHORT)
  x[i] <- ARM_SHORT[x[i]]
  x
}
# blank repeats in leading grouping columns; a column is blanked only within the group to its left
blank_rep <- function(d, cols) {
  key <- rep("", nrow(d))
  d <- as.data.frame(d, check.names = FALSE)
  for (cl in cols) {
    key <- paste(key, d[[cl]], sep = "\r")
    v <- as.character(d[[cl]])
    v[duplicated(key)] <- ""
    d[[cl]] <- v
  }
  d
}
vtab <- function(
  d,
  caption,
  long = FALSE,
  font = 8,
  digits = 3,
  col.names = NA,
  header = NULL,
  group = NULL
) {
  d <- as.data.frame(d, check.names = FALSE)
  if (!is.null(group)) d <- blank_rep(d, group)
  caption <- gsub("(?<!\\\\)%", "\\\\%", caption, perl = TRUE)
  k <- knitr::kable(
    d,
    format = "latex",
    booktabs = TRUE,
    longtable = long,
    digits = digits,
    caption = caption,
    linesep = "",
    col.names = if (identical(col.names, NA)) names(d) else col.names
  )
  k <- kableExtra::kable_styling(
    k,
    font_size = font,
    latex_options = if (long) "repeat_header" else "HOLD_position"
  )
  if (!is.null(header)) k <- kableExtra::add_header_above(k, header)
  k
}

# Appendix B.1: correlation kernels of the error processes and sample latent curves, from
# the study's own kernel code, on the default grid of D = 61 points.
fig_kernels <- function(n_curves = 2, seed = 1) {
  dgp <- new.env()
  sys.source("R/constants.R", envir = dgp)
  sys.source("R/dgp.R", envir = dgp)
  grid <- dgp$t_grid(61)
  processes <- c(
    "independent" = "iid",
    "rough" = "ou",
    "smooth" = "fpc",
    "AR(1)" = "ou0",
    "sign-changing" = "osc"
  )
  kernels <- lapply(processes, dgp$correlation_kernel, t = grid)
  to_factor <- \(x) factor(x, names(processes))
  surface <- purrr::imap_dfr(kernels, \(R, process) {
    tibble::tibble(
      process = to_factor(process),
      t = rep(grid, times = length(grid)),
      t2 = rep(grid, each = length(grid)),
      R = as.vector(R)
    )
  })
  noise <- withr::with_seed(
    seed,
    matrix(stats::rnorm(n_curves * length(grid)), n_curves, length(grid))
  )
  sample_curves <- purrr::imap_dfr(kernels, \(R, process) {
    Z <- noise %*% chol(R)
    tibble::tibble(
      process = to_factor(process),
      curve = factor(rep(seq_len(n_curves), times = length(grid))),
      t = rep(grid, each = n_curves),
      Z = as.vector(Z)
    )
  })
  unit_axis <- list(
    scale_x_continuous(breaks = c(0, 1), expand = c(0, 0)),
    facet_wrap(vars(process), nrow = 1)
  )
  p_kernel <- ggplot(surface, aes(t, t2, fill = R)) +
    geom_raster() +
    unit_axis +
    scale_y_continuous(breaks = c(0, 1), expand = c(0, 0)) +
    scale_fill_gradient2(
      name = "R(t, t')",
      low = "#2166ac",
      mid = "grey96",
      high = "#b2182b",
      midpoint = 0,
      limits = c(-1, 1),
      breaks = c(-1, 0, 1)
    ) +
    labs(x = NULL, y = "t'") +
    theme_doc +
    theme(legend.position = "right", legend.key.height = unit(16, "pt"))
  p_curves <- ggplot(sample_curves, aes(t, Z, colour = curve)) +
    geom_line(linewidth = 0.3) +
    unit_axis +
    scale_colour_manual(values = c("grey20", "#d95f02"), guide = "none") +
    labs(x = "t", y = "sample curves") +
    theme_doc +
    theme(strip.background = element_blank(), strip.text = element_blank())
  p_kernel / p_curves + plot_layout(guides = "collect", heights = c(1.15, 1))
}
