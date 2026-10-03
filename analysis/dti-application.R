# The DTI application: rcst ~ ff(cca), baseline visit, status quo (REML + model-based
# intervals) vs REML + CL2 and NCV + CL2 (+ bias-aware). Writes tables and figures to
# summaries/dti/. Displays are pointwise; no region-wide significance statements.
# Run from the repository root: Rscript analysis/dti-application.R
source("scripts/setup.R")
suppressPackageStartupMessages({
  library(refund)
  library(ggplot2)
  library(patchwork)
})
out_dir <- "summaries/dti"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Data: baseline visit, complete curves (public: refund::DTI).
data("DTI", package = "refund")
base <- DTI[
  DTI$visit == 1 & complete.cases(DTI$cca) & complete.cases(DTI$rcst),
]
dti <- data.frame(Y = I(unclass(base$rcst)), X = I(unclass(base$cca)))
sim <- list(
  data = dti,
  yind = seq(0, 1, length.out = ncol(base$rcst)),
  xind = seq(0, 1, length.out = ncol(base$cca))
)

fit_dti <- function(method) {
  ff <- refund::ff
  xind <- sim$xind
  refund::pffr(
    Y ~ ff(X, xind = xind, splinepars = BASIS_FF),
    data = sim$data,
    yind = sim$yind,
    bs.yindex = BASIS_T,
    bs.int = BASIS_T,
    sandwich = "none",
    method = method
  )
}
fits <- list(REML = fit_dti("REML"), NCV = fit_dti("NCV"))

# beta(s, t) on the evaluation grid, with model-based (Vp for NCV, see R/fit.R) and CL2 SEs.
beta_grid <- expand.grid(s = EVAL_S, t = EVAL_T)
beta_rows <- lapply(names(fits), \(m) {
  fit <- fits[[m]]
  sm <- fit$smooth[[which(vapply(fit$smooth, \(x) "X.smat" %in% x$term, TRUE))]]
  grid <- setNames(beta_grid, sm$term)
  grid[[sm$by]] <- 1
  L <- matrix(0, nrow(grid), length(fit$coefficients))
  L[, sm$first.para:sm$last.para] <- mgcv::PredictMat(sm, grid)
  V_model <- if (m == "NCV") fit$Vp else (fit$Vc %||% fit$Vp)
  V_cl2 <- refund:::pffr_vcov(fit, sandwich = "cl2", cl2_adjustment = "exact")
  se <- \(V) sqrt(rowSums((L %*% V) * L))
  est <- drop(L %*% fit$coefficients) # computed first: tibble() would see its `fit` column
  tibble::tibble(
    beta_grid,
    fit = m,
    est = est,
    se_model = se(V_model),
    se_cl2 = se(V_cl2)
  )
})
beta <- dplyr::bind_rows(beta_rows)
delta <- beta$est[beta$fit == "NCV"] - beta$est[beta$fit == "REML"]
beta$se_bias_aware <- ifelse(
  beta$fit == "NCV",
  sqrt(beta$se_cl2^2 + rep(delta^2, 2)),
  NA
)

z <- qnorm(1 - (1 - LEVEL) / 2)
recipes <- tibble::tribble(
  ~recipe,
  ~fit,
  ~se,
  "REML + model-based (status quo)",
  "REML",
  "se_model",
  "REML + CL2",
  "REML",
  "se_cl2",
  "NCV + CL2",
  "NCV",
  "se_cl2",
  "NCV + CL2 + bias-aware",
  "NCV",
  "se_bias_aware"
)
summary_rows <- lapply(seq_len(nrow(recipes)), \(i) {
  r <- recipes[i, ]
  b <- beta[beta$fit == r$fit, ]
  se <- b[[r$se]]
  tibble::tibble(
    recipe = r$recipe,
    median_halfwidth = median(z * se),
    share_excluding_zero = mean(abs(b$est) > z * se)
  )
})
edf <- vapply(fits, \(f) sum(f$edf), 1)
recipe_table <- dplyr::bind_rows(summary_rows) |>
  dplyr::mutate(total_edf = unname(edf[recipes$fit]))
readr::write_csv(recipe_table, file.path(out_dir, "recipes.csv"))
readr::write_csv(beta, file.path(out_dir, "beta-grid.csv"))
print(as.data.frame(recipe_table), digits = 3)

# Figures: estimates, and pointwise zero-exclusion maps per recipe.
surface <- function(d, fill, title, limits = NULL) {
  ggplot(d, aes(s, t, fill = .data[[fill]])) +
    geom_raster() +
    scale_fill_gradient2(
      low = "#2166ac",
      mid = "grey96",
      high = "#b2182b",
      midpoint = 0,
      limits = limits
    ) +
    coord_fixed() +
    labs(title = title, x = "s (cca)", y = "t (rcst)", fill = NULL) +
    theme_minimal(base_size = 10)
}
lim <- max(abs(beta$est)) * c(-1, 1)
p_est <- surface(beta[beta$fit == "REML", ], "est", "beta(s,t), REML", lim) +
  surface(beta[beta$fit == "NCV", ], "est", "beta(s,t), NCV", lim)
exclusion <- dplyr::bind_rows(lapply(seq_len(nrow(recipes)), \(i) {
  r <- recipes[i, ]
  b <- beta[beta$fit == r$fit, ]
  tibble::tibble(
    s = b$s,
    t = b$t,
    recipe = r$recipe,
    sign = ifelse(abs(b$est) > z * b[[r$se]], sign(b$est), 0)
  )
}))
p_excl <- ggplot(
  exclusion,
  aes(s, t, fill = factor(sign, levels = c(-1, 0, 1)))
) +
  geom_raster() +
  scale_fill_manual(
    values = c(`-1` = "#2166ac", `0` = "grey92", `1` = "#b2182b"),
    labels = c("interval below 0", "covers 0", "interval above 0"),
    name = NULL,
    drop = FALSE
  ) +
  facet_wrap(~recipe, nrow = 1) +
  coord_fixed() +
  labs(
    title = "Pointwise 95% intervals excluding zero",
    x = "s (cca)",
    y = "t (rcst)"
  ) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom")
ggsave(
  file.path(out_dir, "dti-estimates.png"),
  p_est,
  width = 8,
  height = 4,
  dpi = 150
)
ggsave(
  file.path(out_dir, "dti-exclusion.png"),
  p_excl,
  width = 11,
  height = 4,
  dpi = 150
)
cat("written to", out_dir, "\n")
