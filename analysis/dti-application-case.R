# The DTI application with the plasmode study's DTI model: right corticospinal tract FA on
# corpus-callosum FA and case status (MS = 1), baseline visit, the same data, grids and bases
# as the DTI plasmode dataset. Status quo (REML + model-based) vs REML + CL2; the NCV estimate
# as a sensitivity. Writes summaries/dti-case/. Displays are pointwise.
# Run from the repository root: Rscript analysis/dti-application-case.R
source("scripts/setup.R")
source("plasmode/sources.R")
out_dir <- "summaries/dti-case"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

d <- plasmode_data("dti")
sim <- list(
  data = data.frame(Y = I(d$Y), X = I(d$X), z = d$z),
  yind = d$yind,
  xind = d$xind
)
cell <- list(model = "M1", fit_family = "gaussian")
fits <- list(
  REML = fit_pffr(sim, cell, "REML"),
  NCV = fit_pffr(sim, cell, "NCV")
)

term_rows <- function(fit, pick, grid) {
  sm <- fit$smooth[[which(vapply(fit$smooth, pick, TRUE))]]
  grid <- setNames(grid, sm$term)
  if (nzchar(sm$by) && sm$by != "NA") grid[[sm$by]] <- 1
  L <- matrix(0, nrow(grid), length(fit$coefficients))
  L[, sm$first.para:sm$last.para] <- mgcv::PredictMat(sm, grid)
  L
}
is_ff <- \(sm) "X.smat" %in% sm$term
is_gamma <- \(sm) identical(sm$by, "z")

beta_grid <- expand.grid(s = EVAL_S, t = EVAL_T)
gamma_grid <- data.frame(t = d$yind)
rows <- lapply(names(fits), \(m) {
  fit <- fits[[m]]
  V_model <- if (m == "NCV") fit$Vp else (fit$Vc %||% fit$Vp)
  V_cl2 <- refund:::pffr_vcov(fit, sandwich = "cl2", cl2_adjustment = "exact")
  se <- \(L, V) sqrt(rowSums((L %*% V) * L))
  Lb <- term_rows(fit, is_ff, beta_grid)
  Lg <- term_rows(fit, is_gamma, gamma_grid)
  # estimates computed first: tibble() would see its own `fit` column
  est_b <- drop(Lb %*% fit$coefficients)
  est_g <- drop(Lg %*% fit$coefficients)
  list(
    beta = tibble::tibble(
      beta_grid,
      fit = m,
      est = est_b,
      se_model = se(Lb, V_model),
      se_cl2 = se(Lb, V_cl2)
    ),
    gamma = tibble::tibble(
      gamma_grid,
      fit = m,
      est = est_g,
      se_model = se(Lg, V_model),
      se_cl2 = se(Lg, V_cl2)
    )
  )
})
beta <- dplyr::bind_rows(lapply(rows, `[[`, "beta"))
gamma <- dplyr::bind_rows(lapply(rows, `[[`, "gamma"))

z <- qnorm(1 - (1 - LEVEL) / 2)
summ <- \(b, fit_sel, se_col, estimand) {
  b <- b[b$fit == fit_sel, ]
  tibble::tibble(
    estimand = estimand,
    fit = fit_sel,
    se = se_col,
    median_halfwidth = median(z * b[[se_col]]),
    share_excluding_zero = mean(abs(b$est) > z * b[[se_col]]),
    median_se_ratio_cl2_model = median(b$se_cl2 / b$se_model)
  )
}
recipes <- dplyr::bind_rows(
  summ(beta, "REML", "se_model", "beta"),
  summ(beta, "REML", "se_cl2", "beta"),
  summ(gamma, "REML", "se_model", "gamma"),
  summ(gamma, "REML", "se_cl2", "gamma")
)
edf_term <- \(fit, pick) {
  i <- which(vapply(fit$smooth, pick, TRUE))
  sm <- fit$smooth[[i]]
  sum(fit$edf[sm$first.para:sm$last.para])
}
fit_info <- tibble::tibble(
  fit = names(fits),
  G = nrow(d$Y),
  D = length(d$yind),
  total_edf = vapply(fits, \(f) sum(f$edf), 1),
  ff_edf = vapply(fits, edf_term, 1, pick = is_ff),
  gamma_edf = vapply(fits, edf_term, 1, pick = is_gamma),
  r2 = vapply(fits, \(f) summary(f)$r.sq, 1)
)
# Relative difference of the NCV and REML surface estimates, in units of REML's CL2 SE.
rb <- beta[beta$fit == "REML", ]
nb <- beta[beta$fit == "NCV", ]
fit_info$median_abs_diff_in_se <- c(
  NA,
  median(abs(nb$est - rb$est) / rb$se_cl2)
)

readr::write_csv(beta, file.path(out_dir, "beta-grid.csv"))
readr::write_csv(gamma, file.path(out_dir, "gamma.csv"))
readr::write_csv(recipes, file.path(out_dir, "recipes.csv"))
readr::write_csv(fit_info, file.path(out_dir, "fits.csv"))
print(as.data.frame(recipes), digits = 3)
print(as.data.frame(fit_info), digits = 3)
