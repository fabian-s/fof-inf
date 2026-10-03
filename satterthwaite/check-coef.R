# Equivalence check: the df that run-satt.R computes on the study's contrasts equals the
# df column of coef.pffr(..., ci = "pointwise", crit = "satterthwaite", df_gram = ...)
# for the same points. One replicate of one core cell, fit as in the study.
#   Rscript satterthwaite/check-coef.R [cell] [rep]

source("scripts/setup.R")
source("ablation/ablation-lib.R")
source("satterthwaite/satt-lib.R")

args <- commandArgs(trailingOnly = TRUE)
cell_id <- if (length(args) >= 1) as.integer(args[1]) else 7L
rep_id <- if (length(args) >= 2) as.integer(args[2]) else 1L
cell <- ablation_cells()
cell <- cell[cell$cell == cell_id, ]
sim <- simulate(cell, rep_id)
fit <- fit_model(sim, cell, "REML")
V_f <- refund:::pffr_vcov(
  fit,
  sandwich = "cl2",
  freq = TRUE,
  cl2_adjustment = "exact"
)
dc <- df_context(fit, V_f)

for (gram in c("full", "diagonal")) {
  cf <- coef(
    fit,
    sandwich = "cl2",
    cl2_adjustment = "exact",
    ci = "pointwise",
    crit = "satterthwaite",
    df_gram = gram,
    seWithMean = FALSE,
    n1 = 31,
    n2 = 26
  )
  for (estimand in c("beta", "gamma")) {
    sm <- estimand_smooth(fit, estimand)
    term <- cf$smterms[[which(vapply(fit$smooth, \(s) identical(s, sm), TRUE))]]
    grid <- as.data.frame(term$coef[, sm$term, drop = FALSE])
    if (sm$by != "NA") grid[[sm$by]] <- 1
    L <- matrix(0, nrow(grid), length(fit$coefficients))
    L[, sm$first.para:sm$last.para] <- mgcv::PredictMat(sm, grid)
    ours <- refund:::pffr_influence_df(dc$ctx$core, L, df_gram = gram)$df
    est_diff <- max(abs(drop(L %*% fit$coefficients) - term$coef[, "value"]))
    cat(sprintf(
      "%-8s %-5s n = %4d  max|df_ours - df_coef| = %.2e  max|est diff| = %.2e  df range %.1f-%.1f\n",
      gram,
      estimand,
      nrow(L),
      max(abs(ours - term$coef[, "df"])),
      est_diff,
      min(ours),
      max(ours)
    ))
  }
}
