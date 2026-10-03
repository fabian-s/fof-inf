# Re-check of mgcv's neighbourhood-aware covariance for curve-blocked NCV fits
# (nei$jackknife = TRUE) on the current patched mgcv 1.9-5 build. Run from the repository root:
#   Rscript ncv-jack-recheck/jack.R <cell> <n_reps> <out.csv>
args <- commandArgs(trailingOnly = TRUE)
cell_id <- as.integer(args[1])
reps <- seq_len(as.integer(args[2]))
out <- args[3]
source("scripts/setup.R")
source("ablation/ablation-lib.R")
suppressMessages(library(refund))
cat(
  "mgcv",
  format(packageVersion("mgcv")),
  packageDescription("mgcv")$Patched,
  "\n"
)
cell <- design_table()[design_table()$cell == cell_id, ]
linkinv <- FAMILIES[[cell$fit_family]]$fit()$linkinv
z <- qnorm(1 - (1 - LEVEL) / 2)
rows <- list()
for (r in reps) {
  sim <- simulate(cell, r)
  fit <- fit_model(sim, cell, "NCV")
  nei <- fit$pffr$ncv$nei
  stopifnot(!is.null(nei))
  nei$jackknife <- TRUE
  # same call as fit_pffr(), plus the user-supplied nei with jackknife = TRUE
  ff <- refund::ff
  xind <- sim$xind
  yind <- sim$yind
  bases <- cell_bases(cell)
  fj <- refund::pffr(
    Y ~ ff(X, xind = xind, splinepars = bases$ff) + z,
    data = sim$data,
    yind = yind,
    bs.yindex = bases$t,
    bs.int = bases$int %||% bases$t,
    family = FAMILIES[[cell$fit_family]]$fit(),
    sandwich = "none",
    algorithm = "gam",
    method = "NCV",
    nei = nei
  )
  same <- isTRUE(all.equal(
    unname(fj$coefficients),
    unname(fit$coefficients),
    tolerance = 1e-6
  ))
  Vcl2 <- refund:::pffr_vcov(fit, sandwich = "cl2", cl2_adjustment = "exact")
  covs <- list(bayes_Vp = fit$Vp, jack_Vp = fj$Vp, jack_Vc = fj$Vc, cl2 = Vcl2)
  for (e in c("mean", "beta")) {
    L <- evaluation_matrix(fit, sim, e)
    est <- drop(L %*% fit$coefficients)
    for (a in names(covs)) {
      if (is.null(covs[[a]])) next
      se <- sqrt(rowSums((L %*% covs[[a]]) * L))
      sc <- interval_scores(est, se, sim$truth[[e]], z, e, linkinv)
      rows[[length(rows) + 1]] <- data.frame(
        cell = cell_id,
        family = cell$family,
        rep = r,
        estimand = e,
        arm = a,
        coverage = mean(sc$covered),
        width = mean(sc$width),
        same_fit = same,
        se_ratio_to_bayes = median(se / sqrt(rowSums((L %*% fit$Vp) * L)))
      )
    }
  }
  cat("rep", r, "done\n")
}
write.csv(do.call(rbind, rows), out, row.names = FALSE)
