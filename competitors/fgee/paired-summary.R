# Paired MSE ratios (ratio of aggregates over replicates) on the same data sets, with
# percentile bootstrap intervals over replicates. Input: pair-ncv.R's CSV.
#   Rscript competitors/fgee/paired-summary.R <paired csv> <out csv> <cell metadata csv>
# (metadata: the experiment's analyse.R summary, e.g. B-summary.csv or P-summary.csv).
# "fgee_best_vs_ncv" uses the fastFGEE variant with the lowest cell-mean MSE, chosen after
# the fact on the full data: a descriptive sensitivity contrast favourable to fastFGEE whose
# interval does not account for the selection.
args <- commandArgs(TRUE)
x <- read.csv(args[1])
stopifnot(
  !anyDuplicated(x[, c("cell", "rep", "estimand")]),
  max(x$crn_max_diff) < 1e-8, # same data and same REML fit in both runs
  max(x$truth_max_diff) < 1e-10
)
meta <- read.csv(args[3])
cells <- unique(meta[, intersect(
  c("cell", "family", "error", "G", "app"),
  names(meta)
)])
set.seed(1)
ratio_ci <- function(num, den, B = 2000) {
  est <- mean(num) / mean(den)
  boot <- replicate(B, {
    i <- sample.int(length(num), replace = TRUE)
    mean(num[i]) / mean(den[i])
  })
  c(
    est = est,
    lower = unname(quantile(boot, 0.025)),
    upper = unname(quantile(boot, 0.975))
  )
}
out <- do.call(
  rbind,
  lapply(split(x, list(x$cell, x$estimand), drop = TRUE), \(d) {
    r <- rbind(
      ncv_vs_reml = ratio_ci(d$se_ncv, d$se_reml),
      fgee_ind_vs_reml = ratio_ci(d$se_fgee_ind, d$se_reml),
      fgee_ind_vs_ncv = ratio_ci(d$se_fgee_ind, d$se_ncv),
      fgee_best_vs_ncv = ratio_ci(
        pmin(d$se_fgee_ind, d$se_fgee_ar1, d$se_fgee_fpca) *
          0 +
          d[[c(
            "se_fgee_ind",
            "se_fgee_ar1",
            "se_fgee_fpca"
          )[which.min(colMeans(d[, c(
            "se_fgee_ind",
            "se_fgee_ar1",
            "se_fgee_fpca"
          )]))]]],
        d$se_ncv
      )
    )
    data.frame(
      cell = d$cell[1],
      estimand = d$estimand[1],
      R = nrow(d),
      contrast = rownames(r),
      r,
      row.names = NULL
    )
  })
)
out <- merge(out, cells, by = "cell")
write.csv(out, args[2], row.names = FALSE)
for (est in c("beta", "alpha", "gamma")) {
  cat("\n== MSE ratios (paired), ", est, "\n")
  y <- out[out$estimand == est, ]
  y$lab <- if ("app" %in% names(y)) paste(y$app, y$G) else
    paste(y$family, y$error, y$G)
  y$txt <- sprintf("%.2f [%.2f,%.2f]", y$est, y$lower, y$upper)
  print(
    reshape(
      y[, c("lab", "contrast", "txt")],
      idvar = "lab",
      timevar = "contrast",
      direction = "wide"
    ),
    row.names = FALSE
  )
}
