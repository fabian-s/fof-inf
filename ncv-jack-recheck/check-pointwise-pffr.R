# Study-scale check of jackfix_vcov(): one replicate refitted by pffr with per-point
# neighbourhoods (point i predicted with its whole curve dropped; mgcv's
# neighbourhood-aware covariance, jackknife unset) vs jackfix_vcov() of the study's block fit.
#   Rscript ncv-jack-recheck/check-pointwise-pffr.R <cell> <rep>   (from the repository root)
source("scripts/setup.R")
source("ncv-jack-recheck/jackfix-lib.R")
suppressMessages(library(refund))
args <- commandArgs(TRUE)
cid <- as.integer(args[1])
r <- as.integer(args[2])
cat(
  "mgcv",
  format(packageVersion("mgcv")),
  packageDescription("mgcv")$Patched,
  "\n"
)
cell <- design_table()[design_table()$cell == cid, ]
sim <- simulate(cell, r)
fit <- fit_model(sim, cell, "NCV")
jf <- jackfix_vcov(fit, ncv_cluster(fit))
nei_pt <- nei_pointwise(ncv_cluster(fit))
ff <- refund::ff
xind <- sim$xind
yind <- sim$yind
bases <- cell_bases(cell)
t0 <- proc.time()[["elapsed"]]
fp <- refund::pffr(
  Y ~ ff(X, xind = xind, splinepars = bases$ff) + z,
  data = sim$data,
  yind = yind,
  bs.yindex = bases$t,
  bs.int = bases$int %||% bases$t,
  family = FAMILIES[[cell$fit_family]]$fit(),
  sandwich = "none",
  algorithm = "gam",
  method = "NCV",
  nei = nei_pt
)
t_pt <- proc.time()[["elapsed"]] - t0
rd <- \(A, B) max(abs(A - B)) / max(abs(B))
out <- data.frame(
  cell = cid,
  rep = r,
  time_block_fit = fit$restart$time,
  time_pointwise_fit = t_pt,
  max_abs_diff_log_sp = max(abs(log(fp$sp) - log(fit$sp))),
  reldiff_coef = rd(fp$coefficients, fit$coefficients),
  inv_learn_pointwise = attr(fp$Ve, "inv.learn"),
  alpha_jackfix = jf$nei$alpha,
  reldiff_Vp = rd(jf$nei$Vp, fp$Vp),
  reldiff_Vc = rd(jf$nei$Vc, fp$Vc)
)
print(out)
dir.create("ncv-jack-recheck/summaries", showWarnings = FALSE)
write.csv(
  out,
  sprintf(
    "ncv-jack-recheck/summaries/jackfix-pointwise-check-c%d-r%d.csv",
    cid,
    r
  ),
  row.names = FALSE
)
