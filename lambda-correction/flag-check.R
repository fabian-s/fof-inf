# Does a large per-fit SD(rho-hat) flag replicates with high beta relative error or CL2 misses?
# Run from the repository root: Rscript lambda-correction/flag-check.R <lc_results_dir> <out.csv>
args <- commandArgs(trailingOnly = TRUE)
source("scripts/setup.R")
source("ablation/ablation-lib.R")
cells <- ablation_cells()
q <- qnorm(0.975)
auc <- function(score, y) {
  if (length(unique(y)) < 2) return(NA)
  r <- rank(score)
  n1 <- sum(y)
  (sum(r[y]) - n1 * (n1 + 1) / 2) / (n1 * sum(!y))
}
rows <- list()
for (f in list.files(args[1], "rds$", full.names = TRUE)) {
  r <- readRDS(f)
  m <- r$records$beta
  om <- r$omega
  k <- grep("L.X1$", names(om$hat))
  rows[[length(rows) + 1]] <- data.frame(
    cell = r$cell,
    rep = r$rep,
    rel = sqrt(mean((m[, "est_REML"] - m[, "truth"])^2) / mean(m[, "truth"]^2)),
    miss = mean(abs(m[, "est_REML"] - m[, "truth"]) > q * m[, "se_cl2"]),
    rho = om$hat[k],
    sd_ij = om$ij[k],
    sd_os = om$os[k],
    sd_naive = om$naive[k],
    sd_ij_max = max(om$ij[seq_along(r$omega$hat)[
      !grepl("scale", names(om$hat))
    ]])
  )
}
d <- merge(do.call(rbind, rows), cells[, c("cell", "family", "error", "G")])
d$dep <- ifelse(d$error == "iid", "iid", "dependent")
preds <- c("sd_ij", "sd_os", "sd_naive", "sd_ij_max", "rho")
res <- do.call(
  rbind,
  lapply(split(d, d$cell), \(x) {
    hi_rel <- x$rel > quantile(x$rel, 0.8)
    hi_miss <- x$miss > quantile(x$miss, 0.8)
    do.call(
      rbind,
      lapply(preds, \(p) {
        s <- if (p == "rho") -x[[p]] else x[[p]]
        data.frame(
          cell = x$cell[1],
          family = x$family[1],
          dep = x$dep[1],
          G = x$G[1],
          pred = p,
          rho_rel = cor(s, x$rel, method = "spearman"),
          rho_miss = cor(s, x$miss, method = "spearman"),
          auc_rel = auc(s, hi_rel),
          auc_miss = auc(s, hi_miss)
        )
      })
    )
  })
)
options(width = 160)
cat(
  "Within-cell, mean over cells (rho predictor sign-flipped: small lambda-hat = flag)\n"
)
print(
  aggregate(
    cbind(rho_rel, rho_miss, auc_rel, auc_miss) ~ dep + G + pred,
    res,
    mean
  ),
  digits = 2
)
cat("\nDependent cells, per family (sd_ij):\n")
print(
  aggregate(
    cbind(rho_rel, rho_miss, auc_rel, auc_miss) ~ family + G,
    res[res$dep == "dependent" & res$pred == "sd_ij", ],
    mean
  ),
  digits = 2
)
# top-decile flag: miss rate and rel error in flagged vs rest, dependent cells
d2 <- d[d$dep == "dependent", ]
d2 <- do.call(
  rbind,
  lapply(split(d2, d2$cell), \(x) {
    x$flag <- x$sd_ij > quantile(x$sd_ij, 0.9)
    x
  })
)
cat("\nDependent cells: top-decile sd_ij within cell vs rest\n")
print(aggregate(cbind(miss, rel) ~ flag + G, d2, mean), digits = 3)
write.csv(res, args[2], row.names = FALSE)
