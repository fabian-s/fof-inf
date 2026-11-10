# Detection rate of REML + CL2 with the normal and the Satterthwaite critical value, and the
# position of the NCV estimate relative to the REML + CL2 interval, in the 18 main factorial
# settings at the default signal. Re-scores stored task records; fits nothing.
# Run from the repository root:
#   Rscript satterthwaite/analyse-detect.R raw/satt-results raw/results satterthwaite/summaries
# Writes satt-detect.csv to the output directory.
#
# Detection: among grid points with |truth| > 25% of max |truth|, the share whose interval
# excludes zero on the side of the truth (the definition of the main study).
# ncv_outside: share of grid points at which the NCV estimate lies outside the REML + CL2 interval.

args <- commandArgs(trailingOnly = TRUE)
satt_dir <- args[1]
main_dir <- args[2]
out_dir <- args[3]
stopifnot(dir.exists(satt_dir), dir.exists(main_dir), dir.exists(out_dir))

DETECT_SHARE <- 0.25
ESTIMANDS <- c("beta", "gamma")

satt_files <- list.files(satt_dir, pattern = "\\.rds$", full.names = TRUE)
satt <- lapply(satt_files, readRDS)
ok <- vapply(satt, \(r) identical(r$status, "success"), logical(1))
stopifnot(all(ok))
key <- \(cell, rep) paste(cell, rep)
names(satt) <- vapply(satt, \(r) key(r$cell, r$rep), character(1))
cells <- sort(unique(vapply(satt, \(r) r$cell, integer(1))))

# NCV estimates of the same replicates from the main study records
main_files <- list.files(
  main_dir,
  pattern = "\\.rds$",
  full.names = TRUE,
  recursive = TRUE
)
ncv <- list()
for (f in main_files) {
  r <- readRDS(f)
  if (!r$cell %in% cells) next
  ncv[[key(r$cell, r$rep)]] <- lapply(
    r$records[ESTIMANDS],
    \(m) m[, c("truth", "est_REML", "est_NCV")]
  )
}
stopifnot(setequal(names(ncv), names(satt)))

score_one <- function(id, estimand) {
  s <- satt[[id]]$records[[estimand]]
  n <- ncv[[id]][[estimand]]
  # same replicate, same fit: truths identical, REML estimates equal to numerical precision
  stopifnot(
    isTRUE(all.equal(s[, "truth"], n[, "truth"], tolerance = 1e-10)),
    isTRUE(all.equal(s[, "est_REML"], n[, "est_REML"], tolerance = 1e-6))
  )
  truth <- s[, "truth"]
  est <- s[, "est_REML"]
  se <- s[, "se_cl2e_b"]
  crit <- list(
    z = rep(qnorm(0.975), length(est)),
    satt = qt(0.975, s[, "df_full"])
  )
  strong <- abs(truth) > DETECT_SHARE * max(abs(truth))
  out <- lapply(names(crit), \(cr) {
    lo <- est - crit[[cr]] * se
    hi <- est + crit[[cr]] * se
    detected <- ifelse(truth > 0, lo > 0, hi < 0)
    data.frame(
      crit = cr,
      detect = mean(detected[strong]),
      ncv_outside = mean(n[, "est_NCV"] < lo | n[, "est_NCV"] > hi),
      coverage = mean(truth >= lo & truth <= hi)
    )
  })
  cbind(
    cell = satt[[id]]$cell,
    rep = satt[[id]]$rep,
    estimand = estimand,
    do.call(rbind, out)
  )
}

reps <- do.call(
  rbind,
  lapply(names(satt), \(id) {
    do.call(rbind, lapply(ESTIMANDS, \(e) score_one(id, e)))
  })
)
se_mean <- \(x) sd(x) / sqrt(length(x))
cells_out <- do.call(
  rbind,
  lapply(
    split(reps, list(reps$cell, reps$estimand, reps$crit), drop = TRUE),
    \(d)
      data.frame(
        cell = d$cell[1],
        estimand = d$estimand[1],
        crit = d$crit[1],
        n_rep = nrow(d),
        detect_se = se_mean(d$detect),
        detect = mean(d$detect),
        ncv_outside_se = se_mean(d$ncv_outside),
        ncv_outside = mean(d$ncv_outside),
        coverage = mean(d$coverage)
      )
  )
)
cells_out <- cells_out[
  order(cells_out$cell, cells_out$estimand, cells_out$crit),
]
rownames(cells_out) <- NULL
write.csv(cells_out, file.path(out_dir, "satt-detect.csv"), row.names = FALSE)
cat(
  "cells:",
  length(cells),
  " replicates:",
  length(satt),
  " rows:",
  nrow(cells_out),
  "\n"
)
