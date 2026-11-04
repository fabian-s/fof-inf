# Critical values of the calibrated intervals for beta and the mean (run.R, methods
# <method>_wildL and <method>_wildL_noT), summarised per setting; aggregates only.
#   Rscript competitors/fgee/crit-summary.R <experiment results dir> <out csv>
# at_floor: share of replicates in which the t adjustment sits at fastFGEE's floor of one
# degree of freedom (with-t value / without-t value = qt(.975, 1) / qnorm(.975)).
FLOOR_RATIO <- qt(0.975, 1) / qnorm(0.975)
FLOOR_TOL <- 1e-6

crit_of <- function(records, method, estimand) {
  r <- records[records$method == method & records$estimand == estimand, ]
  stopifnot(nrow(r) > 0)
  median((r$upper - r$lower) / (2 * r$se))
}

task_crits <- function(file) {
  x <- readRDS(file)
  grid <- expand.grid(
    base = paste0("fgee_", c("independent", "ar1", "fpca")),
    estimand = c("beta", "mean"),
    stringsAsFactors = FALSE
  )
  grid$cell <- x$cell$cell
  grid$crit_t <- mapply(
    \(b, e) crit_of(x$records, paste0(b, "_wildL"), e),
    grid$base,
    grid$estimand
  )
  grid$crit_no_t <- mapply(
    \(b, e) crit_of(x$records, paste0(b, "_wildL_noT"), e),
    grid$base,
    grid$estimand
  )
  grid
}

args <- commandArgs(TRUE)
files <- list.files(args[1], pattern = "\\.rds$", full.names = TRUE)
crits <- do.call(rbind, lapply(files, task_crits))
groups <- split(
  crits,
  list(crits$cell, crits$base, crits$estimand),
  drop = TRUE
)
out <- do.call(
  rbind,
  lapply(groups, \(d) {
    data.frame(
      cell = d$cell[1],
      base = d$base[1],
      estimand = d$estimand[1],
      R = nrow(d),
      crit_t_median = median(d$crit_t),
      crit_t_max = max(d$crit_t),
      crit_no_t_median = median(d$crit_no_t),
      crit_no_t_max = max(d$crit_no_t),
      at_floor = mean(abs(d$crit_t / d$crit_no_t - FLOOR_RATIO) < FLOOR_TOL)
    )
  })
)
out <- out[order(out$cell, out$base, out$estimand), ]
write.csv(out, args[2], row.names = FALSE)
cat(args[2], ":", length(files), "tasks,", nrow(out), "rows\n")
