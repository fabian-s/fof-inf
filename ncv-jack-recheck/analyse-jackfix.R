# Summaries of the corrected NCV covariance run (run-jackfix.R task files), per cell x
# estimand x arm: grid-averaged coverage of nominal 95% pointwise intervals (MC SE over
# replicates), 5% quantile of pointwise coverage (q05), se_sd = median over grid points of
# RMS SE / replicate SD of the estimate (link scale), mean width relative to NCV + CL2 and
# to REML + CL2 (response scale for the mean, as its intervals), interval score.
# Optionally checks the NCV estimates and NCV/REML CL2 SEs against the study's stored records.
#   Rscript ncv-jack-recheck/analyse-jackfix.R <results_dir> <out_prefix> [records_dir | plasmode]
suppressPackageStartupMessages(library(dplyr))
source("scripts/setup.R")
args <- commandArgs(trailingOnly = TRUE)
dir <- args[1]
prefix <- args[2]
records_dir <- if (length(args) >= 3) args[3] else NA
plasmode <- identical(records_dir, "plasmode")
if (plasmode) {
  source("plasmode/sources.R")
  source("plasmode/plasmode.R")
  records_dir <- NA
}
cells_all <- if (plasmode) plasmode_cells() else design_table()
files <- list.files(dir, "^c[0-9]+-r[0-9]+\\.rds$", full.names = TRUE)
cell_of <- as.integer(sub("^c([0-9]+)-.*", "\\1", basename(files)))
ALPHA <- 1 - LEVEL
z <- qnorm(1 - ALPHA / 2)
linkinv_of <- function(cell_id) {
  FAMILIES[[cells_all$fit_family[cells_all$cell == cell_id]]]$linkinv
}

read_cell <- function(fs) {
  res <- lapply(fs, readRDS)
  ok <- vapply(res, \(r) is.null(r$error), TRUE)
  list(
    points = bind_rows(lapply(res[ok], `[[`, "points")),
    info = bind_rows(lapply(res[ok], `[[`, "info")),
    n_fail = sum(!ok),
    errors = unique(unlist(lapply(res[!ok], `[[`, "error")))
  )
}

summarise_cell <- function(pts, linkinv) {
  pts <- pts |>
    mutate(
      lo = est - z * se,
      hi = est + z * se,
      lo = if_else(estimand == "mean", linkinv(lo), lo),
      hi = if_else(estimand == "mean", linkinv(hi), hi),
      covered = lo <= truth & truth <= hi,
      width = hi - lo,
      score = width +
        2 / ALPHA * ((lo - truth) * (truth < lo) + (truth - hi) * (truth > hi))
    )
  ref <- pts |>
    filter(arm %in% c("ncv_cl2", "reml_cl2")) |>
    group_by(estimand, arm) |>
    summarise(w = mean(width), s = mean(score), .groups = "drop")
  w_ncv <- setNames(
    ref$w[ref$arm == "ncv_cl2"],
    ref$estimand[ref$arm == "ncv_cl2"]
  )
  w_reml <- setNames(
    ref$w[ref$arm == "reml_cl2"],
    ref$estimand[ref$arm == "reml_cl2"]
  )
  s_reml <- setNames(
    ref$s[ref$arm == "reml_cl2"],
    ref$estimand[ref$arm == "reml_cl2"]
  )
  pts |>
    group_by(cell, estimand, arm) |>
    group_modify(\(d, k) {
      per_rep <- d |> group_by(rep) |> summarise(cov = mean(covered))
      per_point <- d |>
        group_by(point) |>
        summarise(
          pc = mean(covered),
          sd_est = sd(est),
          rms_se = sqrt(mean(se^2))
        )
      e <- k$estimand
      tibble(
        n_rep = nrow(per_rep),
        coverage = mean(per_rep$cov),
        coverage_se = sd(per_rep$cov) / sqrt(nrow(per_rep)),
        q05 = unname(quantile(per_point$pc, 0.05)),
        se_sd = median(per_point$rms_se / per_point$sd_est),
        width_vs_ncv_cl2 = mean(d$width) / w_ncv[[e]],
        width_vs_reml_cl2 = mean(d$width) / w_reml[[e]],
        score_vs_reml_cl2 = mean(d$score) / s_reml[[e]]
      )
    }) |>
    ungroup()
}

crn_check <- function(pts, cell_id) {
  reps <- sort(unique(pts$rep))
  bind_rows(lapply(reps, \(r) {
    # named files (LRZ) or numbered files of runs with 95 or 123 cells per replicate;
    # a candidate counts only if its stored cell and rep match
    f <- c(
      file.path(records_dir, sprintf("task-c%03d-r%03d.rds", cell_id, r)),
      file.path(
        records_dir,
        sprintf("task-%d.rds", (r - 1) * c(95, 123) + cell_id)
      )
    )
    rec <- NULL
    for (fi in f[file.exists(f)]) {
      x <- readRDS(fi)
      if (
        identical(as.integer(x$cell), cell_id) &&
          identical(as.integer(x$rep), r)
      ) {
        rec <- x$records
        break
      }
    }
    if (is.null(rec)) return(NULL)
    bind_rows(lapply(names(rec), \(e) {
      a <- pts[pts$rep == r & pts$estimand == e & pts$arm == "ncv_cl2", ]
      b <- pts[pts$rep == r & pts$estimand == e & pts$arm == "reml_cl2", ]
      a <- a[order(a$point), ]
      b <- b[order(b$point), ]
      tibble(
        cell = cell_id,
        rep = r,
        estimand = e,
        est_NCV = max(abs(a$est - rec[[e]][, "est_NCV"])),
        se_cl2_NCV = max(abs(a$se - rec[[e]][, "se_cl2_NCV"])),
        est_REML = max(abs(b$est - rec[[e]][, "est_REML"])),
        se_cl2_REML = max(abs(b$se - rec[[e]][, "se_cl2_REML"]))
      )
    }))
  }))
}

out <- list()
info <- list()
crn <- list()
fails <- list()
for (cid in sort(unique(cell_of))) {
  x <- read_cell(files[cell_of == cid])
  fails[[length(fails) + 1]] <- tibble(
    cell = cid,
    n_ok = nrow(x$info),
    n_fail = x$n_fail,
    errors = paste(x$errors, collapse = " | ")
  )
  if (!nrow(x$points)) next
  out[[length(out) + 1]] <- summarise_cell(x$points, linkinv_of(cid))
  info[[length(info) + 1]] <- x$info
  if (!is.na(records_dir)) crn[[length(crn) + 1]] <- crn_check(x$points, cid)
}
keep <- intersect(
  c("cell", "app", "family", "error", "G", "truth", "residual", "D"),
  names(cells_all)
)
res <- bind_rows(out) |>
  left_join(cells_all[, keep], by = "cell") |>
  relocate(all_of(keep))
info <- bind_rows(info)
dir.create(dirname(prefix), showWarnings = FALSE, recursive = TRUE)
readr::write_csv(res, paste0(prefix, "-cells.csv"))
readr::write_csv(
  info |>
    group_by(cell) |>
    summarise(
      n = n(),
      across(
        c(
          time_ncv,
          time_reml,
          time_jack,
          time_cl2,
          alpha_nei,
          alpha_raw,
          tr_ratio_nei,
          tr_ratio_raw
        ),
        median
      )
    ),
  paste0(prefix, "-info.csv")
)
readr::write_csv(bind_rows(fails), paste0(prefix, "-failures.csv"))
if (length(crn)) {
  crn <- bind_rows(crn)
  readr::write_csv(
    crn |>
      group_by(cell, estimand) |>
      summarise(
        n_rep = n(),
        across(est_NCV:se_cl2_REML, max),
        .groups = "drop"
      ),
    paste0(prefix, "-crn.csv")
  )
}
cat(sprintf(
  "%d files, %d cells -> %s-*.csv\n",
  length(files),
  n_distinct(res$cell),
  prefix
))
