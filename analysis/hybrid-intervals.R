# Hybrid interval "NCV estimate + REML CL2 SE" and pointwise coverage maps, from the
# stored per-point task records (no refitting).
#
# Analysis 1, per cell x estimand, over replicates (MC SEs, 5% quantile of pointwise
# coverage over the grid):
#   inclusion: share of points where est_NCV lies inside est_REML +- z se_cl2_REML;
#   arms reml_cl2, ncv_cl2, ncv_cl2_bias (R/summarise.R) and `hybrid` = est_NCV +- z
#   se_cl2_REML, scored by the same point_scores(); paired contrasts on the same replicates.
# Analysis 2 (beta, selected cells): per (s, t) coverage of each arm over replicates, the
# truth, its curvature |d2 beta/ds2| + |d2 beta/dt2| (second differences on the eval
# grid), and the replicate mean of (est_NCV - est_REML)^2 / se_cl2_NCV^2; plus per map the
# correlation of the coverage shortfall with curvature and with that ratio.
#
# Run from a copy of the repository (it sources scripts/setup.R), BLAS pinned:
#   Rscript <this file> synthetic <results_dir[,results_dir...]> <out_dir>
#   Rscript <this file> plasmode <plasmode_results_root> <out_dir> <app> [<app> ...]
# Output (aggregates only, nothing per subject): hybrid-summary.csv, hybrid-inclusion.csv,
# hybrid-paired.csv, hybrid-maps.csv, hybrid-map-correlations.csv.

# Setup ------------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
MODE <- args[1]
stopifnot(
  "mode must be synthetic or plasmode" = MODE %in% c("synthetic", "plasmode")
)
source("scripts/setup.R")
if (MODE == "plasmode") {
  source("plasmode/sources.R")
  source("plasmode/plasmode.R")
}
suppressPackageStartupMessages(library(dplyr))

HYBRID_ARMS <- c("reml_cl2", "ncv_cl2", "ncv_cl2_bias", "hybrid")
# Only the arms compared here are scored (task_scores() iterates over ARMS).
ARMS <- rbind(
  ARMS[ARMS$arm %in% HYBRID_ARMS, ],
  data.frame(
    arm = "hybrid",
    est = "est_NCV",
    se = "se_cl2_REML",
    bias_aware = FALSE
  )
)
CONTRASTS <- list(
  c("hybrid", "reml_cl2"),
  c("ncv_cl2_bias", "reml_cl2"),
  c("hybrid", "ncv_cl2_bias")
)
Z <- qnorm(1 - (1 - LEVEL) / 2)
CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))

# Cells and files ----------------------------------------------------------------------

# Map cells, synthetic: rough-truth fpc cells, the smooth-truth core fpc cells at the same
# family, signal and G = 100, and the xlarge-basis fpc cells.
synthetic_map_cells <- function(cells) {
  rough <- cells[cells$block == "rough_truth" & cells$error == "fpc", ]
  smooth <- cells$block == "core" &
    cells$error == "fpc" &
    cells$G == 100L &
    paste(cells$family, cells$signal) %in% paste(rough$family, rough$signal)
  xlarge <- cells$block == "basis_size" &
    cells$basis == "xlarge" &
    cells$error == "fpc"
  cells$cell[cells$cell %in% rough$cell | smooth | xlarge]
}

# Map cells, plasmode: running and DTI, every truth source with REML residuals, all
# subjects, attached errors with curve flips.
plasmode_map_cells <- function(cells) {
  cells$cell[
    cells$app %in%
      c("running", "dti") &
      cells$residual == "REML" &
      cells$flip == "curve" &
      cells$error == "attached" &
      cells$G > PLASMODE_G_SMALL
  ]
}

# Task files grouped by the cell stored in each file (older runs name files by task id,
# which depends on that run's design); failed tasks carry no cell and are counted apart.
files_by_cell <- function(dirs, cores) {
  files <- unlist(lapply(
    dirs,
    list.files,
    pattern = "^task-.*\\.rds$",
    full.names = TRUE
  ))
  cell <- unlist(parallel::mclapply(files, mc.cores = cores, \(f) {
    r <- readRDS(f)
    if (identical(r$status, "success")) r$cell else NA_integer_
  }))
  cat(sum(is.na(cell)), "failed task files\n")
  split(files[!is.na(cell)], cell[!is.na(cell)])
}

# Analysis 1 -----------------------------------------------------------------------------

# Per task and scored estimand: inclusion of est_NCV in the REML CL2 interval (a
# monotone link leaves inclusion unchanged, so it is computed on the link scale).
task_inclusion <- function(result, cell) {
  estimands <- scored_estimands(cell, names(result$records))
  dplyr::bind_rows(lapply(estimands, \(estimand) {
    r <- result$records[[estimand]]
    inside <- abs(r[, "est_NCV"] - r[, "est_REML"]) <= Z * r[, "se_cl2_REML"]
    tibble::tibble(
      cell = cell$cell,
      rep = result$rep,
      estimand = estimand,
      inclusion = mean(inside),
      included = list(inside)
    )
  }))
}

summarise_inclusion <- function(inclusion) {
  inclusion |>
    group_by(cell, estimand) |>
    summarise(
      n_rep = n(),
      inclusion_se = sd(inclusion) / sqrt(n_rep),
      inclusion_mean = mean(inclusion),
      point_stats(
        matrix(0, length(included[[1]]), n_rep),
        do.call(cbind, included)
      )[c("pointwise_q05", "pointwise_q05_se")],
      .groups = "drop"
    )
}

paired_all <- function(scores) {
  metrics <- c("score", "coverage", "width")
  dplyr::bind_rows(lapply(CONTRASTS, \(pair) {
    dplyr::bind_rows(lapply(metrics, \(metric) {
      paired_contrast(scores, pair[1], pair[2], metric) |>
        mutate(
          arm_a = pair[1],
          arm_b = pair[2],
          metric = metric,
          ratio_lower = vapply(ratio_ci, `[`, 1, 1),
          ratio_upper = vapply(ratio_ci, `[`, 1, 2)
        ) |>
        select(-ratio_ci)
    }))
  }))
}

# Analysis 2 -----------------------------------------------------------------------------

# Second differences of a matrix along rows (dim 1) or columns, with spacing h; edge rows
# take the value of their interior neighbour.
second_diff <- function(m, h, along) {
  if (along == 2) return(t(second_diff(t(m), h, 1)))
  n <- nrow(m)
  d <- (m[3:n, , drop = FALSE] -
    2 * m[2:(n - 1), , drop = FALSE] +
    m[1:(n - 2), , drop = FALSE]) /
    h^2
  rbind(d[1, ], d, d[n - 2, ])
}

curvature <- function(beta) {
  m <- matrix(beta, length(EVAL_S), length(EVAL_T)) # s varies fastest
  abs(second_diff(m, diff(EVAL_S)[1], 1)) +
    abs(second_diff(m, diff(EVAL_T)[1], 2))
}

# Pointwise beta coverage per arm over replicates, with truth, curvature and delta ratio.
coverage_map <- function(ok, scores, cell_id) {
  records <- lapply(ok, \(r) r$records$beta)
  truth <- records[[1]][, "truth"]
  drift <- max(vapply(records, \(r) max(abs(r[, "truth"] - truth)), 0))
  if (drift > 1e-8)
    stop("cell ", cell_id, ": beta truth varies over replicates")
  ratio <- rowMeans(do.call(
    cbind,
    lapply(records, \(r) {
      (r[, "est_NCV"] - r[, "est_REML"])^2 / r[, "se_cl2_NCV"]^2
    })
  ))
  points <- eval_points("beta")
  beta_scores <- scores[scores$estimand == "beta", ]
  dplyr::bind_rows(lapply(HYBRID_ARMS, \(a) {
    covered <- do.call(cbind, beta_scores$covered[beta_scores$arm == a])
    tibble::tibble(
      cell = cell_id,
      s = points$s,
      t = points$t,
      method = a,
      coverage = rowMeans(covered),
      n_rep = ncol(covered),
      truth_value = truth,
      curvature = as.vector(curvature(truth)),
      delta_ratio = ratio
    )
  }))
}

map_correlations <- function(maps) {
  maps |>
    mutate(shortfall = LEVEL - coverage) |>
    group_by(cell, method) |>
    summarise(
      coverage_mean = mean(coverage),
      coverage_q05 = quantile(coverage, 0.05),
      share_below_90 = mean(coverage < 0.9),
      cor_curv_spearman = cor(shortfall, curvature, method = "spearman"),
      cor_curv_pearson = cor(shortfall, curvature),
      cor_delta_spearman = cor(shortfall, delta_ratio, method = "spearman"),
      cor_delta_pearson = cor(shortfall, delta_ratio),
      cor_curv_delta_spearman = cor(
        curvature,
        delta_ratio,
        method = "spearman"
      ),
      .groups = "drop"
    )
}

# One cell ---------------------------------------------------------------------------------

analyse_cell <- function(files, cells, map_cells) {
  res <- load_results(files)
  if (!length(res$ok)) return(NULL)
  cell <- as.list(cells[cells$cell == res$ok[[1]]$cell, ])
  stopifnot(
    "files of one cell" = all(vapply(res$ok, \(r) r$cell == cell$cell, TRUE))
  )
  scores <- task_scores(res$ok, cells) # checks each stored cell spec
  inclusion <- dplyr::bind_rows(lapply(res$ok, task_inclusion, cell = cell))
  list(
    summary = summarise_cells(scores),
    inclusion = summarise_inclusion(inclusion),
    paired = paired_all(scores),
    maps = if (cell$cell %in% map_cells) {
      coverage_map(res$ok, scores, cell$cell)
    }
  )
}

# Main -----------------------------------------------------------------------------------

out_dir <- args[3]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
if (MODE == "synthetic") {
  cells <- design_table()
  dirs <- strsplit(args[2], ",")[[1]]
  map_cells <- synthetic_map_cells(cells)
  keys <- setdiff(names(cells), "extra")
} else {
  cells <- plasmode_cells()
  dirs <- file.path(args[2], args[-(1:3)])
  map_cells <- plasmode_map_cells(cells)
  keys <- c("cell", "app", "truth", "residual", "flip", "error", "G")
}
by_cell <- files_by_cell(dirs, CORES)
cat(length(by_cell), "cells,", length(unlist(by_cell)), "task files\n")
per_cell <- parallel::mclapply(
  by_cell,
  analyse_cell,
  cells = cells,
  map_cells = map_cells,
  mc.cores = CORES,
  mc.preschedule = FALSE
)
broken <- vapply(per_cell, inherits, TRUE, "try-error")
if (any(broken)) {
  print(per_cell[broken])
  stop("cells failed: ", toString(names(per_cell)[broken]))
}
part <- \(name) dplyr::bind_rows(lapply(per_cell, `[[`, name))
design <- cells[keys]
write <- \(x, name) {
  readr::write_csv(
    dplyr::inner_join(design, x, by = "cell"),
    file.path(out_dir, paste0("hybrid-", name, ".csv"))
  )
}
write(part("summary"), "summary")
write(part("inclusion"), "inclusion")
write(part("paired"), "paired")
maps <- part("maps")
write(maps, "maps")
write(map_correlations(maps), "map-correlations")
