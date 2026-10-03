# Studentised-error diagnostics for REML + CL2 from the stored per-point task records
# (no refitting): per cell and estimand, over all grid points and replicates,
#   z = (est_REML - truth) / se_cl2_REML
# and its mean, SD, and the share of |z| > 1.96, plus the same for the NCV fit with its
# CL2 SE and for the model-based SE. Run on LRZ from a copy of the repository (it sources
# scripts/setup.R), never locally:
#   Rscript writeup/lrz-zstats.R synthetic <results_dir[,results_dir...]> <out.csv>
#   Rscript writeup/lrz-zstats.R plasmode <plasmode_results_root> <out.csv> <app> [<app> ...]
# Output: one row per cell x estimand x fit (aggregates only, nothing per subject).

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
CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
Z <- qnorm(1 - (1 - LEVEL) / 2)
FITS <- c(REML = "REML", NCV = "NCV")

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
  split(files[!is.na(cell)], cell[!is.na(cell)])
}

# z-scores of one task for every estimand and fit, with the CL2 and the model-based SE.
task_z <- function(result, cell) {
  estimands <- scored_estimands(cell, names(result$records))
  dplyr::bind_rows(lapply(estimands, \(estimand) {
    r <- result$records[[estimand]]
    dplyr::bind_rows(lapply(names(FITS), \(fit) {
      est <- r[, paste0("est_", fit)]
      se_cl2 <- r[, paste0("se_cl2_", fit)]
      se_model <- if (paste0("se_model_", fit) %in% colnames(r))
        r[, paste0("se_model_", fit)] else NA_real_
      z_cl2 <- (est - r[, "truth"]) / se_cl2
      z_model <- (est - r[, "truth"]) / se_model
      tibble::tibble(
        cell = cell$cell,
        rep = result$rep,
        estimand = estimand,
        fit = fit,
        n_points = length(z_cl2),
        sum_z = sum(z_cl2),
        sum_z2 = sum(z_cl2^2),
        n_big = sum(abs(z_cl2) > Z),
        sum_zm = sum(z_model),
        sum_zm2 = sum(z_model^2),
        n_bigm = sum(abs(z_model) > Z)
      )
    }))
  }))
}

analyse_cell <- function(files, cells) {
  res <- load_results(files)
  if (!length(res$ok)) return(NULL)
  cell <- as.list(cells[cells$cell == res$ok[[1]]$cell, ])
  dplyr::bind_rows(lapply(res$ok, task_z, cell = cell)) |>
    dplyr::group_by(cell, estimand, fit) |>
    dplyr::summarise(
      n_rep = dplyr::n(),
      n = sum(n_points),
      z_mean_cl2 = sum(sum_z) / n,
      z_sd_cl2 = sqrt(sum(sum_z2) / n - z_mean_cl2^2),
      share_big_cl2 = sum(n_big) / n,
      z_mean_model = sum(sum_zm) / n,
      z_sd_model = sqrt(sum(sum_zm2) / n - z_mean_model^2),
      share_big_model = sum(n_bigm) / n,
      .groups = "drop"
    )
}

out <- args[3]
if (MODE == "synthetic") {
  cells <- design_table()
  dirs <- strsplit(args[2], ",")[[1]]
  keys <- setdiff(names(cells), "extra")
} else {
  cells <- plasmode_cells()
  dirs <- file.path(args[2], args[-(1:3)])
  keys <- c("cell", "app", "truth", "residual", "flip", "error", "G")
}
by_cell <- files_by_cell(dirs, CORES)
cat(length(by_cell), "cells,", length(unlist(by_cell)), "task files\n")
per_cell <- parallel::mclapply(
  by_cell,
  analyse_cell,
  cells = cells,
  mc.cores = CORES,
  mc.preschedule = FALSE
)
broken <- vapply(per_cell, inherits, TRUE, "try-error")
if (any(broken)) stop("cells failed: ", toString(names(per_cell)[broken]))
zstats <- dplyr::inner_join(
  cells[keys],
  dplyr::bind_rows(per_cell),
  by = "cell"
)
readr::write_csv(zstats, out)
cat("written", out, "\n")
