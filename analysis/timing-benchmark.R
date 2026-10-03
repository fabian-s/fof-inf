# Timing of the complete recipes: fit time and CL2 time per method, on the Gaussian cells
# that include the AR(1) comparator. Run from the repository root:
#   Rscript analysis/timing-benchmark.R <out.csv> [n_rep]
source("scripts/setup.R")
suppressPackageStartupMessages(library(refund))
args <- commandArgs(trailingOnly = TRUE)
n_rep <- if (length(args) >= 2) as.integer(args[2]) else 10L

cells <- design_table()
ar1_cells <- cells[grepl("ar1", cells$extra), ]
rows <- list()
for (i in seq_len(nrow(ar1_cells))) {
  cell <- ar1_cells[i, ]
  for (rep in seq_len(n_rep)) {
    sim <- simulate(cell, rep)
    for (method in c("REML", "NCV", "AR1")) {
      inputs <- method_inputs(fit_model(sim, cell, method), method, sim, cell)
      rows[[length(rows) + 1]] <- tibble::tibble(
        cell_id = cell$cell, # distinct names: tibble() sees its own new columns
        error = cell$error,
        rep = rep,
        method = method,
        fit_time = inputs$fit_row$time,
        cl2_time = inputs$fit_row$cl2_time
      )
    }
  }
}
timing <- dplyr::bind_rows(rows)
readr::write_csv(timing, args[1])
print(
  timing |>
    dplyr::group_by(method) |>
    dplyr::summarise(
      fit = median(fit_time),
      cl2 = median(cl2_time),
      total = median(fit_time + dplyr::coalesce(cl2_time, 0))
    )
)
