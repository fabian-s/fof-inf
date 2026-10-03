# Worst-case concurvity of every smooth term of the REML fit to replicate 1 of every cell
# of the design table (the check run_task makes on CONCURVITY_REPS). Run from the repository
# directory:
#   Rscript analysis/concurvity.R [out.csv] [cores]
source("scripts/setup.R")
suppressPackageStartupMessages(library(refund))
args <- commandArgs(trailingOnly = TRUE)
out_file <- if (length(args) >= 1) args[1] else "summaries/final/concurvity.csv"
cores <- if (length(args) >= 2) as.integer(args[2]) else 3L
REP <- 1

cell_concurvity <- function(cell, rep) {
  sim <- simulate(cell, rep)
  conc <- check_concurvity(fit_model(sim, cell, "REML"), cell)
  data.frame(cell = cell$cell, rep = rep, conc)
}

cells <- design_table()
rows <- parallel::mclapply(
  seq_len(nrow(cells)),
  \(i) cell_concurvity(cells[i, ], REP),
  mc.cores = cores
)
failed <- vapply(rows, inherits, logical(1), "try-error")
if (any(failed)) {
  stop("cells failed: ", paste(cells$cell[failed], collapse = ", "))
}
design_columns <- c("cell", "block", "family", "error", "G", "D", "model")
concurvity <- dplyr::bind_rows(rows) |>
  dplyr::left_join(cells[design_columns], by = "cell")
readr::write_csv(concurvity, out_file)

print(
  concurvity |>
    dplyr::group_by(block, model, G) |>
    dplyr::summarise(
      max = max(concurvity),
      term = term[which.max(concurvity)],
      .groups = "drop"
    ),
  n = Inf
)
