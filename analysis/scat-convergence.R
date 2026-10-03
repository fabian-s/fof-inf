# Are the non-converged scat NCV fits real optimization failures? For a few non-converged
# production tasks, refit NCV from different starting penalties and compare the NCV
# criterion and the estimates. Run from the repository root:
#   Rscript analysis/scat-convergence.R <out.csv> [n_tasks]
source("scripts/setup.R")
suppressPackageStartupMessages(library(refund))
args <- commandArgs(trailingOnly = TRUE)
n_tasks <- if (length(args) >= 2) as.integer(args[2]) else 5L

cells <- design_table()
study <- readRDS("summaries/final/study.rds")
scat_cells <- cells$cell[cells$family == "scat"]
failed <- study$fits |>
  dplyr::filter(cell %in% scat_cells, fit == "NCV", !converged) |>
  dplyr::slice_head(n = n_tasks)

refit <- function(sim, cell, start_sp = NULL) {
  fit <- fit_pffr(sim, cell, "NCV")
  if (is.null(start_sp)) return(fit)
  # same model, NCV started from the given penalties (mgcv's in.out)
  fit_start <- try(
    refund::pffr(
      Y ~ ff(X, xind = sim$xind, splinepars = BASIS_FF) + z,
      data = sim$data,
      yind = sim$yind,
      bs.yindex = BASIS_T,
      bs.int = BASIS_T,
      family = FAMILIES$scat$fit(),
      sandwich = "none",
      method = "NCV",
      in.out = list(sp = start_sp, scale = 1)
    ),
    silent = TRUE
  )
  if (inherits(fit_start, "try-error")) NULL else fit_start
}

rows <- lapply(seq_len(nrow(failed)), \(k) {
  cell <- cells[cells$cell == failed$cell[k], ]
  sim <- simulate(cell, failed$rep[k])
  base <- refit(sim, cell)
  reml <- fit_pffr(sim, cell, "REML")
  starts <- list(
    production = NULL,
    reml_sp = reml$sp,
    reml_sp_x10 = 10 * reml$sp
  )
  L <- evaluation_matrix(base, sim, "beta")
  dplyr::bind_rows(lapply(names(starts), \(s) {
    f <- if (s == "production") base else refit(sim, cell, starts[[s]])
    if (is.null(f))
      return(tibble::tibble(task = failed$task[k], start = s, ok = FALSE))
    tibble::tibble(
      task = failed$task[k],
      cell_id = failed$cell[k],
      start = s,
      ok = TRUE,
      conv = f$outer.info$conv %||% NA_character_,
      criterion = as.numeric(f$gcv.ubre),
      max_grad = max(abs(f$outer.info$grad %||% NA)),
      beta_max_diff = max(abs(L %*% (f$coefficients - base$coefficients)))
    )
  }))
})
out <- dplyr::bind_rows(rows)
readr::write_csv(out, args[1])
print(as.data.frame(out), digits = 4)
