# Are the non-converged scat NCV fits real optimization failures? For a few non-converged
# production tasks, refit NCV from different starting penalties and compare the NCV
# criterion, the estimates and the standard errors. Run from the repository root:
#   Rscript analysis/scat-convergence.R <out.csv> [n_tasks] [cell] [part] [n_parts]
# Without `cell`: the first flagged tasks. With `cell`: a seeded random sample of its flagged
# tasks, split into `n_parts` for parallel workers. STUDY_RDS: the study's summarise_study() output.
source("scripts/setup.R")
suppressPackageStartupMessages(library(refund))
args <- commandArgs(trailingOnly = TRUE)
n_tasks <- if (length(args) >= 2) as.integer(args[2]) else 5L
one_cell <- if (length(args) >= 3) as.integer(args[3]) else NA_integer_
part <- if (length(args) >= 4) as.integer(args[4]) else 1L
n_parts <- if (length(args) >= 5) as.integer(args[5]) else 1L

cells <- design_table()
study <- readRDS(Sys.getenv("STUDY_RDS", "raw/study.rds"))
scat_cells <- cells$cell[cells$family == "scat"]
flagged <- study$fits |>
  dplyr::filter(cell %in% scat_cells, fit == "NCV", !converged)
failed <- if (is.na(one_cell)) {
  dplyr::slice_head(flagged, n = n_tasks)
} else {
  sampled <- dplyr::filter(flagged, cell == one_cell) |>
    dplyr::arrange(rep)
  sampled <- withr::with_seed(
    SEED_BOOT,
    dplyr::slice_sample(sampled, n = n_tasks)
  )
  sampled[seq_len(nrow(sampled)) %% n_parts == part %% n_parts, ]
}

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

se_cl2 <- function(fit, L) {
  V <- refund:::pffr_vcov(
    fit,
    sandwich = "cl2",
    freq = FALSE,
    cl2_adjustment = "exact"
  )
  sqrt(rowSums((L %*% V) * L))
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
      beta_max_diff = max(abs(L %*% (f$coefficients - base$coefficients))),
      se_cl2_rel_diff = max(abs(se_cl2(f, L) / se_cl2(base, L) - 1))
    )
  }))
})
out <- dplyr::bind_rows(rows)
readr::write_csv(out, args[1])
print(as.data.frame(out), digits = 4)
