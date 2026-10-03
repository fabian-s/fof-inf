# Gaussian NCV lambda correction on the plasmode cells: aggregate summaries from the task
# records of run-lc-ncv-plasmode.R (run on LRZ: the cardiology task files stay there).
#   Rscript lambda-correction/analyse-lc-ncv-plasmode.R <results_dir> <out_dir>
LC_NCV_LIB_ONLY <- TRUE
source("lambda-correction/analyse-lc-ncv.R")
source("plasmode/sources.R")
source("plasmode/plasmode.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2)
  stop("usage: analyse-lc-ncv-plasmode.R results_dir out_dir")
dir.create(args[2], recursive = TRUE, showWarnings = FALSE)
cells <- plasmode_cells()
cells <- cells[cells$app != "ecg_gls4", ]
cells$family <- cells$app
cells$error <- paste0(
  cells$truth,
  ":",
  cells$residual,
  ifelse(cells$flip == "block", " block", "")
)
files_all <- list.files(args[1], "^task-c\\d+-r\\d+\\.rds$", full.names = TRUE)
out <- lapply(seq_len(nrow(cells)), \(i) {
  f <- grep(sprintf("task-c%03d-", cells$cell[i]), files_all, value = TRUE)
  if (length(f)) ncv_cell(f, cells[i, ])
})
out <- out[!vapply(out, is.null, TRUE)]
bind <- \(nm) dplyr::bind_rows(lapply(out, `[[`, nm))
res <- list(
  cells = bind("cells"),
  contrasts = bind("contrasts"),
  calib = bind("calib"),
  targeting = bind("targeting"),
  crn = bind("crn"),
  failures = bind("failures"),
  times = bind("times")
)
for (nm in names(res))
  write.csv(
    res[[nm]],
    file.path(args[2], paste0("ncvp-", nm, ".csv")),
    row.names = FALSE
  )
cat("wrote", length(out), "cells to", args[2], "\n")
