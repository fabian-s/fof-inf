# LRZ worker for the bias-correction × λ-correction combination (Gaussian NCV), synthetic or
# plasmode cells, rep-major tasks claimed via lock dirs. Plasmode outputs stay on LRZ.
#   Rscript ncv-rbc/worker-combo.R synthetic <out_dir> <cell,cell,...> <rep_from> <rep_to>
#   Rscript ncv-rbc/worker-combo.R plasmode <out_dir> <app,app,...> <rep_from> <rep_to>
#   Rscript ncv-rbc/worker-combo.R misspec <out_dir> <misspec cell,cell,...> <rep_from> <rep_to>
source("scripts/setup.R")
source("lambda-correction/lc-lib.R")
source("lambda-correction/ncv-lib.R")
source("ncv-rbc/combo-lib.R")
args <- commandArgs(trailingOnly = TRUE)
mode <- args[1]
out <- args[2]
reps <- seq(as.integer(args[4]), as.integer(args[5]))
if (mode == "synthetic") {
  cells <- design_table()
  cells <- cells[cells$cell %in% as.integer(strsplit(args[3], ",")[[1]]), ]
  sim_fun <- simulate
} else if (mode == "misspec") {
  source("misspec/misspec-lib.R") # masks make_truth(): simulate() builds the misspecified truths
  cells <- misspec_cells()
  cells <- cells[cells$cell %in% as.integer(strsplit(args[3], ",")[[1]]), ]
  sim_fun <- simulate
} else {
  source("plasmode/sources.R")
  source("plasmode/plasmode.R")
  apps <- strsplit(args[3], ",")[[1]]
  if (any(apps %in% PLASMODE_DEID_APPS) && !identical(PLASMODE_CARDIO, "deid"))
    stop("cardiology apps only de-identified")
  cells <- plasmode_cells()
  cells <- cells[
    cells$app %in%
      apps &
      cells$flip == "curve" &
      cells$G > 40 &
      cells$error == "attached",
  ]
  sim_fun <- plasmode_simulate
}
dir.create(file.path(out, "locks"), recursive = TRUE, showWarnings = FALSE)
tasks <- expand.grid(cell = cells$cell, rep = reps)
tasks <- tasks[order(tasks$rep, tasks$cell), ]
for (i in seq_len(nrow(tasks))) {
  id <- tasks$cell[i]
  rep <- tasks$rep[i]
  f <- file.path(out, sprintf("c%03d-r%03d.rds", id, rep))
  if (file.exists(f)) next
  if (!dir.create(file.path(out, "locks", basename(f)), showWarnings = FALSE))
    next
  res <- tryCatch(
    combo_task(cells[cells$cell == id, ], rep, sim_fun),
    error = \(e) list(error = conditionMessage(e))
  )
  tmp <- tempfile(tmpdir = out)
  saveRDS(res, tmp)
  file.rename(tmp, f)
  cat(sprintf(
    "cell %d rep %d %s\n",
    id,
    rep,
    if (is.null(res$error)) sprintf("%.1fs", res$time) else res$error
  ))
}
