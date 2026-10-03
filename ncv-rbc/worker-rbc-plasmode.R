# LRZ worker for the bias-corrected NCV run on the plasmode cells (all subjects, curve
# flips, attached errors). Claims rep-major tasks via lock dirs. Output stays on LRZ
# (task files hold values derived from the real covariates); only aggregates leave it.
#   Rscript ncv-rbc/worker-rbc-plasmode.R <out_dir> <app,app,...> <n_rep>
source("scripts/setup.R")
source("plasmode/sources.R")
source("plasmode/plasmode.R")
source("ncv-rbc/rbc-lib.R")
args <- commandArgs(trailingOnly = TRUE)
out <- args[1]
apps <- strsplit(args[2], ",")[[1]]
n_rep <- as.integer(args[3])
if (
  any(apps %in% PLASMODE_LOCAL_ONLY & !apps %in% PLASMODE_DEID_APPS) ||
    (any(apps %in% PLASMODE_DEID_APPS) && !identical(PLASMODE_CARDIO, "deid"))
)
  stop(
    "cardiology apps on a cluster only from the de-identified input (PLASMODE_CARDIO=deid)"
  )
cells <- plasmode_cells()
cells <- cells[
  cells$app %in%
    apps &
    cells$flip == "curve" &
    cells$G > 40 &
    cells$error == "attached",
]
dir.create(file.path(out, "locks"), recursive = TRUE, showWarnings = FALSE)
tasks <- expand.grid(cell = cells$cell, rep = seq_len(n_rep))
tasks <- tasks[order(tasks$rep, tasks$cell), ]
for (i in seq_len(nrow(tasks))) {
  id <- tasks$cell[i]
  rep <- tasks$rep[i]
  f <- file.path(out, sprintf("c%03d-r%03d.rds", id, rep))
  if (file.exists(f)) next
  if (!dir.create(file.path(out, "locks", basename(f)), showWarnings = FALSE))
    next
  res <- tryCatch(
    rbc_task(
      cells[cells$cell == id, ],
      rep,
      plasmode_simulate,
      with_reml = TRUE
    ),
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
