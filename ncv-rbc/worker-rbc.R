# LRZ worker for the bias-corrected NCV run: walks the rep-major task table
# (cells x reps) and claims tasks through atomic lock directories (<out>/locks).
#   Rscript ncv-rbc/worker-rbc.R <out_dir> <cell,cell,...> <rep_from> <rep_to> [with_reml]
source("scripts/setup.R")
source("ncv-rbc/rbc-lib.R")
args <- commandArgs(trailingOnly = TRUE)
out <- args[1]
cell_ids <- as.integer(strsplit(args[2], ",")[[1]])
reps <- seq(as.integer(args[3]), as.integer(args[4]))
with_reml <- identical(args[5], "with_reml")
cells <- design_table()
dir.create(file.path(out, "locks"), recursive = TRUE, showWarnings = FALSE)
tasks <- expand.grid(cell = cell_ids, rep = reps)[, c("cell", "rep")]
tasks <- tasks[order(tasks$rep, tasks$cell), ]
for (i in seq_len(nrow(tasks))) {
  id <- tasks$cell[i]
  rep <- tasks$rep[i]
  f <- file.path(out, sprintf("c%03d-r%03d.rds", id, rep))
  if (file.exists(f)) next
  if (!dir.create(file.path(out, "locks", basename(f)), showWarnings = FALSE))
    next
  res <- tryCatch(
    rbc_task(cells[cells$cell == id, ], rep, with_reml = with_reml),
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
