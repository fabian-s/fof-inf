# Run from the repository root. Arguments: chunk_size run out_dir.
source("scripts/setup.R")
args <- commandArgs(trailingOnly = TRUE)
chunk_size <- as.numeric(
  if (length(args) >= 1) args[1] else Sys.getenv("RESTART_CHUNK_SIZE", "1")
)
run <- as.numeric(
  if (length(args) >= 2) args[2] else Sys.getenv("RESTART_RUN", "1")
)
out_dir <- if (length(args) >= 3) args[3] else
  Sys.getenv("RESTART_OUT_DIR", "results")
array_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
cells <- design_table()
tasks <- task_table(cells)
# Optional subset of cells, e.g. RESTART_CELLS="96:99" for an extension run; chunks are
# formed within the subset (ids stay those of the full task table).
if (nzchar(Sys.getenv("RESTART_CELLS"))) {
  subset <- eval(parse(text = paste0("c(", Sys.getenv("RESTART_CELLS"), ")")))
  tasks <- tasks[tasks$cell %in% subset, ]
}
ids <- array_task_ids(tasks, array_id, chunk_size, run)
for (id in ids) run_task(id, out_dir, cells, tasks)
