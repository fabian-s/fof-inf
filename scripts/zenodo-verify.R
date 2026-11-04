# Check one packed archive: file list equals the allowlist, no dropped record is in it, and a random
# sample of records reads back identical to the originals. Run on a compute node:
#   Rscript zenodo-verify.R <results_dir> <keep.txt> <drop.csv> <archive.tar.xz> [n_sample]
N_SAMPLE_DEFAULT <- 20L
SEED <- 20261103L
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 4)
results_dir <- sub("/+$", "", args[1])
keep <- readLines(args[2])
drop <- read.csv(args[3])$file
archive <- args[4]
n_sample <- if (length(args) >= 5) as.integer(args[5]) else N_SAMPLE_DEFAULT
name <- basename(results_dir)
listed <- system2("tar", c("-tJf", shQuote(archive)), stdout = TRUE)
listed <- listed[!grepl("/$", listed)]
inside <- sub(paste0("^", name, "/"), "", listed)
stopifnot(
  "archive differs from allowlist" = setequal(inside, keep),
  "dropped record in archive" = !any(drop %in% inside)
)
set.seed(SEED)
sampled <- sample(keep, min(n_sample, length(keep)))
tmp <- tempfile()
dir.create(tmp)
system2(
  "tar",
  c(
    "-xJf",
    shQuote(archive),
    "-C",
    shQuote(tmp),
    shQuote(file.path(name, sampled))
  )
)
same <- vapply(
  sampled,
  \(f)
    identical(
      readRDS(file.path(results_dir, f)),
      readRDS(file.path(tmp, name, f))
    ),
  logical(1)
)
stopifnot("record changed by packing" = all(same))
unlink(tmp, recursive = TRUE)
orig <- sum(file.size(file.path(results_dir, keep)))
cat(sprintf(
  "%s OK: %d records, %d dropped, %d/%d sampled identical; %.2f GB -> %.2f GB (ratio %.3f)\n",
  name,
  length(keep),
  length(drop),
  sum(same),
  length(same),
  orig / 1e9,
  file.size(archive) / 1e9,
  file.size(archive) / orig
))
