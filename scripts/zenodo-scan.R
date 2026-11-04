# Classify every task record of one results folder for the public deposit. Reads each .rds and walks
# all of it (values, names, attributes, environments) for dataset names and path-like strings.
# Writes one row per file; cardiology records must never be published. From the project directory:
#   Rscript zenodo-scan.R <results_dir> <out.csv> [strings.txt]
# strings.txt: distinct non-numeric strings seen in records that are NOT cardiology, for review.
CARDIO_APPS <- c("af", "ecg_gls4", "ecg_gls4b")
PUBLIC_APPS <- c(
  "running",
  "dti",
  "gait_vc",
  "ecg_roahd",
  "ocean",
  "aemet",
  "adelaide"
)
EXCLUDED_DIRS <- c("locks", "failed-first-attempt")
MAX_STRING_CHARS <- 200
MAX_DISTINCT_PER_FILE <- 2000
PATH_PATTERN <- "^(/|~/)|/home/|/dss/|fabians|ri68suz2|bouchahda"

collect_strings <- function(x) {
  found <- character()
  seen_envs <- list()
  walk <- function(obj) {
    if (is.environment(obj)) {
      if (any(vapply(seen_envs, identical, logical(1), obj))) {
        return(invisible())
      }
      if (
        identical(obj, globalenv()) ||
          identical(obj, emptyenv()) ||
          isNamespace(obj) ||
          identical(obj, baseenv())
      ) {
        return(invisible())
      }
      seen_envs[[length(seen_envs) + 1]] <<- obj
      walk(as.list(obj, all.names = TRUE))
    } else if (is.character(obj)) {
      found <<- c(found, unique(obj))
    } else if (is.factor(obj)) {
      found <<- c(found, levels(obj))
    } else if (is.list(obj) || is.expression(obj)) {
      for (el in as.list(obj)) walk(el)
    } else if (is.function(obj)) {
      walk(environment(obj))
    } else if (is.language(obj)) {
      found <<- c(found, paste(deparse(obj), collapse = " "))
    }
    att <- attributes(obj)
    if (!is.null(att)) {
      found <<- c(found, names(att))
      for (el in att) walk(el)
    }
    invisible()
  }
  walk(x)
  unique(found)
}

classify_values <- function(strings) {
  cardio <- intersect(strings, CARDIO_APPS)
  public <- intersect(strings, PUBLIC_APPS)
  plasmode <- any(strings == "plasmode")
  status <- if (length(cardio) > 0) {
    "cardio"
  } else if (length(public) > 0) {
    "public"
  } else if (plasmode) {
    "unidentified_plasmode"
  } else {
    "no_app_found"
  }
  list(status = status, cardio = cardio, public = public)
}

scan_file <- function(file, results_dir) {
  rel <- substring(file, nchar(results_dir) + 2)
  bytes <- file.size(file)
  obj <- tryCatch(readRDS(file), error = \(e) e)
  if (inherits(obj, "error")) {
    return(list(
      row = data.frame(
        file = rel,
        bytes = bytes,
        status = "unreadable",
        apps = "",
        n_paths = 0L
      ),
      strings = character()
    ))
  }
  strings <- collect_strings(obj)
  cl <- classify_values(strings)
  is_cardio <- cl$status == "cardio"
  paths <- grep(PATH_PATTERN, strings, value = TRUE, ignore.case = TRUE)
  review <- character()
  # Strings of cardiology records stay out of every output: only their status is reported.
  if (!is_cardio) {
    review <- strings[nchar(strings) <= MAX_STRING_CHARS]
    review <- utils::head(review, MAX_DISTINCT_PER_FILE)
  }
  list(
    row = data.frame(
      file = rel,
      bytes = bytes,
      status = cl$status,
      apps = paste(sort(c(cl$cardio, cl$public)), collapse = ";"),
      n_paths = if (is_cardio) NA_integer_ else length(paths)
    ),
    strings = review
  )
}

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 2)
results_dir <- sub("/+$", "", args[1])
out_csv <- args[2]
strings_txt <- if (length(args) >= 3) args[3] else NA_character_
cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))

files <- list.files(results_dir, "\\.rds$", recursive = TRUE, full.names = TRUE)
top <- vapply(
  strsplit(substring(files, nchar(results_dir) + 2), "/"),
  `[`,
  character(1),
  1
)
files <- files[!top %in% EXCLUDED_DIRS]
stopifnot(length(files) > 0)
out <- parallel::mclapply(
  files,
  scan_file,
  results_dir = results_dir,
  mc.cores = cores,
  mc.preschedule = TRUE
)
failed <- !vapply(out, is.list, logical(1))
stopifnot("scan workers failed" = !any(failed))
rows <- do.call(rbind, lapply(out, `[[`, "row"))
write.csv(rows, out_csv, row.names = FALSE)
print(table(rows$status))
print(table(rows$apps))
cat("files with path-like strings:", sum(rows$n_paths > 0, na.rm = TRUE), "\n")
if (!is.na(strings_txt)) {
  writeLines(sort(unique(unlist(lapply(out, `[[`, "strings")))), strings_txt)
}
