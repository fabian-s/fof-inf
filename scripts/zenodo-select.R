# Allowlist for the public deposit: which task records of one results folder may be published.
# Records often do not name their dataset, so plasmode records are identified by cell id through
# plasmode_cells() (cell-app.csv). Cardiology records are never published; a record whose dataset
# cannot be established is held back. Only folders behind the paper and supplement have a rule; any
# other folder stops with an error. From the project directory:
#   Rscript zenodo-select.R <results_dir> <cell-app.csv> <keep.txt> <drop.csv>
CARDIO_APPS <- c("af", "ecg_gls4", "ecg_gls4b")
# 3a93a121c147: the pilot run of the main study (superseded by 373f6fd65cac).
EXCLUDED_DIRS <- c("locks", "failed-first-attempt", "3a93a121c147")
# Folders of plasmode runs whose file names are c<cell>; "c" names in other folders are synthetic.
PLASMODE_C_DIRS <- c(
  "rbc-pm2-results",
  "combo-pm-results",
  "combo2-pm-results",
  "fixlam-pm-results",
  "jackfix-pm-results",
  "lc-plasmode-results",
  "lc-ncvp-results",
  "satt-r4-pm",
  "plasmode-results"
)
SYNTHETIC_DIRS <- c(
  "results",
  "r4-results",
  "ablation-results",
  "satt-results",
  "satt-r4-syn",
  "misspec-results",
  "lc-results",
  "lc-ncv-results",
  "lc-ncvg2-results",
  "rbc-fresh-results",
  "combo-syn-results",
  "combo2-syn-results",
  "combo2-ms-results",
  "jackfix-results",
  "ffpc-results"
)
MIXED_DIRS <- c("comparators-results", "cl2arms-results") # task-p<cell> plasmode, task-s<cell> synthetic
FGEE_DIRS <- "fgee-results" # A, B synthetic; P by cell id; Pc cardiology

plasmode_status <- function(cell, cell_app) {
  app <- cell_app$app[match(cell, cell_app$cell)]
  ifelse(
    is.na(app),
    "unknown_cell",
    ifelse(app %in% CARDIO_APPS, "cardio", "public")
  )
}

classify <- function(dir_name, rel, cell_app) {
  base <- basename(rel)
  sub <- ifelse(grepl("/", rel), sub("/.*", "", rel), "")
  id <- suppressWarnings(as.integer(sub(
    "^.*[cp]([0-9]+)(-r[0-9]+)?\\.rds$",
    "\\1",
    base
  )))
  if (dir_name %in% SYNTHETIC_DIRS) {
    return(rep("synthetic", length(rel)))
  }
  if (dir_name %in% PLASMODE_C_DIRS) {
    return(plasmode_status(id, cell_app))
  }
  if (dir_name %in% MIXED_DIRS) {
    return(ifelse(
      grepl("^task-s[0-9]+-r", base),
      "synthetic",
      plasmode_status(id, cell_app)
    ))
  }
  if (dir_name %in% FGEE_DIRS) {
    return(ifelse(
      sub %in% c("A", "B"),
      "synthetic",
      ifelse(sub == "Pc", "cardio", plasmode_status(id, cell_app))
    ))
  }
  stop("folder has no rule: ", dir_name)
}

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 4)
results_dir <- sub("/+$", "", args[1])
cell_app <- read.csv(args[2])
rel <- list.files(results_dir, "\\.rds$", recursive = TRUE)
rel <- rel[!sub("/.*", "", rel) %in% EXCLUDED_DIRS]
stopifnot(length(rel) > 0)
status <- classify(basename(results_dir), rel, cell_app)
keep <- status %in% c("synthetic", "public")
writeLines(rel[keep], args[3])
write.csv(
  data.frame(file = rel[!keep], status = status[!keep]),
  args[4],
  row.names = FALSE
)
print(table(status))
