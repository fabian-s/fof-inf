# The §5 rule under two calibration criteria, side by side: undercoverage only
# (calibration_error()) and symmetric (over- and undercoverage both count). Applied to the synthetic recommendation pool and to
# the plasmode cells. Point estimates only.
# Run from the repository root: Rscript analysis/rule-variants.R
source("scripts/setup.R")
suppressPackageStartupMessages(library(dplyr))
MAIN_FAMILIES <- c("gaussian", "poisson", "binomial") # as in run-analysis.R
PRIMARY <- c("mean", "beta")

# Symmetric calibration error: every cell's distance from
# nominal counts; MC variance is removed as in calibration_error().
symmetric_error <- function(coverage, coverage_se) {
  sqrt(max(0, mean((coverage - LEVEL)^2 - coverage_se^2)))
}

# The §5 decision given per-cell summaries of the proposal and the fallback (same cells).
decide <- function(prop, fall, criterion) {
  ce <- c(
    proposal = criterion(prop$coverage_mean, prop$coverage_se),
    fallback = criterion(fall$coverage_mean, fall$coverage_se)
  )
  ratio <- exp(mean(log(prop$score_mean / fall$score_mean)))
  adequate <- ce <= ADEQUATE
  decision <- if (sum(adequate) == 1) {
    names(ce)[adequate]
  } else if (!any(adequate)) {
    paste("neither adequate; less bad:", names(which.min(ce)))
  } else if (ratio <= PRACTICAL_RATIO) {
    "proposal"
  } else if (ratio >= 1 / PRACTICAL_RATIO) {
    "fallback"
  } else {
    "fallback (equivalent)"
  }
  tibble(
    ce_proposal = ce[["proposal"]],
    ce_fallback = ce[["fallback"]],
    score_ratio = ratio,
    decision = decision
  )
}

both_criteria <- function(rows) {
  prop <- arrange(filter(rows, arm == PROPOSAL), unit)
  fall <- arrange(filter(rows, arm == FALLBACK), unit)
  stopifnot(identical(prop$unit, fall$unit))
  bind_rows(
    mutate(decide(prop, fall, symmetric_error), criterion = "symmetric"),
    mutate(decide(prop, fall, calibration_error), criterion = "under only")
  ) |>
    mutate(n_cells = nrow(prop), .before = 1)
}

# Synthetic: the §5 pool per family and estimand.
cells <- design_table()
synthetic <- readr::read_csv(
  "summaries/final/cells.csv",
  show_col_types = FALSE
)
synthetic_rules <- tidyr::expand_grid(
  family = MAIN_FAMILIES,
  estimand = PRIMARY
) |>
  purrr::pmap(\(family, estimand) {
    pool <- recommendation_pool(cells, family)
    synthetic |>
      filter(cell %in% pool, estimand == !!estimand) |>
      mutate(unit = cell) |>
      both_criteria() |>
      mutate(
        data = "synthetic",
        group = family,
        estimand = estimand,
        .before = 1
      )
  }) |>
  bind_rows()

# Plasmode: all cells of one application per estimand.
plasmode <- readr::read_csv(
  "summaries/plasmode/coverage.csv",
  show_col_types = FALSE
)
plasmode_rules <- plasmode |>
  mutate(unit = paste(truth, error, G)) |>
  group_by(app, estimand) |>
  group_modify(\(rows, key) both_criteria(rows)) |>
  ungroup() |>
  mutate(data = "plasmode", group = app, .before = 1) |>
  select(-app)

out <- bind_rows(synthetic_rules, plasmode_rules) |>
  relocate(criterion, .after = estimand)
readr::write_csv(out, "summaries/final/rule-variants.csv")
print(as.data.frame(out), digits = 3)
