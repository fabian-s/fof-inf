source(file.path("..", "analysis", "contrasts.R"))

synthetic_tasks <- function(n_rep = 20) {
  cells <- design_table()
  set.seed(3)
  grid <- tidyr::expand_grid(
    cell = cells$cell,
    rep = seq_len(n_rep),
    estimand = c("mean", "beta"),
    arm = c("reml_model", "reml_cl2", "ncv_model", "ncv_cl2", "ncv_cl2_bias")
  )
  target <- c(
    reml_model = 0.8,
    reml_cl2 = 0.93,
    ncv_model = 0.5,
    ncv_cl2 = 0.9,
    ncv_cl2_bias = 0.95
  )
  dplyr::mutate(
    grid,
    task = (rep - 1) * nrow(cells) + cell,
    coverage = pmin(1, target[arm] + rnorm(dplyr::n(), 0, 0.03)),
    width = 1,
    score = ifelse(arm == "ncv_cl2_bias", 0.9, 1) + rnorm(dplyr::n(), 0, 0.01),
    sq_error = ifelse(grepl("ncv", arm), 0.8, 1) * rexp(dplyr::n()),
    deviance = NA_real_
  )
}

test_that("the recommendation with uncertainty reproduces the point decision", {
  cells <- design_table()
  tasks <- synthetic_tasks()
  summaries <- tasks |>
    dplyr::group_by(cell, estimand, arm) |>
    dplyr::summarise(
      coverage_se = sd(coverage) / sqrt(dplyr::n()),
      coverage_mean = mean(coverage),
      score_mean = mean(score),
      .groups = "drop"
    )
  out <- recommend_with_uncertainty(tasks, summaries, cells, "gaussian", "beta")
  expect_identical(out$decision$decision, "proposal")
  ratio <- out$criteria[out$criteria$name == "score_ratio", ]
  expect_equal(ratio$estimate, out$decision$score_ratio, tolerance = 1e-8)
  expect_true(ratio$lower <= ratio$estimate && ratio$estimate <= ratio$upper)
})

test_that("question 1, profiles and matched contrasts run and are ordered", {
  cells <- design_table()
  tasks <- synthetic_tasks()
  q1 <- question1(tasks, cells, "poisson", "mean")
  expect_lt(
    q1$estimate[q1$name == "ce_reml_cl2"],
    q1$estimate[q1$name == "ce_reml_model"]
  )
  prof <- factor_profile(
    tasks,
    cells,
    "binomial",
    "beta",
    "error",
    arm_value("coverage", "reml_model")
  )
  expect_setequal(prof$name, c("iid", "ou", "fpc"))
  expect_equal(prof$estimate, rep(0.8, 3), tolerance = 0.02)
  dense <- cells$cell[cells$block == "dense_grid"][1]
  mc <- matched_contrast(
    tasks,
    dense,
    matching_core_cell(cells, dense),
    "beta",
    log_ratio("sq_error", "ncv_cl2", "reml_cl2")
  )
  expect_equal(nrow(mc), 1)
})
