# Point records ----------------------------------------------------------------

score_task <- function(inputs, sim, cell, task) {
  if (!all(c("REML", "NCV") %in% names(inputs))) {
    stop("Scoring requires REML and NCV inputs.")
  }
  estimands <- model_estimands(cell$model)
  points <- lapply(estimands, \(estimand) {
    truth <- sim$truth[[estimand]]
    n <- nrow(eval_points(estimand, sim$eval_t %||% EVAL_T))
    if (!is.numeric(truth) || length(truth) != n) {
      stop("Truth does not match evaluation points for ", estimand)
    }
    columns <- list(truth = as.vector(truth))
    # every estimate and SE a fit provides (est, se_model, se_cl2, se_cl2f, and the
    # block-clustered se_cl2b, se_cl2bf, se_cl2a, se_cl2af where computed)
    for (method in names(inputs)) {
      value <- inputs[[method]]$estimands[[estimand]]
      columns[paste0(names(value), "_", method)] <- value
    }
    columns <- Filter(Negate(is.null), columns)
    if (is.null(columns$se_cl2_REML) || is.null(columns$se_cl2_NCV)) {
      stop("CL2 standard errors missing for ", estimand)
    }
    if (any(lengths(columns) != n)) stop("Input length mismatch for ", estimand)
    main <- unlist(columns[grepl("_(REML|NCV)$", names(columns))])
    if (!all(is.finite(main)))
      stop("Non-finite REML/NCV estimate or SE for ", estimand)
    do.call(cbind, columns)
  }) |>
    setNames(estimands)
  bind <- \(part) do.call(rbind, lapply(inputs, `[[`, part))
  list(
    task = task$task,
    cell = cell$cell,
    rep = task$rep,
    records = points,
    fits = bind("fit_row"),
    terms = bind("terms"),
    warnings = lapply(inputs, `[[`, "warnings"),
    diag = sim$diag
  )
}
