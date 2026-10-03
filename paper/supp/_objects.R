# Objects for paper/supplement.qmd: run the frozen code of the detailed results write-up
# (paper/supp/src/*.R, purled snapshot) with file output disabled and every figure kept in
# FIGS[[name]]. Sourced after paper/_common.R, the repository root as working directory.
FIGS <- list()
save_fig <- function(plot, name, ...) {
  FIGS[[name]] <<- plot
  invisible(plot)
}
save_tab <- function(df, ...) invisible(df)
show <- function(...) invisible(NULL)
for (src_name in c(
  "RESULTS",
  "_relative-error",
  "_relative-error-plasmode",
  "_revision-synthetic",
  "_revision-misspec",
  "_comparators",
  "_ncv-intervals"
)) {
  sys.source(
    file.path("paper/supp/src", paste0(src_name, ".R")),
    envir = globalenv()
  )
}
rm(src_name)
