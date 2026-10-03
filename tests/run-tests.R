suppressPackageStartupMessages(library(refund))
for (file in sort(list.files("R", pattern = "\\.R$", full.names = TRUE))) {
  source(file)
}
testthat::test_dir("tests")
