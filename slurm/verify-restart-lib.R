# Verify the LRZ restart library in a fresh R session.
lib <- path.expand(Sys.getenv("RESTART_LIB", "~/R/restart-lib"))
.libPaths(c(lib, .libPaths()))
stopifnot(packageVersion("mgcv") == "1.9.5")
suppressPackageStartupMessages(library(refund))
stopifnot(isTRUE(refund:::pffr_ncv_check_blocks()))
stopifnot(requireNamespace("fda", quietly = TRUE)) # pulls in fds/RCurl
cat(R.version.string, "\n")
cat("mgcv", format(packageVersion("mgcv")), find.package("mgcv"), "\n")
cat("refund", format(packageVersion("refund")), find.package("refund"), "\n")
cat("BLAS", extSoftVersion()[["BLAS"]], "\n")
cat("sha", readLines(file.path(lib, "REFUND_SHA")), "\n")
cat("OK\n")
