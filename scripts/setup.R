# Shared setup of the entry scripts (run from the repository root): optional
# study library first (LRZ), then every R/ file.
if (nzchar(Sys.getenv("RESTART_LIB"))) {
  .libPaths(c(Sys.getenv("RESTART_LIB"), .libPaths()))
}
for (file in sort(list.files("R", pattern = "\\.R$", full.names = TRUE))) {
  source(file)
}
