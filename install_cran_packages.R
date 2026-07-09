## install_cran_packages.R
## Installs every CRAN package used across the 20 pipeline scripts plus the
## CLI/web tools. Kept as its own file (rather than an inline Dockerfile
## RUN command) so Docker can cache this layer independently.

pkgs <- c(
  "ggplot2",
  "pheatmap",
  "RColorBrewer",
  "dplyr",
  "tidyr",
  "igraph",
  "ggraph",
  "msigdbr"
)

install.packages(pkgs, repos = "https://cloud.r-project.org", Ncpus = 2)

# Fail the Docker build loudly if anything didn't actually install, rather
# than silently continuing with a broken image.
missing <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(missing) > 0) {
  stop("Failed to install CRAN package(s): ", paste(missing, collapse = ", "))
}
cat("All CRAN packages installed successfully.\n")
