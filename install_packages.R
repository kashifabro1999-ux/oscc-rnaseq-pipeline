cat("Checking/installing required R packages...\n")

lib <- Sys.getenv("R_LIBS_USER")
if (lib == "") {
  lib <- file.path(Sys.getenv("HOME"), "R", "library")
}
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib, .libPaths()))

options(repos = c(CRAN = "https://cloud.r-project.org"))

install_cran <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    cat("Installing CRAN packages:", paste(missing, collapse = ", "), "\n")
    install.packages(missing, lib = lib)
  } else {
    cat("All CRAN packages already installed.\n")
  }
}

install_bioc <- function(pkgs) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", lib = lib)
  }

  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    cat("Installing Bioconductor packages:", paste(missing, collapse = ", "), "\n")
    BiocManager::install(missing, lib = lib, ask = FALSE, update = FALSE)
  } else {
    cat("All Bioconductor packages already installed.\n")
  }
}

cran_pkgs <- c(
  "ggplot2",
  "ggrepel",
  "pheatmap",
  "dplyr",
  "tidyr",
  "readr",
  "stringr",
  "data.table",
  "reshape2",
  "igraph",
  "RColorBrewer",
  "viridis",
  "msigdbr",
  "ggridges",
  "remotes"
)

bioc_pkgs <- c(
  "EnhancedVolcano",
  "org.Hs.eg.db",
  "AnnotationDbi",
  "clusterProfiler",
  "enrichplot",
  "ReactomePA",
  "pathview",
  "STRINGdb",
  "DOSE",
  "fgsea"
)

install_cran(cran_pkgs)
install_bioc(bioc_pkgs)

cat("R package check completed.\n")
