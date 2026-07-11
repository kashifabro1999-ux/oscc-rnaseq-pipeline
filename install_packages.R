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

install_bioc <- function(pkgs, force = FALSE) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", lib = lib)
  }

  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0 || force) {
    target <- if (force) pkgs else missing
    cat("Installing Bioconductor packages:", paste(target, collapse = ", "), "\n")
    BiocManager::install(target, lib = lib, ask = FALSE, update = FALSE, force = force)
  } else {
    cat("All Bioconductor packages already installed.\n")
  }
}

install_github_pkg <- function(repo) {
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes", lib = lib)
  }

  cat("Installing GitHub package:", repo, "\n")
  remotes::install_github(repo, lib = lib, upgrade = "never", force = TRUE)
}

verify_packages <- function(pkgs) {
  pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
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

bioc_basic <- c(
  "EnhancedVolcano",
  "org.Hs.eg.db",
  "AnnotationDbi",
  "pathview",
  "STRINGdb",
  "DOSE",
  "fgsea",
  "GO.db",
  "GOSemSim",
  "reactome.db"
)

bioc_enrichment <- c(
  "treeio",
  "ggtree",
  "enrichplot",
  "clusterProfiler",
  "ReactomePA"
)

critical_pkgs <- c(
  cran_pkgs,
  bioc_basic,
  bioc_enrichment
)

install_cran(cran_pkgs)
install_bioc(bioc_basic)

cat("\nInstalling enrichment dependency chain...\n")
try(install_bioc(bioc_enrichment), silent = FALSE)

missing_after_bioc <- verify_packages(bioc_enrichment)

if (length(missing_after_bioc) > 0) {
  cat("\nBioconductor install left missing packages:", paste(missing_after_bioc, collapse = ", "), "\n")
  cat("Trying YuLab GitHub fallback for treeio/ggtree dependency chain...\n")

  # This fixes common treeio/ggtree/yulab.utils compatibility problems.
  github_repos <- c(
    "YuLab-SMU/yulab.utils",
    "YuLab-SMU/tidytree",
    "YuLab-SMU/treeio",
    "YuLab-SMU/ggfun",
    "YuLab-SMU/ggtree"
  )

  for (repo in github_repos) {
    try(install_github_pkg(repo), silent = FALSE)
  }

  try(install_bioc(c("enrichplot", "clusterProfiler", "ReactomePA"), force = TRUE), silent = FALSE)
}

missing_final <- verify_packages(critical_pkgs)

cat("\nFinal package check:\n")
for (pkg in critical_pkgs) {
  ok <- requireNamespace(pkg, quietly = TRUE)
  cat(sprintf("  %-20s %s\n", pkg, ifelse(ok, "OK", "MISSING")))
}

if (length(missing_final) > 0) {
  cat("\n[ERROR] Missing required R packages after installation:\n")
  cat(paste(" -", missing_final), sep = "\n")
  cat("\nPlease check internet connection and rerun:\n")
  cat("Rscript install_packages.R\n")
  quit(status = 1)
}

cat("\n[OK] All required R packages are installed.\n")
