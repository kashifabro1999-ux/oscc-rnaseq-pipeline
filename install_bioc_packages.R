## install_bioc_packages.R
options(repos = c(CRAN = "https://packagemanager.posit.co/cran/2023-11-15"))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

dir.create("~/.R", showWarnings = FALSE)
writeLines("CXX11STD = -std=gnu++14", con = "~/.R/Makevars")

pkgs <- c(
  "DESeq2",
  "EnhancedVolcano",
  "AnnotationDbi",
  "org.Hs.eg.db",
  "clusterProfiler",
  "enrichplot",
  "pathview",
  "STRINGdb",
  "ReactomePA",
  "reactome.db"
)

BiocManager::install(pkgs, update = FALSE, ask = FALSE)

missing <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(missing) > 0) {
  stop("Failed to install Bioconductor package(s): ", paste(missing, collapse = ", "))
}
cat("All Bioconductor packages installed successfully (pinned to 2023-11-15 CRAN snapshot).\n")
