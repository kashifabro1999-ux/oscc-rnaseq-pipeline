## install_bioc_packages.R
## Installs every Bioconductor package used across the 20 pipeline scripts.
## Kept as its own file (separate Docker layer from the CRAN install) so a
## rebuild triggered by an app-code change doesn't force reinstalling the
## entire Bioconductor stack from scratch.

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

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
cat("All Bioconductor packages installed successfully.\n")
