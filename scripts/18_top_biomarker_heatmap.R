suppressPackageStartupMessages({
  library(pheatmap)
  library(dplyr)
})

dir.create("plots", showWarnings = FALSE, recursive = TRUE)
dir.create("results", showWarnings = FALSE, recursive = TRUE)

vst <- read.delim("data/vst_counts.tsv", check.names = FALSE)
priority <- read.delim("results/biomarker_priority_table.tsv", check.names = FALSE)
sample_info <- read.delim("data/sample_info.tsv", check.names = FALSE)

id_col <- colnames(vst)[1]
vst[[id_col]] <- as.character(vst[[id_col]])
priority$EntrezID <- as.character(priority$EntrezID)

target_genes <- c(
  "PLAU", "MMP14", "COL1A1", "MMP9",
  "LPL", "ACACB", "SLC2A4", "PPARG",
  "GPD1", "PLIN1", "PCK1", "ADH1B", "CYP3A4"
)

priority_sel <- priority %>%
  filter(GeneSymbol %in% target_genes)

if (nrow(priority_sel) < 2) {
  priority_sel <- priority %>%
    arrange(desc(PriorityScore)) %>%
    head(30)
}

match_idx <- match(priority_sel$EntrezID, vst[[id_col]])
valid <- !is.na(match_idx)

priority_sel <- priority_sel[valid, , drop = FALSE]
vst_sel <- vst[match_idx[valid], , drop = FALSE]

if (nrow(vst_sel) < 2) {
  stop("Fewer than 2 biomarker genes found in VST after EntrezID matching.")
}

expr <- as.matrix(vst_sel[, -1])
mode(expr) <- "numeric"
rownames(expr) <- make.unique(priority_sel$GeneSymbol)

sample_col <- colnames(sample_info)[1]
possible_sample_cols <- colnames(sample_info)[tolower(colnames(sample_info)) %in% c("sample", "sample_id", "sampleid", "id")]
if (length(possible_sample_cols) > 0) {
  sample_col <- possible_sample_cols[1]
}

common_samples <- intersect(as.character(sample_info[[sample_col]]), colnames(expr))

annotation_col <- NULL

if (length(common_samples) >= 2) {
  expr <- expr[, common_samples, drop = FALSE]
  sample_info <- sample_info[match(common_samples, sample_info[[sample_col]]), , drop = FALSE]

  condition_cols <- colnames(sample_info)[tolower(colnames(sample_info)) %in% c("condition", "group", "diagnosis", "type")]
  if (length(condition_cols) > 0) {
    annotation_col <- data.frame(Condition = as.factor(sample_info[[condition_cols[1]]]))
    rownames(annotation_col) <- common_samples
  }
}

expr <- expr[apply(expr, 1, sd, na.rm = TRUE) > 0, , drop = FALSE]

if (nrow(expr) < 2) {
  stop("Fewer than 2 variable biomarker genes available for clustering.")
}

expr_z <- t(scale(t(expr)))
expr_z[is.na(expr_z)] <- 0

pheatmap(
  expr_z,
  annotation_col = annotation_col,
  show_colnames = FALSE,
  fontsize_row = 9,
  main = "Prioritized OSCC Candidate Biomarkers",
  filename = "plots/top_biomarker_heatmap.png",
  width = 10,
  height = 7
)

pdf("plots/top_biomarker_heatmap.pdf", width = 10, height = 7)
pheatmap(
  expr_z,
  annotation_col = annotation_col,
  show_colnames = FALSE,
  fontsize_row = 9,
  main = "Prioritized OSCC Candidate Biomarkers"
)
dev.off()

write.table(
  priority_sel,
  "results/top_biomarker_heatmap_genes.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("Top biomarker heatmap completed successfully.\n")
cat("Genes plotted:", nrow(expr_z), "\n")
