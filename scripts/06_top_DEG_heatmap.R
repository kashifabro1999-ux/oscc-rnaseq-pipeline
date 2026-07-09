library(pheatmap)
library(RColorBrewer)

# Read files
vst <- read.delim(
  "data/vst_counts.tsv",
  row.names = 1,
  check.names = FALSE
)

res <- read.delim(
  "results/deseq2_results_annotated.tsv",
  check.names = FALSE
)

sample_info <- read.delim(
  "data/sample_info.tsv",
  check.names = FALSE
)

# Significant genes
res <- res[!is.na(res$padj), ]
sig <- res[res$padj < 0.05 & abs(res$log2FoldChange) >= 2, ]

# Top 50 genes
top50 <- head(sig[order(sig$padj), ], 50)

# Match genes
common_genes <- intersect(
  rownames(vst),
  as.character(top50$EntrezID)
)

heatmap_matrix <- vst[common_genes, ]

# Replace IDs with symbols
gene_symbols <- top50$GeneSymbol[
  match(common_genes, top50$EntrezID)
]

rownames(heatmap_matrix) <- gene_symbols

# Z-score normalization
heatmap_matrix <- t(scale(t(heatmap_matrix)))
heatmap_matrix <- heatmap_matrix[complete.cases(heatmap_matrix), ]

# -------------------------------------------------
# ORDER SAMPLES MANUALLY
# -------------------------------------------------

control_order <- paste0("N", 1:10)
tumor_order <- paste0("T", 1:50)

desired_order <- c(control_order, tumor_order)

desired_order <- desired_order[
  desired_order %in% colnames(heatmap_matrix)
]

heatmap_matrix <- heatmap_matrix[, desired_order]

# Sample annotation
annotation <- data.frame(
  condition = ifelse(
    grepl("^N", desired_order),
    "Control",
    "Tumor"
  )
)

rownames(annotation) <- desired_order

annotation_colors <- list(
  condition = c(
    Control = "#00BFC4",
    Tumor = "#F8766D"
  )
)

# Professional color palette
heat_colors <- colorRampPalette(c(
  "#2166AC",
  "#67A9CF",
  "white",
  "#F4A582",
  "#B2182B"
))(255)

# PNG
png(
  "plots/top50_DEG_heatmap_ordered.png",
  width = 5200,
  height = 4200,
  res = 350
)

pheatmap(
  heatmap_matrix,
  annotation_col = annotation,
  annotation_colors = annotation_colors,
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 10,
  fontsize_col = 8,
  angle_col = 45,
  border_color = NA,
  scale = "none",
  color = heat_colors,
  main = "Top 50 Differentially Expressed Genes"
)

dev.off()

# PDF
pdf(
  "plots/top50_DEG_heatmap_ordered.pdf",
  width = 16,
  height = 13
)

pheatmap(
  heatmap_matrix,
  annotation_col = annotation,
  annotation_colors = annotation_colors,
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 10,
  fontsize_col = 8,
  angle_col = 45,
  border_color = NA,
  scale = "none",
  color = heat_colors,
  main = "Top 50 Differentially Expressed Genes"
)

dev.off()
