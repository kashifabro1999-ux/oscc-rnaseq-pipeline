library(pheatmap)
library(RColorBrewer)

expr <- read.delim(
  "data/vst_counts.tsv",
  row.names = 1,
  check.names = FALSE
)

hub <- read.delim(
  "results/top30_hub_genes.tsv",
  check.names = FALSE
)

sample_info <- read.delim(
  "data/sample_info.tsv",
  check.names = FALSE
)

common_genes <- intersect(
  rownames(expr),
  as.character(hub$EntrezID)
)

heatmap_matrix <- expr[common_genes, ]

symbol_map <- hub$GeneSymbol
names(symbol_map) <- as.character(hub$EntrezID)

rownames(heatmap_matrix) <- symbol_map[rownames(heatmap_matrix)]

heatmap_matrix <- heatmap_matrix[complete.cases(heatmap_matrix), ]

annotation <- data.frame(
  Condition = sample_info$condition
)

rownames(annotation) <- sample_info$sample

annotation <- annotation[colnames(heatmap_matrix), , drop = FALSE]

annotation_colors <- list(
  Condition = c(
    Tumor = "#F8766D",
    Control = "#00BFC4"
  )
)

png("plots/Hub_Genes_Heatmap.png",
    width = 2400,
    height = 3000,
    res = 300)

pheatmap(
  heatmap_matrix,
  scale = "row",
  color = colorRampPalette(
    c("navy", "black", "orange", "yellow")
  )(100),
  annotation_col = annotation,
  annotation_colors = annotation_colors,
  fontsize_row = 9,
  fontsize_col = 7,
  clustering_method = "ward.D2",
  border_color = NA,
  main = "Top 30 Hub Gene Expression Heatmap"
)

dev.off()

pdf("plots/Hub_Genes_Heatmap.pdf",
    width = 10,
    height = 12)

pheatmap(
  heatmap_matrix,
  scale = "row",
  color = colorRampPalette(
    c("navy", "black", "orange", "yellow")
  )(100),
  annotation_col = annotation,
  annotation_colors = annotation_colors,
  fontsize_row = 9,
  fontsize_col = 7,
  clustering_method = "ward.D2",
  border_color = NA,
  main = "Top 30 Hub Gene Expression Heatmap"
)

dev.off()

cat("Hub gene heatmap completed successfully.\n")
cat("Hub genes plotted:", nrow(heatmap_matrix), "\n")
