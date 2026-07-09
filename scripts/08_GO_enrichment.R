library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

# Read annotated DEGs
res <- read.delim(
  "results/deseq2_results_annotated.tsv",
  check.names = FALSE
)

# Remove NA values
res <- res[!is.na(res$padj), ]

# Significant genes
sig <- res[
  res$padj < 0.05 &
  abs(res$log2FoldChange) >= 1,
]

# Separate upregulated and downregulated
up <- sig[sig$log2FoldChange > 1, ]
down <- sig[sig$log2FoldChange < -1, ]

# Entrez IDs
up_genes <- as.character(up$EntrezID)
down_genes <- as.character(down$EntrezID)

# --------------------------------------------------
# GO enrichment for UPREGULATED genes
# --------------------------------------------------

ego_up <- enrichGO(
  gene = up_genes,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

# Save table
write.table(
  as.data.frame(ego_up),
  file = "results/GO_upregulated.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Dotplot PNG
png(
  "plots/GO_upregulated_dotplot.png",
  width = 3200,
  height = 2400,
  res = 300
)

dotplot(
  ego_up,
  showCategory = 15,
  font.size = 12,
  title = "GO Enrichment - Upregulated Genes"
)

dev.off()

# --------------------------------------------------
# GO enrichment for DOWNREGULATED genes
# --------------------------------------------------

ego_down <- enrichGO(
  gene = down_genes,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

# Save table
write.table(
  as.data.frame(ego_down),
  file = "results/GO_downregulated.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Dotplot PNG
png(
  "plots/GO_downregulated_dotplot.png",
  width = 3200,
  height = 2400,
  res = 300
)

dotplot(
  ego_down,
  showCategory = 15,
  font.size = 12,
  title = "GO Enrichment - Downregulated Genes"
)

dev.off()

# --------------------------------------------------
# PDF outputs
# --------------------------------------------------

pdf(
  "plots/GO_upregulated_dotplot.pdf",
  width = 11,
  height = 8
)

dotplot(
  ego_up,
  showCategory = 15,
  font.size = 12,
  title = "GO Enrichment - Upregulated Genes"
)

dev.off()

pdf(
  "plots/GO_downregulated_dotplot.pdf",
  width = 11,
  height = 8
)

dotplot(
  ego_down,
  showCategory = 15,
  font.size = 12,
  title = "GO Enrichment - Downregulated Genes"
)

dev.off()

# --------------------------------------------------
# Summary
# --------------------------------------------------

cat("Upregulated GO terms:", nrow(as.data.frame(ego_up)), "\n")
cat("Downregulated GO terms:", nrow(as.data.frame(ego_down)), "\n")
cat("GO enrichment analysis completed.\n")
