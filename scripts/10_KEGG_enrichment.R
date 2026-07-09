library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(org.Hs.eg.db)

# Read annotated DEGs
res <- read.delim(
  "results/deseq2_results_annotated.tsv",
  check.names = FALSE
)

# Remove NA
res <- res[!is.na(res$padj), ]

# Significant DEGs
sig <- res[
  res$padj < 0.05 &
  abs(res$log2FoldChange) >= 1,
]

# Separate up/down genes
up <- sig[sig$log2FoldChange > 1, ]
down <- sig[sig$log2FoldChange < -1, ]

# Entrez IDs
up_genes <- as.character(up$EntrezID)
down_genes <- as.character(down$EntrezID)

# ---------------------------------------------------
# KEGG enrichment
# ---------------------------------------------------

kegg_up <- enrichKEGG(
  gene = up_genes,
  organism = "hsa",
  pvalueCutoff = 0.05
)

kegg_down <- enrichKEGG(
  gene = down_genes,
  organism = "hsa",
  pvalueCutoff = 0.05
)

# Make readable
kegg_up <- setReadable(
  kegg_up,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID"
)

kegg_down <- setReadable(
  kegg_down,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID"
)

# ---------------------------------------------------
# Save tables
# ---------------------------------------------------

write.table(
  as.data.frame(kegg_up),
  file = "results/KEGG_upregulated.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  as.data.frame(kegg_down),
  file = "results/KEGG_downregulated.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ---------------------------------------------------
# DOTPLOTS
# ---------------------------------------------------

png(
  "plots/KEGG_upregulated_dotplot.png",
  width = 3200,
  height = 2400,
  res = 300
)

dotplot(
  kegg_up,
  showCategory = 15,
  font.size = 12,
  title = "KEGG - Upregulated Genes"
)

dev.off()

png(
  "plots/KEGG_downregulated_dotplot.png",
  width = 3200,
  height = 2400,
  res = 300
)

dotplot(
  kegg_down,
  showCategory = 15,
  font.size = 12,
  title = "KEGG - Downregulated Genes"
)

dev.off()

# ---------------------------------------------------
# BARPLOTS
# ---------------------------------------------------

png(
  "plots/KEGG_upregulated_barplot.png",
  width = 3200,
  height = 2400,
  res = 300
)

barplot(
  kegg_up,
  showCategory = 15,
  title = "KEGG - Upregulated Genes"
)

dev.off()

png(
  "plots/KEGG_downregulated_barplot.png",
  width = 3200,
  height = 2400,
  res = 300
)

barplot(
  kegg_down,
  showCategory = 15,
  title = "KEGG - Downregulated Genes"
)

dev.off()

# ---------------------------------------------------
# CNETPLOTS
# ---------------------------------------------------

png(
  "plots/KEGG_upregulated_cnetplot.png",
  width = 3800,
  height = 3200,
  res = 300
)

cnetplot(
  kegg_up,
  showCategory = 5
)

dev.off()

png(
  "plots/KEGG_downregulated_cnetplot.png",
  width = 3800,
  height = 3200,
  res = 300
)

cnetplot(
  kegg_down,
  showCategory = 5
)

dev.off()

# ---------------------------------------------------
# ENRICHMENT MAPS
# ---------------------------------------------------

kegg_up_sim <- pairwise_termsim(kegg_up)
kegg_down_sim <- pairwise_termsim(kegg_down)

png(
  "plots/KEGG_upregulated_emapplot.png",
  width = 3600,
  height = 3000,
  res = 300
)

emapplot(kegg_up_sim, showCategory = 20)

dev.off()

png(
  "plots/KEGG_downregulated_emapplot.png",
  width = 3600,
  height = 3000,
  res = 300
)

emapplot(kegg_down_sim, showCategory = 20)

dev.off()

# ---------------------------------------------------
# PDFs
# ---------------------------------------------------

pdf("plots/KEGG_upregulated_dotplot.pdf", width = 11, height = 8)
dotplot(kegg_up, showCategory = 15)
dev.off()

pdf("plots/KEGG_downregulated_dotplot.pdf", width = 11, height = 8)
dotplot(kegg_down, showCategory = 15)
dev.off()

cat("KEGG enrichment completed successfully.\n")
