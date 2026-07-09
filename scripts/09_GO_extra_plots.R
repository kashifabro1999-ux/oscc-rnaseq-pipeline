library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

res <- read.delim("results/deseq2_results_annotated.tsv", check.names = FALSE)
res <- res[!is.na(res$padj), ]

sig <- res[res$padj < 0.05 & abs(res$log2FoldChange) >= 1, ]

up <- sig[sig$log2FoldChange > 1, ]
down <- sig[sig$log2FoldChange < -1, ]

up_genes <- as.character(up$EntrezID)
down_genes <- as.character(down$EntrezID)

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

# Barplots
png("plots/GO_upregulated_barplot.png", width = 3200, height = 2400, res = 300)
barplot(ego_up, showCategory = 15, title = "GO BP - Upregulated Genes")
dev.off()

png("plots/GO_downregulated_barplot.png", width = 3200, height = 2400, res = 300)
barplot(ego_down, showCategory = 15, title = "GO BP - Downregulated Genes")
dev.off()

pdf("plots/GO_upregulated_barplot.pdf", width = 11, height = 8)
barplot(ego_up, showCategory = 15, title = "GO BP - Upregulated Genes")
dev.off()

pdf("plots/GO_downregulated_barplot.pdf", width = 11, height = 8)
barplot(ego_down, showCategory = 15, title = "GO BP - Downregulated Genes")
dev.off()

# Cnetplots
png("plots/GO_upregulated_cnetplot.png", width = 3800, height = 3200, res = 300)
cnetplot(ego_up, showCategory = 5)
dev.off()

png("plots/GO_downregulated_cnetplot.png", width = 3800, height = 3200, res = 300)
cnetplot(ego_down, showCategory = 5)
dev.off()

pdf("plots/GO_upregulated_cnetplot.pdf", width = 12, height = 10)
cnetplot(ego_up, showCategory = 5)
dev.off()

pdf("plots/GO_downregulated_cnetplot.pdf", width = 12, height = 10)
cnetplot(ego_down, showCategory = 5)
dev.off()

# Enrichment maps
ego_up_sim <- pairwise_termsim(ego_up)
ego_down_sim <- pairwise_termsim(ego_down)

png("plots/GO_upregulated_emapplot.png", width = 3600, height = 3000, res = 300)
emapplot(ego_up_sim, showCategory = 20)
dev.off()

png("plots/GO_downregulated_emapplot.png", width = 3600, height = 3000, res = 300)
emapplot(ego_down_sim, showCategory = 20)
dev.off()

pdf("plots/GO_upregulated_emapplot.pdf", width = 12, height = 10)
emapplot(ego_up_sim, showCategory = 20)
dev.off()

pdf("plots/GO_downregulated_emapplot.pdf", width = 12, height = 10)
emapplot(ego_down_sim, showCategory = 20)
dev.off()

cat("Extra GO plots completed successfully.\n")
