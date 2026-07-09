library(ReactomePA)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(ggplot2)

# -----------------------------
# LOAD ANNOTATED DESEQ2 RESULTS
# -----------------------------

res <- read.delim(
  "results/deseq2_results_annotated.tsv",
  check.names = FALSE
)

res <- res[!is.na(res$padj), ]

# -----------------------------
# SIGNIFICANT DEGs
# -----------------------------

sig <- res[
  res$padj < 0.05 &
    abs(res$log2FoldChange) >= 1,
]

up <- sig[sig$log2FoldChange > 1, ]
down <- sig[sig$log2FoldChange < -1, ]

up_genes <- as.character(up$EntrezID)
down_genes <- as.character(down$EntrezID)

# -----------------------------
# REACTOME ENRICHMENT
# -----------------------------

reactome_up <- enrichPathway(
  gene = up_genes,
  organism = "human",
  pvalueCutoff = 0.05,
  readable = TRUE
)

reactome_down <- enrichPathway(
  gene = down_genes,
  organism = "human",
  pvalueCutoff = 0.05,
  readable = TRUE
)

# -----------------------------
# SAVE TABLES
# -----------------------------

write.table(
  as.data.frame(reactome_up),
  "results/Reactome_upregulated.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  as.data.frame(reactome_down),
  "results/Reactome_downregulated.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# -----------------------------
# DOTPLOTS
# -----------------------------

png(
  "plots/Reactome_upregulated_dotplot.png",
  width = 3200,
  height = 2600,
  res = 300
)

dotplot(
  reactome_up,
  showCategory = 20,
  font.size = 11,
  title = "Reactome Pathways - Upregulated Genes"
)

dev.off()

png(
  "plots/Reactome_downregulated_dotplot.png",
  width = 3200,
  height = 2600,
  res = 300
)

dotplot(
  reactome_down,
  showCategory = 20,
  font.size = 11,
  title = "Reactome Pathways - Downregulated Genes"
)

dev.off()

# -----------------------------
# BARPLOTS
# -----------------------------

png(
  "plots/Reactome_upregulated_barplot.png",
  width = 3200,
  height = 2600,
  res = 300
)

barplot(
  reactome_up,
  showCategory = 20,
  title = "Reactome Pathways - Upregulated Genes"
)

dev.off()

png(
  "plots/Reactome_downregulated_barplot.png",
  width = 3200,
  height = 2600,
  res = 300
)

barplot(
  reactome_down,
  showCategory = 20,
  title = "Reactome Pathways - Downregulated Genes"
)

dev.off()

# -----------------------------
# ENRICHMENT MAPS
# -----------------------------

reactome_up_sim <- pairwise_termsim(reactome_up)
reactome_down_sim <- pairwise_termsim(reactome_down)

png(
  "plots/Reactome_upregulated_emapplot.png",
  width = 3600,
  height = 3200,
  res = 300
)

emapplot(
  reactome_up_sim,
  showCategory = 25
)

dev.off()

png(
  "plots/Reactome_downregulated_emapplot.png",
  width = 3600,
  height = 3200,
  res = 300
)

emapplot(
  reactome_down_sim,
  showCategory = 25
)

dev.off()

# -----------------------------
# CNETPLOTS
# -----------------------------

png(
  "plots/Reactome_upregulated_cnetplot.png",
  width = 4000,
  height = 3400,
  res = 300
)

cnetplot(
  reactome_up,
  showCategory = 6
)

dev.off()

png(
  "plots/Reactome_downregulated_cnetplot.png",
  width = 4000,
  height = 3400,
  res = 300
)

cnetplot(
  reactome_down,
  showCategory = 6
)

dev.off()

# -----------------------------
# PDF OUTPUTS
# -----------------------------

pdf("plots/Reactome_upregulated_dotplot.pdf", width = 12, height = 9)
dotplot(reactome_up, showCategory = 20)
dev.off()

pdf("plots/Reactome_downregulated_dotplot.pdf", width = 12, height = 9)
dotplot(reactome_down, showCategory = 20)
dev.off()

# -----------------------------
# SUMMARY
# -----------------------------

cat("Reactome enrichment completed successfully.\n")
cat("Upregulated Reactome pathways:", nrow(as.data.frame(reactome_up)), "\n")
cat("Downregulated Reactome pathways:", nrow(as.data.frame(reactome_down)), "\n")
