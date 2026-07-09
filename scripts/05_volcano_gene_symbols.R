library(EnhancedVolcano)

res <- read.delim(
  "results/deseq2_results_annotated.tsv",
  check.names = FALSE
)

res$log2FoldChange <- as.numeric(res$log2FoldChange)
res$padj <- as.numeric(res$padj)

res <- res[!is.na(res$padj), ]

png("plots/volcano_gene_symbols.png",
    width = 3500,
    height = 3000,
    res = 300)

EnhancedVolcano(
  res,
  lab = res$GeneSymbol,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.05,
  FCcutoff = 1,
  pointSize = 2,
  labSize = 4,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  title = "Tumor vs Control",
  subtitle = "Annotated Differentially Expressed Genes",
  caption = "padj < 0.05, |log2FC| > 1"
)

dev.off()

pdf("plots/volcano_gene_symbols.pdf",
    width = 12,
    height = 10)

EnhancedVolcano(
  res,
  lab = res$GeneSymbol,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.05,
  FCcutoff = 1,
  pointSize = 2,
  labSize = 4,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  title = "Tumor vs Control",
  subtitle = "Annotated Differentially Expressed Genes",
  caption = "padj < 0.05, |log2FC| > 1"
)

dev.off()
