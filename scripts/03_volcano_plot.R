library(EnhancedVolcano)

res <- read.delim("data/deseq2_results.tsv",
                  header = FALSE,
                  check.names = FALSE)

colnames(res) <- c(
  "gene",
  "baseMean",
  "log2FoldChange",
  "lfcSE",
  "stat",
  "pvalue",
  "padj"
)

res$log2FoldChange <- as.numeric(res$log2FoldChange)
res$padj <- as.numeric(res$padj)

res <- res[!is.na(res$padj), ]

png("plots/volcano_plot.png",
    width = 3000,
    height = 2500,
    res = 300)

EnhancedVolcano(
  res,
  lab = res$gene,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.05,
  FCcutoff = 1,
  pointSize = 2.0,
  labSize = 3.0,
  title = "Tumor vs Control",
  subtitle = "Differentially Expressed Genes",
  caption = "Cutoff: padj < 0.05, |log2FC| > 1"
)

dev.off()

pdf("plots/volcano_plot.pdf",
    width = 10,
    height = 8)

EnhancedVolcano(
  res,
  lab = res$gene,
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.05,
  FCcutoff = 1,
  pointSize = 2.0,
  labSize = 3.0,
  title = "Tumor vs Control",
  subtitle = "Differentially Expressed Genes",
  caption = "Cutoff: padj < 0.05, |log2FC| > 1"
)

dev.off()
