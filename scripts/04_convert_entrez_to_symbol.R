library(AnnotationDbi)
library(org.Hs.eg.db)

# Read DESeq2 result without header
res <- read.delim("data/deseq2_results.tsv",
                  header = FALSE,
                  check.names = FALSE)

colnames(res) <- c(
  "EntrezID",
  "baseMean",
  "log2FoldChange",
  "lfcSE",
  "stat",
  "pvalue",
  "padj"
)

# Convert Entrez IDs to gene symbols
res$GeneSymbol <- mapIds(
  org.Hs.eg.db,
  keys = as.character(res$EntrezID),
  column = "SYMBOL",
  keytype = "ENTREZID",
  multiVals = "first"
)

# Convert Entrez IDs to gene names/descriptions
res$GeneName <- mapIds(
  org.Hs.eg.db,
  keys = as.character(res$EntrezID),
  column = "GENENAME",
  keytype = "ENTREZID",
  multiVals = "first"
)

# Reorder columns
res <- res[, c(
  "EntrezID",
  "GeneSymbol",
  "GeneName",
  "baseMean",
  "log2FoldChange",
  "lfcSE",
  "stat",
  "pvalue",
  "padj"
)]

# Save annotated result
write.table(res,
            file = "results/deseq2_results_annotated.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

# Save significant DEGs
sig <- res[!is.na(res$padj) & res$padj < 0.05 & abs(res$log2FoldChange) >= 1, ]

write.table(sig,
            file = "results/significant_DEGs_annotated.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

# Save upregulated and downregulated separately
up <- sig[sig$log2FoldChange >= 1, ]
down <- sig[sig$log2FoldChange <= -1, ]

write.table(up,
            file = "results/upregulated_genes.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

write.table(down,
            file = "results/downregulated_genes.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

cat("Total genes:", nrow(res), "\n")
cat("Significant DEGs:", nrow(sig), "\n")
cat("Upregulated:", nrow(up), "\n")
cat("Downregulated:", nrow(down), "\n")
cat("Annotation completed successfully.\n")
