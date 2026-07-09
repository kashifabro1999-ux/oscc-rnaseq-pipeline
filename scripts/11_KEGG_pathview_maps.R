library(pathview)
library(org.Hs.eg.db)
library(AnnotationDbi)

# Read annotated DESeq2 results
res <- read.delim(
  "results/deseq2_results_annotated.tsv",
  check.names = FALSE
)

# Remove NA
res <- res[!is.na(res$padj), ]

# Prepare named log2FC vector
gene_fc <- res$log2FoldChange
names(gene_fc) <- as.character(res$EntrezID)

# Create output folder
dir.create("plots/pathview", showWarnings = FALSE)

# Move to output folder because pathview saves files in working directory
setwd("plots/pathview")

# KEGG pathway maps
pathview(gene.data = gene_fc,
         pathway.id = "04151",   # PI3K-Akt signaling pathway
         species = "hsa",
         out.suffix = "PI3K_Akt",
         kegg.native = TRUE)

pathview(gene.data = gene_fc,
         pathway.id = "04512",   # ECM-receptor interaction
         species = "hsa",
         out.suffix = "ECM_receptor",
         kegg.native = TRUE)

pathview(gene.data = gene_fc,
         pathway.id = "04510",   # Focal adhesion
         species = "hsa",
         out.suffix = "Focal_adhesion",
         kegg.native = TRUE)

pathview(gene.data = gene_fc,
         pathway.id = "05205",   # Proteoglycans in cancer
         species = "hsa",
         out.suffix = "Proteoglycans_cancer",
         kegg.native = TRUE)

pathview(gene.data = gene_fc,
         pathway.id = "05165",   # Human papillomavirus infection
         species = "hsa",
         out.suffix = "HPV_infection",
         kegg.native = TRUE)

cat("KEGG pathway maps completed successfully.\n")
