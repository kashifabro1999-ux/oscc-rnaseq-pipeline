library(clusterProfiler)
library(msigdbr)
library(enrichplot)
library(ggplot2)
library(dplyr)

# -----------------------------
# LOAD DESEQ2 RESULTS
# -----------------------------

res <- read.delim(
  "results/deseq2_results_annotated.tsv",
  check.names = FALSE
)

# Remove missing values
res <- res[!is.na(res$log2FoldChange), ]
res <- res[!is.na(res$GeneSymbol), ]

# Remove duplicate gene symbols
res <- res[!duplicated(res$GeneSymbol), ]

# -----------------------------
# CREATE RANKED GENE LIST
# -----------------------------

gene_list <- res$log2FoldChange
names(gene_list) <- res$GeneSymbol

# Sort decreasing for GSEA
gene_list <- sort(gene_list, decreasing = TRUE)

# -----------------------------
# LOAD MSIGDB HALLMARK GENE SETS
# -----------------------------

hallmark <- msigdbr(
  species = "Homo sapiens",
  category = "H"
)

hallmark_sets <- hallmark %>%
  dplyr::select(gs_name, gene_symbol)

# -----------------------------
# RUN GSEA
# -----------------------------

gsea_hallmark <- GSEA(
  geneList = gene_list,
  TERM2GENE = hallmark_sets,
  pvalueCutoff = 1,
  verbose = FALSE
)

# -----------------------------
# SAVE RESULTS
# -----------------------------

write.table(
  as.data.frame(gsea_hallmark),
  "results/GSEA_Hallmark_results.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Separate activated and suppressed pathways
gsea_df <- as.data.frame(gsea_hallmark)

activated <- gsea_df[gsea_df$NES > 0, ]
suppressed <- gsea_df[gsea_df$NES < 0, ]

write.table(
  activated,
  "results/GSEA_Hallmark_activated_in_tumor.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  suppressed,
  "results/GSEA_Hallmark_suppressed_in_tumor.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# -----------------------------
# DOTPLOT
# -----------------------------

png(
  "plots/GSEA_Hallmark_dotplot.png",
  width = 3400,
  height = 2800,
  res = 300
)

dotplot(
  gsea_hallmark,
  showCategory = 25,
  split = ".sign"
) +
  facet_grid(. ~ .sign) +
  ggtitle("GSEA Hallmark Pathways")

dev.off()

pdf(
  "plots/GSEA_Hallmark_dotplot.pdf",
  width = 14,
  height = 10
)

dotplot(
  gsea_hallmark,
  showCategory = 25,
  split = ".sign"
) +
  facet_grid(. ~ .sign) +
  ggtitle("GSEA Hallmark Pathways")

dev.off()

# -----------------------------
# RIDGEPLOT
# -----------------------------

png(
  "plots/GSEA_Hallmark_ridgeplot.png",
  width = 3600,
  height = 3000,
  res = 300
)

ridgeplot(
  gsea_hallmark,
  showCategory = 20
) +
  ggtitle("GSEA Hallmark Ridgeplot")

dev.off()

# -----------------------------
# NES BARPLOT
# -----------------------------

top_gsea <- gsea_df %>%
  arrange(p.adjust) %>%
  head(25)

png(
  "plots/GSEA_Hallmark_NES_barplot.png",
  width = 3600,
  height = 3000,
  res = 300
)

ggplot(
  top_gsea,
  aes(
    x = reorder(Description, NES),
    y = NES,
    fill = NES
  )
) +
  geom_col() +
  coord_flip() +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0
  ) +
  theme_bw() +
  labs(
    title = "Top Hallmark GSEA Pathways",
    x = "",
    y = "Normalized Enrichment Score (NES)"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 10)
  )

dev.off()

# -----------------------------
# INDIVIDUAL ENRICHMENT PLOTS
# -----------------------------

dir.create("plots/GSEA_enrichment_plots", showWarnings = FALSE)

important_pathways <- c(
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_HYPOXIA",
  "HALLMARK_ANGIOGENESIS",
  "HALLMARK_PI3K_AKT_MTOR_SIGNALING",
  "HALLMARK_MTORC1_SIGNALING",
  "HALLMARK_GLYCOLYSIS",
  "HALLMARK_ADIPOGENESIS",
  "HALLMARK_FATTY_ACID_METABOLISM",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_APOPTOSIS",
  "HALLMARK_P53_PATHWAY"
)

available_pathways <- important_pathways[
  important_pathways %in% gsea_hallmark@result$ID
]

for (pathway in available_pathways) {
  
  clean_name <- gsub("HALLMARK_", "", pathway)
  clean_name <- gsub("_", "_", clean_name)
  
  png(
    paste0("plots/GSEA_enrichment_plots/", clean_name, ".png"),
    width = 3000,
    height = 2200,
    res = 300
  )
  
  print(
    gseaplot2(
      gsea_hallmark,
      geneSetID = pathway,
      title = clean_name,
      pvalue_table = TRUE
    )
  )
  
  dev.off()
}

# -----------------------------
# SUMMARY
# -----------------------------

cat("GSEA Hallmark analysis completed successfully.\n")
cat("Total significant Hallmark pathways:", nrow(gsea_df), "\n")
cat("Activated in tumor:", nrow(activated), "\n")
cat("Suppressed in tumor:", nrow(suppressed), "\n")
cat("Results saved in results/GSEA_Hallmark_results.tsv\n")
cat("Plots saved in plots/\n")
