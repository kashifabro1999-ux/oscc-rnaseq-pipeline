library(pheatmap)
library(RColorBrewer)

# Load files
vst <- read.delim("data/vst_counts.tsv", row.names = 1, check.names = FALSE)
res <- read.delim("results/deseq2_results_annotated.tsv", check.names = FALSE)

# Ordered samples
control_order <- paste0("N", 1:10)
tumor_order <- paste0("T", 1:50)
sample_order <- c(control_order, tumor_order)
sample_order <- sample_order[sample_order %in% colnames(vst)]

# Annotation
annotation <- data.frame(
  condition = ifelse(grepl("^N", sample_order), "Control", "Tumor")
)
rownames(annotation) <- sample_order

annotation_colors <- list(
  condition = c(
    Control = "#00BFC4",
    Tumor = "#F8766D"
  )
)

# Color palette
heat_colors <- colorRampPalette(c(
  "#2166AC",
  "#67A9CF",
  "white",
  "#F4A582",
  "#B2182B"
))(255)

# Candidate tumor-promoting genes / oncogenes
tumor_promoters <- c(
  "MMP3", "MMP13", "MMP14", "ADAM12", "PLAU", "TNC",
  "PTK7", "PDPN", "GREM1", "MSN", "MAGEA2", "MAGEA12",
  "KRT76", "KRT36", "KRT222", "SCIN", "MYRIP", "PPP1R1B",
  "ATP2B2", "CRISP3", "FAM3B", "SH3BGRL2"
)

# Candidate tumor-suppressive / tumor-opposing genes
tumor_suppressors <- c(
  "DLG2", "RBP1", "STATH", "BPIFA1", "BPIFA2",
  "AQP5", "AQP5-AS1", "LACRT", "PRB3", "PIP",
  "KRT4", "KRT13", "CRNN", "SPRR1A", "SPRR1B",
  "SLC25A4", "SLC25A21", "ADH4", "CYP3A4"
)

# Function to make heatmap
make_gene_heatmap <- function(gene_list, output_name, title_text) {
  
  selected <- res[res$GeneSymbol %in% gene_list, ]
  selected <- selected[!is.na(selected$padj), ]
  selected <- selected[order(selected$padj), ]
  
  genes_found <- as.character(selected$EntrezID)
  genes_found <- genes_found[genes_found %in% rownames(vst)]
  
  mat <- vst[genes_found, sample_order, drop = FALSE]
  
  gene_symbols <- selected$GeneSymbol[
    match(genes_found, selected$EntrezID)
  ]
  
  rownames(mat) <- gene_symbols
  
  mat <- t(scale(t(mat)))
  mat <- mat[complete.cases(mat), ]
  
  png(
    paste0("plots/", output_name, ".png"),
    width = 4200,
    height = 2800,
    res = 350
  )
  
  pheatmap(
    mat,
    annotation_col = annotation,
    annotation_colors = annotation_colors,
    cluster_cols = FALSE,
    cluster_rows = TRUE,
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize_row = 10,
    fontsize_col = 7,
    angle_col = 45,
    border_color = NA,
    scale = "none",
    color = heat_colors,
    main = title_text
  )
  
  dev.off()
  
  pdf(
    paste0("plots/", output_name, ".pdf"),
    width = 14,
    height = 9
  )
  
  pheatmap(
    mat,
    annotation_col = annotation,
    annotation_colors = annotation_colors,
    cluster_cols = FALSE,
    cluster_rows = TRUE,
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize_row = 10,
    fontsize_col = 7,
    angle_col = 45,
    border_color = NA,
    scale = "none",
    color = heat_colors,
    main = title_text
  )
  
  dev.off()
  
  cat(title_text, "\n")
  cat("Genes found:", rownames(mat), "\n\n")
}

# Make heatmaps
make_gene_heatmap(
  tumor_promoters,
  "tumor_promoting_genes_heatmap",
  "Tumor-Promoting / Oncogenic Genes"
)

make_gene_heatmap(
  tumor_suppressors,
  "tumor_suppressor_genes_heatmap",
  "Tumor-Suppressive / Tumor-Opposing Genes"
)
