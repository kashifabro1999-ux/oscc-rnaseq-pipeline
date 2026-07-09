library(pheatmap)
library(RColorBrewer)

# -----------------------------
# LOAD FILES
# -----------------------------

expr <- read.delim(
  "data/vst_counts.tsv",
  row.names = 1,
  check.names = FALSE
)

meta <- read.delim(
  "data/sample_info.tsv",
  check.names = FALSE
)

hub <- read.delim(
  "results/top30_hub_genes.tsv",
  check.names = FALSE
)

# -----------------------------
# BIOMARKER GENES
# -----------------------------

genes <- c(
  "COL1A1",
  "MMP14",
  "PLAU",
  "MMP9",
  "PPARG",
  "ACACB",
  "GPD1",
  "LPL",
  "PCK1",
  "PLIN1",
  "SLC2A4",
  "ADH1B",
  "CYP3A4"
)

# -----------------------------
# MAP GENE SYMBOLS
# -----------------------------

symbol_map <- hub$EntrezID
names(symbol_map) <- hub$GeneSymbol

available_genes <- genes[
  genes %in% names(symbol_map)
]

entrez_ids <- as.character(
  symbol_map[available_genes]
)

entrez_ids <- entrez_ids[
  entrez_ids %in% rownames(expr)
]

# -----------------------------
# EXTRACT MATRIX
# -----------------------------

mat <- expr[entrez_ids, ]

gene_symbols <- names(symbol_map)[
  match(
    rownames(mat),
    as.character(symbol_map)
  )
]

rownames(mat) <- gene_symbols

# -----------------------------
# ORDER SAMPLES
# -----------------------------

control_samples <- meta$sample[
  meta$condition == "Control"
]

tumor_samples <- meta$sample[
  meta$condition == "Tumor"
]

ordered_samples <- c(
  control_samples,
  tumor_samples
)

ordered_samples <- ordered_samples[
  ordered_samples %in% colnames(mat)
]

mat <- mat[, ordered_samples]

# -----------------------------
# SCALE MATRIX
# -----------------------------

mat_scaled <- t(
  scale(
    t(mat)
  )
)

# -----------------------------
# ANNOTATION
# -----------------------------

annotation_col <- data.frame(
  Condition = factor(
    c(
      rep(
        "Control",
        length(control_samples)
      ),
      rep(
        "Tumor",
        length(tumor_samples)
      )
    ),
    levels = c(
      "Control",
      "Tumor"
    )
  )
)

rownames(annotation_col) <- ordered_samples

# -----------------------------
# COLORS
# -----------------------------

ann_colors <- list(
  Condition = c(
    Control = "#00BFC4",
    Tumor = "#F8766D"
  )
)

heat_colors <- colorRampPalette(
  c(
    "#313695",
    "#4575B4",
    "#74ADD1",
    "#ABD9E9",
    "#FFFFBF",
    "#FDAE61",
    "#F46D43",
    "#D73027"
  )
)(100)

# -----------------------------
# PNG OUTPUT
# -----------------------------

png(
  "plots/top_biomarker_heatmap.png",
  width = 6900,
  height = 4100,
  res = 300
)

pheatmap(
  mat_scaled,

  annotation_col = annotation_col,
  annotation_colors = ann_colors,

  cluster_rows = TRUE,
  cluster_cols = FALSE,

  fontsize_row = 22,
  fontsize_col = 16,

  border_color = NA,

  color = heat_colors,

  show_colnames = TRUE,
  show_rownames = TRUE,

  angle_col = 90,

  treeheight_row = 80,
  treeheight_col = 0,

  main = "Top Biomarker Expression Heatmap"
)

dev.off()

# -----------------------------
# PDF OUTPUT
# -----------------------------

pdf(
  "plots/top_biomarker_heatmap.pdf",
  width = 28,
  height = 16
)

pheatmap(
  mat_scaled,

  annotation_col = annotation_col,
  annotation_colors = ann_colors,

  cluster_rows = TRUE,
  cluster_cols = FALSE,

  fontsize_row = 18,
  fontsize_col = 8,

  border_color = NA,

  color = heat_colors,

  show_colnames = TRUE,
  show_rownames = TRUE,

  angle_col = 90,

  treeheight_row = 80,
  treeheight_col = 0,

  main = "Top Biomarker Expression Heatmap"
)

dev.off()

# -----------------------------
# DONE
# -----------------------------

cat(
  "Publication-quality heatmap generated successfully.\n"
)
