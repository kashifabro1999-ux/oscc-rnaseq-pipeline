# OSCC RNA-seq Analysis: Complete Step-by-Step Tutorial

**From raw counts to candidate biomarkers — a reproducible walkthrough of all 20 analysis scripts**

This tutorial documents the full computational workflow used to analyze the integrated
five-region OSCC (oral squamous cell carcinoma) RNA-seq dataset: 60 paired-end samples
(50 tumor, 10 control) from Taiwan, the United States, India, Germany, and China.

Each of the 20 numbered R scripts is presented with:
- **What it does** and why it fits in the pipeline
- **Ubuntu installation commands** for every system library and R/Bioconductor package
  it needs (run these *before* the script, in order)
- **Required input files** and **files it produces**
- **The complete, unmodified script**

> **How to use this document:** work through the sections in numeric order (01 → 20).
> Each script assumes the outputs of earlier scripts already exist on disk. Run everything
> from one consistent project folder (see the directory layout below) so relative paths
> (`data/...`, `results/...`, `plots/...`) resolve correctly.

---

## Table of Contents

0. [Prerequisites & Project Setup](#0-prerequisites--project-setup)
1. [Where These Scripts Fit: Upstream Pipeline (Raw Reads → Count Matrix)](#1-where-these-scripts-fit-upstream-pipeline-raw-reads--count-matrix)
2. [Expected Directory Layout](#2-expected-directory-layout)
3. [One-Shot Install (All Packages, All Scripts)](#3-one-shot-install-all-packages-all-scripts)

1. [Step 1 — Principal Component Analysis (PCA)](#01-pca-plotr)
2. [Step 2 — Sample-to-Sample Distance Heatmap](#02-sample-distance-heatmapr)
3. [Step 3 — Volcano Plot (Entrez ID stage)](#03-volcano-plotr)
4. [Step 4 — Annotate DESeq2 Results (Entrez ID → Gene Symbol)](#04-convert-entrez-to-symbolr)
5. [Step 5 — Volcano Plot (Gene Symbol stage)](#05-volcano-gene-symbolsr)
6. [Step 6 — Heatmap of Top 50 Differentially Expressed Genes](#06-top-deg-heatmapr)
7. [Step 7 — Tumor-Promoting vs Tumor-Opposing Gene Heatmaps](#07-onco-tsg-heatmapsr)
8. [Step 8 — Gene Ontology (Biological Process) Enrichment](#08-go-enrichmentr)
9. [Step 9 — Additional GO Visualizations (Barplots, Cnetplots, Enrichment Maps)](#09-go-extra-plotsr)
10. [Step 10 — KEGG Pathway Enrichment](#10-kegg-enrichmentr)
11. [Step 11 — KEGG Pathview Maps](#11-kegg-pathview-mapsr)
12. [Step 12 — STRING Protein-Protein Interaction (PPI) Mapping & Hub-Gene Ranking](#12-hub-gene-string-analysisr)
13. [Step 13 — PPI Network Visualization](#13-string-network-plotr)
14. [Step 14 — Top 30 Hub Gene Expression Heatmap](#14-hub-gene-heatmapr)
15. [Step 15 — Export Candidate Gene List for Survival Analysis](#15-survival-candidate-genesr)
16. [Step 16 — Boxplots of Top 10 Hub Genes](#16-hub-gene-boxplotsr)
17. [Step 17 — Biomarker Priority Table](#17-biomarker-priority-tabler)
18. [Step 18 — Publication-Quality Top Biomarker Heatmap](#18-top-biomarker-heatmapr)
19. [Step 19 — Reactome Pathway Enrichment](#19-reactome-enrichmentr)
20. [Step 20 — Hallmark Gene Set Enrichment Analysis (GSEA)](#20-gsea-hallmark-analysisr)

---


## 0. Prerequisites & Project Setup

This tutorial assumes a fresh **Ubuntu 22.04 / 24.04** machine (or WSL/Docker container).
Run this once before anything else.

### 0.1 Update the system and install R

```bash
sudo apt-get update && sudo apt-get upgrade -y

# Add the official CRAN apt repository so you get a current R version
# (Ubuntu's default repos often ship an outdated R)
sudo apt-get install -y --no-install-recommends software-properties-common dirmngr
wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | \
  sudo tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
sudo add-apt-repository -y "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"
sudo apt-get update
sudo apt-get install -y r-base r-base-dev
```

### 0.2 Install system libraries that R packages compile against

Almost every Bioconductor/CRAN package used in this pipeline links against these
system libraries. Installing them upfront avoids compilation failures later.

```bash
sudo apt-get install -y \
  build-essential gfortran \
  libcurl4-openssl-dev libssl-dev libxml2-dev \
  libfontconfig1-dev libharfbuzz-dev libfribidi-dev \
  libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev \
  libglpk-dev libgmp-dev \
  libcairo2-dev libxt-dev \
  zlib1g-dev libbz2-dev liblzma-dev \
  pandoc git cmake \
  libmagick++-dev
```

### 0.3 Install BiocManager (one-time, inside R)

Every Bioconductor package below (`EnhancedVolcano`, `clusterProfiler`, `org.Hs.eg.db`,
`STRINGdb`, `pathview`, `ReactomePA`, etc.) is installed through `BiocManager`, not
`install.packages()` directly. Set this up once:

```bash
sudo Rscript -e 'install.packages("BiocManager", repos = "https://cloud.r-project.org")'
sudo Rscript -e 'BiocManager::install(version = "3.19", ask = FALSE, update = FALSE)'
```

> Run R commands with `sudo Rscript -e '...'` (or from an R session with admin rights)
> if you want packages installed system-wide in a shared library path. If you're working
> in a personal user library instead, drop `sudo` and just run `Rscript -e '...'`.

With this done, you're ready to install packages for individual scripts, either
one section at a time (recommended — see each step below) or all at once
(Section 3).

---

## 1. Where These Scripts Fit: Upstream Pipeline (Raw Reads → Count Matrix)

Scripts 01–20 all start from **already-processed** data:
`data/vst_counts.tsv`, `data/sample_info.tsv`, and `data/deseq2_results.tsv`.
Those files are the *output* of an upstream raw-data processing pipeline that is
not part of these 20 scripts but is needed to produce their inputs. For
completeness, here is that upstream stage and how to install it on Ubuntu:

```bash
# Quality control
sudo apt-get install -y fastqc

# fastp (adapter/quality trimming) - precompiled static binary is simplest
wget http://opengene.org/fastp/fastp
chmod a+x ./fastp
sudo mv ./fastp /usr/local/bin/

# MultiQC (aggregates FastQC/fastp/HISAT2 reports)
sudo apt-get install -y python3-pip
pip3 install multiqc

# HISAT2 (spliced alignment to hg38)
sudo apt-get install -y hisat2

# Subread / featureCounts (gene-level quantification)
sudo apt-get install -y subread
```

Typical order of operations upstream of Script 01:

1. `fastqc` on all raw FASTQ files → inspect quality.
2. `fastp` to trim adapters/low-quality bases.
3. `multiqc` to summarize QC reports before/after trimming.
4. `hisat2` to align trimmed reads to the hg38 reference genome.
5. `featureCounts` (Subread) to produce a raw gene-level count matrix.
6. **DESeq2** (in R) to normalize counts, run the tumor-vs-control differential
   expression test, and produce a variance-stabilizing-transformed (VST) matrix.
   This is where `data/vst_counts.tsv` and `data/deseq2_results.tsv` come from,
   and `data/sample_info.tsv` is simply a two/three-column table of
   `sample`, `condition` (Tumor/Control), and optionally `region`.

Install DESeq2 (needed to generate those inputs, and implicitly relied upon
throughout this pipeline):

```bash
sudo Rscript -e 'BiocManager::install("DESeq2", ask = FALSE, update = FALSE)'
```

Once `data/vst_counts.tsv`, `data/sample_info.tsv`, and `data/deseq2_results.tsv`
exist, you're ready to start at Script 01 below.

---

## 2. Expected Directory Layout

Create this structure in your project folder before running any script:

```
project/
├── data/
│   ├── vst_counts.tsv        # VST-normalized expression matrix (genes x samples)
│   ├── sample_info.tsv       # sample, condition (Tumor/Control), region, etc.
│   └── deseq2_results.tsv    # raw DESeq2 output, no header: EntrezID, baseMean,
│                              #   log2FoldChange, lfcSE, stat, pvalue, padj
├── results/                  # all .tsv result tables are written here
├── plots/                    # all .png / .pdf figures are written here
└── scripts/
    ├── 01_pca_plot.R
    ├── 02_sample_distance_heatmap.R
    ├── ... (through 20_GSEA_hallmark_analysis.R)
```

Create the output folders once:

```bash
mkdir -p project/data project/results project/plots project/scripts
```

Each script below is written assuming your **working directory is `project/`**
(i.e., you run `Rscript scripts/01_pca_plot.R` from inside `project/`, not from
inside `scripts/`), since all paths in the scripts are relative (`data/...`,
`results/...`, `plots/...`).

---

## 3. One-Shot Install (All Packages, All Scripts)

If you'd rather install everything up front instead of package-by-package,
run this once (it covers every script from 01 to 20):

```bash
sudo Rscript -e '
install.packages(c(
  "ggplot2", "pheatmap", "RColorBrewer", "dplyr", "tidyr",
  "igraph", "ggraph", "msigdbr"
), repos = "https://cloud.r-project.org")

BiocManager::install(c(
  "DESeq2", "EnhancedVolcano", "AnnotationDbi", "org.Hs.eg.db",
  "clusterProfiler", "enrichplot", "pathview", "STRINGdb",
  "ReactomePA", "reactome.db"
), update = FALSE, ask = FALSE)
'
```

The rest of this document repeats the *specific* subset of these commands
needed before each individual script, in case you prefer to install
incrementally as you go.

---

## Step 1 — Principal Component Analysis (PCA)

**Script file:** `01_pca_plot.R`


**What it does:** Runs PCA on the variance-stabilized (VST) expression matrix to see how tumor and control samples separate globally, and saves a PCA scatter plot colored by condition.


**Required input file(s):**

- `data/vst_counts.tsv`

- `data/sample_info.tsv`


**Output file(s) produced:**

- `plots/PCA_Tumor_vs_Control.png`

- `plots/PCA_Tumor_vs_Control.pdf`



### Install packages needed for this script (Ubuntu)


```bash
sudo Rscript -e '
install.packages(c("ggplot2"), repos = "https://cloud.r-project.org")
'
```


### Full script: `01_pca_plot.R`


```r
library(ggplot2)

vst <- read.delim("data/vst_counts.tsv", row.names = 1, check.names = FALSE)
sample_info <- read.delim("data/sample_info.tsv", check.names = FALSE)

vst <- as.matrix(vst)

pca <- prcomp(t(vst), scale. = FALSE)

pca_data <- as.data.frame(pca$x)
pca_data$sample <- rownames(pca_data)

pca_data <- merge(pca_data, sample_info, by = "sample")

percent_var <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)), 1)

p <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition, label = sample)) +
  geom_point(size = 3) +
  geom_text(vjust = -0.8, size = 2.5) +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  theme_bw() +
  ggtitle("PCA Plot: Tumor vs Control")

ggsave("plots/PCA_Tumor_vs_Control.png", p, width = 8, height = 6, dpi = 300)
ggsave("plots/PCA_Tumor_vs_Control.pdf", p, width = 8, height = 6)

print(p)
```


---

## Step 2 — Sample-to-Sample Distance Heatmap

**Script file:** `02_sample_distance_heatmap.R`


**What it does:** Computes Euclidean distances between all samples on VST-normalized expression and draws a clustered heatmap, letting you spot outliers and confirm tumor/control clustering.


**Required input file(s):**

- `data/vst_counts.tsv`

- `data/sample_info.tsv`


**Output file(s) produced:**

- `plots/sample_distance_heatmap.png`

- `plots/sample_distance_heatmap.pdf`



### Install packages needed for this script (Ubuntu)


```bash
sudo Rscript -e '
install.packages(c("pheatmap", "RColorBrewer"), repos = "https://cloud.r-project.org")
'
```


### Full script: `02_sample_distance_heatmap.R`


```r
library(pheatmap)
library(RColorBrewer)

vst <- read.delim("data/vst_counts.tsv",
                  row.names = 1,
                  check.names = FALSE)

sample_info <- read.delim("data/sample_info.tsv",
                          check.names = FALSE)

vst <- as.matrix(vst)

sample_dists <- dist(t(vst))

sample_dist_matrix <- as.matrix(sample_dists)

rownames(sample_dist_matrix) <- colnames(vst)
colnames(sample_dist_matrix) <- colnames(vst)

annotation <- data.frame(condition = sample_info$condition)
rownames(annotation) <- sample_info$sample

pheatmap(sample_dist_matrix,
         annotation_col = annotation,
         annotation_row = annotation,
         clustering_distance_rows = sample_dists,
         clustering_distance_cols = sample_dists,
         col = colorRampPalette(
           rev(brewer.pal(9, "Blues"))
         )(255),
         main = "Sample-to-Sample Distance Heatmap",
         filename = "plots/sample_distance_heatmap.png",
         width = 10,
         height = 10)

pdf("plots/sample_distance_heatmap.pdf",
    width = 10,
    height = 10)

pheatmap(sample_dist_matrix,
         annotation_col = annotation,
         annotation_row = annotation,
         clustering_distance_rows = sample_dists,
         clustering_distance_cols = sample_dists,
         col = colorRampPalette(
           rev(brewer.pal(9, "Blues"))
         )(255),
         main = "Sample-to-Sample Distance Heatmap")

dev.off()
```


---

## Step 3 — Volcano Plot (Entrez ID stage)

**Script file:** `03_volcano_plot.R`


**What it does:** Builds the first volcano plot straight from the raw DESeq2 results table (still keyed by Entrez ID), using padj < 0.05 and |log2FC| >= 1 as the significance cutoffs.


**Required input file(s):**

- `data/deseq2_results.tsv`


**Output file(s) produced:**

- `plots/volcano_plot.png`

- `plots/volcano_plot.pdf`



### Install packages needed for this script (Ubuntu)


```bash
sudo Rscript -e '
BiocManager::install(c("EnhancedVolcano"), update = FALSE, ask = FALSE)
'
```


### Full script: `03_volcano_plot.R`


```r
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
```


---

## Step 4 — Annotate DESeq2 Results (Entrez ID → Gene Symbol)

**Script file:** `04_convert_entrez_to_symbol.R`


**What it does:** Maps every Entrez Gene ID in the DESeq2 output to its HGNC gene symbol and full gene name using org.Hs.eg.db, then writes out the full annotated table plus separate significant / upregulated / downregulated gene tables. This annotated file is the main input for almost every later script.


**Required input file(s):**

- `data/deseq2_results.tsv`


**Output file(s) produced:**

- `results/deseq2_results_annotated.tsv`

- `results/significant_DEGs_annotated.tsv`

- `results/upregulated_genes.tsv`

- `results/downregulated_genes.tsv`



### Install packages needed for this script (Ubuntu)


```bash
sudo Rscript -e '
BiocManager::install(c("AnnotationDbi", "org.Hs.eg.db"), update = FALSE, ask = FALSE)
'
```


### Full script: `04_convert_entrez_to_symbol.R`


```r
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
```


---

## Step 5 — Volcano Plot (Gene Symbol stage)

**Script file:** `05_volcano_gene_symbols.R`


**What it does:** Re-draws the volcano plot using the annotated table from Step 4, so points are now labeled with readable gene symbols instead of Entrez IDs, with connector lines for label clarity.


**Required input file(s):**

- `results/deseq2_results_annotated.tsv`


**Output file(s) produced:**

- `plots/volcano_gene_symbols.png`

- `plots/volcano_gene_symbols.pdf`



### Install packages needed for this script (Ubuntu)


```bash
sudo Rscript -e '
BiocManager::install(c("EnhancedVolcano"), update = FALSE, ask = FALSE)
'
```


### Full script: `05_volcano_gene_symbols.R`


```r
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
```


---

## Step 6 — Heatmap of Top 50 Differentially Expressed Genes

**Script file:** `06_top_DEG_heatmap.R`


**What it does:** Selects the 50 most significant DEGs (padj < 0.05, |log2FC| >= 2), z-score normalizes their VST expression, manually orders samples as Controls-then-Tumors, and draws an annotated heatmap.


**Required input file(s):**

- `data/vst_counts.tsv`

- `results/deseq2_results_annotated.tsv`


**Output file(s) produced:**

- `plots/top50_DEG_heatmap_ordered.png`

- `plots/top50_DEG_heatmap_ordered.pdf`



### Install packages needed for this script (Ubuntu)


```bash
sudo Rscript -e '
install.packages(c("pheatmap", "RColorBrewer"), repos = "https://cloud.r-project.org")
'
```


### Full script: `06_top_DEG_heatmap.R`


```r
library(pheatmap)
library(RColorBrewer)

# Read files
vst <- read.delim(
  "data/vst_counts.tsv",
  row.names = 1,
  check.names = FALSE
)

res <- read.delim(
  "results/deseq2_results_annotated.tsv",
  check.names = FALSE
)

sample_info <- read.delim(
  "data/sample_info.tsv",
  check.names = FALSE
)

# Significant genes
res <- res[!is.na(res$padj), ]
sig <- res[res$padj < 0.05 & abs(res$log2FoldChange) >= 2, ]

# Top 50 genes
top50 <- head(sig[order(sig$padj), ], 50)

# Match genes
common_genes <- intersect(
  rownames(vst),
  as.character(top50$EntrezID)
)

heatmap_matrix <- vst[common_genes, ]

# Replace IDs with symbols
gene_symbols <- top50$GeneSymbol[
  match(common_genes, top50$EntrezID)
]

rownames(heatmap_matrix) <- gene_symbols

# Z-score normalization
heatmap_matrix <- t(scale(t(heatmap_matrix)))
heatmap_matrix <- heatmap_matrix[complete.cases(heatmap_matrix), ]

# -------------------------------------------------
# ORDER SAMPLES MANUALLY
# -------------------------------------------------

control_order <- paste0("N", 1:10)
tumor_order <- paste0("T", 1:50)

desired_order <- c(control_order, tumor_order)

desired_order <- desired_order[
  desired_order %in% colnames(heatmap_matrix)
]

heatmap_matrix <- heatmap_matrix[, desired_order]

# Sample annotation
annotation <- data.frame(
  condition = ifelse(
    grepl("^N", desired_order),
    "Control",
    "Tumor"
  )
)

rownames(annotation) <- desired_order

annotation_colors <- list(
  condition = c(
    Control = "#00BFC4",
    Tumor = "#F8766D"
  )
)

# Professional color palette
heat_colors <- colorRampPalette(c(
  "#2166AC",
  "#67A9CF",
  "white",
  "#F4A582",
  "#B2182B"
))(255)

# PNG
png(
  "plots/top50_DEG_heatmap_ordered.png",
  width = 5200,
  height = 4200,
  res = 350
)

pheatmap(
  heatmap_matrix,
  annotation_col = annotation,
  annotation_colors = annotation_colors,
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 10,
  fontsize_col = 8,
  angle_col = 45,
  border_color = NA,
  scale = "none",
  color = heat_colors,
  main = "Top 50 Differentially Expressed Genes"
)

dev.off()

# PDF
pdf(
  "plots/top50_DEG_heatmap_ordered.pdf",
  width = 16,
  height = 13
)

pheatmap(
  heatmap_matrix,
  annotation_col = annotation,
  annotation_colors = annotation_colors,
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 10,
  fontsize_col = 8,
  angle_col = 45,
  border_color = NA,
  scale = "none",
  color = heat_colors,
  main = "Top 50 Differentially Expressed Genes"
)

dev.off()
```


---

## Step 7 — Tumor-Promoting vs Tumor-Opposing Gene Heatmaps

**Script file:** `07_onco_tsg_heatmaps.R`


**What it does:** Uses two curated gene lists (invasion/oncogene-like genes vs. normal tissue-identity genes) and draws a dedicated heatmap for each, showing coordinated up- or down-regulation across tumor vs control samples.


**Required input file(s):**

- `data/vst_counts.tsv`

- `results/deseq2_results_annotated.tsv`


**Output file(s) produced:**

- `plots/tumor_promoting_genes_heatmap.png/.pdf`

- `plots/tumor_suppressor_genes_heatmap.png/.pdf`



### Install packages needed for this script (Ubuntu)


```bash
sudo Rscript -e '
install.packages(c("pheatmap", "RColorBrewer"), repos = "https://cloud.r-project.org")
'
```


### Full script: `07_onco_tsg_heatmaps.R`


```r
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
```


---

## Step 8 — Gene Ontology (Biological Process) Enrichment

**Script file:** `08_GO_enrichment.R`


**What it does:** Splits significant DEGs into upregulated and downregulated sets and runs GO Biological Process over-representation analysis (clusterProfiler::enrichGO) on each, saving result tables and dotplots.


**Required input file(s):**

- `results/deseq2_results_annotated.tsv`


**Output file(s) produced:**

- `results/GO_upregulated.tsv`

- `results/GO_downregulated.tsv`

- `plots/GO_upregulated_dotplot.png/.pdf`

- `plots/GO_downregulated_dotplot.png/.pdf`



### Install packages needed for this script (Ubuntu)


```bash
sudo Rscript -e '
install.packages(c("ggplot2"), repos = "https://cloud.r-project.org")
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot"), update = FALSE, ask = FALSE)
'
```


### Full script: `08_GO_enrichment.R`


```r
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

# Read annotated DEGs
res <- read.delim(
  "results/deseq2_results_annotated.tsv",
  check.names = FALSE
)

# Remove NA values
res <- res[!is.na(res$padj), ]

# Significant genes
sig <- res[
  res$padj < 0.05 &
  abs(res$log2FoldChange) >= 1,
]

# Separate upregulated and downregulated
up <- sig[sig$log2FoldChange > 1, ]
down <- sig[sig$log2FoldChange < -1, ]

# Entrez IDs
up_genes <- as.character(up$EntrezID)
down_genes <- as.character(down$EntrezID)

# --------------------------------------------------
# GO enrichment for UPREGULATED genes
# --------------------------------------------------

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

# Save table
write.table(
  as.data.frame(ego_up),
  file = "results/GO_upregulated.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Dotplot PNG
png(
  "plots/GO_upregulated_dotplot.png",
  width = 3200,
  height = 2400,
  res = 300
)

dotplot(
  ego_up,
  showCategory = 15,
  font.size = 12,
  title = "GO Enrichment - Upregulated Genes"
)

dev.off()

# --------------------------------------------------
# GO enrichment for DOWNREGULATED genes
# --------------------------------------------------

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

# Save table
write.table(
  as.data.frame(ego_down),
  file = "results/GO_downregulated.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Dotplot PNG
png(
  "plots/GO_downregulated_dotplot.png",
  width = 3200,
  height = 2400,
  res = 300
)

dotplot(
  ego_down,
  showCategory = 15,
  font.size = 12,
  title = "GO Enrichment - Downregulated Genes"
)

dev.off()

# --------------------------------------------------
# PDF outputs
# --------------------------------------------------

pdf(
  "plots/GO_upregulated_dotplot.pdf",
  width = 11,
  height = 8
)

dotplot(
  ego_up,
  showCategory = 15,
  font.size = 12,
  title = "GO Enrichment - Upregulated Genes"
)

dev.off()

pdf(
  "plots/GO_downregulated_dotplot.pdf",
  width = 11,
  height = 8
)

dotplot(
  ego_down,
  showCategory = 15,
  font.size = 12,
  title = "GO Enrichment - Downregulated Genes"
)

dev.off()

# --------------------------------------------------
# Summary
# --------------------------------------------------

cat("Upregulated GO terms:", nrow(as.data.frame(ego_up)), "\n")
cat("Downregulated GO terms:", nrow(as.data.frame(ego_down)), "\n")
cat("GO enrichment analysis completed.\n")
```


---

## Step 9 — Additional GO Visualizations (Barplots, Cnetplots, Enrichment Maps)

**Script file:** `09_GO_extra_plots.R`


**What it does:** Re-runs the same GO enrichment as Step 8 and produces extra visualization types (barplot, cnetplot, emapplot) for both up- and down-regulated gene sets.


**Required input file(s):**

- `results/deseq2_results_annotated.tsv`


**Output file(s) produced:**

- `plots/GO_upregulated_barplot.png/.pdf`

- `plots/GO_downregulated_barplot.png/.pdf`

- `plots/GO_upregulated_cnetplot.png/.pdf`

- `plots/GO_downregulated_cnetplot.png/.pdf`

- `plots/GO_upregulated_emapplot.png/.pdf`

- `plots/GO_downregulated_emapplot.png/.pdf`



### Install packages needed for this script (Ubuntu)


```bash
sudo Rscript -e '
install.packages(c("ggplot2"), repos = "https://cloud.r-project.org")
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot"), update = FALSE, ask = FALSE)
'
```


### Full script: `09_GO_extra_plots.R`


```r
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
```


---

## Step 10 — KEGG Pathway Enrichment

**Script file:** `10_KEGG_enrichment.R`


**What it does:** Runs KEGG pathway over-representation analysis (clusterProfiler::enrichKEGG) on upregulated and downregulated genes, converts Entrez IDs to readable symbols, and generates tables plus dotplots, barplots, cnetplots, and enrichment maps.


**Required input file(s):**

- `results/deseq2_results_annotated.tsv`


**Output file(s) produced:**

- `results/KEGG_upregulated.tsv`

- `results/KEGG_downregulated.tsv`

- `plots/KEGG_upregulated_*.png/.pdf`

- `plots/KEGG_downregulated_*.png/.pdf`



> **Note:** enrichKEGG() queries the KEGG REST API live over the internet, so an active internet connection is required when this script runs.


### Install packages needed for this script (Ubuntu)


```bash
sudo Rscript -e '
install.packages(c("ggplot2"), repos = "https://cloud.r-project.org")
BiocManager::install(c("clusterProfiler", "enrichplot", "org.Hs.eg.db"), update = FALSE, ask = FALSE)
'
```


### Full script: `10_KEGG_enrichment.R`


```r
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
```


---

## Step 11 — KEGG Pathview Maps

**Script file:** `11_KEGG_pathview_maps.R`


**What it does:** Maps log2 fold-change values directly onto native KEGG pathway diagrams (PI3K-Akt signaling, ECM-receptor interaction, focal adhesion, proteoglycans in cancer, HPV infection) using the pathview package, coloring genes by up/down direction.


**Required input file(s):**

- `results/deseq2_results_annotated.tsv`


**Output file(s) produced:**

- `plots/pathview/*.png (one annotated KEGG diagram per pathway)`



> **Note:** pathview downloads each KEGG pathway diagram from the KEGG server on first use, so an internet connection is required. The script changes the R working directory to plots/pathview/ before running (pathview always writes to the current working directory), so run it from your project root and do not run other scripts in the same R session afterward without resetting the working directory.


### Install packages needed for this script (Ubuntu)


```bash
sudo Rscript -e '
BiocManager::install(c("pathview", "org.Hs.eg.db", "AnnotationDbi"), update = FALSE, ask = FALSE)
'
```


### Full script: `11_KEGG_pathview_maps.R`


```r
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
```


---

## Step 12 — STRING Protein-Protein Interaction (PPI) Mapping & Hub-Gene Ranking

**Script file:** `12_hub_gene_STRING_analysis.R`


**What it does:** Takes the top 300 most significant DEGs, maps them onto the STRING protein interaction database (v12, human, confidence score >= 400), retrieves their interaction network, and ranks genes by node degree (number of connections) to prioritize hub genes.


**Required input file(s):**

- `results/deseq2_results_annotated.tsv`


**Output file(s) produced:**

- `results/PPI_nodes_top300_DEGs.tsv`

- `results/PPI_edges_top300_DEGs.tsv`

- `results/hub_genes_ranked.tsv`

- `results/top30_hub_genes.tsv`



> **Note:** STRINGdb downloads the STRING interaction database for human (~9606) on first run and caches it in stringdb_cache/. This download is several hundred MB and requires internet access; subsequent runs reuse the cache.


### Install packages needed for this script (Ubuntu)


```bash
sudo Rscript -e '
install.packages(c("dplyr"), repos = "https://cloud.r-project.org")
BiocManager::install(c("STRINGdb"), update = FALSE, ask = FALSE)
'
```


### Full script: `12_hub_gene_STRING_analysis.R`


```r
library(STRINGdb)
library(dplyr)

options(timeout = 600)

res <- read.delim(
  "results/deseq2_results_annotated.tsv",
  check.names = FALSE
)

res <- res[!is.na(res$padj), ]
res <- res[!is.na(res$GeneSymbol), ]

sig <- res[
  res$padj < 0.05 &
    abs(res$log2FoldChange) >= 1,
]

top_genes <- sig %>%
  arrange(padj) %>%
  head(300)

dir.create("stringdb_cache", showWarnings = FALSE)

string_db <- STRINGdb$new(
  version = "12",
  species = 9606,
  score_threshold = 400,
  input_directory = "stringdb_cache"
)

mapped <- string_db$map(
  top_genes,
  "GeneSymbol",
  removeUnmappedRows = TRUE
)

interactions <- string_db$get_interactions(mapped$STRING_id)

interactions_filtered <- interactions[
  interactions$from %in% mapped$STRING_id &
    interactions$to %in% mapped$STRING_id,
]

nodes <- mapped[, c(
  "STRING_id",
  "EntrezID",
  "GeneSymbol",
  "GeneName",
  "log2FoldChange",
  "padj"
)]

degree_table <- data.frame(
  STRING_id = c(interactions_filtered$from, interactions_filtered$to)
) %>%
  group_by(STRING_id) %>%
  summarise(Degree = n(), .groups = "drop") %>%
  arrange(desc(Degree))

hub_genes <- merge(
  degree_table,
  nodes,
  by = "STRING_id"
)

hub_genes <- hub_genes[order(-hub_genes$Degree), ]

write.table(
  nodes,
  "results/PPI_nodes_top300_DEGs.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  interactions_filtered,
  "results/PPI_edges_top300_DEGs.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  hub_genes,
  "results/hub_genes_ranked.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  head(hub_genes, 30),
  "results/top30_hub_genes.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("PPI and hub gene analysis completed successfully.\n")
cat("Significant DEGs:", nrow(sig), "\n")
cat("Top genes used:", nrow(top_genes), "\n")
cat("Mapped genes:", nrow(mapped), "\n")
cat("Interactions:", nrow(interactions_filtered), "\n")
cat("Top hub genes saved in results/top30_hub_genes.tsv\n")
```


---

## Step 13 — PPI Network Visualization

**Script file:** `13_STRING_network_plot.R`


**What it does:** Loads the node and edge tables produced in Step 12, builds an igraph network object, and draws a force-directed (Fruchterman-Reingold) network plot where node size reflects connectivity (degree) and node color reflects log2 fold-change.


**Required input file(s):**

- `results/PPI_nodes_top300_DEGs.tsv`

- `results/PPI_edges_top300_DEGs.tsv`


**Output file(s) produced:**

- `plots/STRING_network_top300.png`

- `plots/STRING_network_top300.pdf`



### Install packages needed for this script (Ubuntu)


```bash
# system libraries
sudo apt-get install -y libglpk-dev libgmp-dev
```


```bash
sudo Rscript -e '
install.packages(c("igraph", "ggraph", "ggplot2"), repos = "https://cloud.r-project.org")
'
```


### Full script: `13_STRING_network_plot.R`


```r
library(igraph)
library(ggraph)
library(ggplot2)

nodes <- read.delim("results/PPI_nodes_top300_DEGs.tsv", check.names = FALSE)
edges <- read.delim("results/PPI_edges_top300_DEGs.tsv", check.names = FALSE)

network <- graph_from_data_frame(
  d = edges[, c("from", "to")],
  vertices = nodes,
  directed = FALSE
)

V(network)$degree <- degree(network)
V(network)$logFC <- as.numeric(V(network)$log2FoldChange)
V(network)$gene <- V(network)$GeneSymbol

png("plots/STRING_network_top300.png",
    width = 2600,
    height = 2400,
    res = 300)

ggraph(network, layout = "fr") +
  geom_edge_link(alpha = 0.25, colour = "grey70") +
  geom_node_point(aes(size = degree, color = logFC)) +
  geom_node_text(aes(label = gene), repel = TRUE, size = 3) +
  scale_color_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    name = "log2FC"
  ) +
  scale_size_continuous(name = "Degree") +
  theme_void() +
  ggtitle("Top 300 DEG Protein-Protein Interaction Network")

dev.off()

pdf("plots/STRING_network_top300.pdf",
    width = 12,
    height = 10)

ggraph(network, layout = "fr") +
  geom_edge_link(alpha = 0.25, colour = "grey70") +
  geom_node_point(aes(size = degree, color = logFC)) +
  geom_node_text(aes(label = gene), repel = TRUE, size = 3) +
  scale_color_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    name = "log2FC"
  ) +
  scale_size_continuous(name = "Degree") +
  theme_void() +
  ggtitle("Top 300 DEG Protein-Protein Interaction Network")

dev.off()

cat("STRING network visualization completed successfully.\n")
```


---

## Step 14 — Top 30 Hub Gene Expression Heatmap

**Script file:** `14_hub_gene_heatmap.R`


**What it does:** Pulls VST expression for the top 30 hub genes identified in Step 12 and draws a row-scaled (z-score) heatmap across all samples, annotated by tumor/control condition, using Ward's clustering method.


**Required input file(s):**

- `data/vst_counts.tsv`

- `results/top30_hub_genes.tsv`

- `data/sample_info.tsv`


**Output file(s) produced:**

- `plots/Hub_Genes_Heatmap.png`

- `plots/Hub_Genes_Heatmap.pdf`



### Install packages needed for this script (Ubuntu)


```bash
sudo Rscript -e '
install.packages(c("pheatmap", "RColorBrewer"), repos = "https://cloud.r-project.org")
'
```


### Full script: `14_hub_gene_heatmap.R`


```r
library(pheatmap)
library(RColorBrewer)

expr <- read.delim(
  "data/vst_counts.tsv",
  row.names = 1,
  check.names = FALSE
)

hub <- read.delim(
  "results/top30_hub_genes.tsv",
  check.names = FALSE
)

sample_info <- read.delim(
  "data/sample_info.tsv",
  check.names = FALSE
)

common_genes <- intersect(
  rownames(expr),
  as.character(hub$EntrezID)
)

heatmap_matrix <- expr[common_genes, ]

symbol_map <- hub$GeneSymbol
names(symbol_map) <- as.character(hub$EntrezID)

rownames(heatmap_matrix) <- symbol_map[rownames(heatmap_matrix)]

heatmap_matrix <- heatmap_matrix[complete.cases(heatmap_matrix), ]

annotation <- data.frame(
  Condition = sample_info$condition
)

rownames(annotation) <- sample_info$sample

annotation <- annotation[colnames(heatmap_matrix), , drop = FALSE]

annotation_colors <- list(
  Condition = c(
    Tumor = "#F8766D",
    Control = "#00BFC4"
  )
)

png("plots/Hub_Genes_Heatmap.png",
    width = 2400,
    height = 3000,
    res = 300)

pheatmap(
  heatmap_matrix,
  scale = "row",
  color = colorRampPalette(
    c("navy", "black", "orange", "yellow")
  )(100),
  annotation_col = annotation,
  annotation_colors = annotation_colors,
  fontsize_row = 9,
  fontsize_col = 7,
  clustering_method = "ward.D2",
  border_color = NA,
  main = "Top 30 Hub Gene Expression Heatmap"
)

dev.off()

pdf("plots/Hub_Genes_Heatmap.pdf",
    width = 10,
    height = 12)

pheatmap(
  heatmap_matrix,
  scale = "row",
  color = colorRampPalette(
    c("navy", "black", "orange", "yellow")
  )(100),
  annotation_col = annotation,
  annotation_colors = annotation_colors,
  fontsize_row = 9,
  fontsize_col = 7,
  clustering_method = "ward.D2",
  border_color = NA,
  main = "Top 30 Hub Gene Expression Heatmap"
)

dev.off()

cat("Hub gene heatmap completed successfully.\n")
cat("Hub genes plotted:", nrow(heatmap_matrix), "\n")
```


---

## Step 15 — Export Candidate Gene List for Survival Analysis

**Script file:** `15_survival_candidate_genes.R`


**What it does:** Takes the ranked hub-gene table and writes out a clean, sorted candidate list (symbol, name, degree, log2FC, padj) intended to be looked up manually in GEPIA2 / TCGA-HNSC for overall-survival curves.


**Required input file(s):**

- `results/top30_hub_genes.tsv`


**Output file(s) produced:**

- `results/survival_candidate_genes.tsv`



> **Note:** This script does not run survival analysis itself — GEPIA2 (http://gepia2.cancer-pku.cn) is a web tool, not an R package, so the resulting gene list is meant to be entered there manually or via its API separately.


### Install packages needed for this script (Ubuntu)


```bash
sudo Rscript -e '
install.packages(c("dplyr"), repos = "https://cloud.r-project.org")
'
```


### Full script: `15_survival_candidate_genes.R`


```r
library(dplyr)

hub <- read.delim(
  "results/top30_hub_genes.tsv",
  check.names = FALSE
)

candidate_genes <- hub %>%
  arrange(desc(Degree)) %>%
  select(
    GeneSymbol,
    GeneName,
    Degree,
    log2FoldChange,
    padj
  )

write.table(
  candidate_genes,
  "results/survival_candidate_genes.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("Survival candidate gene list exported successfully.\n")
```


---

## Step 16 — Boxplots of Top 10 Hub Genes

**Script file:** `16_hub_gene_boxplots.R`


**What it does:** For each of the top 10 hub genes, draws a boxplot (with jittered points) of VST expression split by condition (Control vs Tumor), saving one PNG per gene.


**Required input file(s):**

- `data/vst_counts.tsv`

- `data/sample_info.tsv`

- `results/top30_hub_genes.tsv`


**Output file(s) produced:**

- `plots/hub_boxplots/<GeneSymbol>_boxplot.png (x10)`



### Install packages needed for this script (Ubuntu)


```bash
sudo Rscript -e '
install.packages(c("ggplot2", "dplyr", "tidyr"), repos = "https://cloud.r-project.org")
'
```


### Full script: `16_hub_gene_boxplots.R`


```r
library(ggplot2)
library(dplyr)
library(tidyr)

expr <- read.delim(
  "data/vst_counts.tsv",
  row.names = 1,
  check.names = FALSE
)

sample_info <- read.delim(
  "data/sample_info.tsv",
  check.names = FALSE
)

hub <- read.delim(
  "results/top30_hub_genes.tsv",
  check.names = FALSE
)

top_genes <- hub %>%
  arrange(desc(Degree)) %>%
  head(10)

dir.create("plots/hub_boxplots", showWarnings = FALSE)

for (i in 1:nrow(top_genes)) {
  
  entrez <- as.character(top_genes$EntrezID[i])
  symbol <- top_genes$GeneSymbol[i]
  
  if (entrez %in% rownames(expr)) {
    
    df <- data.frame(
      sample = colnames(expr),
      expression = as.numeric(expr[entrez, ])
    )
    
    df <- merge(df, sample_info, by = "sample")
    
    p <- ggplot(df, aes(x = condition, y = expression, fill = condition)) +
      geom_boxplot(outlier.shape = NA, alpha = 0.7) +
      geom_jitter(width = 0.2, size = 2, alpha = 0.8) +
      theme_bw() +
      labs(
        title = paste0(symbol, " Expression"),
        x = "",
        y = "VST normalized expression"
      ) +
      scale_fill_manual(values = c("Control" = "#00BFC4", "Tumor" = "#F8766D")) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "none"
      )
    
    ggsave(
      paste0("plots/hub_boxplots/", symbol, "_boxplot.png"),
      p,
      width = 5,
      height = 5,
      dpi = 300
    )
  }
}

cat("Hub gene boxplots completed successfully.\n")
```


---

## Step 17 — Biomarker Priority Table

**Script file:** `17_biomarker_priority_table.R`


**What it does:** Scores every hub gene by Degree x |log2FoldChange| and manually classifies genes into biological categories (Invasion/ECM remodeling, Metabolic reprogramming, Carcinogen/drug metabolism, Other) to produce a ranked candidate-biomarker table.


**Required input file(s):**

- `results/top30_hub_genes.tsv`


**Output file(s) produced:**

- `results/biomarker_priority_table.tsv`



### Install packages needed for this script (Ubuntu)


```bash
sudo Rscript -e '
install.packages(c("dplyr"), repos = "https://cloud.r-project.org")
'
```


### Full script: `17_biomarker_priority_table.R`


```r
library(dplyr)

hub <- read.delim(
  "results/top30_hub_genes.tsv",
  check.names = FALSE
)

biomarkers <- hub %>%
  mutate(
    Regulation = ifelse(log2FoldChange > 0, "Upregulated in Tumor", "Downregulated in Tumor"),
    PriorityScore = Degree * abs(log2FoldChange),
    Category = case_when(
      GeneSymbol %in% c("MMP9", "MMP14", "PLAU", "COL1A1") ~ "Invasion / ECM remodeling",
      GeneSymbol %in% c("PPARG", "ACACB", "GPD1", "LPL", "PCK1", "PLIN1", "SLC2A4") ~ "Metabolic reprogramming",
      GeneSymbol %in% c("ADH1B", "CYP3A4") ~ "Carcinogen / drug metabolism",
      TRUE ~ "Other cancer-associated hub"
    )
  ) %>%
  arrange(desc(PriorityScore))

write.table(
  biomarkers,
  "results/biomarker_priority_table.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("Biomarker priority table created successfully.\n")
```


---

## Step 18 — Publication-Quality Top Biomarker Heatmap

**Script file:** `18_top_biomarker_heatmap.R`


**What it does:** Builds a polished, presentation-ready heatmap restricted to the 13 named candidate biomarker genes (COL1A1, MMP14, PLAU, MMP9, PPARG, ACACB, GPD1, LPL, PCK1, PLIN1, SLC2A4, ADH1B, CYP3A4), row-scaled and ordered Control-then-Tumor.


**Required input file(s):**

- `data/vst_counts.tsv`

- `data/sample_info.tsv`

- `results/top30_hub_genes.tsv`


**Output file(s) produced:**

- `plots/top_biomarker_heatmap.png`

- `plots/top_biomarker_heatmap.pdf`



### Install packages needed for this script (Ubuntu)


```bash
sudo Rscript -e '
install.packages(c("pheatmap", "RColorBrewer"), repos = "https://cloud.r-project.org")
'
```


### Full script: `18_top_biomarker_heatmap.R`


```r
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
```


---

## Step 19 — Reactome Pathway Enrichment

**Script file:** `19_reactome_enrichment.R`


**What it does:** Runs Reactome pathway over-representation analysis (ReactomePA::enrichPathway) separately on upregulated and downregulated genes, producing result tables plus dotplots, barplots, enrichment maps, and cnetplots.


**Required input file(s):**

- `results/deseq2_results_annotated.tsv`


**Output file(s) produced:**

- `results/Reactome_upregulated.tsv`

- `results/Reactome_downregulated.tsv`

- `plots/Reactome_*_dotplot/barplot/emapplot/cnetplot.png/.pdf`



### Install packages needed for this script (Ubuntu)


```bash
sudo Rscript -e '
install.packages(c("ggplot2"), repos = "https://cloud.r-project.org")
BiocManager::install(c("ReactomePA", "clusterProfiler", "enrichplot", "org.Hs.eg.db", "reactome.db"), update = FALSE, ask = FALSE)
'
```


### Full script: `19_reactome_enrichment.R`


```r
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
```


---

## Step 20 — Hallmark Gene Set Enrichment Analysis (GSEA)

**Script file:** `20_GSEA_hallmark_analysis.R`


**What it does:** Builds a full ranked gene list (all genes by log2FC, not just significant ones) and runs pre-ranked GSEA against the MSigDB Hallmark gene sets, then produces dotplots, a ridgeplot, an NES barplot, and individual enrichment plots for key cancer-relevant hallmark pathways.


**Required input file(s):**

- `results/deseq2_results_annotated.tsv`


**Output file(s) produced:**

- `results/GSEA_Hallmark_results.tsv`

- `results/GSEA_Hallmark_activated_in_tumor.tsv`

- `results/GSEA_Hallmark_suppressed_in_tumor.tsv`

- `plots/GSEA_Hallmark_dotplot.png/.pdf`

- `plots/GSEA_Hallmark_ridgeplot.png`

- `plots/GSEA_Hallmark_NES_barplot.png`

- `plots/GSEA_enrichment_plots/*.png`



> **Note:** msigdbr() downloads MSigDB gene sets on first use and requires internet access.


### Install packages needed for this script (Ubuntu)


```bash
sudo Rscript -e '
install.packages(c("msigdbr", "dplyr", "ggplot2"), repos = "https://cloud.r-project.org")
BiocManager::install(c("clusterProfiler", "enrichplot"), update = FALSE, ask = FALSE)
'
```


### Full script: `20_GSEA_hallmark_analysis.R`


```r
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
  pvalueCutoff = 0.05,
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
```


---


## Quick-Reference: Full Run Order

Once every package above is installed, the entire pipeline can be executed
end-to-end from your `project/` directory with:

```bash
cd project/

Rscript scripts/01_pca_plot.R
Rscript scripts/02_sample_distance_heatmap.R
Rscript scripts/03_volcano_plot.R
Rscript scripts/04_convert_entrez_to_symbol.R
Rscript scripts/05_volcano_gene_symbols.R
Rscript scripts/06_top_DEG_heatmap.R
Rscript scripts/07_onco_tsg_heatmaps.R
Rscript scripts/08_GO_enrichment.R
Rscript scripts/09_GO_extra_plots.R
Rscript scripts/10_KEGG_enrichment.R
Rscript scripts/11_KEGG_pathview_maps.R
Rscript scripts/12_hub_gene_STRING_analysis.R
Rscript scripts/13_STRING_network_plot.R
Rscript scripts/14_hub_gene_heatmap.R
Rscript scripts/15_survival_candidate_genes.R
Rscript scripts/16_hub_gene_boxplots.R
Rscript scripts/17_biomarker_priority_table.R
Rscript scripts/18_top_biomarker_heatmap.R
Rscript scripts/19_reactome_enrichment.R
Rscript scripts/20_GSEA_hallmark_analysis.R
```

**Important ordering notes:**
- Scripts **04 must run before 05–20** — it produces `results/deseq2_results_annotated.tsv`,
  which every later script reads.
- Script **12 must run before 13, 14, 15, 16, 17, 18** — it produces the hub-gene and
  PPI tables those scripts depend on.
- Script **11** changes R's working directory to `plots/pathview/` internally. If you're
  running scripts interactively in a single R session (rather than one `Rscript` call
  per script), restart the session or `setwd()` back to your project root before
  running any script after it.
- Scripts **10, 11, 12, and 20** need an active internet connection (they query the
  KEGG REST API, download KEGG diagrams, download the STRING interaction database,
  and download MSigDB gene sets, respectively).

### After Script 15: manual survival step (not R)

`results/survival_candidate_genes.tsv` is meant to be used outside R. Open
**GEPIA2** (http://gepia2.cancer-pku.cn), choose the **TCGA-HNSC** dataset, and
enter each candidate gene symbol under "Survival Analysis" to generate
overall-survival Kaplan-Meier curves with hazard ratios and log-rank p-values,
as referenced in the manuscript (Figure 5D).
