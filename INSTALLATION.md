# Installation Guide for Beginners

This guide is for new users who want to run the OSCC RNA-seq post-DESeq2 analysis pipeline from GitHub.

The pipeline automatically installs/checks required Python and R analysis packages after cloning, but the system must already have basic tools installed.

## 1. Install basic system requirements

For Ubuntu, Google Cloud Shell, or WSL2 Ubuntu, run this first:

    sudo apt update

    sudo apt install -y git python3 python3-pip r-base r-base-dev build-essential gfortran make cmake libcurl4-openssl-dev libssl-dev libxml2-dev libpng-dev libjpeg-dev libtiff5-dev libcairo2-dev libxt-dev libfontconfig1-dev libfreetype6-dev libharfbuzz-dev libfribidi-dev

These packages install Git, Python, R, compilers, and system libraries needed for R/Bioconductor packages.

## 2. Clone the repository

    git clone https://github.com/kashifabro1999-ux/oscc-rnaseq-pipeline.git
    cd oscc-rnaseq-pipeline

## 3. Start the beginner-friendly menu

    python3 pipeline_menu.py

Inside the menu, type:

    setup
    run all
    organize
    zip
    quit

## 4. Required input files

The pipeline needs four files:

    deseq2_results.tsv
    vst_counts.tsv
    normalized_counts.tsv
    sample_info.tsv

Your original files may have different names. The setup menu will ask for their locations and automatically standardize them into:

    raw_data/
    data/

## 5. Output folders

After analysis, plots are organized into:

    plots_sorted/
    01_exploratory_PCA_distance
    02_volcano_and_DEG_heatmaps
    03_GO_enrichment
    04_KEGG_enrichment
    05_KEGG_pathview_maps
    06_Reactome_enrichment
    07_GSEA_Hallmark
    08_STRING_PPI_network
    09_hub_gene_heatmaps
    10_biomarker_boxplots
    11_top_biomarker_heatmap

Individual biomarker boxplots are saved in:

    plots_sorted/10_biomarker_boxplots/individual_boxplots/

The combined landscape collage is saved as:

    plots_sorted/10_biomarker_boxplots/all_biomarker_boxplots_landscape_collage.png

## 6. Create final ZIP

Inside the menu, type:

    zip

This creates:

    OSCC_real_pipeline_outputs_sorted_final.zip

In Google Cloud Shell, download it with:

    cloudshell download OSCC_real_pipeline_outputs_sorted_final.zip

## 7. Notes

Use Linux, Ubuntu, Google Cloud Shell, or WSL2 Ubuntu.

A stable internet connection is required because some steps download KEGG, STRING, MSigDB, and Bioconductor resources.

Recommended free disk space: 4–6 GB minimum.

Beginners should use:

    python3 pipeline_menu.py

Do not run the R scripts manually unless you know the workflow.
