# OSCC RNA-seq Integrative Analysis Pipeline

## Beginner Installation

New users should first read:

[INSTALLATION.md](INSTALLATION.md)

Quick start after installing system requirements:

    git clone https://github.com/kashifabro1999-ux/oscc-rnaseq-pipeline.git
    cd oscc-rnaseq-pipeline
    python3 pipeline_menu.py

Inside the menu:

    setup
    run all
    organize
    zip
    quit

Reproducible pipeline and tooling for an integrative whole-transcriptome
RNA-seq analysis of oral squamous cell carcinoma (OSCC), combining public
datasets from five geographic regions (Taiwan, USA, India, Germany, China).

This repository contains:

- **`scripts/`** — the 20 numbered R scripts that make up the full analysis:
  QC/PCA, differential expression (DESeq2), volcano plots, heatmaps,
  functional enrichment (GO, KEGG, Reactome, Hallmark GSEA), STRING
  protein-protein interaction analysis, hub-gene prioritization, and
  candidate biomarker visualization.
- **`oscc_pipeline.py`** — a pure-stdlib Python CLI that orchestrates the
  R scripts: checks your environment, installs required R/system
  packages, runs steps in the correct dependency order, and logs everything.
- **`oscc_web.py`** — a Streamlit web control panel (same functionality as
  the CLI, browser-based): live pipeline status, one-click environment
  check, dependency installer, step runner, log viewer, and plot gallery.
- **`Dockerfile`**, **`install_cran_packages.R`**, **`install_bioc_packages.R`**,
  **`requirements.txt`** — everything needed to build a container image with
  R + all required CRAN/Bioconductor packages baked in, ready to deploy
  (tested against Google Cloud Run).
- **`docs/OSCC_RNAseq_Analysis_Tutorial.md`** — a full step-by-step tutorial
  covering Ubuntu installation commands for every package, in front of
  every script, plus the upstream raw-read-processing stage (FastQC → fastp
  → HISAT2 → featureCounts → DESeq2) that produces this pipeline's inputs.
- **`docs/CLI_AND_WEB_UI_GUIDE.md`** — usage guide for `oscc_pipeline.py`
  and `oscc_web.py`.

## Quick start (local)

```bash
git clone https://github.com/kashifabro1999-ux/oscc-rnaseq-pipeline.git
cd oscc-rnaseq-pipeline

python3 oscc_pipeline.py setup
# place your data/vst_counts.tsv, data/sample_info.tsv, data/deseq2_results.tsv
python3 oscc_pipeline.py check-env
python3 oscc_pipeline.py install-deps --all
python3 oscc_pipeline.py run --all
python3 oscc_pipeline.py status
```

See `docs/OSCC_RNAseq_Analysis_Tutorial.md` for full details, or
`docs/CLI_AND_WEB_UI_GUIDE.md` for the web UI (`streamlit run oscc_web.py`).

## Deploying the web UI

The `Dockerfile` in this repo builds a container with R and every
CRAN/Bioconductor package the pipeline needs already installed. It has been
tested deploying to Google Cloud Run:

```bash
gcloud run deploy oscc-rnaseq-pipeline \
  --source . \
  --region us-central1 \
  --memory 4Gi \
  --cpu 2 \
  --timeout 900 \
  --max-instances 2 \
  --allow-unauthenticated
```

## Data availability

All RNA-seq data used with this pipeline is drawn from public repositories
(GEO/SRA). No proprietary or patient-identifiable data is included in this
repository.

## Citation / context

This pipeline supports the manuscript *"Integrative Whole-Transcriptome
Profiling of Oral Squamous Cell Carcinoma Using Public RNA-seq Datasets from
Five Geographic Regions Identifies Candidate Diagnostic and Prognostic
Biomarkers."* Candidate biomarkers identified (PLAU, MMP14, COL1A1, MMP9,
and metabolic/tissue-identity genes LPL, ACACB, SLC2A4, PPARG, GPD1, PLIN1,
PCK1, ADH1B) are discovery-level and require independent experimental
validation.

## License

MIT — see `LICENSE`.
