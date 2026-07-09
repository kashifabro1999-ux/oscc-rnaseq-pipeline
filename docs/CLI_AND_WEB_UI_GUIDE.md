# oscc_pipeline.py — CLI Orchestrator for the OSCC RNA-seq Pipeline

A single-file, dependency-free (pure standard library) Python 3 command-line
tool that drives your existing 20 R scripts (`01_pca_plot.R` ...
`20_GSEA_hallmark_analysis.R`). It does **not** reimplement the analysis —
it manages running it: checking your environment, installing R/system
packages, running scripts in the correct dependency order, logging every
run, and reporting what's done vs. pending.

## Requirements

- Python 3.7+ (already on virtually any Ubuntu system: `python3 --version`)
- R and Rscript installed separately (see the tutorial document, Section 0,
  for the exact Ubuntu install commands) — this tool calls `Rscript` for you,
  it does not install R itself.
- No pip packages needed — the tool only uses the Python standard library.

## Quick Start

```bash
# 1. Put oscc_pipeline.py anywhere, e.g. your project folder
cd my_oscc_project/

# 2. Create the folder structure (data/, results/, plots/, scripts/, logs/)
python3 oscc_pipeline.py setup

# 3. Copy your 20 numbered .R scripts into scripts/
#    and your data/vst_counts.tsv, data/sample_info.tsv, data/deseq2_results.tsv
#    into data/

# 4. Check R + package + file status
python3 oscc_pipeline.py check-env

# 5. Install anything missing (all steps at once)
python3 oscc_pipeline.py install-deps --all

# 6. Run the whole pipeline
python3 oscc_pipeline.py run --all

# 7. Check what finished / what's pending
python3 oscc_pipeline.py status
```

## Commands

| Command | Purpose |
|---|---|
| `setup` | Creates `data/`, `results/`, `plots/`, `scripts/`, `logs/` folders |
| `list` | Lists all 20 steps, their titles, dependencies, and internet requirements |
| `check-env` | Verifies `Rscript` is on PATH, checks every CRAN/Bioconductor package used across all 20 scripts is installed, and checks core data files exist |
| `install-deps` | Installs `apt` system libraries + CRAN + Bioconductor packages for selected step(s) |
| `run` | Runs selected step(s) via `Rscript`, checking inputs/step-dependencies first, logging output to `logs/` |
| `status` | Shows DONE / PENDING / BLOCKED for every step, based on whether each step's primary output file exists |

## Selecting which steps to act on

`install-deps` and `run` both accept the same selection flags:

```bash
--all                 # every step, 1 through 20
--step 12              # a single step
--from 8 --to 11       # an inclusive range
```

(default with no flags: all steps)

## Useful flags for `run`

- `--force` — re-run a step even if its output file already exists
  (by default, already-completed steps are skipped, so re-running the whole
  pipeline after a partial failure won't redo finished work)
- `--keep-going` — don't stop the whole pipeline if one step fails or is
  blocked by a missing dependency; log the failure and continue to the next step

## Useful flag for `install-deps`

- `--dry-run` — print the exact `apt-get` / `Rscript -e` commands that would
  run, without executing them (handy for reviewing before you `sudo` anything)

## How step dependencies work

The tool knows, from the pipeline structure, that:
- Steps 5–11, 19, 20 all require Step 4's output
  (`results/deseq2_results_annotated.tsv`) to exist first
- Steps 13–18 all require Step 12's output
  (`results/top30_hub_genes.tsv`) to exist first

If you try to `run` a step whose dependency hasn't produced its output yet,
the tool stops and tells you exactly which step to run first, rather than
launching an R script that would just fail on a missing file.

## Logs

Every run of every step writes a timestamped log file to `logs/`, e.g.
`logs/step12_20260709_143022.log`, containing the full stdout/stderr from
that R process — useful for debugging a failed step without having to
re-run it.

## What's next

This CLI tool is step one of a two-part plan. **Step two is done too** — see
`oscc_web.py` below for a local browser UI built on the exact same step
metadata, so the CLI and web tool never drift out of sync.

---

# oscc_web.py — Local Web Control Panel (Streamlit)

A browser-based control panel for the same pipeline, built with
[Streamlit](https://streamlit.io). It reuses `oscc_pipeline.py`'s `STEPS`
list directly (by importing it from the same folder), so both tools always
agree on inputs, outputs, dependencies, and package requirements.

## Setup

```bash
# 1. Keep oscc_web.py in the SAME folder as oscc_pipeline.py — it loads
#    oscc_pipeline.py directly from its own directory.

# 2. Install Streamlit (one-time)
pip install streamlit --break-system-packages

# 3. Launch
streamlit run oscc_web.py
```

Streamlit will print a local URL (typically `http://localhost:8501`) —
open it in your browser.

## What's in the UI

- **📊 Overview & Status** — a live table of all 20 steps with DONE / PENDING
  / BLOCKED status and a progress bar, based on whether each step's output
  file exists on disk.
- **🔎 Environment Check** — one-click check of your R installation, every
  CRAN/Bioconductor package the pipeline needs, and whether your core data
  files (`vst_counts.tsv`, `sample_info.tsv`, `deseq2_results.tsv`) exist.
- **📦 Install Dependencies** — select step(s) (or leave blank for all),
  and get the exact `apt-get` and `Rscript -e '...'` install commands.
  `apt` commands are shown for you to run yourself in a terminal (a web
  app can't supply your `sudo` password); R package installs can be
  triggered directly from the button if your R library path is
  user-writable.
- **▶️ Run Steps** — run all steps, a single step, or a numeric range,
  with the same dependency/skip/force logic as the CLI's `run` command.
  Each step's output streams into its own expandable panel, and a log
  file is saved to `logs/` exactly like the CLI does.
- **📜 Logs** — browse and download any past run's log file.
- **🖼️ Plot Gallery** — a searchable grid of every PNG the pipeline has
  produced so far in `plots/` (including subfolders like
  `plots/pathview/` and `plots/hub_boxplots/`), so you can review figures
  without leaving the browser.

## Notes

- The web UI and the CLI can be used interchangeably on the same project
  folder — they read and write the exact same `data/`, `results/`,
  `plots/`, and `logs/` structure, so switching between them mid-pipeline
  is safe.
- Like the CLI, this tool orchestrates your existing R scripts; it does
  not reimplement the analysis in Python.
- `sudo apt-get` commands are intentionally *not* auto-executed from the
  web UI (Streamlit has no way to prompt for your password) — copy them
  into a terminal yourself.

