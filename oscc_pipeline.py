#!/usr/bin/env python3
"""
oscc_pipeline.py — CLI orchestrator for the OSCC RNA-seq analysis pipeline.

Runs, checks, and installs dependencies for the 20 numbered R scripts
(01_pca_plot.R ... 20_GSEA_hallmark_analysis.R) that make up the OSCC
integrative transcriptomics workflow. This tool does not reimplement the
analysis — it drives the existing R scripts: verifying inputs exist,
installing required R/system packages, running scripts in dependency
order, logging output, and reporting pipeline status.

Usage:
    python3 oscc_pipeline.py setup
    python3 oscc_pipeline.py check-env
    python3 oscc_pipeline.py install-deps --all
    python3 oscc_pipeline.py install-deps --step 12
    python3 oscc_pipeline.py run --all
    python3 oscc_pipeline.py run --step 5
    python3 oscc_pipeline.py run --from 8 --to 11
    python3 oscc_pipeline.py status
    python3 oscc_pipeline.py list

Run `python3 oscc_pipeline.py --help` or `<command> --help` for full options.
"""

import argparse
import datetime
import json
import shutil
import subprocess
import sys
import textwrap
from pathlib import Path

# ---------------------------------------------------------------------------
# ANSI colors (no external deps — falls back gracefully on non-tty output)
# ---------------------------------------------------------------------------
class C:
    GREEN = "\033[92m"
    YELLOW = "\033[93m"
    RED = "\033[91m"
    CYAN = "\033[96m"
    BOLD = "\033[1m"
    END = "\033[0m"

def _c(text, color):
    if not sys.stdout.isatty():
        return text
    return f"{color}{text}{C.END}"

def ok(text):    print(_c("[OK]   ", C.GREEN) + text)
def warn(text):  print(_c("[WARN] ", C.YELLOW) + text)
def err(text):   print(_c("[FAIL] ", C.RED) + text)
def info(text):  print(_c("[INFO] ", C.CYAN) + text)
def head(text):  print("\n" + _c(text, C.BOLD))

# ---------------------------------------------------------------------------
# Pipeline definition — one entry per script. This is the single source of
# truth for filenames, required inputs/outputs, package dependencies, step
# ordering, and internet requirements.
# ---------------------------------------------------------------------------
STEPS = [
    dict(num=1, file="01_pca_plot.R", title="PCA plot",
         inputs=["data/vst_counts.tsv", "data/sample_info.tsv"],
         primary_output="plots/PCA_Tumor_vs_Control.png",
         requires=[], cran=["ggplot2"], bioc=[], apt=[], needs_internet=False),

    dict(num=2, file="02_sample_distance_heatmap.R", title="Sample distance heatmap",
         inputs=["data/vst_counts.tsv", "data/sample_info.tsv"],
         primary_output="plots/sample_distance_heatmap.png",
         requires=[], cran=["pheatmap", "RColorBrewer"], bioc=[], apt=[], needs_internet=False),

    dict(num=3, file="03_volcano_plot.R", title="Volcano plot (Entrez ID stage)",
         inputs=["data/deseq2_results.tsv"],
         primary_output="plots/volcano_plot.png",
         requires=[], cran=[], bioc=["EnhancedVolcano"], apt=[], needs_internet=False),

    dict(num=4, file="04_convert_entrez_to_symbol.R", title="Annotate DESeq2 results (Entrez -> Symbol)",
         inputs=["data/deseq2_results.tsv"],
         primary_output="results/deseq2_results_annotated.tsv",
         requires=[], cran=[], bioc=["AnnotationDbi", "org.Hs.eg.db"], apt=[], needs_internet=False),

    dict(num=5, file="05_volcano_gene_symbols.R", title="Volcano plot (gene symbol stage)",
         inputs=["results/deseq2_results_annotated.tsv"],
         primary_output="plots/volcano_gene_symbols.png",
         requires=[4], cran=[], bioc=["EnhancedVolcano"], apt=[], needs_internet=False),

    dict(num=6, file="06_top_DEG_heatmap.R", title="Top 50 DEG heatmap",
         inputs=["data/vst_counts.tsv", "results/deseq2_results_annotated.tsv"],
         primary_output="plots/top50_DEG_heatmap_ordered.png",
         requires=[4], cran=["pheatmap", "RColorBrewer"], bioc=[], apt=[], needs_internet=False),

    dict(num=7, file="07_onco_tsg_heatmaps.R", title="Tumor-promoting vs tumor-opposing heatmaps",
         inputs=["data/vst_counts.tsv", "results/deseq2_results_annotated.tsv"],
         primary_output="plots/tumor_promoting_genes_heatmap.png",
         requires=[4], cran=["pheatmap", "RColorBrewer"], bioc=[], apt=[], needs_internet=False),

    dict(num=8, file="08_GO_enrichment.R", title="GO Biological Process enrichment",
         inputs=["results/deseq2_results_annotated.tsv"],
         primary_output="results/GO_upregulated.tsv",
         requires=[4], cran=["ggplot2"], bioc=["clusterProfiler", "org.Hs.eg.db", "enrichplot"],
         apt=[], needs_internet=False),

    dict(num=9, file="09_GO_extra_plots.R", title="GO extra visualizations (bar/cnet/emap)",
         inputs=["results/deseq2_results_annotated.tsv"],
         primary_output="plots/GO_upregulated_barplot.png",
         requires=[4], cran=["ggplot2"], bioc=["clusterProfiler", "org.Hs.eg.db", "enrichplot"],
         apt=[], needs_internet=False),

    dict(num=10, file="10_KEGG_enrichment.R", title="KEGG pathway enrichment",
         inputs=["results/deseq2_results_annotated.tsv"],
         primary_output="results/KEGG_upregulated.tsv",
         requires=[4], cran=["ggplot2"], bioc=["clusterProfiler", "enrichplot", "org.Hs.eg.db"],
         apt=[], needs_internet=True),

    dict(num=11, file="11_KEGG_pathview_maps.R", title="KEGG Pathview maps",
         inputs=["results/deseq2_results_annotated.tsv"],
         primary_output="plots/pathview",
         requires=[4], cran=[], bioc=["pathview", "org.Hs.eg.db", "AnnotationDbi"],
         apt=[], needs_internet=True,
         note="Writes to plots/pathview/ and changes R's working directory internally "
              "(only within that R subprocess — safe to run from this tool)."),

    dict(num=12, file="12_hub_gene_STRING_analysis.R", title="STRING PPI mapping & hub-gene ranking",
         inputs=["results/deseq2_results_annotated.tsv"],
         primary_output="results/top30_hub_genes.tsv",
         requires=[4], cran=["dplyr"], bioc=["STRINGdb"], apt=[], needs_internet=True,
         note="Downloads and caches the STRING human interaction database (~hundreds of MB) on first run."),

    dict(num=13, file="13_STRING_network_plot.R", title="PPI network visualization",
         inputs=["results/PPI_nodes_top300_DEGs.tsv", "results/PPI_edges_top300_DEGs.tsv"],
         primary_output="plots/STRING_network_top300.png",
         requires=[12], cran=["igraph", "ggraph", "ggplot2"], bioc=[],
         apt=["libglpk-dev", "libgmp-dev"], needs_internet=False),

    dict(num=14, file="14_hub_gene_heatmap.R", title="Top 30 hub-gene heatmap",
         inputs=["data/vst_counts.tsv", "results/top30_hub_genes.tsv", "data/sample_info.tsv"],
         primary_output="plots/Hub_Genes_Heatmap.png",
         requires=[12], cran=["pheatmap", "RColorBrewer"], bioc=[], apt=[], needs_internet=False),

    dict(num=15, file="15_survival_candidate_genes.R", title="Export survival candidate gene list",
         inputs=["results/top30_hub_genes.tsv"],
         primary_output="results/survival_candidate_genes.tsv",
         requires=[12], cran=["dplyr"], bioc=[], apt=[], needs_internet=False,
         note="Output is meant for manual entry into GEPIA2 (http://gepia2.cancer-pku.cn) — a web tool, not run here."),

    dict(num=16, file="16_hub_gene_boxplots.R", title="Top 10 hub-gene boxplots",
         inputs=["data/vst_counts.tsv", "data/sample_info.tsv", "results/top30_hub_genes.tsv"],
         primary_output="plots/hub_boxplots",
         requires=[12], cran=["ggplot2", "dplyr", "tidyr"], bioc=[], apt=[], needs_internet=False),

    dict(num=17, file="17_biomarker_priority_table.R", title="Biomarker priority table",
         inputs=["results/top30_hub_genes.tsv"],
         primary_output="results/biomarker_priority_table.tsv",
         requires=[12], cran=["dplyr"], bioc=[], apt=[], needs_internet=False),

    dict(num=18, file="18_top_biomarker_heatmap.R", title="Publication-quality top biomarker heatmap",
         inputs=["data/vst_counts.tsv", "data/sample_info.tsv", "results/top30_hub_genes.tsv"],
         primary_output="plots/top_biomarker_heatmap.png",
         requires=[12], cran=["pheatmap", "RColorBrewer"], bioc=[], apt=[], needs_internet=False),

    dict(num=19, file="19_reactome_enrichment.R", title="Reactome pathway enrichment",
         inputs=["results/deseq2_results_annotated.tsv"],
         primary_output="results/Reactome_upregulated.tsv",
         requires=[4], cran=["ggplot2"],
         bioc=["ReactomePA", "clusterProfiler", "enrichplot", "org.Hs.eg.db", "reactome.db"],
         apt=[], needs_internet=False),

    dict(num=20, file="20_GSEA_hallmark_analysis.R", title="Hallmark GSEA",
         inputs=["results/deseq2_results_annotated.tsv"],
         primary_output="results/GSEA_Hallmark_results.tsv",
         requires=[4], cran=["msigdbr", "dplyr", "ggplot2"], bioc=["clusterProfiler", "enrichplot"],
         apt=[], needs_internet=True),
]

STEP_BY_NUM = {s["num"]: s for s in STEPS}

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def project_paths(project_dir: Path):
    return {
        "root": project_dir,
        "data": project_dir / "data",
        "results": project_dir / "results",
        "plots": project_dir / "plots",
        "scripts": project_dir / "scripts",
        "logs": project_dir / "logs",
    }

def resolve_steps(args):
    """Return a sorted list of step dicts based on --step / --from/--to / --all."""
    if args.all:
        nums = [s["num"] for s in STEPS]
    elif args.step:
        nums = [args.step]
    elif args.from_ or args.to:
        start = args.from_ or 1
        end = args.to or 20
        nums = [n for n in range(start, end + 1)]
    else:
        nums = [s["num"] for s in STEPS]  # default: everything
    return [STEP_BY_NUM[n] for n in nums if n in STEP_BY_NUM]

def run_shell(cmd, cwd=None):
    """Run a shell command, streaming output live. Returns exit code."""
    proc = subprocess.Popen(cmd, shell=True, cwd=cwd, stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT, text=True, bufsize=1)
    lines = []
    for line in proc.stdout:
        print(line, end="")
        lines.append(line)
    proc.wait()
    return proc.returncode, "".join(lines)

def rscript_available():
    return shutil.which("Rscript") is not None

# ---------------------------------------------------------------------------
# Commands
# ---------------------------------------------------------------------------

def cmd_setup(args):
    paths = project_paths(Path(args.project_dir))
    head(f"Setting up project directory at {paths['root'].resolve()}")
    for key in ["data", "results", "plots", "scripts", "logs"]:
        paths[key].mkdir(parents=True, exist_ok=True)
        ok(f"Created/confirmed: {paths[key]}")
    info("Next: place your 20 numbered .R scripts into the 'scripts/' folder,")
    info("and your data/vst_counts.tsv, data/sample_info.tsv, data/deseq2_results.tsv")
    info("into the 'data/' folder. Then run: python3 oscc_pipeline.py check-env")

def cmd_list(args):
    head("Pipeline steps (in required run order)")
    for s in STEPS:
        deps = f"requires step(s) {s['requires']}" if s["requires"] else "no step dependencies"
        net = " [needs internet]" if s["needs_internet"] else ""
        print(f"  {s['num']:>2}. {s['file']:<32} {s['title']:<48} ({deps}){net}")

def cmd_check_env(args):
    head("Checking environment")
    if not rscript_available():
        err("Rscript not found on PATH. Install R first (see tutorial Section 0).")
        return 1
    ok("Rscript found: " + shutil.which("Rscript"))

    version_code, version_out = run_shell("Rscript -e 'cat(R.version.string)'")
    if version_code == 0:
        ok(f"R version: {version_out.strip()}")

    all_cran = sorted({p for s in STEPS for p in s["cran"]})
    all_bioc = sorted({p for s in STEPS for p in s["bioc"]})

    head("Checking installed R packages")
    missing_cran, missing_bioc = [], []
    for pkg in all_cran:
        rc, _ = run_shell(f'Rscript -e \'if (!requireNamespace("{pkg}", quietly=TRUE)) quit(status=1)\'')
        (ok if rc == 0 else warn)(f"CRAN pkg '{pkg}': {'installed' if rc == 0 else 'MISSING'}")
        if rc != 0:
            missing_cran.append(pkg)
    for pkg in all_bioc:
        rc, _ = run_shell(f'Rscript -e \'if (!requireNamespace("{pkg}", quietly=TRUE)) quit(status=1)\'')
        (ok if rc == 0 else warn)(f"Bioc pkg '{pkg}': {'installed' if rc == 0 else 'MISSING'}")
        if rc != 0:
            missing_bioc.append(pkg)

    if missing_cran or missing_bioc:
        warn("Missing packages found. Run: python3 oscc_pipeline.py install-deps --all")
    else:
        ok("All required R packages are installed.")

    paths = project_paths(Path(args.project_dir))
    head("Checking project directory")
    for key in ["data", "results", "plots", "scripts"]:
        exists = paths[key].exists()
        (ok if exists else warn)(f"{paths[key]}: {'found' if exists else 'missing (run: setup)'}")

    head("Checking core data files")
    for f in ["data/vst_counts.tsv", "data/sample_info.tsv", "data/deseq2_results.tsv"]:
        p = paths["root"] / f
        (ok if p.exists() else err)(f"{f}: {'found' if p.exists() else 'MISSING'}")
    return 0

def cmd_install_deps(args):
    steps = resolve_steps(args)
    apt_pkgs = sorted({p for s in steps for p in s["apt"]})
    cran_pkgs = sorted({p for s in steps for p in s["cran"]})
    bioc_pkgs = sorted({p for s in steps for p in s["bioc"]})

    head(f"Installing dependencies for step(s): {[s['num'] for s in steps]}")

    if apt_pkgs:
        cmd = "sudo apt-get install -y " + " ".join(apt_pkgs)
        info("System packages: " + cmd)
        if not args.dry_run:
            rc, _ = run_shell(cmd)
            (ok if rc == 0 else err)("apt-get install " + ("succeeded" if rc == 0 else "failed"))
        else:
            print("  (dry-run, not executed)")

    r_lines = []
    if cran_pkgs:
        pkg_list = ", ".join(f'"{p}"' for p in cran_pkgs)
        r_lines.append(f'install.packages(c({pkg_list}), repos = "https://cloud.r-project.org")')
    if bioc_pkgs:
        pkg_list = ", ".join(f'"{p}"' for p in bioc_pkgs)
        r_lines.append(f'if (!requireNamespace("BiocManager", quietly=TRUE)) '
                        f'install.packages("BiocManager", repos="https://cloud.r-project.org"); '
                        f'BiocManager::install(c({pkg_list}), update = FALSE, ask = FALSE)')

    if r_lines:
        r_expr = "; ".join(r_lines)
        info("R packages: Rscript -e '" + r_expr + "'")
        if not args.dry_run:
            rc, _ = run_shell(f"Rscript -e '{r_expr}'")
            (ok if rc == 0 else err)("R package install " + ("succeeded" if rc == 0 else "failed"))
        else:
            print("  (dry-run, not executed)")

    if not apt_pkgs and not r_lines:
        ok("No new packages required for the selected step(s).")

def cmd_status(args):
    paths = project_paths(Path(args.project_dir))
    head("Pipeline status")
    for s in STEPS:
        out_path = paths["root"] / s["primary_output"]
        done = out_path.exists()
        missing_inputs = [i for i in s["inputs"] if not (paths["root"] / i).exists()]
        if done:
            ok(f"Step {s['num']:>2} ({s['file']}): DONE  -> {s['primary_output']}")
        elif missing_inputs:
            warn(f"Step {s['num']:>2} ({s['file']}): BLOCKED - missing input(s) {missing_inputs}")
        else:
            info(f"Step {s['num']:>2} ({s['file']}): PENDING")

def cmd_run(args):
    paths = project_paths(Path(args.project_dir))
    steps = resolve_steps(args)
    paths["logs"].mkdir(parents=True, exist_ok=True)

    head(f"Running step(s): {[s['num'] for s in steps]}")
    for s in steps:
        script_path = paths["scripts"] / s["file"]
        primary_out = paths["root"] / s["primary_output"]

        # Skip if already done and --force not given
        if primary_out.exists() and not args.force:
            info(f"Step {s['num']} ({s['file']}) already has output {s['primary_output']} — skipping "
                 f"(use --force to re-run).")
            continue

        # Check step dependencies
        unmet = [r for r in s["requires"] if not (paths["root"] / STEP_BY_NUM[r]["primary_output"]).exists()]
        if unmet:
            err(f"Step {s['num']} ({s['file']}) requires step(s) {unmet} to run first. Skipping.")
            if not args.keep_going:
                return 1
            continue

        # Check input files
        missing = [i for i in s["inputs"] if not (paths["root"] / i).exists()]
        if missing:
            err(f"Step {s['num']} ({s['file']}) missing input file(s): {missing}. Skipping.")
            if not args.keep_going:
                return 1
            continue

        if not script_path.exists():
            err(f"Script not found: {script_path}. Place it in the scripts/ folder. Skipping.")
            if not args.keep_going:
                return 1
            continue

        if s["needs_internet"]:
            warn(f"Step {s['num']} ({s['file']}) requires an internet connection.")
        if s.get("note"):
            info(f"Note: {s['note']}")

        head(f"Running Step {s['num']}: {s['title']}  [{s['file']}]")
        timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        log_path = paths["logs"] / f"step{s['num']:02d}_{timestamp}.log"

        cmd = f'Rscript "{script_path}"'
        rc, output = run_shell(cmd, cwd=paths["root"])
        log_path.write_text(output)

        if rc == 0:
            ok(f"Step {s['num']} completed successfully. Log: {log_path}")
        else:
            err(f"Step {s['num']} FAILED (exit code {rc}). See log: {log_path}")
            if not args.keep_going:
                return 1
    ok("Run finished.")
    return 0

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def build_parser():
    parser = argparse.ArgumentParser(
        prog="oscc_pipeline.py",
        description="CLI orchestrator for the OSCC RNA-seq 20-script R pipeline.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent(__doc__),
    )
    parser.add_argument("--project-dir", default=".", help="Project root directory (default: current dir)")
    sub = parser.add_subparsers(dest="command", required=True)

    sub.add_parser("setup", help="Create data/results/plots/scripts/logs folders").set_defaults(func=cmd_setup)
    sub.add_parser("list", help="List all 20 pipeline steps and their dependencies").set_defaults(func=cmd_list)
    sub.add_parser("check-env", help="Check R, Rscript, and required packages are installed").set_defaults(func=cmd_check_env)
    sub.add_parser("status", help="Show which steps have been completed").set_defaults(func=cmd_status)

    def add_step_selection(p):
        p.add_argument("--all", action="store_true", help="Apply to all 20 steps")
        p.add_argument("--step", type=int, help="Apply to a single step number (1-20)")
        p.add_argument("--from", dest="from_", type=int, help="Start of an inclusive step range")
        p.add_argument("--to", type=int, help="End of an inclusive step range")

    p_install = sub.add_parser("install-deps", help="Install apt/CRAN/Bioconductor packages for selected step(s)")
    add_step_selection(p_install)
    p_install.add_argument("--dry-run", action="store_true", help="Print install commands without executing them")
    p_install.set_defaults(func=cmd_install_deps)

    p_run = sub.add_parser("run", help="Run selected step(s) with the Rscript interpreter")
    add_step_selection(p_run)
    p_run.add_argument("--force", action="store_true", help="Re-run steps even if their output already exists")
    p_run.add_argument("--keep-going", action="store_true", help="Continue to next step even if one fails/is blocked")
    p_run.set_defaults(func=cmd_run)

    return parser

def main():
    parser = build_parser()
    args = parser.parse_args()
    rc = args.func(args)
    sys.exit(rc if isinstance(rc, int) else 0)

# OSCC_PIPELINE_MANAGED_MAIN_START
def _oscc_preflight_install_if_needed():
    try:
        import os
        import sys
        import subprocess
        from pathlib import Path

        if "run" not in sys.argv:
            return

        if os.environ.get("OSCC_SKIP_AUTO_INSTALL") == "1":
            return

        root = Path(__file__).resolve().parent
        marker = root / ".oscc_setup_complete"
        setup_script = root / "setup_pipeline.py"

        watched = [
            root / "requirements.txt",
            root / "install_packages.R",
            setup_script
        ]

        needs_setup = not marker.exists()
        if marker.exists():
            marker_time = marker.stat().st_mtime
            for f in watched:
                if f.exists() and f.stat().st_mtime > marker_time:
                    needs_setup = True
                    break

        if needs_setup and setup_script.exists():
            print("[INFO] First run or setup files changed. Installing/checking dependencies...")
            code = subprocess.call(
                [sys.executable, str(setup_script), "--install-only", "--non-interactive"],
                cwd=str(root)
            )
            if code != 0:
                print("[WARN] Dependency installer reported a problem.")
                print("[INFO] You can manually run: python3 setup_pipeline.py --install-only")

    except Exception as e:
        print(f"[WARN] Dependency preflight skipped: {e}")


def _oscc_organize_after_run():
    try:
        import sys
        import subprocess
        from pathlib import Path

        if "run" not in sys.argv:
            return

        root = Path(__file__).resolve().parent
        organizer = root / "organize_outputs.py"

        if organizer.exists():
            print("[INFO] Organizing plots into sorted folders...")
            subprocess.run(
                [sys.executable, str(organizer)],
                cwd=str(root),
                check=False
            )

    except Exception as e:
        print(f"[WARN] Output organization skipped: {e}")


if __name__ == "__main__":
    _oscc_preflight_install_if_needed()
    main()
    _oscc_organize_after_run()
# OSCC_PIPELINE_MANAGED_MAIN_END
