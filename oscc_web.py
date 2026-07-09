#!/usr/bin/env python3
"""
oscc_web.py — Local web control panel (Streamlit) for the OSCC RNA-seq pipeline.

Reuses the exact same step metadata (inputs, outputs, dependencies, package
requirements) defined in oscc_pipeline.py, so the CLI and the web UI never
drift out of sync. This file must live in the same folder as
oscc_pipeline.py.

Run with:
    pip install streamlit --break-system-packages   # one-time
    streamlit run oscc_web.py

Then open the URL Streamlit prints (usually http://localhost:8501).
"""

import datetime
import importlib.util
import shutil
import subprocess
import sys
from pathlib import Path

import streamlit as st

# ---------------------------------------------------------------------------
# Load oscc_pipeline.py from the same directory as this file, regardless of
# the working directory Streamlit was launched from, so STEPS/paths logic
# is defined in exactly one place.
# ---------------------------------------------------------------------------
THIS_DIR = Path(__file__).resolve().parent
_pipeline_path = THIS_DIR / "oscc_pipeline.py"
if not _pipeline_path.exists():
    st.error(f"Could not find oscc_pipeline.py next to this file ({THIS_DIR}). "
             "Place oscc_web.py in the same folder as oscc_pipeline.py.")
    st.stop()

_spec = importlib.util.spec_from_file_location("oscc_pipeline", _pipeline_path)
oscc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(oscc)

STEPS = oscc.STEPS
STEP_BY_NUM = oscc.STEP_BY_NUM
project_paths = oscc.project_paths

# ---------------------------------------------------------------------------
# Helpers (web-specific: capture output as text instead of streaming to a
# terminal, since that's what the CLI's run_shell() does)
# ---------------------------------------------------------------------------

def run_capture(cmd, cwd=None):
    """Run a shell command, capturing all stdout+stderr as text. Returns (rc, text)."""
    proc = subprocess.run(cmd, shell=True, cwd=cwd, capture_output=True, text=True)
    return proc.returncode, (proc.stdout or "") + (proc.stderr or "")

def step_status(paths, s):
    out_path = paths["root"] / s["primary_output"]
    if out_path.exists():
        return "DONE"
    missing_inputs = [i for i in s["inputs"] if not (paths["root"] / i).exists()]
    if missing_inputs:
        return "BLOCKED"
    return "PENDING"

def unmet_step_deps(paths, s):
    return [r for r in s["requires"] if not (paths["root"] / STEP_BY_NUM[r]["primary_output"]).exists()]

def missing_inputs(paths, s):
    return [i for i in s["inputs"] if not (paths["root"] / i).exists()]

STATUS_ICON = {"DONE": "✅", "PENDING": "🕒", "BLOCKED": "⛔"}

# ---------------------------------------------------------------------------
# Page setup
# ---------------------------------------------------------------------------
st.set_page_config(page_title="OSCC RNA-seq Pipeline", layout="wide")
st.title("🧬 OSCC RNA-seq Pipeline — Control Panel")
st.caption("Local web UI for the 20-script OSCC integrative transcriptomics pipeline. "
           "This UI runs the same Rscript calls as oscc_pipeline.py on the command line.")

# --- Sidebar: project directory + quick actions ---------------------------
if "project_dir" not in st.session_state:
    st.session_state.project_dir = str(Path.cwd())

st.sidebar.header("Project")
project_dir_input = st.sidebar.text_input("Project directory", st.session_state.project_dir)
st.session_state.project_dir = project_dir_input
paths = project_paths(Path(project_dir_input))

st.sidebar.markdown("---")
if st.sidebar.button("📁 Create folder structure (setup)"):
    for key in ["data", "results", "plots", "scripts", "logs"]:
        paths[key].mkdir(parents=True, exist_ok=True)
    st.sidebar.success("data/, results/, plots/, scripts/, logs/ created or confirmed.")

st.sidebar.markdown("---")
st.sidebar.caption(
    "Before running steps, place your 20 numbered `.R` scripts in `scripts/` "
    "and `vst_counts.tsv`, `sample_info.tsv`, `deseq2_results.tsv` in `data/`."
)

rscript_bin = shutil.which("Rscript")
if rscript_bin:
    st.sidebar.success(f"Rscript found: {rscript_bin}")
else:
    st.sidebar.error("Rscript not found on PATH. Install R before running any step.")

tabs = st.tabs([
    "📊 Overview & Status",
    "🔎 Environment Check",
    "📦 Install Dependencies",
    "▶️ Run Steps",
    "📜 Logs",
    "🖼️ Plot Gallery",
])

# ===========================================================================
# TAB 1: Overview & Status
# ===========================================================================
with tabs[0]:
    st.subheader("Pipeline Status")
    st.caption("Status is determined by whether each step's primary output file already exists.")

    rows = []
    for s in STEPS:
        status = step_status(paths, s)
        rows.append({
            "Step": s["num"],
            "Script": s["file"],
            "Title": s["title"],
            "Status": f"{STATUS_ICON[status]} {status}",
            "Requires step(s)": ", ".join(map(str, s["requires"])) if s["requires"] else "—",
            "Needs internet": "🌐" if s["needs_internet"] else "",
        })
    st.dataframe(rows, use_container_width=True, hide_index=True)

    done_count = sum(1 for s in STEPS if step_status(paths, s) == "DONE")
    st.progress(done_count / len(STEPS))
    st.caption(f"{done_count} / {len(STEPS)} steps completed.")

    if st.button("🔄 Refresh"):
        st.rerun()

# ===========================================================================
# TAB 2: Environment Check
# ===========================================================================
with tabs[1]:
    st.subheader("Environment Check")

    if st.button("Run environment check"):
        with st.spinner("Checking R, packages, and data files..."):
            if not rscript_bin:
                st.error("Rscript not found on PATH. Install R first.")
            else:
                rc, out = run_capture("Rscript -e 'cat(R.version.string)'")
                st.success(f"R version: {out.strip()}" if rc == 0 else "Could not read R version.")

                all_cran = sorted({p for s in STEPS for p in s["cran"]})
                all_bioc = sorted({p for s in STEPS for p in s["bioc"]})

                col1, col2 = st.columns(2)
                with col1:
                    st.markdown("**CRAN packages**")
                    for pkg in all_cran:
                        rc, _ = run_capture(
                            f'Rscript -e \'if (!requireNamespace("{pkg}", quietly=TRUE)) quit(status=1)\''
                        )
                        st.write(("✅ " if rc == 0 else "⚠️ MISSING: ") + pkg)
                with col2:
                    st.markdown("**Bioconductor packages**")
                    for pkg in all_bioc:
                        rc, _ = run_capture(
                            f'Rscript -e \'if (!requireNamespace("{pkg}", quietly=TRUE)) quit(status=1)\''
                        )
                        st.write(("✅ " if rc == 0 else "⚠️ MISSING: ") + pkg)

                st.markdown("**Core data files**")
                for f in ["data/vst_counts.tsv", "data/sample_info.tsv", "data/deseq2_results.tsv"]:
                    p = paths["root"] / f
                    st.write(("✅ " if p.exists() else "❌ MISSING: ") + f)

    st.info("Missing packages? Go to the **Install Dependencies** tab.")

# ===========================================================================
# TAB 3: Install Dependencies
# ===========================================================================
with tabs[2]:
    st.subheader("Install Dependencies")

    step_options = {f"{s['num']:>2} — {s['file']}": s["num"] for s in STEPS}
    selected_labels = st.multiselect(
        "Select step(s) to install packages for (leave empty = all steps)",
        list(step_options.keys()),
    )
    selected_nums = [step_options[l] for l in selected_labels] or [s["num"] for s in STEPS]
    selected_steps = [STEP_BY_NUM[n] for n in selected_nums]

    apt_pkgs = sorted({p for s in selected_steps for p in s["apt"]})
    cran_pkgs = sorted({p for s in selected_steps for p in s["cran"]})
    bioc_pkgs = sorted({p for s in selected_steps for p in s["bioc"]})

    if apt_pkgs:
        apt_cmd = "sudo apt-get install -y " + " ".join(apt_pkgs)
        st.markdown("**System (apt) packages** — run this yourself in a terminal "
                     "(a web app cannot supply your `sudo` password):")
        st.code(apt_cmd, language="bash")

    r_lines = []
    if cran_pkgs:
        pkg_list = ", ".join(f'"{p}"' for p in cran_pkgs)
        r_lines.append(f'install.packages(c({pkg_list}), repos = "https://cloud.r-project.org")')
    if bioc_pkgs:
        pkg_list = ", ".join(f'"{p}"' for p in bioc_pkgs)
        r_lines.append(
            'if (!requireNamespace("BiocManager", quietly=TRUE)) '
            'install.packages("BiocManager", repos="https://cloud.r-project.org"); '
            f'BiocManager::install(c({pkg_list}), update = FALSE, ask = FALSE)'
        )

    if r_lines:
        r_expr = "; ".join(r_lines)
        st.markdown("**R packages (CRAN + Bioconductor)**")
        st.code(f"Rscript -e '{r_expr}'", language="bash")

        st.warning(
            "Clicking below installs into whatever R library path the Streamlit "
            "process's user account can write to. If R was installed system-wide "
            "as root, run the command above yourself in a terminal with `sudo` instead."
        )
        if st.button("▶️ Run R package install now (no sudo)"):
            with st.spinner("Installing R packages — this can take several minutes..."):
                rc, out = run_capture(f"Rscript -e '{r_expr}'")
            if rc == 0:
                st.success("R package installation finished.")
            else:
                st.error(f"Installation exited with code {rc}.")
            st.code(out or "(no output captured)", language="text")
    else:
        st.info("No R packages required for the current selection.")

# ===========================================================================
# TAB 4: Run Steps
# ===========================================================================
with tabs[3]:
    st.subheader("Run Steps")

    run_mode = st.radio("Selection mode", ["All steps", "Single step", "Range"], horizontal=True)

    if run_mode == "Single step":
        chosen = st.selectbox("Step", [f"{s['num']:>2} — {s['file']} — {s['title']}" for s in STEPS])
        nums_to_run = [int(chosen.split("—")[0].strip())]
    elif run_mode == "Range":
        c1, c2 = st.columns(2)
        start = c1.number_input("From step", min_value=1, max_value=20, value=1)
        end = c2.number_input("To step", min_value=1, max_value=20, value=20)
        nums_to_run = list(range(int(start), int(end) + 1))
    else:
        nums_to_run = [s["num"] for s in STEPS]

    col_a, col_b = st.columns(2)
    force = col_a.checkbox("Force re-run (ignore existing output)", value=False)
    keep_going = col_b.checkbox("Keep going past failures/blocked steps", value=True)

    if st.button("▶️ Run selected step(s)", type="primary"):
        paths["logs"].mkdir(parents=True, exist_ok=True)
        for n in nums_to_run:
            s = STEP_BY_NUM.get(n)
            if not s:
                continue

            primary_out = paths["root"] / s["primary_output"]
            with st.expander(f"Step {s['num']} — {s['title']} [{s['file']}]", expanded=True):

                if primary_out.exists() and not force:
                    st.info(f"Output already exists ({s['primary_output']}) — skipped. "
                            "Check 'Force re-run' to redo it.")
                    continue

                unmet = unmet_step_deps(paths, s)
                if unmet:
                    st.error(f"Requires step(s) {unmet} to run first.")
                    if not keep_going:
                        break
                    continue

                miss = missing_inputs(paths, s)
                if miss:
                    st.error(f"Missing input file(s): {miss}")
                    if not keep_going:
                        break
                    continue

                script_path = paths["scripts"] / s["file"]
                if not script_path.exists():
                    st.error(f"Script not found: {script_path}")
                    if not keep_going:
                        break
                    continue

                if s["needs_internet"]:
                    st.warning("This step requires an internet connection.")
                if s.get("note"):
                    st.caption(f"Note: {s['note']}")

                with st.spinner(f"Running {s['file']}..."):
                    rc, output = run_capture(f'Rscript "{script_path}"', cwd=paths["root"])

                timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
                log_path = paths["logs"] / f"step{s['num']:02d}_{timestamp}.log"
                log_path.write_text(output)

                if rc == 0:
                    st.success(f"Completed successfully. Log saved to {log_path.name}")
                else:
                    st.error(f"Failed (exit code {rc}). Log saved to {log_path.name}")
                    if not keep_going:
                        st.code(output, language="text")
                        break

                st.code(output or "(no output)", language="text")

# ===========================================================================
# TAB 5: Logs
# ===========================================================================
with tabs[4]:
    st.subheader("Run Logs")
    log_dir = paths["logs"]
    if not log_dir.exists():
        st.info("No logs/ folder yet. Run a step first, or click 'setup' in the sidebar.")
    else:
        log_files = sorted(log_dir.glob("*.log"), key=lambda p: p.stat().st_mtime, reverse=True)
        if not log_files:
            st.info("No log files yet.")
        else:
            chosen_log = st.selectbox("Select a log file", [p.name for p in log_files])
            content = (log_dir / chosen_log).read_text()
            st.code(content, language="text")
            st.download_button("Download this log", content, file_name=chosen_log)

# ===========================================================================
# TAB 6: Plot Gallery
# ===========================================================================
with tabs[5]:
    st.subheader("Plot Gallery")
    plots_dir = paths["plots"]
    if not plots_dir.exists():
        st.info("No plots/ folder yet.")
    else:
        image_files = sorted(plots_dir.rglob("*.png"))
        if not image_files:
            st.info("No PNG plots found yet. Run some steps first.")
        else:
            search = st.text_input("Filter by filename contains…", "")
            filtered = [p for p in image_files if search.lower() in p.name.lower()] if search else image_files
            st.caption(f"Showing {len(filtered)} of {len(image_files)} plot(s).")

            cols = st.columns(3)
            for i, img_path in enumerate(filtered):
                with cols[i % 3]:
                    st.image(str(img_path), caption=str(img_path.relative_to(plots_dir)), use_container_width=True)
