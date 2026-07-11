#!/usr/bin/env python3
from pathlib import Path
import argparse
import os
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
RAW_DIR = ROOT / "raw_data"
DATA_DIR = ROOT / "data"
SETUP_MARKER = ROOT / ".oscc_setup_complete"

REQUIRED = {
    "deseq2": {
        "canonical": "deseq2_results.tsv",
        "label": "DESeq2 results table",
        "examples": [
            "~/deseq2_results.tsv",
            "/home/username/Downloads/deseq2_results.csv",
            "my_DESeq2_output.txt"
        ]
    },
    "vst": {
        "canonical": "vst_counts.tsv",
        "label": "VST counts table",
        "examples": [
            "~/vst_counts.tsv",
            "/home/username/Downloads/vst_matrix.csv",
            "vst_normalized_counts.txt"
        ]
    },
    "normalized": {
        "canonical": "normalized_counts.tsv",
        "label": "Normalized counts table",
        "examples": [
            "~/normalized_counts.tsv",
            "/home/username/Downloads/normalized_counts.csv",
            "norm_counts.txt"
        ]
    },
    "sample_info": {
        "canonical": "sample_info.tsv",
        "label": "Sample metadata table",
        "examples": [
            "~/sample_info.tsv",
            "/home/username/Downloads/sample_metadata.csv",
            "metadata.txt"
        ]
    }
}


def run(cmd, check=False):
    print("\n$ " + " ".join(str(x) for x in cmd))
    return subprocess.run(cmd, cwd=ROOT, check=check)


def ensure_r_environment():
    r_lib = Path.home() / "R" / "library"
    r_lib.mkdir(parents=True, exist_ok=True)

    renviron = Path.home() / ".Renviron"
    line = f"R_LIBS_USER={r_lib}\n"
    old = renviron.read_text() if renviron.exists() else ""
    if "R_LIBS_USER=" not in old:
        renviron.write_text(old.rstrip() + "\n" + line)

    r_dir = Path.home() / ".R"
    r_dir.mkdir(exist_ok=True)
    makevars = r_dir / "Makevars"
    needed = [
        "CXX11STD = -std=gnu++14",
        "CXX14STD = -std=gnu++14",
        "CXX17STD = -std=gnu++17"
    ]
    old = makevars.read_text() if makevars.exists() else ""
    new_lines = []
    for line in needed:
        if line not in old:
            new_lines.append(line)
    if new_lines:
        makevars.write_text(old.rstrip() + "\n" + "\n".join(new_lines) + "\n")

    print("[OK] R environment prepared.")


def install_python_packages():
    req = ROOT / "requirements.txt"
    if not req.exists():
        print("[WARN] requirements.txt not found. Skipping Python package install.")
        return 0

    return run([sys.executable, "-m", "pip", "install", "--user", "-r", str(req)]).returncode


def install_r_packages():
    installer = ROOT / "install_packages.R"
    if not installer.exists():
        print("[WARN] install_packages.R not found. Skipping R package install.")
        return 0

    code = run(["Rscript", str(installer)]).returncode

    if code != 0:
        print("[ERROR] R package installer failed.")
        return code

    check_code = subprocess.run(
        [
            "Rscript",
            "-e",
            "pkgs <- c('clusterProfiler','ReactomePA','enrichplot','treeio','ggtree','pathview','STRINGdb','fgsea','ggridges','org.Hs.eg.db'); missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly=TRUE)]; if(length(missing)>0){cat('Missing:', paste(missing, collapse=', '), '\\n'); quit(status=1)} else {cat('[OK] Critical R packages verified.\\n')}"
        ],
        cwd=ROOT
    ).returncode

    return check_code


def install_all():
    print("\n=== OSCC pipeline dependency setup ===")
    ensure_r_environment()
    py_code = install_python_packages()
    r_code = install_r_packages()

    if py_code == 0 and r_code == 0:
        SETUP_MARKER.write_text("setup completed\n")
        print("\n[OK] Dependency setup completed.")
        return 0

    print("\n[ERROR] Dependency setup did not complete successfully.")
    print("The analysis will not run safely until dependencies are fixed.")
    print("Run again:")
    print("python3 setup_pipeline.py --install-only")
    return 1


def resolve_path(user_text):
    raw = os.path.expandvars(user_text.strip())
    p = Path(raw).expanduser()

    candidates = [
        p,
        ROOT / raw,
        Path.home() / raw,
        ROOT / "data" / raw,
        ROOT / "raw_data" / raw
    ]

    for c in candidates:
        if c.exists() and c.is_file():
            return c.resolve()

    return None


def copy_input_file(src, canonical):
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    DATA_DIR.mkdir(parents=True, exist_ok=True)

    raw_dst = RAW_DIR / canonical
    data_dst = DATA_DIR / canonical

    try:
        # Remove old standardized destination files/links only.
        if data_dst.exists() or data_dst.is_symlink():
            data_dst.unlink()

        if raw_dst.exists() or raw_dst.is_symlink():
            raw_dst.unlink()

        # Store canonical raw file without wasting space when possible.
        # Hardlink = no duplicate storage on same filesystem.
        try:
            os.link(src, raw_dst)
            raw_mode = "hardlink"
        except Exception:
            shutil.copy2(src, raw_dst)
            raw_mode = "copy"

        # Make data/ point to raw_data/ without duplicating large files.
        try:
            os.link(raw_dst, data_dst)
            data_mode = "hardlink"
        except Exception:
            try:
                data_dst.symlink_to(raw_dst)
                data_mode = "symlink"
            except Exception:
                shutil.copy2(raw_dst, data_dst)
                data_mode = "copy"

        print(f"[OK] Added and standardized:")
        print(f"     raw_data/{canonical} ({raw_mode})")
        print(f"     data/{canonical} ({data_mode} to raw_data file)")
        return True

    except OSError as e:
        if getattr(e, "errno", None) == 28:
            print("\n[STORAGE ERROR] No space left on device.")
            print("The program will not close.")
            print("Please free space and try this file again.")
            print("Suggested cleanup commands:")
            print("  rm -rf ~/.cache/pip")
            print("  rm -rf results plots logs plots_sorted")
            print("  find data raw_data -type f ! -name .gitkeep -delete")
            try:
                if data_dst.exists() or data_dst.is_symlink():
                    data_dst.unlink()
                if raw_dst.exists() or raw_dst.is_symlink():
                    raw_dst.unlink()
            except Exception:
                pass
            return False
        raise


def prompt_file(key, cfg):
    canonical = cfg["canonical"]
    existing = DATA_DIR / canonical

    while True:
        print(f"\nRequired file: {cfg['label']}")
        print(f"Pipeline name: data/{canonical}")
        print("Example file paths:")
        for ex in cfg["examples"]:
            print(f"  - {ex}")

        if existing.exists():
            print(f"[INFO] Existing file found: data/{canonical}")
            answer = input("Press Enter to keep it, or paste a new file path: ").strip()
            if answer == "":
                print(f"[OK] Keeping existing data/{canonical}")
                return
        else:
            answer = input("Paste file path here, or type help/quit: ").strip()

        if answer.lower() in {"quit", "q", "exit"}:
            print("Setup cancelled.")
            sys.exit(0)

        if answer.lower() in {"help", "h", "?"}:
            print("\nSyntax examples:")
            print(f"  ~/Downloads/{canonical}")
            print(f"  /home/yourname/{canonical}")
            print(f"  {canonical}")
            continue

        src = resolve_path(answer)
        if src is None:
            print("\n[SYNTAX ERROR] File not found.")
            print("Please paste the full path, for example:")
            print(f"  /home/{Path.home().name}/{canonical}")
            print("The program will not close; try again.")
            continue

        ok = copy_input_file(src, canonical)
        if ok:
            return
        continue


def prepare_data():
    print("\n=== Prepare raw input files ===")
    print("This will copy your original files into raw_data/ and standardized files into data/.")
    print("The pipeline needs these exact standardized names:")

    for cfg in REQUIRED.values():
        print(f"  - data/{cfg['canonical']}")

    for key, cfg in REQUIRED.items():
        prompt_file(key, cfg)

    return verify_data()


def verify_data():
    print("\n=== Checking required input files ===")
    missing = []
    for cfg in REQUIRED.values():
        f = DATA_DIR / cfg["canonical"]
        if f.exists() and f.stat().st_size > 0:
            print(f"[OK] data/{cfg['canonical']}")
        else:
            print(f"[MISSING] data/{cfg['canonical']}")
            missing.append(cfg["canonical"])

    if missing:
        print("\n[WARN] Missing required files:")
        for m in missing:
            print(f"  - {m}")
        print("\nRun:")
        print("python3 setup_pipeline.py --prepare-data")
        return 1

    print("[OK] All required input files are ready.")
    return 0


def run_pipeline():
    code = verify_data()
    if code != 0:
        return code

    return run([
        sys.executable,
        "oscc_pipeline.py",
        "run",
        "--from",
        "1",
        "--to",
        "20",
        "--force",
        "--keep-going"
    ]).returncode


def interactive_setup():
    print("\nOSCC RNA-seq Pipeline Setup")
    print("===========================")

    while True:
        print("\nChoose an option:")
        print("1. Install/check all Python and R packages")
        print("2. Prepare raw input files")
        print("3. Verify input files")
        print("4. Run full analysis")
        print("5. Exit")

        choice = input("Type 1, 2, 3, 4, or 5: ").strip()

        if choice == "1":
            install_all()
        elif choice == "2":
            prepare_data()
        elif choice == "3":
            verify_data()
        elif choice == "4":
            run_pipeline()
        elif choice == "5":
            print("Goodbye.")
            break
        else:
            print("\n[SYNTAX ERROR] Wrong option.")
            print("Example: type 1 to install packages, or type 2 to prepare raw files.")
            print("The program will not close; try again.")


def main():
    parser = argparse.ArgumentParser(
        description="Interactive setup for the OSCC post-DESeq2 RNA-seq pipeline."
    )
    parser.add_argument("--install-only", action="store_true", help="Install/check Python and R packages only.")
    parser.add_argument("--prepare-data", action="store_true", help="Interactively copy and rename raw input files.")
    parser.add_argument("--verify", action="store_true", help="Verify required input files.")
    parser.add_argument("--run", action="store_true", help="Run the full pipeline after verification.")
    parser.add_argument("--non-interactive", action="store_true", help="Do not open the menu.")

    args = parser.parse_args()

    if args.install_only:
        sys.exit(install_all())

    if args.prepare_data:
        sys.exit(prepare_data())

    if args.verify:
        sys.exit(verify_data())

    if args.run:
        sys.exit(run_pipeline())

    if args.non_interactive:
        sys.exit(0)

    interactive_setup()


if __name__ == "__main__":
    main()
