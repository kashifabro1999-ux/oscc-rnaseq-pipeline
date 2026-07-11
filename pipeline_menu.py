#!/usr/bin/env python3
from pathlib import Path
import subprocess
import sys
import zipfile
import re

ROOT = Path(__file__).resolve().parent


def run(cmd):
    print("\n$ " + " ".join(str(x) for x in cmd))
    return subprocess.run(cmd, cwd=ROOT).returncode


def help_text():
    print("""
OSCC Pipeline Interactive Menu

Correct commands:
  setup              install packages + prepare input files
  install            install/check Python and R packages
  prepare            add raw files and auto-rename them
  verify             check required input files
  status             show pipeline status
  run all            run steps 1 to 20
  run 7              run only step 7
  run 8-12           run steps 8 to 12
  organize           create plots_sorted folders and biomarker collage
  zip                create final output ZIP
  help               show this help
  quit               exit

Examples:
  run all
  run 20
  run 1-20
""")


def make_zip():
    out = ROOT / "OSCC_real_pipeline_outputs_sorted_final.zip"
    include_dirs = ["results", "plots_sorted", "logs"]

    if out.exists():
        out.unlink()

    with zipfile.ZipFile(out, "w", compression=zipfile.ZIP_DEFLATED) as z:
        for d in include_dirs:
            folder = ROOT / d
            if not folder.exists():
                continue
            for f in folder.rglob("*"):
                if f.is_file():
                    z.write(f, f.relative_to(ROOT))

    print(f"[OK] ZIP created: {out.name}")
    print("Download in Cloud Shell with:")
    print(f"cloudshell download {out.name}")


def parse_run_command(text):
    text = text.strip().lower()
    parts = text.split()

    if len(parts) < 2:
        print("[SYNTAX ERROR] Missing run target.")
        print("Examples: run all | run 20 | run 8-12")
        return

    target = parts[1]

    if target == "all":
        run([sys.executable, "oscc_pipeline.py", "run", "--from", "1", "--to", "20", "--force", "--keep-going"])
        return

    if re.fullmatch(r"\d+", target):
        run([sys.executable, "oscc_pipeline.py", "run", "--step", target, "--force", "--keep-going"])
        return

    m = re.fullmatch(r"(\d+)-(\d+)", target)
    if m:
        start, end = m.group(1), m.group(2)
        run([sys.executable, "oscc_pipeline.py", "run", "--from", start, "--to", end, "--force", "--keep-going"])
        return

    print("[SYNTAX ERROR] Wrong run command.")
    print("Examples: run all | run 20 | run 8-12")


def main():
    print("OSCC RNA-seq Pipeline Interactive Menu")
    print("Type help to see commands.")

    while True:
        try:
            text = input("\noscc> ").strip()
        except KeyboardInterrupt:
            print("\nType quit to exit.")
            continue

        if not text:
            continue

        cmd = text.lower()

        if cmd in {"quit", "exit", "q"}:
            print("Goodbye.")
            break

        elif cmd in {"help", "h", "?"}:
            help_text()

        elif cmd == "setup":
            run([sys.executable, "setup_pipeline.py"])

        elif cmd == "install":
            run([sys.executable, "setup_pipeline.py", "--install-only"])

        elif cmd == "prepare":
            run([sys.executable, "setup_pipeline.py", "--prepare-data"])

        elif cmd == "verify":
            run([sys.executable, "setup_pipeline.py", "--verify"])

        elif cmd == "status":
            run([sys.executable, "oscc_pipeline.py", "status"])

        elif cmd.startswith("run"):
            parse_run_command(cmd)

        elif cmd == "organize":
            run([sys.executable, "organize_outputs.py"])

        elif cmd == "zip":
            make_zip()

        else:
            print("[SYNTAX ERROR] Unknown command.")
            print("Examples: setup | run all | run 20 | status | organize | zip | quit")
            print("The program will not close; type help for all commands.")


if __name__ == "__main__":
    main()
