#!/usr/bin/env python3
"""
Master pipeline for mitochondrial longevity analysis using GenBank-based workflow.

Execution order:
  1. anage_LQ.py                  -> builds LQ groups and target lists in data/LQ/
  2. download_genbank.py          -> downloads mitochondrial GenBank files to data/genbank_files/
  3. extract_genes_genbank.py     -> extracts genes into data/extracted_genes/
  4. align_genes.py               -> builds alignments in data/alignments/
  5. caas_discovery_from_lists.py -> performs CAAS discovery/validation from alignments

Before each pipeline run, the following data subdirectories are removed so you
always start from a clean state:
  - data/alignments
  - data/LQ
  - data/extracted_genes
  - data/genbank_files
"""

import shutil
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent

DATA_DIRS_TO_CLEAN = [
    ROOT / "data" / "alignments",
    ROOT / "data" / "LQ",
    ROOT / "data" / "extracted_genes",
    ROOT / "data" / "genbank_files",
]

SCRIPTS_IN_ORDER = [
    ("anage_LQ.py", []),
    ("download_genbank.py", []),
    ("extract_genes_genbank.py", []),
    ("align_genes.py", []),
    ("caas_discovery_from_lists.py", []),
]


def clean_output_dirs() -> None:
    """Remove all output/data directories so each run starts fresh."""
    print("=" * 70)
    print("Cleaning data output directories")
    print("=" * 70)
    for d in DATA_DIRS_TO_CLEAN:
        if d.exists():
            print(f"  Removing {d}")
            shutil.rmtree(d, ignore_errors=True)
        else:
            print(f"  (skip) {d} does not exist")
    print()


def run_step(script_name: str, args: list[str]) -> None:
    """Run a single pipeline step as a subprocess."""
    script_path = ROOT / script_name
    if not script_path.exists():
        raise FileNotFoundError(f"Script not found: {script_path}")

    cmd = [sys.executable, str(script_path), *args]
    print("=" * 70)
    print(f"Running: {' '.join(cmd)}")
    print("=" * 70)

    # Run from project root so all scripts' relative paths (data/, results/) work
    subprocess.run(cmd, check=True, cwd=ROOT)
    print()


def main() -> None:
    clean_output_dirs()

    for script_name, args in SCRIPTS_IN_ORDER:
        run_step(script_name, args)

    print("=" * 70)
    print("Pipeline completed successfully.")
    print("=" * 70)


if __name__ == "__main__":
    main()


