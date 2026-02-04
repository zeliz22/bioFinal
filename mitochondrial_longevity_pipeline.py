#!/usr/bin/env python3
"""
Mitochondrial Genome Longevity Analysis Pipeline
Based on Farré et al. 2021 methodology adapted for mitochondrial genomes

This pipeline identifies convergent amino acid substitutions (CAAS) and genes
associated with lifespan variation across species using mitochondrial genomes.
"""

import os
import sys
from pathlib import Path

# Pipeline configuration
PIPELINE_DIR = Path(__file__).parent
SCRIPTS_DIR = PIPELINE_DIR / "scripts"
DATA_DIR = PIPELINE_DIR / "data"
RESULTS_DIR = PIPELINE_DIR / "results"

# Create directory structure
for directory in [DATA_DIR, RESULTS_DIR, SCRIPTS_DIR]:
    directory.mkdir(exist_ok=True)
    
(DATA_DIR / "fasta_files").mkdir(exist_ok=True)
(DATA_DIR / "alignments").mkdir(exist_ok=True)
(DATA_DIR / "phenotypes").mkdir(exist_ok=True)
(RESULTS_DIR / "caas_discovery").mkdir(exist_ok=True)
(RESULTS_DIR / "caas_validation").mkdir(exist_ok=True)
(RESULTS_DIR / "dnds_analysis").mkdir(exist_ok=True)
(RESULTS_DIR / "functional_analysis").mkdir(exist_ok=True)
(RESULTS_DIR / "figures").mkdir(exist_ok=True)

print("=" * 80)
print("MITOCHONDRIAL GENOME LONGEVITY ANALYSIS PIPELINE")
print("=" * 80)
print("\nPipeline structure created:")
print(f"  Data directory: {DATA_DIR}")
print(f"  Scripts directory: {SCRIPTS_DIR}")
print(f"  Results directory: {RESULTS_DIR}")
print("\n" + "=" * 80)
print("\nPIPELINE STEPS:")
print("=" * 80)
print("""
Step 1: Data Preparation
   - Place FASTA files in data/fasta_files/
   - Create phenotype file (CSV with species, max_lifespan, body_mass)
   
Step 2: Extract & Align Mitochondrial Genes
   - Extract 13 protein-coding genes from mitochondrial genomes
   - Perform multiple sequence alignment
   
Step 3: Calculate Longevity Quotient (LQ)
   - Compute LQ = observed_lifespan / expected_lifespan(body_mass)
   - Classify species into extreme and intermediate groups
   
Step 4: CAAS Discovery Phase
   - Identify convergent amino acid substitutions in extreme species
   - Statistical validation against random expectations
   
Step 5: CAAS Validation Phase
   - Phylogenetic ANOVA on intermediate species
   - Keep only validated positions
   
Step 6: dN/dS Analysis
   - Calculate evolutionary rates per gene
   - PGLS regression with longevity
   
Step 7: Functional Analysis
   - Protein stability predictions
   - Pathway enrichment
   - Visualization and reporting

Run individual scripts in order or use this master script.
""")
print("=" * 80)

def run_pipeline():
    """Execute the full pipeline"""
    print("\nTo run the pipeline, execute scripts in order:")
    print("1. python scripts/01_extract_mt_genes.py")
    print("2. python scripts/02_align_genes.py")
    print("3. python scripts/03_calculate_lq.py")
    print("4. python scripts/04_caas_discovery.py")
    print("5. python scripts/05_caas_validation.py")
    print("6. python scripts/06_dnds_analysis.py")
    print("7. python scripts/07_functional_analysis.py")
    print("\nOr use: python mitochondrial_longevity_pipeline.py --run-all")

if __name__ == "__main__":
    if "--run-all" in sys.argv:
        print("\nFull pipeline execution not yet implemented.")
        print("Please run scripts individually first to ensure each step works.")
    else:
        run_pipeline()
