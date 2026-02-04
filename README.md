# Mitochondrial Genome Longevity Analysis Pipeline

A computational pipeline for identifying convergent amino acid substitutions (CAAS) associated with longevity across species using mitochondrial genomes. Based on the methodology from Farré et al. 2021 (Mol Biol Evol).

## Overview

This pipeline analyzes the 13 protein-coding genes in mitochondrial genomes to discover amino acid changes that distinguish long-lived from short-lived species, correcting for body mass effects.

### Mitochondrial Protein-Coding Genes Analyzed
- **ATP synthase**: ATP6, ATP8
- **Cytochrome c oxidase**: COX1, COX2, COX3
- **Cytochrome b**: CYTB
- **NADH dehydrogenase**: ND1, ND2, ND3, ND4, ND4L, ND5, ND6

## Installation

### Requirements
```bash
# Python packages
pip install biopython pandas numpy scipy matplotlib seaborn statsmodels

# Alignment tools (choose one)
# MAFFT (recommended)
sudo apt-get install mafft  # Ubuntu/Debian
brew install mafft           # macOS
conda install -c bioconda mafft  # Conda

# OR MUSCLE
sudo apt-get install muscle  # Ubuntu/Debian
```

## Directory Structure

```
mitochondrial_longevity_pipeline/
├── data/
│   ├── fasta_files/          # Place your mitochondrial genome FASTA files here
│   ├── extracted_genes/      # Extracted genes (auto-generated)
│   ├── alignments/           # Gene alignments (auto-generated)
│   └── phenotypes/
│       └── phenotype_data.csv  # Species metadata (YOU MUST CREATE THIS)
├── scripts/
│   ├── 01_extract_mt_genes.py
│   ├── 02_align_genes.py
│   ├── 03_calculate_lq.py
│   ├── 04_caas_discovery.py
│   └── 05_caas_validation.py
├── results/
│   ├── caas_discovery/       # Discovery results (auto-generated)
│   └── caas_validation/      # Validation results (auto-generated)
└── README.md
```

## Step-by-Step Usage

### Step 1: Prepare Your Data

#### 1.1 Mitochondrial Genome Files
Place your mitochondrial genome FASTA files in `data/fasta_files/`.

**File naming convention**: Use species identifiers as filenames (e.g., `homo_sapiens.fasta`, `mus_musculus.fasta`)

**Supported formats**:
- Complete mitochondrial genomes (FASTA)
- GenBank format with CDS annotations (preferred for accurate gene extraction)

#### 1.2 Phenotype Data
Create `data/phenotypes/phenotype_data.csv` with the following columns:

| Column | Description | Example |
|--------|-------------|---------|
| species_id | Unique identifier (match FASTA filename) | homo_sapiens |
| species_name | Full scientific name | Homo sapiens |
| max_lifespan_years | Maximum documented lifespan | 122.5 |
| body_mass_kg | Adult body mass in kilograms | 62.0 |
| notes | Optional notes | Reference: Smith et al. 2020 |

**Example:**
```csv
species_id,species_name,max_lifespan_years,body_mass_kg,notes
homo_sapiens,Homo sapiens,122.5,62.0,Human
mus_musculus,Mus musculus,4.0,0.02,House mouse
myotis_lucifugus,Myotis lucifugus,34,0.007,Little brown bat
```

**Data sources**:
- AnAge Database: https://genomics.senescence.info/species/
- NCBI Taxonomy: https://www.ncbi.nlm.nih.gov/taxonomy
- Primary literature

**Recommended minimum**: 20-30 species for robust analysis

### Step 2: Extract Mitochondrial Genes

```bash
python scripts/01_extract_mt_genes.py -i data/fasta_files -o data/extracted_genes -f fasta
```

**Options**:
- `-i, --input`: Input directory with FASTA/GenBank files
- `-o, --output`: Output directory for extracted genes
- `-f, --format`: File format (`fasta` or `genbank`)

**Output**:
- Creates separate directories for each gene
- Extracts both nucleotide and protein sequences
- Translation uses vertebrate mitochondrial genetic code (Table 2)

**Note**: For raw FASTA files, this script uses approximate gene coordinates. For best results, use GenBank format files with CDS annotations.

### Step 3: Align Genes

```bash
python scripts/02_align_genes.py -i data/extracted_genes -o data/alignments -a mafft -t protein
```

**Options**:
- `-i, --input`: Input directory with extracted genes
- `-o, --output`: Output directory for alignments
- `-a, --aligner`: Alignment tool (`mafft` or `muscle`)
- `-t, --type`: Sequence type (`protein` or `nucleotide`)

**Output**:
- Multiple sequence alignments for each gene
- Alignment statistics (length, gaps percentage)

**Recommendation**: Use protein alignments for CAAS analysis

### Step 4: Calculate Longevity Quotient (LQ)

```bash
python scripts/03_calculate_lq.py -i data/phenotypes/phenotype_data.csv -o data/phenotypes/lq_data.csv --plot
```

**Options**:
- `-i, --input`: Input phenotype CSV file
- `-o, --output`: Output CSV with LQ values
- `--top-decile`: Top decile fraction (default: 0.1 = top 10%)
- `--bottom-decile`: Bottom decile fraction (default: 0.1)
- `--plot`: Generate visualization plots

**Output**:
- `lq_data.csv`: Species with LQ values and classifications
- `long-lived_species.txt`: List of long-lived species
- `short-lived_species.txt`: List of short-lived species
- `intermediate_species.txt`: List of intermediate species
- `lq_distribution.png`: Visualization (if --plot used)

**LQ Formula**: LQ = Observed MLS / Expected MLS
- Expected MLS = 4.88 × body_mass^0.19 (de Magalhães et al. 2007)
- LQ > 1: Lives longer than expected for body size
- LQ < 1: Lives shorter than expected

### Step 5: CAAS Discovery

```bash
python scripts/04_caas_discovery.py -a data/alignments -l data/phenotypes/lq_data.csv -o results/caas_discovery
```

**Options**:
- `-a, --alignments`: Directory with aligned sequences
- `-l, --lq-data`: CSV with LQ data
- `-o, --output`: Output directory
- `--permutations`: Number of permutations for significance (default: 1000)
- `--skip-permutations`: Skip permutation testing

**Three Scenarios Detected**:

1. **Scenario 1**: All long-lived have same AA, all short-lived have different fixed AA
   - Example: Position 100: Long-lived = A, Short-lived = T
   
2. **Scenario 2**: All long-lived have same AA, short-lived have variable AAs
   - Example: Position 200: Long-lived = G, Short-lived = {C, T, A}
   
3. **Scenario 3**: Short-lived have fixed AA, long-lived have variable AAs
   - Example: Position 300: Short-lived = F, Long-lived = {Y, W, H}

**Output**:
- `caas_discovered.csv`: All discovered CAAS
- `gene_summary.json`: Per-gene statistics
- `discovery_summary.txt`: Overall summary

**Statistical Test**: Permutation test comparing observed vs. random expectations (p < 0.05)

### Step 6: CAAS Validation

```bash
python scripts/05_caas_validation.py -a data/alignments -l data/phenotypes/lq_data.csv -c results/caas_discovery/caas_discovered.csv -o results/caas_validation
```

**Options**:
- `-a, --alignments`: Directory with alignments
- `-l, --lq-data`: CSV with LQ data
- `-c, --caas`: CSV with discovered CAAS
- `-o, --output`: Output directory
- `--fdr`: FDR threshold (default: 0.05)

**Validation Method**:
1. For each discovered CAAS position
2. Group intermediate species by AA (long-lived AA vs short-lived AA)
3. Compare LQ values between groups (t-test)
4. Apply FDR correction (Benjamini-Hochberg)
5. Keep only validated positions (FDR < 0.05)

**Output**:
- `caas_validation_results.csv`: All validation results
- `caas_validated.csv`: Only validated CAAS
- `validation_summary.txt`: Summary statistics

**Expected**: ~40% validation rate (similar to Farré et al. 2021)

## Interpreting Results

### Key Output Files

1. **caas_validated.csv**: Your main results
   - Columns: gene, position, scenario, p_value, direction
   - Focus on validated CAAS for biological interpretation

2. **gene_summary.json**: Which genes have most CAAS
   - Identifies genes under strongest selection for longevity

3. **validation_summary.txt**: Overall statistics
   - Compare your results to expected ~42% validation rate

### What to Look For

1. **Genes with most validated CAAS**: Strongest candidates for longevity
2. **Specific positions**: Can be mapped to protein structures
3. **Amino acid changes**: Functional implications (charge, size, hydrophobicity)
4. **Validation rates**: Should be significantly higher than 5% (random expectation)

### Follow-Up Analyses (Not Included Yet)

1. **Functional annotation**: Map CAAS to protein domains
2. **Protein stability**: Compare stability between long/short-lived proteins
3. **dN/dS analysis**: Gene-level evolutionary rates
4. **Pathway enrichment**: Which biological pathways are enriched

## Troubleshooting

### Common Issues

**Problem**: No genes extracted
- **Solution**: Check FASTA file format. Use GenBank format if possible.

**Problem**: Low alignment quality (>50% gaps)
- **Solution**: Remove divergent species or problematic sequences

**Problem**: No CAAS discovered
- **Solution**: Check species classification. May need more species or different decile thresholds.

**Problem**: Low validation rate (<10%)
- **Solution**: May indicate insufficient power. Try:
  - Add more species
  - Adjust decile thresholds
  - Check data quality

### Getting Help

1. Check output logs for error messages
2. Verify input file formats match specifications
3. Ensure sufficient species in each category (min 5-6 per group)

## Example Workflow

```bash
# 1. Setup (one time)
mkdir -p data/fasta_files data/phenotypes

# 2. Add your data
# - Copy FASTA files to data/fasta_files/
# - Create data/phenotypes/phenotype_data.csv

# 3. Run pipeline
python scripts/01_extract_mt_genes.py
python scripts/02_align_genes.py
python scripts/03_calculate_lq.py -i data/phenotypes/phenotype_data.csv --plot
python scripts/04_caas_discovery.py
python scripts/05_caas_validation.py

# 4. Check results
cat results/caas_validation/validation_summary.txt
```

## Citation

If you use this pipeline, please cite:

**Original methodology**:
Farré X, et al. (2021) Comparative Analysis of Mammal Genomes Unveils Key Genomic Variability for Human Life Span. Mol Biol Evol. 38(11):4948-4961.

## License

This pipeline is provided as-is for academic and research purposes.

## Contact

For questions about this pipeline, check the troubleshooting section or review the original paper's methodology.
