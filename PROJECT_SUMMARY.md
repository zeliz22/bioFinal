# Mitochondrial Longevity Analysis - Project Summary

## What I've Created For You

A complete bioinformatics pipeline to replicate the Farré et al. (2021) analysis using **mitochondrial genomes** instead of nuclear genomes. This will allow you to write a similar research paper comparing longevity across species.

## Files Provided

### Documentation
1. **README.md** - Complete documentation with detailed instructions
2. **QUICKSTART.md** - Quick reference guide to get started fast
3. **PROJECT_SUMMARY.md** (this file) - Overview of the project

### Main Pipeline
- **mitochondrial_longevity_pipeline.py** - Master pipeline orchestrator

### Analysis Scripts (in `scripts/` directory)

1. **01_extract_mt_genes.py**
   - Extracts 13 mitochondrial protein-coding genes from your FASTA files
   - Handles both FASTA and GenBank formats
   - Translates to protein sequences automatically

2. **02_align_genes.py**
   - Aligns sequences using MAFFT or MUSCLE
   - Creates multiple sequence alignments for each gene
   - Quality metrics and gap statistics

3. **03_calculate_lq.py**
   - Calculates Longevity Quotient (LQ) for each species
   - Corrects for body mass effects
   - Classifies species into long-lived, short-lived, and intermediate
   - Generates visualization plots

4. **04_caas_discovery.py**
   - Discovers Convergent Amino Acid Substitutions (CAAS)
   - Implements 3 scenarios from the paper
   - Statistical validation via permutation testing
   - Identifies positions distinguishing long vs short-lived species

5. **05_caas_validation.py**
   - Validates discovered CAAS using intermediate species
   - Phylogenetic ANOVA testing
   - FDR correction (Benjamini-Hochberg)
   - Expected ~40% validation rate

### Template Data
- **data/phenotypes/phenotype_template.csv** - Template for your species data

## How the Analysis Works

### Step 1: Data Preparation
- You provide mitochondrial genome sequences (FASTA)
- You create a phenotype CSV with lifespan and body mass data

### Step 2: Gene Extraction
- Extracts 13 protein-coding genes:
  - ATP6, ATP8 (ATP synthase)
  - COX1, COX2, COX3 (Cytochrome c oxidase)
  - CYTB (Cytochrome b)
  - ND1-ND6, ND4L (NADH dehydrogenase)

### Step 3: Alignment
- Multiple sequence alignment for each gene
- Uses MAFFT (recommended) or MUSCLE

### Step 4: LQ Calculation
- Formula: LQ = Observed Lifespan / Expected Lifespan
- Expected = 4.88 × body_mass^0.19
- Species classified by LQ deciles

### Step 5: CAAS Discovery
Three scenarios detected:
- **Scenario 1**: Fixed AA difference (long vs short)
- **Scenario 2**: Long-lived fixed, short-lived variable
- **Scenario 3**: Short-lived fixed, long-lived variable

### Step 6: CAAS Validation
- Uses intermediate species as independent test
- Phylogenetic ANOVA
- Only keeps statistically significant positions (FDR < 0.05)

## Expected Results

Based on the original paper adapted for mitochondria:

- **Input**: 20-30 species minimum
- **Discovery**: 50-500 CAAS positions
- **Validation**: ~40% (20-200 validated CAAS)
- **Top genes**: Usually ND genes (largest, most variable)

## Key Differences from Original Paper

| Aspect | Original Paper | Your Analysis |
|--------|---------------|---------------|
| Genome | Nuclear (~19,000 genes) | Mitochondrial (13 genes) |
| Genes analyzed | 13,035 protein-coding genes | 13 mitochondrial genes |
| Species | 57 mammals | Your dataset (20-30+ recommended) |
| CAAS discovered | 2,737 | Expected: 50-500 |
| Validation rate | 42.3% | Expected: ~40% |

## What You Need to Do

### Before Running

1. **Collect Data**:
   - Download mitochondrial genomes from NCBI
   - Get lifespan data from AnAge database
   - Get body mass data from literature

2. **Prepare Files**:
   - Place FASTA files in `data/fasta_files/`
   - Create `data/phenotypes/phenotype_data.csv`

3. **Install Requirements**:
   ```bash
   pip install biopython pandas numpy scipy matplotlib seaborn statsmodels
   sudo apt-get install mafft  # or brew install mafft on macOS
   ```

### Running the Pipeline

```bash
# Run scripts in order
python scripts/01_extract_mt_genes.py
python scripts/02_align_genes.py
python scripts/03_calculate_lq.py -i data/phenotypes/phenotype_data.csv --plot
python scripts/04_caas_discovery.py
python scripts/05_caas_validation.py

# Check results
cat results/caas_validation/validation_summary.txt
```

### After Running

1. **Analyze Results**:
   - Which genes have most CAAS?
   - What amino acid changes occurred?
   - Are they functionally important?

2. **Biological Interpretation**:
   - Map CAAS to protein structures
   - Check if positions are in functional domains
   - Literature search for gene functions

3. **Write Paper**:
   - Use methods template in QUICKSTART.md
   - Report discovery and validation statistics
   - Discuss functional implications

## For Your Paper

### Title Suggestions
- "Comparative Analysis of Mitochondrial Genomes Reveals Amino Acid Substitutions Associated with Longevity"
- "Convergent Evolution of Mitochondrial Proteins in Long-Lived Species"
- "Mitochondrial Genomic Signatures of Extended Lifespan Across Mammals"

### Key Sections

**Introduction**:
- Mitochondria and aging connection
- Variation in mammalian lifespan
- Comparative genomics approach

**Methods**:
- Use template in QUICKSTART.md
- Cite Farré et al. 2021 for methodology

**Results**:
- Report discovery statistics
- Validation rates
- Per-gene breakdown
- Example CAAS positions

**Discussion**:
- Compare to nuclear genome findings
- Mitochondrial-specific mechanisms
- Functional implications
- Clinical relevance

### Figures to Include

1. LQ distribution with species classification
2. CAAS discovery and validation flowchart
3. CAAS per gene bar chart
4. Example protein structure with CAAS highlighted
5. Validation statistics comparison

## Potential Extensions

After completing the basic analysis:

1. **Protein Modeling**:
   - Use AlphaFold to predict structures
   - Use FoldX to assess stability changes

2. **Functional Analysis**:
   - Map to mitochondrial pathways
   - Check expression patterns
   - Link to oxidative stress response

3. **Population Genetics**:
   - Check human variation at CAAS positions
   - Compare to human longevity GWAS

4. **Evolutionary Analysis**:
   - dN/dS ratios per gene
   - Ancestral state reconstruction
   - Selection tests

5. **Clinical Relevance**:
   - Link to mitochondrial diseases
   - Aging-related conditions
   - Therapeutic targets

## Troubleshooting

Send me:
- Error messages
- Summary statistics from each step
- Number of species and genes
- Alignment quality metrics

I'll help you:
- Debug issues
- Interpret results
- Optimize parameters
- Refine analysis

## Expected Timeline

- **Data collection**: 1-2 weeks
- **Pipeline execution**: 1-2 hours
- **Results analysis**: 1 week
- **Paper writing**: 2-4 weeks
- **Total**: 1-2 months

## Success Criteria

✓ Successfully extract genes from all species
✓ Good quality alignments (<30% gaps)
✓ Clear LQ classification (distinct groups)
✓ Significant CAAS enrichment (p < 0.05)
✓ Reasonable validation rate (30-50%)
✓ Biologically meaningful results

## Resources

- **AnAge Database**: https://genomics.senescence.info/species/
- **NCBI Genome**: https://www.ncbi.nlm.nih.gov/genome
- **Original Paper**: Farré et al. 2021, Mol Biol Evol
- **Mitochondrial Genes**: MitoMap - https://www.mitomap.org

## Next Steps

1. Read QUICKSTART.md for immediate actions
2. Collect your data (FASTA + phenotypes)
3. Run the pipeline step by step
4. Send me the outputs from each step
5. I'll help you interpret and analyze results
6. Together we'll write the paper!

## Questions?

When you run into issues or get results, send me:
- Output files from `results/` directory
- Error messages if any
- Your phenotype CSV
- Summary statistics

I'm here to help you write a great paper! 🚀
