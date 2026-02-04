# QUICK START GUIDE

## What You Need Before Starting

1. **Mitochondrial genome sequences** (FASTA or GenBank format)
   - Download from NCBI: https://www.ncbi.nlm.nih.gov/nuccore
   - Search: "mitochondrion complete genome [organism name]"
   - Recommended: 20-30 species minimum

2. **Phenotype data** (lifespan and body mass)
   - AnAge Database: https://genomics.senescence.info/species/
   - Create CSV file with: species_id, species_name, max_lifespan_years, body_mass_kg

## 5-Minute Setup

```bash
# 1. Create directory structure
mkdir -p data/fasta_files data/phenotypes results

# 2. Place your FASTA files in data/fasta_files/
# Example: homo_sapiens.fasta, mus_musculus.fasta, etc.

# 3. Create data/phenotypes/phenotype_data.csv
# Use the template in data/phenotypes/phenotype_template.csv
```

## Run All Analysis (One Command at a Time)

```bash
# Step 1: Extract genes (5-10 min)
python scripts/01_extract_mt_genes.py

# Step 2: Align genes (10-20 min)  
python scripts/02_align_genes.py

# Step 3: Calculate LQ and classify species (< 1 min)
python scripts/03_calculate_lq.py \
    -i data/phenotypes/phenotype_data.csv \
    --plot

# Step 4: Discover CAAS (10-30 min depending on permutations)
python scripts/04_caas_discovery.py

# Step 5: Validate CAAS (5-10 min)
python scripts/05_caas_validation.py
```

## Check Your Results

```bash
# Overall summary
cat results/caas_validation/validation_summary.txt

# View validated CAAS
head results/caas_validation/caas_validated.csv

# View plots
open data/phenotypes/lq_distribution.png
```

## Expected Output

- **Discovery**: 50-500 CAAS positions (varies by dataset)
- **Validation**: ~40% of discovered CAAS (20-200 validated)
- **Per gene**: 0-50 validated CAAS per gene
- **Most CAAS**: Usually in ND genes (largest, most variable)

## What Each File Contains

| File | Content |
|------|---------|
| `data/phenotypes/lq_data.csv` | Species with LQ values and classifications |
| `results/caas_discovery/caas_discovered.csv` | All discovered CAAS positions |
| `results/caas_validation/caas_validated.csv` | **MAIN RESULTS** - Validated CAAS only |
| `results/caas_validation/validation_summary.txt` | **READ THIS FIRST** - Summary statistics |

## Common Problems & Solutions

### "No genes extracted"
- Use GenBank format instead of FASTA
- Or manually annotate gene boundaries

### "Alignment failed"
- Install MAFFT: `sudo apt-get install mafft`
- Or use MUSCLE: add `-a muscle` flag

### "No CAAS found"
- Need more species (aim for 20-30)
- Check if species are diverse enough in lifespan
- Try different decile thresholds: `--top-decile 0.15 --bottom-decile 0.15`

### "Low validation rate (<10%)"
- Normal if dataset is small
- Try adding more intermediate species
- Check data quality (alignments, phenotypes)

## For Your Paper

### Key Results to Report

1. **Number of species analyzed**
   - Long-lived: X species
   - Short-lived: X species  
   - Intermediate: X species

2. **CAAS discovery**
   - Total discovered: X positions
   - Scenario 1: X, Scenario 2: X, Scenario 3: X
   - Permutation p-value: X

3. **CAAS validation**
   - Validated: X/X (X%)
   - Per gene breakdown
   - Top genes with most CAAS

4. **Example CAAS positions**
   - Pick 3-5 interesting positions
   - Describe AA change and functional implications

### Figures to Make

1. **LQ distribution** (already generated if you used --plot)
2. **CAAS per gene bar chart**
3. **Validation rate comparison**
4. **Example protein structure** with CAAS highlighted (use PyMOL)

### Methods Section Template

```
We analyzed mitochondrial genomes from X species to identify 
convergent amino acid substitutions (CAAS) associated with longevity.
Species were classified based on longevity quotient (LQ = observed 
lifespan / expected lifespan based on body mass) into long-lived 
(top 10%, n=X), short-lived (bottom 10%, n=X), and intermediate 
groups (n=X).

We identified CAAS using the three-scenario approach of Farré et al. 
(2021). Discovered CAAS were validated using phylogenetic ANOVA on 
intermediate species with FDR correction (α=0.05). Statistical 
significance was assessed via permutation testing (n=1000 permutations).

We discovered X CAAS positions across 13 mitochondrial protein-coding 
genes, of which X (X%) were validated. Genes with the most validated 
CAAS were [list top 3 genes].
```

## Next Steps After This Pipeline

1. **Functional annotation**: Map CAAS to protein domains
2. **Protein modeling**: Predict structural effects (FoldX, AlphaFold)
3. **Literature search**: Link genes to known longevity pathways
4. **Pathway analysis**: Gene ontology enrichment (if enough genes)
5. **Population genetics**: Check if AAs vary in human populations

## Need More Help?

1. Read the full README.md
2. Check the original paper: Farré et al. 2021 (Mol Biol Evol)
3. Review script comments for detailed explanations
