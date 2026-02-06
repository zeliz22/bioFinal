## Mitochondrial Genomic Signatures of Mammalian Longevity

### Abstract

Comparative genomics has revealed that patterns of genomic variability are associated with maximum lifespan across mammals, as elegantly demonstrated by Farré et al. in their nuclear-genome study of convergent amino-acid substitutions (CAAS) and longevity \[[Farré et al., 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC8557403/)\]. Inspired by that work, this project focuses specifically on **mitochondrial DNA**, using complete mitochondrial genomes to search for amino-acid changes associated with extended lifespan.  
Here, we (1) compute **Longevity Quotient (LQ)** for mammalian species from the **AnAge** database \[[AnAge database](https://genomics.senescence.info/species/)\], (2) classify species into long‑lived, short‑lived, and validation groups, (3) download their complete mitochondrial genomes from **NCBI** in GenBank format, (4) extract the 13 canonical mitochondrial protein‑coding genes, (5) align nucleotide and protein sequences, and (6) identify CAAS distinguishing long‑lived and short‑lived species, with statistical validation on intermediate‑LQ species. The goal is to obtain a small but interpretable set of mitochondrial amino‑acid substitutions that are enriched in long‑lived mammals and could represent **mitochondrial genomic signatures of longevity**.

---

## Overview of the Repository

- **anage_LQ.py**  
  **Step 1 – Longevity Quotient calculation and target selection.**  
  - Reads the AnAge tab‑delimited file (e.g. `data/anage_data.txt`).  
  - Filters species by taxonomic class (default: **Mammalia**).  
  - Computes Longevity Quotient (LQ) and writes `data/LQ/lq_results.csv`.  
  - Creates species lists for downstream steps:
    - `data/LQ/fetch_targets.txt` – all species used for genome download  
    - `data/LQ/long_lived_targets.txt` – long‑lived species  
    - `data/LQ/short_lived_targets.txt` – short‑lived species  
    - `data/LQ/validation_targets.txt` – intermediate‑LQ species (for validation)

- **download_genbank.py**  
  **Step 2 – Download complete mitochondrial genomes from NCBI.**  
  - Reads species IDs from `data/LQ/fetch_targets.txt`.  
  - Queries NCBI `nuccore` (RefSeq preferred) for each species’ **complete mitochondrial genome**.  
  - Saves GenBank files to `data/genbank_files/*.gb`.  

- **extract_genes_genbank.py**  
  **Step 3 – Extract mitochondrial genes from GenBank.**  
  - Reads each `.gb` file in `data/genbank_files/`.  
  - Uses CDS annotations to extract the 13 protein‑coding mitochondrial genes:
    - **ATP6, ATP8, COX1, COX2, COX3, CYTB, ND1, ND2, ND3, ND4, ND4L, ND5, ND6**  
  - Handles multiple gene name variants in GenBank annotations (e.g. COI, COXI, NAD1, COB, etc.).  
  - Writes:
    - Nucleotide FASTA: `data/extracted_genes/nucleotides/GENE/species.fasta`  
    - Protein FASTA: `data/extracted_genes/proteins/GENE/species.fasta`

- **align_genes.py**  
  **Step 4 – Multiple sequence alignment of mitochondrial genes.**  
  - Combines per‑species FASTA files per gene into multi‑FASTA.  
  - Aligns sequences using **MAFFT** (default) or **MUSCLE**.  
  - Produces alignments for both:
    - `data/alignments/proteins/GENE_aligned.fasta`  
    - `data/alignments/nucleotides/GENE_aligned.fasta`  
  - Reports basic alignment statistics (number of sequences, alignment length, gap percentage).

- **caas_discovery_from_lists.py**  
  **Step 5 – CAAS discovery and validation using predefined species lists.**  
  - Uses protein and/or nucleotide alignments in `data/alignments/`.  
  - Reads long‑lived and short‑lived species from:
    - `data/LQ/long_lived_targets.txt`  
    - `data/LQ/short_lived_targets.txt`  
  - Optionally reads `data/LQ/lq_results.csv` to use LQ values in validation.  
  - Detects positions where amino‑acid patterns differ systematically between groups (CAAS).  
  - Performs a Mann–Whitney U‑test–based validation and a permutation test to assess enrichment.

- **mitochondrial_longevity_pipeline.py**  
  **Master pipeline runner.**  
  - Cleans all intermediate data directories:
    - `data/alignments/`  
    - `data/LQ/`  
    - `data/extracted_genes/`  
    - `data/genbank_files/`  
  - Then runs the full pipeline in order:
    1. `anage_LQ.py`  
    2. `download_genbank.py`  
    3. `extract_genes_genbank.py`  
    4. `align_genes.py`  
    5. `caas_discovery_from_lists.py`

---

## Installation and Requirements

- **Python**: 3.x  
- **Python packages** (install via `pip`):
  - **biopython**
  - **pandas**
  - **numpy**
  - **scipy**
  - (optionally) **matplotlib**, **seaborn**, **statsmodels** if you extend the analyses
- **External tools**:
  - **MAFFT** (recommended) or **MUSCLE** for multiple sequence alignment  
  - Internet access to query NCBI and download GenBank records

Example installation:

```bash
pip install biopython pandas numpy scipy
sudo apt-get install mafft   # or: sudo apt-get install muscle
```

---

## Data Inputs

- **AnAge longevity and body mass data**
  - Tab‑delimited file such as `data/anage_data.txt` obtained from the **AnAge** database \[[AnAge database](https://genomics.senescence.info/species/)\].
  - Must contain (at least) columns similar to:
    - `Class`, `Genus`, `Species`  
    - `Maximum longevity (yrs)`  
    - `Body mass (g)`  

- **NCBI mitochondrial genomes**
  - Downloaded automatically by `download_genbank.py` from NCBI’s `nuccore` database into `data/genbank_files/` as `.gb` files.  
  - Searches are constrained to **mitochondrial complete genomes**, preferring curated RefSeq accessions.

---

## Pipeline Usage

- **Run everything from scratch (recommended)**:

```bash
python mitochondrial_longevity_pipeline.py
```

This will:
- Wipe `data/alignments/`, `data/LQ/`, `data/extracted_genes/`, and `data/genbank_files/`.  
- Recompute LQ and target species.  
- Re‑download mitochondrial genomes.  
- Re‑extract genes and realign them.  
- Run CAAS discovery and validation.

- **Run individual steps manually** (if you want more control):

```bash
python anage_LQ.py
python download_genbank.py
python extract_genes_genbank.py
python align_genes.py
python caas_discovery_from_lists.py
```

---

## Methodological Details

### Longevity Quotient (LQ) Calculation

- **Objective**: Correct maximum lifespan for body mass to identify species that live longer or shorter than expected for their size, following the comparative framework of \[[Farré et al., 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC8557403/)\].  
- **Input**: AnAge mammalian dataset (`data/anage_data.txt`).  
- **Steps (implemented in `anage_LQ.py`)**:
  - Filter to **Mammalia** (or another class via `--class_filter`).  
  - Compute:
    - \( \text{body\_mass\_kg} = \frac{\text{Body mass (g)}}{1000} \)  
    - \( \text{expected\_mls} = 4.88 \times \text{body\_mass\_kg}^{0.19} \)  
    - \( \text{LQ} = \dfrac{\text{max\_lifespan\_years}}{\text{expected\_mls}} \)  
  - Define species ID:
    - `species_id = Genus + '_' + Species` (e.g. `Homo_sapiens`).  
  - Save processed table (`lq_results.csv`) with `species_id`, `longevity_quotient`, lifespan, body mass, and taxonomic information.  
  - Select groups:
    - **Long‑lived**: highest‑LQ species (top tail).  
    - **Short‑lived**: lowest‑LQ species (bottom tail).  
    - **Validation**: species from the mid‑LQ range.  
  - Write text files listing species per group for downstream scripts.

### Mitochondrial Gene Extraction

- **Objective**: Obtain consistent sets of the 13 mitochondrial protein‑coding genes across species directly from annotated GenBank records (mitochondrial genomes).  
- **Implementation (`extract_genes_genbank.py`)**:
  - For each `*.gb` file:
    - Parse GenBank annotations using **Biopython**.  
    - Iterate over CDS features; try to infer a standard gene name from `gene` or `product` qualifiers.
    - Normalize various annotation styles (e.g. “Cytochrome c oxidase subunit I” → `COX1`).
    - Extract the nucleotide sequence using feature locations.
    - Translate using the **vertebrate mitochondrial genetic code (table 2)**:
      - Remove trailing stop codon (`*`) when present.
  - Save per‑gene nucleotide and protein FASTA files organized by gene folder.

### Multiple Sequence Alignment

- **Objective**: Produce per‑gene multiple sequence alignments for both nucleotide and amino‑acid sequences.  
- **Implementation (`align_genes.py`)**:
  - For each gene:
    - Combine individual FASTA files into a gene‑specific multi‑FASTA.  
    - Align with MAFFT (`--auto`) or MUSCLE.  
    - Save:
      - `data/alignments/proteins/GENE_aligned.fasta`  
      - `data/alignments/nucleotides/GENE_aligned.fasta`  
  - Compute basic alignment statistics:
    - Number of sequences.  
    - Alignment length.  
    - Overall gap percentage.

### CAAS Discovery and Validation

- **Objective**: Identify **Convergent Amino‑Acid Substitutions (CAAS)** that distinguish long‑lived from short‑lived mammals, and test whether these patterns are associated with LQ in intermediate‑LQ (validation) species.  
- **Implementation (`caas_discovery_from_lists.py`)**:
  - Uses predefined sets of species:
    - `discovery_long`: intersection of `long_lived_targets.txt` with alignment species IDs.  
    - `discovery_short`: intersection of `short_lived_targets.txt` with alignment species IDs.  
    - `validation_species`: remaining species present in alignments.  
  - For each aligned position in each gene:
    - Skip sites where any discovery species is missing or carries a gap.  
    - Evaluate three classical CAAS scenarios (adapted from \[[Farré et al., 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC8557403/)\]):
      - **Scenario 1**: long‑lived fixed for amino acid A, short‑lived fixed for amino acid B, \(A \neq B\).  
      - **Scenario 2**: long‑lived fixed (A), short‑lived variable for amino acids not containing A.  
      - **Scenario 3**: short‑lived fixed (B), long‑lived variable for amino acids not containing B.  
    - Record gene, alignment position, scenario, and group‑specific residues.
  - **Validation using intermediate‑LQ species**:
    - For each CAAS position, collect LQ values (from `lq_results.csv`) for validation species grouped by their amino acid (long‑like vs short‑like configuration).  
    - Apply **Mann–Whitney U test** (one‑sided; typically testing whether long‑like configuration is associated with higher LQ).  
    - Mark CAAS as validated if the P‑value is below a relaxed threshold (e.g. 0.2, as implemented in the script).  
  - **Permutation test**:
    - Randomly shuffle species labels between long‑lived and short‑lived groups many times.  
    - Recount S1 + S2 CAAS under random groupings.  
    - Compare observed CAAS count to this null distribution to obtain an empirical P‑value for enrichment.
  - **Outputs** (per alignment type – protein/nucleotide):
    - `caas_discovered_*.csv` – all discovered CAAS.  
    - `caas_validated_*.csv` – subset passing validation criteria.  
    - `summary_*.json` – summary statistics (group sizes, counts per scenario, validation counts, permutation P‑value, etc.).

---

## Interpretation and Extensions

- **Biological interpretation**:
  - Positions where long‑lived species share one amino acid and short‑lived species another may indicate **convergent adaptation** in mitochondrial proteins affecting lifespan.  
  - Enrichment of CAAS in particular genes or complexes (e.g. NADH dehydrogenase vs. cytochrome c oxidase) may point to specific mitochondrial pathways linked to longevity.

- **Possible extensions**:
  - Structural mapping of validated CAAS on mitochondrial protein 3D models.  
  - Integration with nuclear genomic data (as in \[[Farré et al., 2021](https://pmc.ncbi.nlm.nih.gov/articles/PMC8557403/)\]) to compare nuclear and mitochondrial signatures.  
  - More formal phylogenetic methods (e.g. PAML \[[Yang, 1997](https://pmc.ncbi.nlm.nih.gov/articles/PMC8557403/)\] cited in the same paper) or dN/dS analyses for selection tests on mitochondrial genes.

---

## Reference

- **Comparative Analysis of Mammal Genomes Unveils Key Genomic Variability for Human Life Span** – Farré et al., *Molecular Biology and Evolution* (2021) \[[link](https://pmc.ncbi.nlm.nih.gov/articles/PMC8557403/)\].  
  - This project adapts the CAAS‑and‑longevity framework from that nuclear‑genome study to focus exclusively on **mitochondrial DNA**.

