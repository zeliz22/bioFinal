#!/usr/bin/env python3
"""
Script: Extract Mitochondrial Protein-Coding Genes from GenBank Files

Extracts the 13 mitochondrial protein-coding genes from GenBank format files
using the actual CDS annotations (not human coordinates).

Mitochondrial protein-coding genes:
- ATP6, ATP8: ATP synthase subunits
- COX1, COX2, COX3: Cytochrome c oxidase subunits
- CYTB: Cytochrome b
- ND1, ND2, ND3, ND4, ND4L, ND5, ND6: NADH dehydrogenase subunits
"""

import os
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse

# Mitochondrial protein-coding genes
MT_GENES = [
    'ATP6', 'ATP8',
    'COX1', 'COX2', 'COX3',
    'CYTB',
    'ND1', 'ND2', 'ND3', 'ND4', 'ND4L', 'ND5', 'ND6'
]

# Gene name variations to handle different annotation formats
GENE_NAME_VARIANTS = {
    'COX1': ['COX1', 'CO1', 'COI', 'COXI'],
    'COX2': ['COX2', 'CO2', 'COII', 'COXII'],
    'COX3': ['COX3', 'CO3', 'COIII', 'COXIII'],
    'CYTB': ['CYTB', 'CYB', 'COB'],
    'ND1': ['ND1', 'NAD1'],
    'ND2': ['ND2', 'NAD2'],
    'ND3': ['ND3', 'NAD3'],
    'ND4': ['ND4', 'NAD4'],
    'ND4L': ['ND4L', 'NAD4L'],
    'ND5': ['ND5', 'NAD5'],
    'ND6': ['ND6', 'NAD6'],
    'ATP6': ['ATP6', 'ATPASE6'],
    'ATP8': ['ATP8', 'ATPASE8'],
}


def normalize_gene_name(gene_text):
    """Convert various gene name formats to standard name"""
    gene_upper = gene_text.upper().strip()

    # Remove common prefixes/suffixes
    gene_upper = gene_upper.replace('NADH DEHYDROGENASE SUBUNIT', 'ND')
    gene_upper = gene_upper.replace('CYTOCHROME C OXIDASE SUBUNIT', 'COX')
    gene_upper = gene_upper.replace('CYTOCHROME B', 'CYTB')
    gene_upper = gene_upper.replace('ATP SYNTHASE', 'ATP')
    gene_upper = gene_upper.replace('SUBUNIT', '').strip()

    # Check against variants
    for standard_name, variants in GENE_NAME_VARIANTS.items():
        if any(variant in gene_upper for variant in variants):
            return standard_name

    return None


def extract_gene_from_genbank(genbank_file, gene_name):
    """
    Extract gene from GenBank format file using CDS annotations

    Args:
        genbank_file: Path to GenBank file
        gene_name: Standard gene name (e.g., 'COX1', 'ND1')

    Returns:
        Nucleotide sequence or None if not found
    """
    try:
        record = SeqIO.read(genbank_file, "genbank")

        for feature in record.features:
            if feature.type == "CDS":
                # Get gene name from feature
                feature_gene = None

                # Try 'gene' qualifier
                if 'gene' in feature.qualifiers:
                    feature_gene = normalize_gene_name(feature.qualifiers['gene'][0])

                # Try 'product' qualifier if gene not found
                if not feature_gene and 'product' in feature.qualifiers:
                    feature_gene = normalize_gene_name(feature.qualifiers['product'][0])

                # Check if this is our target gene
                if feature_gene == gene_name:
                    # Extract sequence using the location in the annotation
                    gene_seq = feature.location.extract(record.seq)
                    return gene_seq

        return None

    except Exception as e:
        print(f"    Error reading GenBank file: {e}")
        return None


def translate_sequence(seq):
    """
    Translate DNA sequence to protein using mitochondrial genetic code

    Args:
        seq: DNA sequence

    Returns:
        Protein sequence or None if translation fails
    """
    try:
        # Vertebrate mitochondrial genetic code is table 2
        protein = seq.translate(table=2, cds=False)

        # Remove stop codon if present at the end
        protein_str = str(protein)
        if protein_str.endswith('*'):
            protein = protein[:-1]

        return protein

    except Exception as e:
        print(f"    Translation error: {e}")
        return None


def process_mitochondrial_genomes(input_dir, output_dir):
    """
    Process all GenBank files and extract genes

    Args:
        input_dir: Directory containing GenBank files
        output_dir: Directory for output gene files
    """
    input_path = Path(input_dir)
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True, parents=True)

    # Create output directories for nucleotides and proteins
    nucleotide_dir = output_path / "nucleotides"
    protein_dir = output_path / "proteins"
    nucleotide_dir.mkdir(exist_ok=True)
    protein_dir.mkdir(exist_ok=True)

    # Create subdirectories for each gene
    for gene in MT_GENES:
        (nucleotide_dir / gene).mkdir(exist_ok=True)
        (protein_dir / gene).mkdir(exist_ok=True)

    # Find all GenBank files
    gb_files = list(input_path.glob("*.gb")) + list(input_path.glob("*.gbk"))

    if not gb_files:
        print("ERROR: No GenBank files found!")
        print(f"Looking in: {input_path}")
        return

    print(f"\nFound {len(gb_files)} GenBank files to process")
    print("=" * 70)

    # Track extraction statistics
    extraction_summary = {gene: {'success': 0, 'failed': 0} for gene in MT_GENES}

    # Process each file
    for gb_file in gb_files:
        species_id = gb_file.stem
        species_name = species_id.replace('_', ' ')

        print(f"\nProcessing: {species_name}")
        print("-" * 70)

        for gene in MT_GENES:
            # Extract nucleotide sequence using GenBank annotations
            gene_seq = extract_gene_from_genbank(gb_file, gene)

            if gene_seq:
                # Save nucleotide sequence
                nuc_file = nucleotide_dir / gene / f"{species_id}.fasta"
                nuc_record = SeqRecord(
                    gene_seq,
                    id=species_id,
                    description=f"{species_name} {gene} nucleotide"
                )
                SeqIO.write(nuc_record, nuc_file, "fasta")

                # Translate to protein
                protein_seq = translate_sequence(gene_seq)

                if protein_seq:
                    # Save protein sequence
                    prot_file = protein_dir / gene / f"{species_id}.fasta"
                    prot_record = SeqRecord(
                        protein_seq,
                        id=species_id,
                        description=f"{species_name} {gene} protein"
                    )
                    SeqIO.write(prot_record, prot_file, "fasta")

                    extraction_summary[gene]['success'] += 1
                    print(f"  ✓ {gene:6s}: {len(gene_seq):4d} bp → {len(protein_seq):3d} aa")
                else:
                    extraction_summary[gene]['failed'] += 1
                    print(f"  ✗ {gene:6s}: Extracted but translation failed")
            else:
                extraction_summary[gene]['failed'] += 1
                print(f"  ✗ {gene:6s}: Not found in annotations")

    # Print summary
    print("\n" + "=" * 70)
    print("EXTRACTION SUMMARY")
    print("=" * 70)
    print(f"{'Gene':<8} {'Success':<10} {'Failed':<10} {'Total':<10} {'Success %'}")
    print("-" * 70)

    for gene in MT_GENES:
        success = extraction_summary[gene]['success']
        failed = extraction_summary[gene]['failed']
        total = success + failed
        pct = (success / total * 100) if total > 0 else 0
        print(f"{gene:<8} {success:<10} {failed:<10} {total:<10} {pct:6.1f}%")

    print("=" * 70)

    # Overall statistics
    total_success = sum(s['success'] for s in extraction_summary.values())
    total_attempts = sum(s['success'] + s['failed'] for s in extraction_summary.values())
    overall_pct = (total_success / total_attempts * 100) if total_attempts > 0 else 0

    print(f"\nOverall: {total_success}/{total_attempts} ({overall_pct:.1f}%) successful extractions")

    return extraction_summary


def main():
    parser = argparse.ArgumentParser(
        description='Extract mitochondrial protein-coding genes from GenBank files'
    )
    parser.add_argument('-i', '--input', default='data/genbank_files',
                        help='Input directory with GenBank files')
    parser.add_argument('-o', '--output', default='data/extracted_genes',
                        help='Output directory for extracted genes')

    args = parser.parse_args()

    print("=" * 70)
    print("MITOCHONDRIAL GENE EXTRACTION FROM GENBANK FILES")
    print("=" * 70)
    print(f"Input directory:  {args.input}")
    print(f"Output directory: {args.output}")

    process_mitochondrial_genomes(args.input, args.output)

    print("\n" + "=" * 70)
    print("✓ Gene extraction complete!")
    print("=" * 70)
    print(f"\nOutput structure:")
    print(f"  {args.output}/nucleotides/GENE/species.fasta  ← Nucleotide sequences")
    print(f"  {args.output}/proteins/GENE/species.fasta     ← Protein sequences")
    print("\nNext steps:")
    print("  1. Combine sequences per gene")
    print("  2. Align nucleotides with MAFFT")
    print("  3. Translate aligned nucleotides to proteins")
    print("=" * 70)


if __name__ == "__main__":
    main()
