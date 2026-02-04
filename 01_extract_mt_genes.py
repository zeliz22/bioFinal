#!/usr/bin/env python3
"""
Script 1: Extract Mitochondrial Protein-Coding Genes

Extracts the 13 mitochondrial protein-coding genes from complete mitochondrial
genome FASTA files. Handles different annotation formats.

Mitochondrial protein-coding genes:
- ATP6, ATP8: ATP synthase subunits
- COX1, COX2, COX3: Cytochrome c oxidase subunits
- CYTB: Cytochrome b
- ND1, ND2, ND3, ND4, ND4L, ND5, ND6: NADH dehydrogenase subunits
"""

import os
import re
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

# Standard gene coordinates for human mitochondrial genome (as reference)
# These are approximate - actual extraction should use annotations
HUMAN_MT_COORDS = {
    'ND1': (3307, 4262),
    'ND2': (4470, 5511),
    'COX1': (5904, 7445),
    'COX2': (7586, 8269),
    'ATP8': (8366, 8572),
    'ATP6': (8527, 9207),
    'COX3': (9207, 9990),
    'ND3': (10059, 10404),
    'ND4L': (10470, 10766),
    'ND4': (10760, 12137),
    'ND5': (12337, 14148),
    'ND6': (14149, 14673),
    'CYTB': (14747, 15887)
}

def extract_gene_from_genbank(genbank_file, gene_name):
    """Extract gene from GenBank format file"""
    try:
        record = SeqIO.read(genbank_file, "genbank")
        for feature in record.features:
            if feature.type == "CDS":
                # Check various gene name formats
                gene = feature.qualifiers.get('gene', [''])[0].upper()
                product = feature.qualifiers.get('product', [''])[0].upper()
                
                if gene_name in gene or gene_name in product:
                    seq = feature.location.extract(record.seq)
                    return seq
    except Exception as e:
        print(f"Error reading GenBank file: {e}")
    return None

def extract_gene_from_fasta(fasta_file, gene_name, coords=None):
    """
    Extract gene from FASTA file
    If coords provided, use coordinates; otherwise try to find based on similarity
    """
    try:
        record = SeqIO.read(fasta_file, "fasta")
        
        if coords:
            start, end = coords
            # Convert to 0-based indexing
            gene_seq = record.seq[start-1:end]
            return gene_seq
        
        # If no coordinates, return None - will need manual annotation
        return None
        
    except Exception as e:
        print(f"Error reading FASTA file: {e}")
    return None

def translate_sequence(seq):
    """Translate DNA sequence to protein"""
    # Mitochondrial genetic code (vertebrate)
    # Table 2 in NCBI
    seq_obj = Seq(str(seq))
    try:
        protein = seq_obj.translate(table=2)  # Vertebrate mitochondrial code
        return protein
    except Exception as e:
        print(f"Translation error: {e}")
        return None

def process_mitochondrial_genomes(input_dir, output_dir, file_format='fasta'):
    """
    Process all mitochondrial genome files and extract genes
    
    Args:
        input_dir: Directory containing input files
        output_dir: Directory for output gene files
        file_format: 'fasta' or 'genbank'
    """
    input_path = Path(input_dir)
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True, parents=True)
    
    # Create output directories for each gene
    for gene in MT_GENES:
        (output_path / gene).mkdir(exist_ok=True)
        (output_path / f"{gene}_protein").mkdir(exist_ok=True)
    
    # Process each file
    file_pattern = "*.fasta" if file_format == 'fasta' else "*.gb"
    files = list(input_path.glob(file_pattern))
    
    if not files:
        file_pattern = "*.fa" if file_format == 'fasta' else "*.gbk"
        files = list(input_path.glob(file_pattern))
    
    print(f"\nFound {len(files)} {file_format} files to process")
    
    extraction_summary = {gene: {'success': 0, 'failed': 0} for gene in MT_GENES}
    
    for file in files:
        # Extract species name from filename
        species_name = file.stem.replace('_', ' ').replace('.', ' ')
        species_id = file.stem
        
        print(f"\nProcessing: {species_name}")
        
        for gene in MT_GENES:
            gene_seq = None
            
            if file_format == 'genbank':
                gene_seq = extract_gene_from_genbank(file, gene)
            else:
                # For FASTA, you may need to provide coordinates
                # This is a simplified version - may need manual curation
                coords = HUMAN_MT_COORDS.get(gene)
                gene_seq = extract_gene_from_fasta(file, gene, coords)
            
            if gene_seq:
                # Save nucleotide sequence
                nuc_file = output_path / gene / f"{species_id}.fasta"
                record = SeqRecord(gene_seq, id=species_id, description=f"{species_name} {gene}")
                SeqIO.write(record, nuc_file, "fasta")
                
                # Translate and save protein sequence
                protein_seq = translate_sequence(gene_seq)
                if protein_seq:
                    prot_file = output_path / f"{gene}_protein" / f"{species_id}.fasta"
                    prot_record = SeqRecord(protein_seq, id=species_id, 
                                          description=f"{species_name} {gene} protein")
                    SeqIO.write(prot_record, prot_file, "fasta")
                    extraction_summary[gene]['success'] += 1
                    print(f"  ✓ {gene}: Extracted and translated")
                else:
                    extraction_summary[gene]['failed'] += 1
                    print(f"  ✗ {gene}: Translation failed")
            else:
                extraction_summary[gene]['failed'] += 1
                print(f"  ✗ {gene}: Extraction failed")
    
    # Print summary
    print("\n" + "="*60)
    print("EXTRACTION SUMMARY")
    print("="*60)
    for gene in MT_GENES:
        success = extraction_summary[gene]['success']
        failed = extraction_summary[gene]['failed']
        total = success + failed
        pct = (success/total*100) if total > 0 else 0
        print(f"{gene:6s}: {success:3d}/{total:3d} ({pct:5.1f}%) successful")
    print("="*60)
    
    return extraction_summary

def main():
    parser = argparse.ArgumentParser(
        description='Extract mitochondrial protein-coding genes from genome files'
    )
    parser.add_argument('-i', '--input', default='data/fasta_files',
                       help='Input directory with FASTA/GenBank files')
    parser.add_argument('-o', '--output', default='data/extracted_genes',
                       help='Output directory for extracted genes')
    parser.add_argument('-f', '--format', choices=['fasta', 'genbank'], default='fasta',
                       help='Input file format')
    
    args = parser.parse_args()
    
    print("="*60)
    print("MITOCHONDRIAL GENE EXTRACTION")
    print("="*60)
    print(f"Input directory: {args.input}")
    print(f"Output directory: {args.output}")
    print(f"File format: {args.format}")
    
    process_mitochondrial_genomes(args.input, args.output, args.format)
    
    print("\n✓ Gene extraction complete!")
    print(f"Output saved to: {args.output}")
    print("\nNOTE: If using raw FASTA files without annotations,")
    print("you may need to manually annotate gene boundaries or")
    print("use GenBank format files with CDS features.")

if __name__ == "__main__":
    main()
