#!/usr/bin/env python3
"""
Script 2: Align Mitochondrial Genes

Updated to support the new directory structure:
data/extracted_genes/nucleotides/GENE/
data/extracted_genes/proteins/GENE/
"""

import os
import subprocess
from pathlib import Path
import argparse
from Bio import SeqIO, AlignIO

MT_GENES = ['ATP6', 'ATP8', 'COX1', 'COX2', 'COX3', 'CYTB',
            'ND1', 'ND2', 'ND3', 'ND4', 'ND4L', 'ND5', 'ND6']


def check_aligner_available(aligner='mafft'):
    """Check if alignment tool is available"""
    try:
        if aligner == 'mafft':
            subprocess.run(['mafft', '--version'],
                           capture_output=True, check=True)
        elif aligner == 'muscle':
            subprocess.run(['muscle', '-version'],
                           capture_output=True, check=True)
        return True
    except (subprocess.CalledProcessError, FileNotFoundError):
        return False


def combine_sequences(gene_dir, output_file):
    """Combine individual FASTA files into a single multi-FASTA"""
    gene_path = Path(gene_dir)
    # Search for .fasta, .fa, or .fas extensions
    fasta_files = list(gene_path.glob("*.fasta")) + list(gene_path.glob("*.fa")) + list(gene_path.glob("*.fas"))

    if not fasta_files:
        print(f"  Warning: No FASTA files found in {gene_dir}")
        return False

    sequences = []
    for fasta_file in fasta_files:
        try:
            for record in SeqIO.parse(fasta_file, "fasta"):
                # Use filename as part of ID if sequences don't have unique IDs
                sequences.append(record)
        except Exception as e:
            print(f"  Error reading {fasta_file}: {e}")

    if sequences:
        SeqIO.write(sequences, output_file, "fasta")
        print(f"  Combined {len(sequences)} sequences into {output_file.name}")
        return True
    return False


def align_with_mafft(input_file, output_file):
    """Align sequences using MAFFT"""
    try:
        cmd = ['mafft', '--auto', '--thread', '-1', str(input_file)]
        with open(output_file, 'w') as out:
            subprocess.run(cmd, stdout=out, stderr=subprocess.DEVNULL, check=True)
        return True
    except Exception as e:
        print(f"  MAFFT alignment failed: {e}")
        return False


def align_with_muscle(input_file, output_file):
    """Align sequences using MUSCLE"""
    try:
        # Muscle v5 syntax
        cmd = ['muscle', '-align', str(input_file), '-output', str(output_file)]
        subprocess.run(cmd, stderr=subprocess.DEVNULL, check=True)
        return True
    except Exception as e:
        print(f"  MUSCLE alignment failed: {e}")
        return False


def calculate_alignment_stats(alignment_file):
    """Calculate basic statistics for an alignment"""
    try:
        alignment = AlignIO.read(alignment_file, "fasta")
        num_sequences = len(alignment)
        alignment_length = alignment.get_alignment_length()
        total_positions = num_sequences * alignment_length
        gap_count = sum(str(record.seq).count('-') for record in alignment)
        gap_percentage = (gap_count / total_positions) * 100 if total_positions > 0 else 0

        return {
            'num_sequences': num_sequences,
            'length': alignment_length,
            'gap_percentage': gap_percentage
        }
    except Exception as e:
        print(f"  Error calculating stats: {e}")
        return None


def align_genes(input_dir, output_dir, aligner='mafft', sequence_type='proteins'):
    """
    Align all mitochondrial genes based on new folder structure
    """
    input_path = Path(input_dir)
    # The crucial fix: point to the correct subfolder (nucleotides/ or proteins/)
    base_gene_path = input_path / sequence_type

    output_path = Path(output_dir) / sequence_type
    output_path.mkdir(exist_ok=True, parents=True)

    temp_dir = output_path / "temp"
    temp_dir.mkdir(exist_ok=True)

    alignment_summary = {}

    for gene in MT_GENES:
        print(f"\n{'=' * 60}")
        print(f"Aligning {gene} ({sequence_type})")
        print(f"{'=' * 60}")

        # Look inside the subfolder for the specific gene
        gene_dir = base_gene_path / gene

        if not gene_dir.exists():
            print(f"  ✗ Directory not found: {gene_dir}")
            alignment_summary[gene] = {'status': 'failed', 'reason': 'directory_not_found'}
            continue

        combined_file = temp_dir / f"{gene}_combined.fasta"
        if not combine_sequences(gene_dir, combined_file):
            alignment_summary[gene] = {'status': 'failed', 'reason': 'no_sequences'}
            continue

        aligned_file = output_path / f"{gene}_aligned.fasta"

        if aligner == 'mafft':
            success = align_with_mafft(combined_file, aligned_file)
        else:
            success = align_with_muscle(combined_file, aligned_file)

        if success:
            stats = calculate_alignment_stats(aligned_file)
            if stats:
                print(f"  ✓ Alignment successful: {stats['num_sequences']} seqs, {stats['length']} positions")
                alignment_summary[gene] = {'status': 'success', **stats}
        else:
            alignment_summary[gene] = {'status': 'failed', 'reason': 'alignment_error'}

    # Clean up
    for file in temp_dir.glob("*"):
        file.unlink()
    temp_dir.rmdir()

    return alignment_summary


def main():
    parser = argparse.ArgumentParser(description='Align mitochondrial genes')
    parser.add_argument('-i', '--input', default='data/extracted_genes', help='Root input dir')
    parser.add_argument('-o', '--output', default='data/alignments', help='Output dir')
    parser.add_argument('-a', '--aligner', choices=['mafft', 'muscle'], default='mafft')
    # Changed default/choices to match your folder names: 'nucleotides' or 'proteins'
    parser.add_argument('-t', '--type', choices=['nucleotides', 'proteins'], default='proteins')

    args = parser.parse_args()

    if not check_aligner_available(args.aligner):
        print(f"ERROR: {args.aligner} not found.")
        return

    # Change this in your main() function to do both automatically:
    for seq_type in ['proteins', 'nucleotides']:
        print(f"Starting alignment for: {seq_type}")
        align_genes(args.input, args.output, args.aligner, seq_type)

    # summary = align_genes(args.input, args.output, args.aligner, args.type)
    print(f"\n✓ Finished! Aligned files are in {args.output}/{args.type}")


if __name__ == "__main__":
    main()