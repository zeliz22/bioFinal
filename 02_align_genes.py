#!/usr/bin/env python3
"""
Script 2: Align Mitochondrial Genes

Performs multiple sequence alignment for each of the 13 mitochondrial genes
using MAFFT (preferred) or MUSCLE. Aligns both nucleotide and protein sequences.
"""

import os
import subprocess
from pathlib import Path
import argparse
from Bio import SeqIO, AlignIO
#from Bio.Align.Applications import MafftCommandline, MuscleCommandline

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
    fasta_files = list(gene_path.glob("*.fasta")) + list(gene_path.glob("*.fa"))
    
    if not fasta_files:
        print(f"  Warning: No FASTA files found in {gene_dir}")
        return False
    
    sequences = []
    for fasta_file in fasta_files:
        for record in SeqIO.parse(fasta_file, "fasta"):
            sequences.append(record)
    
    if sequences:
        SeqIO.write(sequences, output_file, "fasta")
        print(f"  Combined {len(sequences)} sequences into {output_file.name}")
        return True
    return False

def align_with_mafft(input_file, output_file, sequence_type='protein'):
    """Align sequences using MAFFT"""
    try:
        # MAFFT parameters
        # --auto: automatically selects strategy
        # --thread -1: use all cores
        cmd = [
            'mafft',
            '--auto',
            '--thread', '-1',
            str(input_file)
        ]
        
        with open(output_file, 'w') as out:
            subprocess.run(cmd, stdout=out, stderr=subprocess.DEVNULL, check=True)
        
        return True
    except Exception as e:
        print(f"  MAFFT alignment failed: {e}")
        return False

def align_with_muscle(input_file, output_file, sequence_type='protein'):
    """Align sequences using MUSCLE"""
    try:
        cmd = [
            'muscle',
            '-align', str(input_file),
            '-output', str(output_file)
        ]
        
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
        
        # Calculate % gaps
        total_positions = num_sequences * alignment_length
        gap_count = sum(str(record.seq).count('-') for record in alignment)
        gap_percentage = (gap_count / total_positions) * 100
        
        return {
            'num_sequences': num_sequences,
            'length': alignment_length,
            'gap_percentage': gap_percentage
        }
    except Exception as e:
        print(f"  Error calculating stats: {e}")
        return None

def align_genes(input_dir, output_dir, aligner='mafft', sequence_type='protein'):
    """
    Align all mitochondrial genes
    
    Args:
        input_dir: Directory with extracted genes
        output_dir: Directory for alignments
        aligner: 'mafft' or 'muscle'
        sequence_type: 'protein' or 'nucleotide'
    """
    input_path = Path(input_dir)
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True, parents=True)
    
    # Create temporary directory for combined sequences
    temp_dir = output_path / "temp"
    temp_dir.mkdir(exist_ok=True)
    
    alignment_summary = {}
    
    for gene in MT_GENES:
        print(f"\n{'='*60}")
        print(f"Aligning {gene} ({sequence_type})")
        print(f"{'='*60}")
        
        # Determine input directory based on sequence type
        if sequence_type == 'protein':
            gene_dir = input_path / f"{gene}_protein"
        else:
            gene_dir = input_path / gene
        
        if not gene_dir.exists():
            print(f"  ✗ Directory not found: {gene_dir}")
            alignment_summary[gene] = {'status': 'failed', 'reason': 'directory_not_found'}
            continue
        
        # Combine sequences
        combined_file = temp_dir / f"{gene}_{sequence_type}_combined.fasta"
        if not combine_sequences(gene_dir, combined_file):
            alignment_summary[gene] = {'status': 'failed', 'reason': 'no_sequences'}
            continue
        
        # Align sequences
        aligned_file = output_path / f"{gene}_{sequence_type}_aligned.fasta"
        
        if aligner == 'mafft':
            success = align_with_mafft(combined_file, aligned_file, sequence_type)
        elif aligner == 'muscle':
            success = align_with_muscle(combined_file, aligned_file, sequence_type)
        else:
            print(f"  Unknown aligner: {aligner}")
            success = False
        
        if success:
            # Calculate statistics
            stats = calculate_alignment_stats(aligned_file)
            if stats:
                print(f"  ✓ Alignment successful")
                print(f"    Sequences: {stats['num_sequences']}")
                print(f"    Length: {stats['length']} positions")
                print(f"    Gaps: {stats['gap_percentage']:.2f}%")
                
                alignment_summary[gene] = {
                    'status': 'success',
                    'num_sequences': stats['num_sequences'],
                    'length': stats['length'],
                    'gap_percentage': stats['gap_percentage']
                }
            else:
                alignment_summary[gene] = {'status': 'success', 'stats': 'unavailable'}
        else:
            alignment_summary[gene] = {'status': 'failed', 'reason': 'alignment_error'}
    
    # Clean up temp directory
    for file in temp_dir.glob("*"):
        file.unlink()
    temp_dir.rmdir()
    
    # Print summary
    print("\n" + "="*60)
    print("ALIGNMENT SUMMARY")
    print("="*60)
    successful = sum(1 for g in alignment_summary.values() if g['status'] == 'success')
    print(f"Successfully aligned: {successful}/{len(MT_GENES)} genes")
    print("\nDetailed results:")
    for gene, info in alignment_summary.items():
        if info['status'] == 'success':
            if 'num_sequences' in info:
                print(f"  {gene:6s}: {info['num_sequences']} sequences, "
                      f"{info['length']} bp, {info['gap_percentage']:.1f}% gaps")
            else:
                print(f"  {gene:6s}: Success (stats unavailable)")
        else:
            reason = info.get('reason', 'unknown')
            print(f"  {gene:6s}: Failed ({reason})")
    print("="*60)
    
    return alignment_summary

def main():
    parser = argparse.ArgumentParser(
        description='Align mitochondrial genes'
    )
    parser.add_argument('-i', '--input', default='data/extracted_genes',
                       help='Input directory with extracted genes')
    parser.add_argument('-o', '--output', default='data/alignments',
                       help='Output directory for alignments')
    parser.add_argument('-a', '--aligner', choices=['mafft', 'muscle'],
                       default='mafft', help='Alignment tool to use')
    parser.add_argument('-t', '--type', choices=['protein', 'nucleotide'],
                       default='protein', help='Sequence type to align')
    
    args = parser.parse_args()
    
    # Check if aligner is available
    if not check_aligner_available(args.aligner):
        print(f"ERROR: {args.aligner} is not installed or not in PATH")
        print(f"\nTo install {args.aligner}:")
        if args.aligner == 'mafft':
            print("  Ubuntu/Debian: sudo apt-get install mafft")
            print("  macOS: brew install mafft")
            print("  Conda: conda install -c bioconda mafft")
        else:
            print("  Ubuntu/Debian: sudo apt-get install muscle")
            print("  macOS: brew install muscle")
            print("  Conda: conda install -c bioconda muscle")
        return
    
    print("="*60)
    print("MITOCHONDRIAL GENE ALIGNMENT")
    print("="*60)
    print(f"Input directory: {args.input}")
    print(f"Output directory: {args.output}")
    print(f"Aligner: {args.aligner}")
    print(f"Sequence type: {args.type}")
    
    align_genes(args.input, args.output, args.aligner, args.type)
    
    print(f"\n✓ Alignment complete!")
    print(f"Aligned sequences saved to: {args.output}")

if __name__ == "__main__":
    main()
