#!/usr/bin/env python3
"""
CAAS Discovery using PRE-DEFINED species lists from text files.
Reads long-lived and short-lived species from your target files.
"""

import pandas as pd
import numpy as np
from pathlib import Path
from Bio import AlignIO
import argparse
import json
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

MT_GENES = ['ATP6', 'ATP8', 'COX1', 'COX2', 'COX3', 'CYTB',
            'ND1', 'ND2', 'ND3', 'ND4', 'ND4L', 'ND5', 'ND6']

class CAASDiscoveryFromLists:
    """
    CAAS Discovery using pre-defined species lists.
    Reads species from long_lived_targets.txt and short_lived_targets.txt
    """
    
    def __init__(self, alignment_dir, long_lived_file, short_lived_file, 
                 alignment_type='protein', lq_data_file=None):
        """
        Parameters:
        -----------
        alignment_dir : str
            Path to directory containing alignments
        long_lived_file : str
            Path to file with long-lived species IDs (one per line)
        short_lived_file : str
            Path to file with short-lived species IDs (one per line)
        alignment_type : str
            'protein' or 'nucleotide'
        lq_data_file : str, optional
            Path to CSV with longevity quotient data (for validation if available)
        """
        self.base_dir = Path(alignment_dir)
        self.alignment_type = alignment_type
        
        if alignment_type == 'protein':
            self.align_dir = self.base_dir / "proteins"
        else:
            self.align_dir = self.base_dir / "nucleotides"
        
        # Read species lists from files
        self._load_species_lists(long_lived_file, short_lived_file)
        
        # Optionally load LQ data for validation
        self.lq_data = None
        self.lq_dict = {}
        if lq_data_file and Path(lq_data_file).exists():
            print(f"\nLoading LQ data from: {lq_data_file}")
            self.lq_data = pd.read_csv(lq_data_file)
            print(f"  Columns in LQ file: {list(self.lq_data.columns)}")
            print(f"  Shape: {self.lq_data.shape}")
            
            # Try to find the right columns
            if 'species_id' in self.lq_data.columns and 'longevity_quotient' in self.lq_data.columns:
                self.lq_dict = dict(zip(
                    self.lq_data['species_id'],
                    self.lq_data['longevity_quotient']
                ))
                print(f"  ✓ Successfully loaded LQ data for {len(self.lq_dict)} species")
            else:
                # Show what columns exist to help debug
                print(f"  ✗ Expected columns 'species_id' and 'longevity_quotient' not found!")
                print(f"  Available columns: {list(self.lq_data.columns)}")
                print(f"  First few rows:")
                print(self.lq_data.head())
        elif lq_data_file:
            print(f"\nWARNING: LQ data file not found: {lq_data_file}")
        
        # Detect which species are in alignments
        self._detect_alignment_species()
        
        # Create groups
        self._create_groups()
        
        self.caas_results = []
        self.permutation_pvalue = None
        
    def _load_species_lists(self, long_lived_file, short_lived_file):
        """Load species IDs from text files."""
        # Read long-lived species
        with open(long_lived_file, 'r') as f:
            self.target_long = set(line.strip() for line in f if line.strip())
        
        # Read short-lived species
        with open(short_lived_file, 'r') as f:
            self.target_short = set(line.strip() for line in f if line.strip())
        
        print(f"\n{'='*70}")
        print(f"LOADING SPECIES LISTS")
        print(f"{'='*70}")
        print(f"Long-lived file:  {long_lived_file}")
        print(f"  → {len(self.target_long)} species")
        print(f"Short-lived file: {short_lived_file}")
        print(f"  → {len(self.target_short)} species")
        
        # Check for overlap (should be none)
        overlap = self.target_long & self.target_short
        if overlap:
            print(f"WARNING: {len(overlap)} species in both lists: {overlap}")
        
    def _detect_alignment_species(self):
        """Detect which target species are in alignment files."""
        sample_files = list(self.align_dir.glob("*_aligned.fasta"))
        if not sample_files:
            raise FileNotFoundError(f"No alignment files in {self.align_dir}")
        
        # Read first alignment to see which species exist
        sample_aln = AlignIO.read(sample_files[0], "fasta")
        self.alignment_species = {rec.id for rec in sample_aln}
        
        print(f"\nAlignment file check: {sample_files[0].name}")
        print(f"  → {len(self.alignment_species)} species in alignments")
        
    def _create_groups(self):
        """
        Create discovery groups from the target lists.
        Only use species that are in both the target lists AND alignments.
        """
        # Long-lived: intersection of target list and alignment species
        self.discovery_long = self.target_long & self.alignment_species
        
        # Short-lived: intersection of target list and alignment species
        self.discovery_short = self.target_short & self.alignment_species
        
        # Validation species: any other species in alignments (not in target lists)
        self.validation_species = (
            self.alignment_species - self.target_long - self.target_short
        )
        
        print(f"\n{'='*70}")
        print(f"CAAS DISCOVERY - {self.alignment_type.upper()}")
        print(f"{'='*70}")
        print(f"Alignment Directory: {self.align_dir}")
        print(f"")
        print(f"Discovery groups (from target files):")
        print(f"  Long-lived:     {len(self.discovery_long)} species")
        for sp in sorted(self.discovery_long):
            print(f"    - {sp}")
        print(f"  Short-lived:    {len(self.discovery_short)} species")
        for sp in sorted(self.discovery_short):
            print(f"    - {sp}")
        
        if self.validation_species:
            print(f"  Validation:     {len(self.validation_species)} species")
            for sp in sorted(self.validation_species):
                print(f"    - {sp}")
        else:
            print(f"  Validation:     0 species (all species used for discovery)")
        
        # Check if any target species are missing from alignments
        missing_long = self.target_long - self.alignment_species
        missing_short = self.target_short - self.alignment_species
        
        if missing_long:
            print(f"\nWARNING: {len(missing_long)} long-lived targets NOT in alignments:")
            for sp in sorted(missing_long):
                print(f"    - {sp}")
        
        if missing_short:
            print(f"\nWARNING: {len(missing_short)} short-lived targets NOT in alignments:")
            for sp in sorted(missing_short):
                print(f"    - {sp}")
        
        print(f"{'='*70}\n")
        
        # Validation check
        if len(self.discovery_long) < 2:
            raise ValueError(
                f"Only {len(self.discovery_long)} long-lived species found in alignments. "
                f"Need at least 2."
            )
        if len(self.discovery_short) < 2:
            raise ValueError(
                f"Only {len(self.discovery_short)} short-lived species found in alignments. "
                f"Need at least 2."
            )
    
    def load_alignment(self, gene):
        """Load alignment for a gene."""
        path = self.align_dir / f"{gene}_aligned.fasta"
        if not path.exists():
            return None
        return AlignIO.read(path, "fasta")
    
    def check_scenario_1(self, long_chars, short_chars):
        """Scenario 1: Fixed in both, different."""
        long_set = set(long_chars.values())
        short_set = set(short_chars.values())
        
        if len(long_set) == 1 and len(short_set) == 1:
            long_char = list(long_set)[0]
            short_char = list(short_set)[0]
            if long_char != short_char:
                return ('scenario_1', long_char, short_char)
        return None
    
    def check_scenario_2(self, long_chars, short_chars):
        """Scenario 2: Fixed long, variable short, non-overlapping."""
        long_set = set(long_chars.values())
        short_set = set(short_chars.values())
        
        if len(long_set) == 1 and len(short_set) > 1:
            long_char = list(long_set)[0]
            if long_char not in short_set:
                short_str = "".join(sorted(list(short_set)))
                return ('scenario_2', long_char, short_str)
        return None
    
    def check_scenario_3(self, long_chars, short_chars):
        """Scenario 3: Variable long, fixed short, non-overlapping."""
        long_set = set(long_chars.values())
        short_set = set(short_chars.values())
        
        if len(short_set) == 1 and len(long_set) > 1:
            short_char = list(short_set)[0]
            if short_char not in long_set:
                long_str = "".join(sorted(list(long_set)))
                return ('scenario_3', long_str, short_char)
        return None
    
    def discover_caas(self):
        """Discovery phase."""
        print("Starting CAAS discovery...")
        
        for gene in MT_GENES:
            aln = self.load_alignment(gene)
            if aln is None:
                print(f"  [SKIP] {gene} - file not found")
                continue
            
            gene_caas = {'s1': 0, 's2': 0, 's3': 0}
            
            for pos in range(aln.get_alignment_length()):
                long_chars = {
                    r.id: str(r.seq[pos]).upper() 
                    for r in aln if r.id in self.discovery_long
                }
                short_chars = {
                    r.id: str(r.seq[pos]).upper() 
                    for r in aln if r.id in self.discovery_short
                }
                
                # Skip if any discovery species missing
                if (len(long_chars) != len(self.discovery_long) or 
                    len(short_chars) != len(self.discovery_short)):
                    continue
                
                # Skip gaps
                if '-' in long_chars.values() or '-' in short_chars.values():
                    continue
                
                # Check scenarios
                result = None
                for check_func in [self.check_scenario_1, 
                                  self.check_scenario_2, 
                                  self.check_scenario_3]:
                    result = check_func(long_chars, short_chars)
                    if result:
                        break
                
                if result:
                    scenario, long_val, short_val = result
                    
                    self.caas_results.append({
                        'gene': gene,
                        'position': pos + 1,
                        'scenario': scenario,
                        'long_lived_char': long_val,
                        'short_lived_char': short_val,
                        'alignment_type': self.alignment_type
                    })
                    
                    if scenario == 'scenario_1':
                        gene_caas['s1'] += 1
                    elif scenario == 'scenario_2':
                        gene_caas['s2'] += 1
                    else:
                        gene_caas['s3'] += 1
            
            total = sum(gene_caas.values())
            print(f"  [DONE] {gene}: {total} CAAS "
                  f"(S1={gene_caas['s1']}, S2={gene_caas['s2']}, S3={gene_caas['s3']})")
        
        # Summary
        s1 = sum(1 for c in self.caas_results if c['scenario'] == 'scenario_1')
        s2 = sum(1 for c in self.caas_results if c['scenario'] == 'scenario_2')
        s3 = sum(1 for c in self.caas_results if c['scenario'] == 'scenario_3')
        
        print(f"\nDiscovery Summary:")
        print(f"  Total CAAS:      {len(self.caas_results)}")
        print(f"  Scenario 1:      {s1}")
        print(f"  Scenario 2:      {s2}")
        print(f"  Scenario 3:      {s3}")
        if s3 > 0:
            print(f"  S2:S3 ratio:     {s2/s3:.2f}")
            print(f"  (Farré et al. found S2:S3 = 4.6:1, suggesting evolutionary trend to longevity)")
    
    def validate_caas(self):
        """Validation phase - only if we have validation species and LQ data."""
        if not self.validation_species:
            print("\nNo validation species (all species used for discovery)")
            for caas in self.caas_results:
                caas['validated'] = None
                caas['validation_note'] = "no_validation_species"
            return
        
        if not self.lq_dict:
            print("\nNo LQ data available - skipping validation")
            for caas in self.caas_results:
                caas['validated'] = None
                caas['validation_note'] = "no_lq_data"
            return
        
        print(f"\nStarting validation with {len(self.validation_species)} species...")
        
        validated = 0
        
        for caas in self.caas_results:
            gene = caas['gene']
            pos = caas['position'] - 1
            scenario = caas['scenario']
            
            aln = self.load_alignment(gene)
            if aln is None:
                continue
            
            validation_chars = {
                r.id: str(r.seq[pos]).upper()
                for r in aln if r.id in self.validation_species
            }
            validation_chars = {k: v for k, v in validation_chars.items() if v != '-'}
            
            if len(validation_chars) < 2:
                caas['validated'] = False
                caas['validation_pvalue'] = np.nan
                continue
            
            long_char = caas['long_lived_char']
            short_char = caas['short_lived_char']
            
            if scenario in ['scenario_1', 'scenario_2']:
                group_long = [
                    self.lq_dict[sp] for sp, char in validation_chars.items()
                    if char == long_char and sp in self.lq_dict
                ]
                group_short = [
                    self.lq_dict[sp] for sp, char in validation_chars.items()
                    if char in short_char and sp in self.lq_dict
                ]
            else:  # scenario_3
                group_long = [
                    self.lq_dict[sp] for sp, char in validation_chars.items()
                    if char in long_char and sp in self.lq_dict
                ]
                group_short = [
                    self.lq_dict[sp] for sp, char in validation_chars.items()
                    if char == short_char and sp in self.lq_dict
                ]
            
            if len(group_long) < 1 or len(group_short) < 1:
                caas['validated'] = False
                caas['validation_pvalue'] = np.nan
                continue
            
            try:
                statistic, pvalue = stats.mannwhitneyu(
                    group_long, group_short, 
                    alternative='greater'
                )
                
                caas['validated'] = pvalue < 0.2
                caas['validation_pvalue'] = pvalue
                caas['n_long_val'] = len(group_long)
                caas['n_short_val'] = len(group_short)
                
                if caas['validated']:
                    validated += 1
                    
            except Exception as e:
                caas['validated'] = False
                caas['validation_pvalue'] = np.nan
        
        total_tested = sum(1 for c in self.caas_results if 'validated' in c and c['validated'] is not None)
        if total_tested > 0:
            print(f"  Validated:       {validated}/{total_tested} ({validated/total_tested*100:.1f}%)")
    
    def permutation_test(self, n_permutations=1000):
        """Permutation test."""
        print(f"\nRunning permutation test ({n_permutations} permutations)...")
        
        all_species = list(self.discovery_long | self.discovery_short)
        n_long = len(self.discovery_long)
        
        alignments = {}
        for gene in MT_GENES:
            aln = self.load_alignment(gene)
            if aln:
                alignments[gene] = aln
        
        observed = len([
            c for c in self.caas_results 
            if c['scenario'] in ['scenario_1', 'scenario_2']
        ])
        
        null_counts = []
        
        for i in range(n_permutations):
            if (i + 1) % 100 == 0:
                print(f"  Permutation {i + 1}/{n_permutations}...", end='\r')
            
            shuffled = np.random.permutation(all_species)
            rand_long = set(shuffled[:n_long])
            rand_short = set(shuffled[n_long:])
            
            perm_count = 0
            
            for gene, aln in alignments.items():
                for pos in range(aln.get_alignment_length()):
                    long_chars = {
                        r.id: str(r.seq[pos]).upper()
                        for r in aln if r.id in rand_long
                    }
                    short_chars = {
                        r.id: str(r.seq[pos]).upper()
                        for r in aln if r.id in rand_short
                    }
                    
                    if (len(long_chars) != n_long or 
                        len(short_chars) != len(rand_short)):
                        continue
                    if '-' in long_chars.values() or '-' in short_chars.values():
                        continue
                    
                    if (self.check_scenario_1(long_chars, short_chars) or
                        self.check_scenario_2(long_chars, short_chars)):
                        perm_count += 1
            
            null_counts.append(perm_count)
        
        print()
        
        self.permutation_pvalue = sum(
            1 for x in null_counts if x >= observed
        ) / n_permutations
        
        print(f"  Observed (S1+S2):    {observed}")
        print(f"  Mean null:           {np.mean(null_counts):.1f} ± {np.std(null_counts):.1f}")
        print(f"  P-value:             {self.permutation_pvalue:.4f}")
        
        if self.permutation_pvalue < 0.15:
            print(f"  *** SIGNIFICANT: More CAAS than expected by chance ***")
        else:
            print(f"  Not significant at p < 0.15")
        
        return null_counts
    
    def save_results(self, output_dir):
        """Save results."""
        out = Path(output_dir)
        out.mkdir(parents=True, exist_ok=True)
        
        if self.caas_results:
            df_all = pd.DataFrame(self.caas_results)
            df_all.to_csv(
                out / f"caas_discovered_{self.alignment_type}.csv",
                index=False
            )
            print(f"\n[SAVED] {out}/caas_discovered_{self.alignment_type}.csv")
            
            validated = [c for c in self.caas_results if c.get('validated', False)]
            if validated:
                df_val = pd.DataFrame(validated)
                df_val.to_csv(
                    out / f"caas_validated_{self.alignment_type}.csv",
                    index=False
                )
                print(f"[SAVED] {out}/caas_validated_{self.alignment_type}.csv")
        
        summary = {
            'alignment_type': self.alignment_type,
            'n_total_species': len(self.alignment_species),
            'n_long_lived': len(self.discovery_long),
            'n_short_lived': len(self.discovery_short),
            'n_validation': len(self.validation_species),
            'total_caas': len(self.caas_results),
            'scenario_1': len([c for c in self.caas_results if c['scenario'] == 'scenario_1']),
            'scenario_2': len([c for c in self.caas_results if c['scenario'] == 'scenario_2']),
            'scenario_3': len([c for c in self.caas_results if c['scenario'] == 'scenario_3']),
            'validated_caas': len([c for c in self.caas_results if c.get('validated', False)]),
            'permutation_pvalue': self.permutation_pvalue,
            'long_lived_species': sorted(list(self.discovery_long)),
            'short_lived_species': sorted(list(self.discovery_short))
        }
        
        with open(out / f"summary_{self.alignment_type}.json", 'w') as f:
            json.dump(summary, f, indent=2)
        print(f"[SAVED] {out}/summary_{self.alignment_type}.json")


def main():
    parser = argparse.ArgumentParser(
        description='CAAS Discovery using pre-defined species lists from text files'
    )
    parser.add_argument(
        '-a', '--alignments',
        default='data/alignments',
        help='Directory containing protein and nucleotide alignment folders'
    )
    parser.add_argument(
        '-ll', '--long-lived',
        default='data/LQ/long_lived_targets.txt',
        help='Text file with long-lived species IDs (one per line)'
    )
    parser.add_argument(
        '-sl', '--short-lived',
        default='data/LQ/short_lived_targets.txt',
        help='Text file with short-lived species IDs (one per line)'
    )
    parser.add_argument(
        '-l', '--lq',
        default='data/LQ/lq_results.csv',
        help='CSV with longevity quotient data for validation (columns: species_id, longevity_quotient)'
    )
    parser.add_argument(
        '-o', '--output',
        default='results/caas_discovery',
        help='Output directory'
    )
    parser.add_argument(
        '-p', '--permutations',
        type=int,
        default=1000,
        help='Number of permutations for significance testing'
    )
    parser.add_argument(
        '--protein-only',
        action='store_true',
        help='Only analyze protein alignments'
    )
    parser.add_argument(
        '--nucleotide-only',
        action='store_true',
        help='Only analyze nucleotide alignments'
    )
    
    args = parser.parse_args()
    
    # Check that list files exist
    if not Path(args.long_lived).exists():
        raise FileNotFoundError(f"Long-lived species file not found: {args.long_lived}")
    if not Path(args.short_lived).exists():
        raise FileNotFoundError(f"Short-lived species file not found: {args.short_lived}")
    
    run_protein = not args.nucleotide_only
    run_nucleotide = not args.protein_only
    
    if run_protein:
        print("\n" + "#"*70)
        print("# PROTEIN ANALYSIS")
        print("#"*70)
        
        prot = CAASDiscoveryFromLists(
            args.alignments, 
            args.long_lived, 
            args.short_lived,
            'protein',
            args.lq
        )
        prot.discover_caas()
        prot.validate_caas()
        prot.permutation_test(args.permutations)
        prot.save_results(args.output)
    
    if run_nucleotide:
        print("\n" + "#"*70)
        print("# NUCLEOTIDE ANALYSIS")
        print("#"*70)
        
        nuc = CAASDiscoveryFromLists(
            args.alignments, 
            args.long_lived, 
            args.short_lived,
            'nucleotide',
            args.lq
        )
        nuc.discover_caas()
        nuc.validate_caas()
        nuc.permutation_test(args.permutations)
        nuc.save_results(args.output)
    
    print("\n" + "="*70)
    print("ANALYSIS COMPLETE")
    print("="*70)


if __name__ == "__main__":
    main()