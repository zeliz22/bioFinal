#!/usr/bin/env python3
"""
Script 4: CAAS Discovery Phase

Identifies Convergent Amino Acid Substitutions (CAAS) shared by long-lived mammals.
Implements the three scenarios from Farré et al. 2021:

Scenario 1: All long-lived have same AA, all short-lived have different fixed AA
Scenario 2: All long-lived have same AA, short-lived have variable AAs (all different from long-lived)
Scenario 3: Short-lived have fixed AA, long-lived have variable AAs
"""

import pandas as pd
import numpy as np
from pathlib import Path
from Bio import AlignIO
from collections import Counter, defaultdict
import argparse
import json
from scipy import stats

MT_GENES = ['ATP6', 'ATP8', 'COX1', 'COX2', 'COX3', 'CYTB',
            'ND1', 'ND2', 'ND3', 'ND4', 'ND4L', 'ND5', 'ND6']

class CAASDiscovery:
    def __init__(self, alignment_dir, lq_data_file):
        """
        Initialize CAAS discovery
        
        Args:
            alignment_dir: Directory with aligned protein sequences
            lq_data_file: CSV file with LQ data and species classification
        """
        self.alignment_dir = Path(alignment_dir)
        self.lq_data = pd.read_csv(lq_data_file)
        
        # Get species lists for each category
        self.long_lived = set(self.lq_data[self.lq_data['longevity_class'] == 'long-lived']['species_id'])
        self.short_lived = set(self.lq_data[self.lq_data['longevity_class'] == 'short-lived']['species_id'])
        self.intermediate = set(self.lq_data[self.lq_data['longevity_class'] == 'intermediate']['species_id'])
        
        print(f"Loaded {len(self.long_lived)} long-lived species")
        print(f"Loaded {len(self.short_lived)} short-lived species")
        print(f"Loaded {len(self.intermediate)} intermediate species")
        
        self.caas_results = []
        self.gene_summary = {}
    
    def load_alignment(self, gene):
        """Load protein alignment for a gene"""
        alignment_file = self.alignment_dir / f"{gene}_protein_aligned.fasta"
        
        if not alignment_file.exists():
            print(f"  Warning: Alignment file not found: {alignment_file}")
            return None
        
        try:
            alignment = AlignIO.read(alignment_file, "fasta")
            return alignment
        except Exception as e:
            print(f"  Error loading alignment: {e}")
            return None
    
    def get_aa_at_position(self, alignment, position, species_set):
        """
        Get amino acids at a specific position for a set of species
        
        Returns:
            dict: {species_id: amino_acid}
        """
        aa_dict = {}
        
        for record in alignment:
            species_id = record.id
            if species_id in species_set:
                aa = str(record.seq[position]).upper()
                aa_dict[species_id] = aa
        
        return aa_dict
    
    def check_scenario_1(self, long_aa, short_aa):
        """
        Scenario 1: All long-lived have same AA, all short-lived have different fixed AA
        
        Returns:
            bool: True if scenario 1 is satisfied
        """
        # Check no gaps
        if '-' in long_aa.values() or '-' in short_aa.values():
            return False
        
        # All long-lived must have same AA
        long_aas = set(long_aa.values())
        if len(long_aas) != 1:
            return False
        
        # All short-lived must have same AA
        short_aas = set(short_aa.values())
        if len(short_aas) != 1:
            return False
        
        # The AAs must be different
        if long_aas == short_aas:
            return False
        
        return True
    
    def check_scenario_2(self, long_aa, short_aa):
        """
        Scenario 2: All long-lived have same AA, short-lived have variable AAs 
                    (all different from long-lived)
        
        Returns:
            bool: True if scenario 2 is satisfied
        """
        # Check no gaps
        if '-' in long_aa.values() or '-' in short_aa.values():
            return False
        
        # All long-lived must have same AA
        long_aas = set(long_aa.values())
        if len(long_aas) != 1:
            return False
        
        # Short-lived must have variable AAs
        short_aas = set(short_aa.values())
        if len(short_aas) <= 1:
            return False
        
        # None of the short-lived AAs should be in long-lived
        long_aa_value = list(long_aas)[0]
        if long_aa_value in short_aas:
            return False
        
        return True
    
    def check_scenario_3(self, long_aa, short_aa):
        """
        Scenario 3: Short-lived have fixed AA, long-lived have variable AAs
        
        Returns:
            bool: True if scenario 3 is satisfied
        """
        # Check no gaps
        if '-' in long_aa.values() or '-' in short_aa.values():
            return False
        
        # All short-lived must have same AA
        short_aas = set(short_aa.values())
        if len(short_aas) != 1:
            return False
        
        # Long-lived must have variable AAs
        long_aas = set(long_aa.values())
        if len(long_aas) <= 1:
            return False
        
        # None of the long-lived AAs should be in short-lived
        short_aa_value = list(short_aas)[0]
        if short_aa_value in long_aas:
            return False
        
        return True
    
    def discover_caas_in_gene(self, gene):
        """
        Discover CAAS in a specific gene
        
        Returns:
            list: List of CAAS dictionaries
        """
        print(f"\n{'='*60}")
        print(f"Analyzing {gene}")
        print(f"{'='*60}")
        
        alignment = self.load_alignment(gene)
        if alignment is None:
            return []
        
        alignment_length = alignment.get_alignment_length()
        print(f"Alignment length: {alignment_length} positions")
        
        gene_caas = []
        scenario_counts = {'scenario_1': 0, 'scenario_2': 0, 'scenario_3': 0}
        
        # Analyze each position
        for pos in range(alignment_length):
            # Get AAs for long-lived and short-lived species
            long_aa = self.get_aa_at_position(alignment, pos, self.long_lived)
            short_aa = self.get_aa_at_position(alignment, pos, self.short_lived)
            
            # Skip if not all species present
            if len(long_aa) != len(self.long_lived) or len(short_aa) != len(self.short_lived):
                continue
            
            # Check scenarios
            scenario = None
            if self.check_scenario_1(long_aa, short_aa):
                scenario = 'scenario_1'
            elif self.check_scenario_2(long_aa, short_aa):
                scenario = 'scenario_2'
            elif self.check_scenario_3(long_aa, short_aa):
                scenario = 'scenario_3'
            
            if scenario:
                scenario_counts[scenario] += 1
                
                caas = {
                    'gene': gene,
                    'position': pos + 1,  # 1-based
                    'scenario': scenario,
                    'long_lived_aa': list(set(long_aa.values())),
                    'short_lived_aa': list(set(short_aa.values())),
                    'long_lived_species': list(long_aa.keys()),
                    'short_lived_species': list(short_aa.keys())
                }
                
                gene_caas.append(caas)
        
        # Summary for this gene
        total_caas = len(gene_caas)
        print(f"Total CAAS found: {total_caas}")
        print(f"  Scenario 1: {scenario_counts['scenario_1']}")
        print(f"  Scenario 2: {scenario_counts['scenario_2']}")
        print(f"  Scenario 3: {scenario_counts['scenario_3']}")
        
        self.gene_summary[gene] = {
            'total_caas': total_caas,
            'scenario_1': scenario_counts['scenario_1'],
            'scenario_2': scenario_counts['scenario_2'],
            'scenario_3': scenario_counts['scenario_3'],
            'alignment_length': alignment_length
        }
        
        return gene_caas
    
    def run_discovery(self):
        """Run CAAS discovery for all genes"""
        print("\n" + "="*60)
        print("STARTING CAAS DISCOVERY")
        print("="*60)
        
        all_caas = []
        
        for gene in MT_GENES:
            gene_caas = self.discover_caas_in_gene(gene)
            all_caas.extend(gene_caas)
        
        self.caas_results = all_caas
        
        return all_caas
    
    def calculate_random_expectations(self, n_permutations=1000):
        """
        Calculate random expectations by permuting species labels
        
        This tests whether the number of CAAS is higher than expected by chance
        """
        print("\n" + "="*60)
        print("CALCULATING RANDOM EXPECTATIONS")
        print("="*60)
        print(f"Running {n_permutations} permutations...")
        
        all_species = list(self.long_lived | self.short_lived)
        n_long = len(self.long_lived)
        n_short = len(self.short_lived)
        
        random_counts = []
        
        for i in range(n_permutations):
            if (i + 1) % 100 == 0:
                print(f"  Permutation {i + 1}/{n_permutations}")
            
            # Randomly shuffle species
            shuffled = np.random.permutation(all_species)
            random_long = set(shuffled[:n_long])
            random_short = set(shuffled[n_long:n_long + n_short])
            
            # Count CAAS with random grouping
            total_random_caas = 0
            
            for gene in MT_GENES:
                alignment = self.load_alignment(gene)
                if alignment is None:
                    continue
                
                for pos in range(alignment.get_alignment_length()):
                    long_aa = self.get_aa_at_position(alignment, pos, random_long)
                    short_aa = self.get_aa_at_position(alignment, pos, random_short)
                    
                    if len(long_aa) != n_long or len(short_aa) != n_short:
                        continue
                    
                    if (self.check_scenario_1(long_aa, short_aa) or 
                        self.check_scenario_2(long_aa, short_aa)):
                        total_random_caas += 1
            
            random_counts.append(total_random_caas)
        
        # Calculate p-value
        observed_caas = sum(1 for c in self.caas_results 
                          if c['scenario'] in ['scenario_1', 'scenario_2'])
        p_value = sum(1 for count in random_counts if count >= observed_caas) / n_permutations
        
        print(f"\nObserved CAAS (Scenarios 1+2): {observed_caas}")
        print(f"Mean random CAAS: {np.mean(random_counts):.1f}")
        print(f"SD random CAAS: {np.std(random_counts):.1f}")
        print(f"P-value: {p_value:.4f}")
        
        if p_value < 0.05:
            print("✓ Significant enrichment of CAAS!")
        else:
            print("⚠ No significant enrichment detected")
        
        return {
            'observed': observed_caas,
            'random_mean': np.mean(random_counts),
            'random_sd': np.std(random_counts),
            'p_value': p_value,
            'random_counts': random_counts
        }
    
    def save_results(self, output_dir):
        """Save CAAS discovery results"""
        output_path = Path(output_dir)
        output_path.mkdir(exist_ok=True, parents=True)
        
        # Save detailed CAAS results
        caas_file = output_path / "caas_discovered.csv"
        
        # Convert to DataFrame
        caas_df = pd.DataFrame([
            {
                'gene': c['gene'],
                'position': c['position'],
                'scenario': c['scenario'],
                'long_lived_aa': ','.join(c['long_lived_aa']),
                'short_lived_aa': ','.join(c['short_lived_aa']),
                'n_long_lived': len(c['long_lived_species']),
                'n_short_lived': len(c['short_lived_species'])
            }
            for c in self.caas_results
        ])
        
        caas_df.to_csv(caas_file, index=False)
        print(f"\n✓ CAAS results saved to: {caas_file}")
        
        # Save gene summary
        summary_file = output_path / "gene_summary.json"
        with open(summary_file, 'w') as f:
            json.dump(self.gene_summary, f, indent=2)
        print(f"✓ Gene summary saved to: {summary_file}")
        
        # Save overall summary
        overall_summary = {
            'total_caas': len(self.caas_results),
            'scenario_1_caas': sum(1 for c in self.caas_results if c['scenario'] == 'scenario_1'),
            'scenario_2_caas': sum(1 for c in self.caas_results if c['scenario'] == 'scenario_2'),
            'scenario_3_caas': sum(1 for c in self.caas_results if c['scenario'] == 'scenario_3'),
            'genes_with_caas': len([g for g in self.gene_summary if self.gene_summary[g]['total_caas'] > 0]),
            'n_long_lived_species': len(self.long_lived),
            'n_short_lived_species': len(self.short_lived)
        }
        
        summary_text_file = output_path / "discovery_summary.txt"
        with open(summary_text_file, 'w') as f:
            f.write("CAAS DISCOVERY SUMMARY\n")
            f.write("="*60 + "\n\n")
            for key, value in overall_summary.items():
                f.write(f"{key}: {value}\n")
            f.write("\n" + "="*60 + "\n")
            f.write("PER-GENE SUMMARY\n")
            f.write("="*60 + "\n")
            for gene in MT_GENES:
                if gene in self.gene_summary:
                    s = self.gene_summary[gene]
                    f.write(f"\n{gene}:\n")
                    f.write(f"  Total CAAS: {s['total_caas']}\n")
                    f.write(f"  Scenario 1: {s['scenario_1']}\n")
                    f.write(f"  Scenario 2: {s['scenario_2']}\n")
                    f.write(f"  Scenario 3: {s['scenario_3']}\n")
                    f.write(f"  Alignment length: {s['alignment_length']}\n")
        
        print(f"✓ Summary saved to: {summary_text_file}")
        
        return overall_summary

def main():
    parser = argparse.ArgumentParser(
        description='Discover convergent amino acid substitutions (CAAS)'
    )
    parser.add_argument('-a', '--alignments', default='data/alignments',
                       help='Directory with aligned protein sequences')
    parser.add_argument('-l', '--lq-data', default='data/phenotypes/lq_data.csv',
                       help='CSV file with LQ data and species classification')
    parser.add_argument('-o', '--output', default='results/caas_discovery',
                       help='Output directory for results')
    parser.add_argument('--permutations', type=int, default=1000,
                       help='Number of permutations for significance testing')
    parser.add_argument('--skip-permutations', action='store_true',
                       help='Skip permutation testing (faster)')
    
    args = parser.parse_args()
    
    print("="*60)
    print("CAAS DISCOVERY PHASE")
    print("="*60)
    
    # Initialize discovery
    discovery = CAASDiscovery(args.alignments, args.lq_data)
    
    # Run discovery
    caas_results = discovery.run_discovery()
    
    # Calculate random expectations (unless skipped)
    if not args.skip_permutations:
        random_stats = discovery.calculate_random_expectations(args.permutations)
    
    # Save results
    summary = discovery.save_results(args.output)
    
    # Print final summary
    print("\n" + "="*60)
    print("DISCOVERY COMPLETE")
    print("="*60)
    print(f"Total CAAS discovered: {summary['total_caas']}")
    print(f"  Scenario 1: {summary['scenario_1_caas']}")
    print(f"  Scenario 2: {summary['scenario_2_caas']}")
    print(f"  Scenario 3: {summary['scenario_3_caas']}")
    print(f"Genes with CAAS: {summary['genes_with_caas']}/{len(MT_GENES)}")
    print("="*60)

if __name__ == "__main__":
    main()
