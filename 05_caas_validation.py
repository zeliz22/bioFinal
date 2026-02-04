#!/usr/bin/env python3
"""
Script 5: CAAS Validation Phase

Validates discovered CAAS using intermediate species with phylogenetic ANOVA.
For each discovered CAAS position:
1. Group intermediate species by whether they have long-lived or short-lived AA
2. Test if the group with long-lived AA has significantly higher LQ
3. Keep only validated CAAS (FDR < 0.05)
"""

import pandas as pd
import numpy as np
from pathlib import Path
from Bio import AlignIO
from scipy import stats
from statsmodels.stats.multitest import fdrcorrection
import argparse
import json

class CAASValidator:
    def __init__(self, alignments_dir, lq_data_file, caas_file):
        """
        Initialize CAAS validator
        
        Args:
            alignments_dir: Directory with aligned protein sequences
            lq_data_file: CSV file with LQ data
            caas_file: CSV file with discovered CAAS
        """
        self.alignments_dir = Path(alignments_dir)
        self.lq_data = pd.read_csv(lq_data_file)
        self.caas_data = pd.read_csv(caas_file)
        
        # Get intermediate species
        self.intermediate = self.lq_data[
            self.lq_data['longevity_class'] == 'intermediate'
        ].copy()
        
        # Create LQ lookup
        self.lq_lookup = dict(zip(self.lq_data['species_id'], self.lq_data['LQ']))
        
        print(f"Loaded {len(self.intermediate)} intermediate species for validation")
        print(f"Loaded {len(self.caas_data)} CAAS to validate")
        
        self.validation_results = []
    
    def load_alignment(self, gene):
        """Load protein alignment for a gene"""
        alignment_file = self.alignments_dir / f"{gene}_protein_aligned.fasta"
        
        if not alignment_file.exists():
            return None
        
        try:
            return AlignIO.read(alignment_file, "fasta")
        except Exception as e:
            print(f"Error loading {gene}: {e}")
            return None
    
    def get_aa_at_position(self, alignment, position):
        """
        Get amino acid at position for all species in alignment
        
        Returns:
            dict: {species_id: amino_acid}
        """
        aa_dict = {}
        for record in alignment:
            species_id = record.id
            # Position is 1-based in CAAS data, convert to 0-based
            aa = str(record.seq[position - 1]).upper()
            aa_dict[species_id] = aa
        return aa_dict
    
    def validate_caas_position(self, caas_row):
        """
        Validate a single CAAS position using intermediate species
        
        Returns:
            dict: Validation results including p-value
        """
        gene = caas_row['gene']
        position = caas_row['position']
        scenario = caas_row['scenario']
        
        # Get long-lived and short-lived AAs
        long_lived_aas = set(caas_row['long_lived_aa'].split(','))
        short_lived_aas = set(caas_row['short_lived_aa'].split(','))
        
        # Load alignment
        alignment = self.load_alignment(gene)
        if alignment is None:
            return None
        
        # Get AAs for intermediate species
        all_aas = self.get_aa_at_position(alignment, position)
        
        # Classify intermediate species
        has_long_aa = []
        has_short_aa = []
        has_other_aa = []
        
        for species_id in self.intermediate['species_id']:
            if species_id not in all_aas:
                continue
            
            aa = all_aas[species_id]
            
            if aa == '-':  # Skip gaps
                continue
            
            if aa in long_lived_aas:
                has_long_aa.append(species_id)
            elif aa in short_lived_aas:
                has_short_aa.append(species_id)
            else:
                has_other_aa.append(species_id)
        
        # Need at least 2 species in each group for meaningful test
        if len(has_long_aa) < 2 or len(has_short_aa) < 2:
            return {
                'gene': gene,
                'position': position,
                'scenario': scenario,
                'validated': False,
                'reason': 'insufficient_data',
                'n_long_aa': len(has_long_aa),
                'n_short_aa': len(has_short_aa),
                'n_other_aa': len(has_other_aa),
                'p_value': np.nan,
                'direction': np.nan
            }
        
        # Get LQ values for each group
        long_aa_lqs = [self.lq_lookup[s] for s in has_long_aa if s in self.lq_lookup]
        short_aa_lqs = [self.lq_lookup[s] for s in has_short_aa if s in self.lq_lookup]
        
        # Perform t-test
        # We expect long_aa_lqs > short_aa_lqs
        try:
            # One-tailed t-test: long-lived AA should have higher LQ
            t_stat, p_value_two_tailed = stats.ttest_ind(long_aa_lqs, short_aa_lqs)
            
            # Convert to one-tailed
            if np.mean(long_aa_lqs) > np.mean(short_aa_lqs):
                p_value = p_value_two_tailed / 2
                direction = 'correct'
            else:
                p_value = 1 - (p_value_two_tailed / 2)
                direction = 'wrong'
            
            return {
                'gene': gene,
                'position': position,
                'scenario': scenario,
                'validated': False,  # Will be set after FDR correction
                'reason': None,
                'n_long_aa': len(has_long_aa),
                'n_short_aa': len(has_short_aa),
                'n_other_aa': len(has_other_aa),
                'mean_lq_long_aa': np.mean(long_aa_lqs),
                'mean_lq_short_aa': np.mean(short_aa_lqs),
                'p_value': p_value,
                'direction': direction,
                't_statistic': t_stat
            }
        
        except Exception as e:
            print(f"Error in t-test for {gene} position {position}: {e}")
            return None
    
    def validate_all_caas(self):
        """Validate all discovered CAAS"""
        print("\n" + "="*60)
        print("VALIDATING CAAS")
        print("="*60)
        
        validation_results = []
        
        for idx, row in self.caas_data.iterrows():
            if (idx + 1) % 50 == 0:
                print(f"Validating {idx + 1}/{len(self.caas_data)}...")
            
            result = self.validate_caas_position(row)
            if result is not None:
                validation_results.append(result)
        
        self.validation_results = validation_results
        
        # Apply FDR correction
        self.apply_fdr_correction()
        
        return validation_results
    
    def apply_fdr_correction(self, alpha=0.05):
        """Apply FDR correction to p-values"""
        print("\nApplying FDR correction...")
        
        # Get p-values (only for those with correct direction)
        testable = [r for r in self.validation_results 
                   if not np.isnan(r['p_value']) and r['direction'] == 'correct']
        
        if not testable:
            print("No testable CAAS positions found!")
            return
        
        p_values = [r['p_value'] for r in testable]
        
        # FDR correction
        rejected, corrected_p = fdrcorrection(p_values, alpha=alpha)
        
        # Update validation results
        testable_idx = 0
        for r in self.validation_results:
            if not np.isnan(r['p_value']) and r['direction'] == 'correct':
                r['fdr_corrected_p'] = corrected_p[testable_idx]
                r['validated'] = rejected[testable_idx]
                testable_idx += 1
            else:
                r['fdr_corrected_p'] = np.nan
                r['validated'] = False
        
        # Summary
        n_validated = sum(1 for r in self.validation_results if r['validated'])
        n_tested = len(testable)
        n_total = len(self.validation_results)
        
        print(f"\nValidation summary:")
        print(f"  Total CAAS: {n_total}")
        print(f"  Testable (correct direction): {n_tested}")
        print(f"  Validated (FDR < {alpha}): {n_validated} ({n_validated/n_total*100:.1f}%)")
        
        if n_tested > 0:
            print(f"  Validation rate: {n_validated/n_tested*100:.1f}%")
    
    def save_results(self, output_dir):
        """Save validation results"""
        output_path = Path(output_dir)
        output_path.mkdir(exist_ok=True, parents=True)
        
        # Save detailed validation results
        results_df = pd.DataFrame(self.validation_results)
        results_file = output_path / "caas_validation_results.csv"
        results_df.to_csv(results_file, index=False)
        print(f"\n✓ Validation results saved to: {results_file}")
        
        # Save only validated CAAS
        validated_df = results_df[results_df['validated'] == True].copy()
        validated_file = output_path / "caas_validated.csv"
        validated_df.to_csv(validated_file, index=False)
        print(f"✓ Validated CAAS saved to: {validated_file}")
        
        # Generate summary by gene
        gene_summary = {}
        for gene in results_df['gene'].unique():
            gene_data = results_df[results_df['gene'] == gene]
            gene_summary[gene] = {
                'total_caas': len(gene_data),
                'validated_caas': sum(gene_data['validated']),
                'validation_rate': sum(gene_data['validated']) / len(gene_data) * 100
            }
        
        summary_file = output_path / "validation_summary.json"
        with open(summary_file, 'w') as f:
            json.dump(gene_summary, f, indent=2)
        print(f"✓ Summary saved to: {summary_file}")
        
        # Generate summary statistics
        summary_stats = {
            'total_discovered': len(self.caas_data),
            'total_tested': len(self.validation_results),
            'total_validated': sum(r['validated'] for r in self.validation_results),
            'validation_rate': sum(r['validated'] for r in self.validation_results) / len(self.validation_results) * 100,
            'by_scenario': {}
        }
        
        for scenario in ['scenario_1', 'scenario_2', 'scenario_3']:
            scenario_data = [r for r in self.validation_results if r['scenario'] == scenario]
            if scenario_data:
                summary_stats['by_scenario'][scenario] = {
                    'tested': len(scenario_data),
                    'validated': sum(r['validated'] for r in scenario_data),
                    'rate': sum(r['validated'] for r in scenario_data) / len(scenario_data) * 100
                }
        
        # Save text summary
        summary_text = output_path / "validation_summary.txt"
        with open(summary_text, 'w') as f:
            f.write("CAAS VALIDATION SUMMARY\n")
            f.write("="*60 + "\n\n")
            f.write(f"Total CAAS discovered: {summary_stats['total_discovered']}\n")
            f.write(f"Total CAAS tested: {summary_stats['total_tested']}\n")
            f.write(f"Total CAAS validated: {summary_stats['total_validated']}\n")
            f.write(f"Overall validation rate: {summary_stats['validation_rate']:.1f}%\n\n")
            
            f.write("BY SCENARIO:\n")
            f.write("-"*60 + "\n")
            for scenario, stats in summary_stats['by_scenario'].items():
                f.write(f"{scenario}:\n")
                f.write(f"  Tested: {stats['tested']}\n")
                f.write(f"  Validated: {stats['validated']}\n")
                f.write(f"  Rate: {stats['rate']:.1f}%\n\n")
            
            f.write("\nBY GENE:\n")
            f.write("-"*60 + "\n")
            for gene, stats in sorted(gene_summary.items()):
                f.write(f"{gene}: {stats['validated_caas']}/{stats['total_caas']} "
                       f"({stats['validation_rate']:.1f}%)\n")
        
        print(f"✓ Text summary saved to: {summary_text}")
        
        return summary_stats

def main():
    parser = argparse.ArgumentParser(
        description='Validate discovered CAAS using intermediate species'
    )
    parser.add_argument('-a', '--alignments', default='data/alignments',
                       help='Directory with aligned protein sequences')
    parser.add_argument('-l', '--lq-data', default='data/phenotypes/lq_data.csv',
                       help='CSV file with LQ data')
    parser.add_argument('-c', '--caas', default='results/caas_discovery/caas_discovered.csv',
                       help='CSV file with discovered CAAS')
    parser.add_argument('-o', '--output', default='results/caas_validation',
                       help='Output directory')
    parser.add_argument('--fdr', type=float, default=0.05,
                       help='FDR threshold (default: 0.05)')
    
    args = parser.parse_args()
    
    print("="*60)
    print("CAAS VALIDATION PHASE")
    print("="*60)
    
    # Initialize validator
    validator = CAASValidator(args.alignments, args.lq_data, args.caas)
    
    # Validate all CAAS
    validator.validate_all_caas()
    
    # Save results
    summary = validator.save_results(args.output)
    
    print("\n" + "="*60)
    print("VALIDATION COMPLETE")
    print("="*60)
    print(f"Validated CAAS: {summary['total_validated']}/{summary['total_tested']}")
    print(f"Validation rate: {summary['validation_rate']:.1f}%")
    print("="*60)

if __name__ == "__main__":
    main()
