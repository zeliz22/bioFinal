#!/usr/bin/env python3
"""
Script 3: Calculate Longevity Quotient (LQ)

Calculates the Longevity Quotient for each species and classifies them into:
- Long-lived (top decile)
- Short-lived (bottom decile)
- Intermediate (middle 80%)

LQ = Observed Maximum Lifespan / Expected Maximum Lifespan (based on body mass)

Expected lifespan formula from de Magalhães et al. 2007:
MLS = 4.88 * body_mass^0.19 (for mammals)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import argparse
import json

def calculate_expected_lifespan(body_mass_kg):
    """
    Calculate expected maximum lifespan based on body mass
    
    Formula from de Magalhães et al. 2007 for mammals:
    MLS (years) = 4.88 * body_mass(kg)^0.19
    
    Args:
        body_mass_kg: Body mass in kilograms
    
    Returns:
        Expected maximum lifespan in years
    """
    return 4.88 * (body_mass_kg ** 0.19)

def calculate_longevity_quotient(observed_mls, body_mass_kg):
    """
    Calculate Longevity Quotient (LQ)
    
    LQ = Observed MLS / Expected MLS
    
    LQ > 1: Species lives longer than expected
    LQ < 1: Species lives shorter than expected
    LQ = 1: Species lives as expected for its body size
    """
    expected_mls = calculate_expected_lifespan(body_mass_kg)
    return observed_mls / expected_mls

def load_phenotype_data(phenotype_file):
    """
    Load phenotype data from CSV file
    
    Expected columns:
    - species_id: Unique identifier (matching FASTA file names)
    - species_name: Full species name
    - max_lifespan_years: Maximum lifespan in years
    - body_mass_kg: Adult body mass in kilograms
    """
    df = pd.read_csv(phenotype_file)
    
    # Validate required columns
    required_cols = ['species_id', 'species_name', 'max_lifespan_years', 'body_mass_kg']
    missing_cols = [col for col in required_cols if col not in df.columns]
    
    if missing_cols:
        raise ValueError(f"Missing required columns: {missing_cols}")
    
    # Remove rows with missing data
    initial_rows = len(df)
    df = df.dropna(subset=['max_lifespan_years', 'body_mass_kg'])
    
    if len(df) < initial_rows:
        print(f"Warning: Removed {initial_rows - len(df)} rows with missing data")
    
    return df

def classify_species(df, top_decile=0.1, bottom_decile=0.1):
    """
    Classify species into long-lived, short-lived, and intermediate groups
    
    Args:
        df: DataFrame with LQ values
        top_decile: Fraction for top decile (default 0.1 = top 10%)
        bottom_decile: Fraction for bottom decile (default 0.1 = bottom 10%)
    """
    # Calculate decile thresholds
    top_threshold = df['LQ'].quantile(1 - top_decile)
    bottom_threshold = df['LQ'].quantile(bottom_decile)
    
    # Classify species
    def classify(lq):
        if lq >= top_threshold:
            return 'long-lived'
        elif lq <= bottom_threshold:
            return 'short-lived'
        else:
            return 'intermediate'
    
    df['longevity_class'] = df['LQ'].apply(classify)
    
    return df, top_threshold, bottom_threshold

def generate_summary_statistics(df):
    """Generate summary statistics for the dataset"""
    stats = {
        'total_species': len(df),
        'mean_lifespan': df['max_lifespan_years'].mean(),
        'median_lifespan': df['max_lifespan_years'].median(),
        'min_lifespan': df['max_lifespan_years'].min(),
        'max_lifespan': df['max_lifespan_years'].max(),
        'mean_body_mass': df['body_mass_kg'].mean(),
        'median_body_mass': df['body_mass_kg'].median(),
        'mean_lq': df['LQ'].mean(),
        'median_lq': df['LQ'].median(),
        'n_long_lived': len(df[df['longevity_class'] == 'long-lived']),
        'n_short_lived': len(df[df['longevity_class'] == 'short-lived']),
        'n_intermediate': len(df[df['longevity_class'] == 'intermediate'])
    }
    return stats

def plot_lq_distribution(df, output_file):
    """Plot LQ distribution with classification thresholds"""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # 1. LQ distribution histogram
    ax1 = axes[0, 0]
    ax1.hist(df['LQ'], bins=30, edgecolor='black', alpha=0.7)
    ax1.axvline(df[df['longevity_class'] == 'long-lived']['LQ'].min(),
                color='green', linestyle='--', label='Long-lived threshold')
    ax1.axvline(df[df['longevity_class'] == 'short-lived']['LQ'].max(),
                color='red', linestyle='--', label='Short-lived threshold')
    ax1.set_xlabel('Longevity Quotient (LQ)')
    ax1.set_ylabel('Frequency')
    ax1.set_title('Distribution of Longevity Quotient')
    ax1.legend()
    ax1.grid(alpha=0.3)
    
    # 2. Body mass vs Lifespan
    ax2 = axes[0, 1]
    colors = {'long-lived': 'green', 'short-lived': 'red', 'intermediate': 'gray'}
    for category in ['long-lived', 'short-lived', 'intermediate']:
        subset = df[df['longevity_class'] == category]
        ax2.scatter(subset['body_mass_kg'], subset['max_lifespan_years'],
                   c=colors[category], label=category, alpha=0.6, s=50)
    
    # Add expected lifespan curve
    mass_range = np.logspace(np.log10(df['body_mass_kg'].min()),
                            np.log10(df['body_mass_kg'].max()), 100)
    expected = calculate_expected_lifespan(mass_range)
    ax2.plot(mass_range, expected, 'k--', label='Expected MLS', linewidth=2)
    
    ax2.set_xlabel('Body Mass (kg)')
    ax2.set_ylabel('Maximum Lifespan (years)')
    ax2.set_title('Body Mass vs Maximum Lifespan')
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.legend()
    ax2.grid(alpha=0.3)
    
    # 3. LQ by species (sorted)
    ax3 = axes[1, 0]
    df_sorted = df.sort_values('LQ')
    colors_list = [colors[c] for c in df_sorted['longevity_class']]
    ax3.barh(range(len(df_sorted)), df_sorted['LQ'], color=colors_list, alpha=0.7)
    ax3.axvline(1.0, color='black', linestyle='-', linewidth=1, label='LQ = 1')
    ax3.set_ylabel('Species (sorted by LQ)')
    ax3.set_xlabel('Longevity Quotient (LQ)')
    ax3.set_title('Longevity Quotient by Species')
    ax3.set_yticks([])
    ax3.legend()
    ax3.grid(alpha=0.3, axis='x')
    
    # 4. Box plot by category
    ax4 = axes[1, 1]
    data_to_plot = [df[df['longevity_class'] == cat]['LQ'].values 
                    for cat in ['short-lived', 'intermediate', 'long-lived']]
    bp = ax4.boxplot(data_to_plot, labels=['Short-lived', 'Intermediate', 'Long-lived'],
                     patch_artist=True)
    for patch, color in zip(bp['boxes'], ['red', 'gray', 'green']):
        patch.set_facecolor(color)
        patch.set_alpha(0.6)
    ax4.set_ylabel('Longevity Quotient (LQ)')
    ax4.set_title('LQ Distribution by Category')
    ax4.grid(alpha=0.3, axis='y')
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"  Plot saved to: {output_file}")
    plt.close()

def main():
    parser = argparse.ArgumentParser(
        description='Calculate Longevity Quotient and classify species'
    )
    parser.add_argument('-i', '--input', required=True,
                       help='Input CSV file with phenotype data')
    parser.add_argument('-o', '--output', default='data/phenotypes/lq_data.csv',
                       help='Output CSV file with LQ values')
    parser.add_argument('--top-decile', type=float, default=0.1,
                       help='Top decile fraction (default: 0.1)')
    parser.add_argument('--bottom-decile', type=float, default=0.1,
                       help='Bottom decile fraction (default: 0.1)')
    parser.add_argument('--plot', action='store_true',
                       help='Generate plots')
    
    args = parser.parse_args()
    
    print("="*60)
    print("LONGEVITY QUOTIENT CALCULATION")
    print("="*60)
    
    # Load data
    print(f"\nLoading phenotype data from: {args.input}")
    df = load_phenotype_data(args.input)
    print(f"Loaded {len(df)} species")
    
    # Calculate LQ
    print("\nCalculating Longevity Quotient...")
    df['expected_lifespan_years'] = df['body_mass_kg'].apply(calculate_expected_lifespan)
    df['LQ'] = calculate_longevity_quotient(df['max_lifespan_years'], df['body_mass_kg'])
    
    # Classify species
    print(f"\nClassifying species (top {args.top_decile*100}%, bottom {args.bottom_decile*100}%)...")
    df, top_thresh, bottom_thresh = classify_species(df, args.top_decile, args.bottom_decile)
    
    # Generate statistics
    stats = generate_summary_statistics(df)
    
    # Print summary
    print("\n" + "="*60)
    print("SUMMARY STATISTICS")
    print("="*60)
    print(f"Total species: {stats['total_species']}")
    print(f"\nLifespan (years):")
    print(f"  Mean: {stats['mean_lifespan']:.2f}")
    print(f"  Median: {stats['median_lifespan']:.2f}")
    print(f"  Range: {stats['min_lifespan']:.2f} - {stats['max_lifespan']:.2f}")
    print(f"\nBody mass (kg):")
    print(f"  Mean: {stats['mean_body_mass']:.2f}")
    print(f"  Median: {stats['median_body_mass']:.2f}")
    print(f"\nLongevity Quotient:")
    print(f"  Mean: {stats['mean_lq']:.3f}")
    print(f"  Median: {stats['median_lq']:.3f}")
    print(f"\nClassification:")
    print(f"  Long-lived: {stats['n_long_lived']} species (LQ ≥ {top_thresh:.3f})")
    print(f"  Short-lived: {stats['n_short_lived']} species (LQ ≤ {bottom_thresh:.3f})")
    print(f"  Intermediate: {stats['n_intermediate']} species")
    print("="*60)
    
    # Show top and bottom species
    print("\nTOP 5 LONG-LIVED SPECIES:")
    top5 = df.nlargest(5, 'LQ')[['species_name', 'max_lifespan_years', 'body_mass_kg', 'LQ']]
    print(top5.to_string(index=False))
    
    print("\nTOP 5 SHORT-LIVED SPECIES:")
    bottom5 = df.nsmallest(5, 'LQ')[['species_name', 'max_lifespan_years', 'body_mass_kg', 'LQ']]
    print(bottom5.to_string(index=False))
    
    # Save results
    output_path = Path(args.output)
    output_path.parent.mkdir(exist_ok=True, parents=True)
    df.to_csv(output_path, index=False)
    print(f"\n✓ LQ data saved to: {output_path}")
    
    # Save species lists for each category
    for category in ['long-lived', 'short-lived', 'intermediate']:
        category_file = output_path.parent / f"{category}_species.txt"
        species_list = df[df['longevity_class'] == category]['species_id'].tolist()
        with open(category_file, 'w') as f:
            f.write('\n'.join(species_list))
        print(f"  {category.capitalize()} species list: {category_file}")
    
    # Save summary statistics
    stats_file = output_path.parent / "summary_statistics.json"
    with open(stats_file, 'w') as f:
        json.dump(stats, f, indent=2)
    print(f"  Summary statistics: {stats_file}")
    
    # Generate plots
    if args.plot:
        plot_file = output_path.parent / "lq_distribution.png"
        print(f"\nGenerating plots...")
        plot_lq_distribution(df, plot_file)
    
    print("\n✓ Analysis complete!")

if __name__ == "__main__":
    main()
