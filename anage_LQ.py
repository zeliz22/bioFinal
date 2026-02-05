#!/usr/bin/env python3
"""
Script 3: AnAge to LQ Pipeline

1. Reads raw anage_data.txt (Tab-delimited)
2. Cleans and converts units (Grams to KG)
3. Calculates Longevity Quotient (LQ)
4. Saves top 4 and bottom 4 species to fetch_targets.txt
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import argparse

def calculate_expected_lifespan(body_mass_kg, taxon_class):
    """
    Calculates expected lifespan. 
    Note: The constants 4.88 and 0.19 are standard for Mammals.
    """
    return 4.88 * (body_mass_kg ** 0.19)

def main():
    parser = argparse.ArgumentParser(description='Convert AnAge TXT to LQ Targets')
    parser.add_argument('-i', '--input', default='data/anage_data.txt', help='Path to anage_data.txt')
    parser.add_argument('-c', '--class_filter', default='Mammalia', help='Taxonomic class (e.g., Mammalia, Aves)')
    parser.add_argument('-o', '--output_dir', default='data/LQ', help='Output directory')
    args = parser.parse_args()

    # 1. READ AND CLEAN
    # AnAge is tab-separated. We use low_memory=False to handle mixed types.
    df_raw = pd.read_csv(args.input, sep='\t', low_memory=False)

    # Filter by Class if specified
    if args.class_filter:
        df = df_raw[df_raw['Class'] == args.class_filter].copy()
    else:
        df = df_raw.copy()

    # Create IDs and handle units
    # We combine Genus and Species to match NCBI format (e.g., Homo_sapiens)
    df['species_id'] = df['Genus'] + '_' + df['Species']
    df['species_name'] = df['Common name']
    
    # Map AnAge columns to our variables
    # We divide Body mass (g) by 1000 to get KG for the formula
    df['max_lifespan_years'] = pd.to_numeric(df['Maximum longevity (yrs)'], errors='coerce')
    df['body_mass_kg'] = pd.to_numeric(df['Body mass (g)'], errors='coerce') / 1000

    # Remove rows with missing critical data
    df = df.dropna(subset=['max_lifespan_years', 'body_mass_kg'])

    # 2. CALCULATE LQ
    df['expected_mls'] = df.apply(lambda x: calculate_expected_lifespan(x['body_mass_kg'], x['Class']), axis=1)
    df['LQ'] = df['max_lifespan_years'] / df['expected_mls']

    # 3. CLASSIFY
    top_t = df['LQ'].quantile(0.9)
    bot_t = df['LQ'].quantile(0.1)

    def classify(lq):
        if lq >= top_t: return 'long-lived'
        elif lq <= bot_t: return 'short-lived'
        return 'intermediate'
    df['longevity_class'] = df['LQ'].apply(classify)

    # 4. SAVE FILES
    out_dir = Path(args.output_dir)
    out_dir.mkdir(exist_ok=True, parents=True)

    # 1. Get the lists
    top4 = df.nlargest(4, 'LQ')['species_id'].tolist()
    bot4 = df.nsmallest(4, 'LQ')['species_id'].tolist()
    
    # 2. Save the combined list for the downloader
    with open(out_dir / "fetch_targets.txt", 'w') as f:
        f.write('\n'.join(top4 + bot4))

    # 3. Save the separate lists for the comparison script
    with open(out_dir / "long_lived_targets.txt", 'w') as f:
        f.write('\n'.join(top4))
        
    with open(out_dir / "short_lived_targets.txt", 'w') as f:
        f.write('\n'.join(bot4))

    # # Save categorical lists
    # for cat in ['long-lived', 'short-lived']:
    #     ids = df[df['longevity_class'] == cat]['species_id'].tolist()
    #     with open(out_dir / f"{cat}_species.txt", 'w') as f:
    #         f.write('\n'.join(ids))

    # # 5. PLOT
    # plt.figure(figsize=(10, 6))
    # plt.scatter(df['body_mass_kg'], df['max_lifespan_years'], c='gray', alpha=0.3)
    # plt.scatter(df.nlargest(4, 'LQ')['body_mass_kg'], df.nlargest(4, 'LQ')['max_lifespan_years'], c='green', label='Top 4 LQ')
    # plt.scatter(df.nsmallest(4, 'LQ')['body_mass_kg'], df.nsmallest(4, 'LQ')['max_lifespan_years'], c='red', label='Bottom 4 LQ')
    # plt.xscale('log'); plt.yscale('log')
    # plt.xlabel('Mass (kg)'); plt.ylabel('Lifespan (years)')
    # plt.title(f'Longevity Outliers in {args.class_filter}')
    # plt.legend()
    # plt.savefig(out_dir / 'lq_chart.png')

    # print(f"Done! Targets saved to {out_dir}/fetch_targets.txt")

if __name__ == "__main__":
    main()