#!/usr/bin/env python3
"""
Script 3: AnAge to LQ Pipeline (With Validation Support)
Creates lq_results.csv with proper column names for CAAS validation
"""

import pandas as pd
import numpy as np
from pathlib import Path
import argparse

def calculate_expected_lifespan(body_mass_kg, taxon_class):
    return 4.88 * (body_mass_kg ** 0.19)

def main():
    parser = argparse.ArgumentParser(description='Convert AnAge TXT to LQ Targets')
    parser.add_argument('-i', '--input', default='data/anage_data.txt', help='Path to anage_data.txt')
    parser.add_argument('-c', '--class_filter', default='Mammalia', help='Taxonomic class')
    parser.add_argument('-o', '--output_dir', default='data/LQ', help='Output directory')
    args = parser.parse_args()

    df_raw = pd.read_csv(args.input, sep='\t', low_memory=False)

    if args.class_filter:
        df = df_raw[df_raw['Class'] == args.class_filter].copy()
    else:
        df = df_raw.copy()

    df['species_id'] = df['Genus'] + '_' + df['Species']
    df['max_lifespan_years'] = pd.to_numeric(df['Maximum longevity (yrs)'], errors='coerce')
    df['body_mass_kg'] = pd.to_numeric(df['Body mass (g)'], errors='coerce') / 1000

    df = df.dropna(subset=['max_lifespan_years', 'body_mass_kg'])

    df['expected_mls'] = df.apply(lambda x: calculate_expected_lifespan(x['body_mass_kg'], x['Class']), axis=1)
    df['LQ'] = df['max_lifespan_years'] / df['expected_mls']
    
    df['longevity_quotient'] = df['LQ']

    
    df_sorted = df.sort_values('LQ', ascending=False)


    amoundOfSpeciesInEachGroup = 6

    top4 = df_sorted.head(amoundOfSpeciesInEachGroup)['species_id'].tolist()
    
    bot4 = df_sorted.tail(amoundOfSpeciesInEachGroup)['species_id'].tolist()
    
 
    mid_index = len(df_sorted) // 2
    val4 = df_sorted.iloc[mid_index-4*amoundOfSpeciesInEachGroup : mid_index+4*amoundOfSpeciesInEachGroup]['species_id'].tolist()

    out_dir = Path(args.output_dir)
    out_dir.mkdir(exist_ok=True, parents=True)

   
    output_cols = ['species_id', 'longevity_quotient', 'max_lifespan_years', 
                   'body_mass_kg', 'expected_mls', 'LQ', 'Genus', 'Species', 
                   'Common name', 'Class', 'Order', 'Family']
    
    output_cols = [col for col in output_cols if col in df.columns]
    
    df[output_cols].to_csv(out_dir / 'lq_results.csv', index=False)

   
    all_12 = top4 + bot4 + val4
    with open(out_dir / "fetch_targets.txt", 'w') as f:
        f.write('\n'.join(all_12))

    with open(out_dir / "long_lived_targets.txt", 'w') as f:
        f.write('\n'.join(top4))
        
    with open(out_dir / "short_lived_targets.txt", 'w') as f:
        f.write('\n'.join(bot4))

    with open(out_dir / "validation_targets.txt", 'w') as f:
        f.write('\n'.join(val4))

    print(f"Success! Created target files in {out_dir}")
    print(f"Total species to download: {len(all_12)}")
    print(f"Discovery: 4 Long vs 4 Short")
    print(f"Validation: 4 Middle-lived species")
    print(f"\nCreated lq_results.csv with columns:")
    print(f"  - species_id: {df['species_id'].nunique()} species")
    print(f"  - longevity_quotient: range {df['longevity_quotient'].min():.2f} to {df['longevity_quotient'].max():.2f}")
    print(f"\nLong-lived species (High LQ):")
    for sp in top4:
        lq_val = df[df['species_id'] == sp]['longevity_quotient'].values[0]
        print(f"  {sp}: LQ = {lq_val:.2f}")
    print(f"\nShort-lived species (Low LQ):")
    for sp in bot4:
        lq_val = df[df['species_id'] == sp]['longevity_quotient'].values[0]
        print(f"  {sp}: LQ = {lq_val:.2f}")
    print(f"\nValidation species (Middle LQ):")
    for sp in val4:
        lq_val = df[df['species_id'] == sp]['longevity_quotient'].values[0]
        print(f"  {sp}: LQ = {lq_val:.2f}")

if __name__ == "__main__":
    main()