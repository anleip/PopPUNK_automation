#!/usr/bin/env python3
"""
Generic script to extract and randomly subsample genomes by species from the 2kk dataset.

Usage:
python extract_species_subset.py "Escherichia coli" 5000
python extract_species_subset.py "Salmonella enterica" 2000
"""

import pandas as pd
import glob
import random
import argparse
import os
import re
from collections import defaultdict, Counter

def sanitize_filename(species_name):
    """Convert species name to a safe filename"""
    # Replace spaces with underscores and remove special characters
    safe_name = re.sub(r'[^\w\s-]', '', species_name)
    safe_name = re.sub(r'[-\s]+', '_', safe_name)
    return safe_name.lower()

def load_data():
    """Load cluster and annotation data"""
    print("Loading cluster assignments...")
    clusters_df = pd.read_csv('/hps/nobackup/rdf/metagenomics/research-team/johanna/projects/PhD_main/analysis/2025-06-18-2kk+space/hq/2kk_hq_clusters.csv')
    print(f"Loaded {len(clusters_df)} genome-cluster assignments")

    print("\nLoading species annotations...")
    annotation_files = glob.glob('/hps/nobackup/rdf/metagenomics/research-team/johanna/projects/PhD_main/analysis/2025-06-18-2kk+space/hq/2kk_hq_annotations*.csv')
    all_annotations = []
    for file in annotation_files:
        print(f"Reading {os.path.basename(file)}...")
        df = pd.read_csv(file)
        all_annotations.append(df)

    annotations_df = pd.concat(all_annotations, ignore_index=True)
    print(f"Loaded {len(annotations_df)} genome annotations")

    print("\nMerging cluster and species data...")
    merged_df = pd.merge(clusters_df, annotations_df, left_on='sample_id', right_on='ID', how='left')
    print(f"Merged dataset has {len(merged_df)} entries")
    
    return merged_df

def find_species_clusters(merged_df, target_species, min_fraction=0.5):
    """Find clusters that are predominantly the target species"""
    print(f"\nAnalyzing clusters for {target_species}...")
    
    cluster_species = defaultdict(list)
    for _, row in merged_df.iterrows():
        cluster_id = row['cluster']
        species = row['GTDB_species']
        if pd.notna(species):
            cluster_species[cluster_id].append(species)

    target_clusters = []
    cluster_stats = {}

    for cluster_id, species_list in cluster_species.items():
        species_counts = Counter(species_list)
        total_with_species = len(species_list)
        
        if total_with_species == 0:
            continue
            
        target_count = species_counts.get(target_species, 0)
        target_fraction = target_count / total_with_species
        
        cluster_stats[cluster_id] = {
            'total_genomes_in_cluster': len(merged_df[merged_df['cluster'] == cluster_id]),
            'genomes_with_species': total_with_species,
            'target_count': target_count,
            'target_fraction': target_fraction,
            'most_common_species': species_counts.most_common(3)
        }
        
        # Consider a cluster as target species if >50% of annotated genomes are the target
        if target_fraction > min_fraction:
            target_clusters.append(cluster_id)

    print(f"Found {len(target_clusters)} {target_species} predominant clusters")
    
    # Show some statistics
    if target_clusters:
        print(f"\n{target_species} cluster statistics (showing first 10):")
        for cluster_id in sorted(target_clusters, key=lambda x: cluster_stats[x]['total_genomes_in_cluster'], reverse=True)[:10]:
            stats = cluster_stats[cluster_id]
            print(f"Cluster {cluster_id}: {stats['target_count']}/{stats['genomes_with_species']} "
                  f"({stats['target_fraction']:.1%}) {target_species}, "
                  f"total genomes: {stats['total_genomes_in_cluster']}")
    
    return target_clusters, cluster_stats

def extract_and_subsample(merged_df, target_clusters, target_species, sample_size, seed=42):
    """Extract all genomes from target clusters and randomly subsample"""
    print(f"\nExtracting all genome IDs from {target_species} clusters...")
    target_genomes = merged_df[merged_df['cluster'].isin(target_clusters)]['sample_id'].unique()
    print(f"Total genomes available: {len(target_genomes)}")
    
    # Random sampling
    random.seed(seed)
    if len(target_genomes) >= sample_size:
        selected_genomes = random.sample(target_genomes.tolist(), sample_size)
        print(f"Randomly selected {len(selected_genomes)} genomes")
    else:
        print(f"Warning: Only {len(target_genomes)} genomes available for {target_species}")
        print(f"Using all {len(target_genomes)} genomes instead of requested {sample_size}")
        selected_genomes = target_genomes.tolist()
    
    return selected_genomes

def create_output_files(selected_genomes, target_species, output_dir):
    """Create the output files: genome IDs and rfile"""
    # Load the full rfile
    print("\nLoading full rfile to get paths...")
    rfile_path = '/hps/nobackup/rdf/metagenomics/research-team/johanna/projects/PhD_main/analysis/2025-06-18-2kk+space/2kk.rfile'
    rfile_df = pd.read_csv(rfile_path, sep='\t', header=None, names=['genome_id', 'path'])
    print(f"Loaded {len(rfile_df)} genome paths from rfile")

    # Filter for selected genomes
    subset_rfile = rfile_df[rfile_df['genome_id'].isin(selected_genomes)]
    print(f"Found paths for {len(subset_rfile)} out of {len(selected_genomes)} selected genomes")

    if len(subset_rfile) < len(selected_genomes):
        missing = set(selected_genomes) - set(subset_rfile['genome_id'])
        print(f"Warning: {len(missing)} genomes not found in rfile")
        if len(missing) <= 5:
            print(f"Missing genomes: {list(missing)}")
        else:
            print(f"First 5 missing: {list(missing)[:5]}")

    # Create safe filename
    safe_species_name = sanitize_filename(target_species)
    
    # Save files
    print(f"\nSaving files to {output_dir}...")
    
    # Save genome IDs (without header for compatibility)
    ids_filename = f"{safe_species_name}.id"
    with open(os.path.join(output_dir, ids_filename), 'w') as f:
        for genome_id in subset_rfile['genome_id']:
            f.write(f"{genome_id}\n")
    
    # Save rfile
    rfile_filename = f"{safe_species_name}.rfile"
    subset_rfile.to_csv(os.path.join(output_dir, rfile_filename), sep='\t', header=False, index=False)
    
    print(f"Files saved:")
    print(f"- {ids_filename}: {len(subset_rfile)} genome IDs")
    print(f"- {rfile_filename}: {len(subset_rfile)} genome ID-path pairs")
    
    return ids_filename, rfile_filename

def main():
    parser = argparse.ArgumentParser(description='Extract and randomly subsample genomes by species from 2kk dataset')
    parser.add_argument('species', type=str, help='Species name (e.g., "Escherichia coli")')
    parser.add_argument('sample_size', type=int, help='Number of genomes to randomly sample')
    parser.add_argument('--output_dir', '-o', type=str, 
                       default='/hps/nobackup/rdf/metagenomics/research-team/johanna/projects/PhD_main/analysis/2025-09-01-Anlei-data',
                       help='Output directory (default: current analysis directory)')
    parser.add_argument('--min_fraction', '-f', type=float, default=0.5,
                       help='Minimum fraction of target species in cluster to include cluster (default: 0.5)')
    parser.add_argument('--seed', '-s', type=int, default=42,
                       help='Random seed for reproducibility (default: 42)')
    
    args = parser.parse_args()
    
    print(f"Extracting and subsampling {args.species}")
    print(f"Target sample size: {args.sample_size}")
    print(f"Minimum species fraction in cluster: {args.min_fraction}")
    print(f"Random seed: {args.seed}")
    print(f"Output directory: {args.output_dir}")
    
    # Create output directory if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Load data
    merged_df = load_data()
    
    # Check if species exists in dataset
    unique_species = merged_df['GTDB_species'].dropna().unique()
    if args.species not in unique_species:
        print(f"\nError: Species '{args.species}' not found in dataset!")
        print(f"Available species with '{args.species.split()[0]}' genus:")
        genus = args.species.split()[0]
        matching_species = [s for s in unique_species if s.startswith(genus)]
        for species in sorted(matching_species)[:10]:  # Show first 10 matches
            print(f"  - {species}")
        if len(matching_species) > 10:
            print(f"  ... and {len(matching_species) - 10} more")
        return
    
    # Find clusters
    target_clusters, cluster_stats = find_species_clusters(merged_df, args.species, args.min_fraction)
    
    if not target_clusters:
        print(f"\nNo clusters found with >{args.min_fraction*100}% {args.species}")
        return
    
    # Extract and subsample
    selected_genomes = extract_and_subsample(merged_df, target_clusters, args.species, args.sample_size, args.seed)
    
    # Create output files
    create_output_files(selected_genomes, args.species, args.output_dir)
    
    print(f"\nDone! Use these files for PopPUNK analysis of {args.species}")

if __name__ == "__main__":
    main()