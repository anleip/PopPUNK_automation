#!/usr/bin/env python3

import argparse
import itertools
import math
import random
import csv
from collections import defaultdict
from pathlib import Path
import numpy as np
from ete3 import Tree
import matplotlib.pyplot as plt


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compute per-cluster and tree-level polyphyly from a Newick tree and cluster CSV."
    )

    parser.add_argument(
        "--tree",
        required=True,
        help="Path to Newick tree file"
    )

    parser.add_argument(
        "--clusters_root",
        required=True,
        help="Path to folder that contains fitted thresholds with their own clusters CSV file"
    )

    parser.add_argument(
        "--max-pairs",
        type=int,
        default=1000,
        help="Maximum number of leaf pairs sampled per cluster (default: 1000)"
    )

    parser.add_argument(
        "--plot",
        default=None,
        help="Plot polyphyly distribution for each threhosld."
    )

    return parser.parse_args()


def sample_pairs(items, max_pairs=1500):
    n_pairs = len(items) * (len(items) - 1) // 2

    if n_pairs <= max_pairs:
        return list(itertools.combinations(items, 2))

    pairs = set()
    while len(pairs) < max_pairs:
        i, j = random.sample(range(len(items)), 2)
        if i > j:
            i, j = j, i
        pairs.add((items[i], items[j]))

    return list(pairs)


def main():
    args = parse_args()

    # Load tree
    tree = Tree(args.tree, format=1)

    # Load clusters
    root = Path(args.clusters_root)
    
    # Output record
    tree_level_polyphyly = []

    for subdir in sorted(root.iterdir()):
        if not subdir.is_dir():
            continue
        csv_file = subdir / f"{subdir.name}_clusters.csv"
        print(csv_file)
        clusters = defaultdict(list)
        cluster_map = {}
        with open(csv_file, newline="") as f:
            reader = csv.DictReader(f)
            for row in reader:
                clusters[row["Cluster"]].append(row["Taxon"])
                cluster_map[row["Taxon"]] = row["Cluster"]
        clusters = dict(clusters)
        K = len(clusters)

        print(f"Number of clusters: {K}")

        if K != 1:
            # Cache MRCA cluster counts
            mrca_cluster_count = {}                
            # Compute per-cluster polyphyly
            P_C = {}

            def polyphyly(tree,leaf1, leaf2):
                mrca = tree.get_common_ancestor(leaf1, leaf2)
                key = id(mrca)
                if key not in mrca_cluster_count:
                    descendant_clusters = {
                        cluster_map[leaf.name]
                        for leaf in mrca.get_leaves()
                    }
                    mrca_cluster_count[key] = len(descendant_clusters)

                k = mrca_cluster_count[key]
                return (k - 1) / (K - 1)

            for cid, leaves in clusters.items():
                if len(leaves) < 2:
                    continue

                pairs = sample_pairs(leaves, max_pairs=args.max_pairs)
                values = [polyphyly(tree,a, b) for a, b in pairs]
                P_C[cid] = sum(values) / len(values)

            # Tree-level statistics
            tree_mean_unweighted = sum(P_C.values()) / len(P_C)

            num = 0.0
            den = 0.0
            for cid, pc in P_C.items():
                w = math.log(len(clusters[cid]))
                num += w * pc
                den += w

            tree_mean_log_weighted = num / den
        
        else:
            tree_mean_unweighted = 1
            tree_mean_log_weighted = 1

        print("Unweighted mean tree-level polyphyly:\t", tree_mean_unweighted)
        print("Log weighted mean tree-level polyphyly:\t", tree_mean_log_weighted)

        tree_level_polyphyly.append([tree_mean_unweighted,tree_mean_log_weighted])
        

        # Plot
        if args.plot:
            plt.figure(figsize=(7, 4))
            plt.hist(P_C.values(), bins=30)
            plt.xlabel("Per-cluster mean polyphyly (P_C)")
            plt.ylabel("Number of clusters")
            plt.title("Distribution of cluster polyphyly")
            plt.tight_layout()
        
            plt.savefig(args.plot, dpi=300)
            print(f"Plot saved to {args.plot}")
        

    out_file = Path(f"{args.clusters_root}_polyphyly.tsv")
    np.savetxt(out_file, tree_level_polyphyly, delimiter="\t")
    print("Polyphyly scores saved.")

if __name__ == "__main__":
    main()
