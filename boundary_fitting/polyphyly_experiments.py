# %%
from ete3 import Tree
import itertools
import math
import random
import csv
from collections import defaultdict
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
# %%
# Use Newick tree from grapetree visualisation
tree = Tree("/nfs/research/jlees/anlei/analysis/e.coli/e.coli_20k_qc_fit0.005_viz_gt/e.coli_20k_qc_fit0.005_viz_gt_core_NJ.nwk", format=1)
#tree = Tree("/nfs/research/jlees/anlei/analysis/e.coli/e.coli_5k_num3_x1.1_db_fit0.005_viz/e.coli_5k_num3_x1.1_db_fit0.005_viz_core_NJ.nwk",format=1)
# %%
clusters = defaultdict(list)
cluster_map = {}
# Use cluster csv from fit output
with open("/nfs/research/jlees/anlei/analysis/e.coli/e.coli_20k_qc_fit/e.coli_20k_qc_fit0.005000000000000001/e.coli_20k_qc_fit0.005000000000000001_clusters.csv", newline="") as f:
#with open("/nfs/research/jlees/anlei/analysis/e.coli/e.coli_5k_num3_x1.1_db_fit/e.coli_5k_num3_x1.1_db_fit0.005/e.coli_5k_num3_x1.1_db_fit0.005_clusters.csv",newline="") as f:
    reader = csv.DictReader(f)
    for row in reader:
        clusters[row["Cluster"]].append(row["Taxon"])
        cluster_map[row["Taxon"]] = row["Cluster"]
# convert to regular dict
clusters = dict(clusters)
K = len(clusters)
#print(clusters)
print("number of clusters: ",K)
#print(cluster_map)
# %%
mrca_cluster_count = {}

def polyphyly(tree, leaf1, leaf2):
    mrca = tree.get_common_ancestor(leaf1, leaf2)

    if mrca not in mrca_cluster_count:
        descendant_clusters = {
            cluster_map[leaf.name]
            for leaf in mrca.get_leaves()
        }
        mrca_cluster_count[mrca] = len(descendant_clusters)

    k = mrca_cluster_count[mrca]
    #print(k)
    # normalized polyphyly
    return (k - 1) / (K - 1)  

def sample_pairs(items, max_pairs=1000):
    n_pairs = len(items) * (len(items)-1) //2
    if n_pairs <= max_pairs:
        pairs = list(itertools.combinations(items, 2))
    
    else:
        #print("n pairs > 1000")
        pairs = set()
        while len(pairs) < max_pairs:
            i, j = random.sample(range(len(items)), 2)
            if i > j:
                i, j = j, i
            pairs.add((items[i], items[j]))
            
    return pairs
# %%
P_C = {}  # cluster_id ->  polyphyly

for cid, leaves in clusters.items():
    #print("#\n#\n#\n")
    #print(cid)
    
    if len(leaves) < 2:
        continue  # skip singletons

    pairs = sample_pairs(leaves, max_pairs=1000)
    #print(pairs)
    
    values = [
        polyphyly(tree, a, b)
        for a, b in pairs
    ]
    
    P_C[cid] = sum(values) / len(values)
# %%
print(P_C)

# %%
# Tree-level unweighted mean of P_C
tree_mean_unweighted = sum(P_C.values()) / len(P_C)

# Tree-level log-weighted mean of P_C, avoid domination by huge clusters
num = 0.0
den = 0.0

for cid, pc in P_C.items():
    w = math.log(len(clusters[cid]))
    num += w * pc
    den += w

tree_mean_log_weighted = num / den
print("Unweighted mean tree-level polyphyly:\t",tree_mean_unweighted)
print("Log weighted mean tree-level polyphyly:\t",tree_mean_log_weighted)
# %%
plt.figure(figsize=(7, 4))
plt.hist(P_C.values(), bins=30)
plt.xlabel("Per-cluster median polyphyly (P_C)")
plt.ylabel("Number of clusters")
plt.title("Distribution of cluster polyphyly")
plt.tight_layout()
plt.show()

# %%
root = Path("/nfs/research/jlees/anlei/analysis/e.coli/e.coli_5k_num3_x1.1_db_fit")
for subdir in sorted(root.iterdir()):
    print(subdir)
# %%
# Plot polyphyly scores as threshold changes
thresholds = np.linspace(0.00025,0.0085,34)
print(thresholds)
# %%
scores = np.loadtxt("/nfs/research/jlees/anlei/analysis/s.pneumoniae/s.pneumonia1_5k_x0.2_db_fit_polyphyly.tsv")
print(scores)
print(len(scores))
# %%
plt.scatter(thresholds,scores[:,0],label="unweighted")
plt.scatter(thresholds,scores[:,1],label="log-weighted")
plt.title("S. pneumoniae, N = 5000")
plt.xlabel("Core boundary")
plt.ylabel("Polyphyly scores")
plt.legend()
plt.show()

# %%
# Rank
unweighted_rank = scores[:,0].argsort().argsort()
logweighted_rank = scores[:,1].argsort().argsort()
print(unweighted_rank)
print(logweighted_rank)
# %%
network_score = np.loadtxt("/nfs/research/jlees/anlei/analysis/s.pneumoniae/streptococcus_pneumoniae_5000_scores_double.tsv")
print(network_score)
print(len(network_score))
# %%

# %%
def rank_highest_first(x):
    rank_lowest_first = x.argsort().argsort()
    rank_hightest_first = (rank_lowest_first - max(rank_lowest_first))*-1
    return rank_hightest_first
# %%
ns_default_rank = rank_highest_first(network_score[:,0])
ns_btw_rank = rank_highest_first(network_score[:,1])
ns_weibtw_rank = rank_highest_first(network_score[:,2])
print(ns_default)
print(ns_btw)
print(ns_weibtw)
# %%
# Smaller the rank, better the threshold
ns_default_log_weighted = ns_default_rank + logweighted_rank
ns_weibtw_log_weighted = ns_weibtw + logweighted_rank
plt.plot(thresholds,ns_default_log_weighted,label="default NS + log polyphyly")
plt.plot(thresholds,ns_weibtw_log_weighted,label="w/ weighted btw + log polyphyly")
plt.xlabel("Core boundary")
plt.ylabel("Rank")
plt.legend()
plt.show()
# %%
