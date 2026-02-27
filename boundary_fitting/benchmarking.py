# %%
import numpy as np
import pypickle as ppk 
import pandas as pd
from sklearn.metrics.cluster import adjusted_rand_score 
# %%
# converting pkl file of 20k qc sample ids into tsv file
readpkl = ppk.load("/nfs/research/jlees/anlei/analysis/n.gonorrhoeae/neisseria_gonorrhoeae_20000_qc/neisseria_gonorrhoeae_20000_qc.dists.pkl")
# print(readpkl)
# %%
print(len(readpkl[0]))
qcid = np.array(readpkl[0])
# %%
np.savetxt("/nfs/research/jlees/anlei/analysis/n.gonorrhoeae/n.gonorrhoeae_20k_qc.tsv", qcid,delimiter="\t",fmt='%s')
# %%
# Remove non-relavent lines from download (file no longer exist)
#isolate_df = pd.read_csv("/nfs/research/jlees/anlei/bin/boundary_fitting/benchmarking_data/e.coli_isolates.tsv",delimiter="\t")
# %%
print(isolate_df.head)
isolate_cleaned = pd.concat([isolate_df["BioSample"],isolate_df["SNP cluster"]],axis=1)
print(isolate_cleaned.shape)

# %%
# Export to tsv
isolate_cleaned.to_csv("/nfs/research/jlees/anlei/analysis/e.coli/e.coli_isolates_SNP_all.tsv")
# %%
print(isolate_cleaned.head)
# %%
isolate_20k = pd.read_csv("/nfs/research/jlees/anlei/analysis/e.coli/e.coli_isolates_SNP_20k.tsv",delimiter="\t")
print(isolate_20k.shape)
# %%
#
#
# Join the columns
ec_ppc00175 = pd.read_csv("/nfs/research/jlees/anlei/analysis/e.coli/e.coli_20k_qc_fit/e.coli_20k_qc_fit0.00175/e.coli_20k_qc_fit0.00175_clusters.csv")
ec_ppc005 = pd.read_csv("/nfs/research/jlees/anlei/analysis/e.coli/e.coli_20k_qc_fit/e.coli_20k_qc_fit0.005000000000000001/e.coli_20k_qc_fit0.005000000000000001_clusters.csv")
# %%
ec_stcc = pd.read_csv("/nfs/research/jlees/anlei/analysis/e.coli/e.coli_MLST/20k_stcc_qc.csv",na_values=["-","new_alleles","new_ST"])
print(ec_stcc.head)
print(len(ec_stcc))
# %%
ec_ppc005.rename(columns={"Taxon":"Genome","Cluster":"ppc_0.005"},inplace=True)
ec_ppc00175.rename(columns={"Taxon":"Genome","Cluster":"ppc_0.00175"},inplace=True)
ec_ppc = pd.merge(ec_ppc005,ec_ppc00175,"left")
print(ec_ppc.head)
# %%
ec_merged = pd.merge(ec_ppc,ec_stcc,"left")
print(ec_merged)
# %%
ec_merged.to_csv("/nfs/research/jlees/anlei/analysis/e.coli/e.coli_20k_all_clusters.csv")
# %%
# Adjusted Rand index
#
#
a = [0,0,1,1]; b = [1,1,0,0]
print(adjusted_rand_score(a,b))
# %%
# E. coli 0.00175 as an example. ST as ground truth
print("Not typed genomes",ec_merged["ST"].isna().sum())
print("No clonal complex genomes",ec_merged["clonal_complex"].isna().sum())
# %%
# Clonal complex as ground truth
ec_noNA = ec_merged.dropna()
print(len(ec_noNA))
# %%
pp005_clusters = ec_noNA["ppc_0.005"].values.tolist()
pp00175_clusters = ec_noNA["ppc_0.00175"].values.tolist()
st_clusters = ec_noNA["ST"].values.tolist()
cc_clusters = ec_noNA["clonal_complex"].values.tolist()
print("PopPUNK 0.005 against ST:\t",adjusted_rand_score(st_clusters,pp005_clusters))
print("PopPUNK 0.005 against clonal complex:\t",adjusted_rand_score(cc_clusters,pp005_clusters))
print("PopPUNK 0.00175 against ST:\t",adjusted_rand_score(st_clusters,pp00175_clusters))
print("PopPUNK 0.00175 against clonal complex:\t",adjusted_rand_score(cc_clusters,pp00175_clusters))
# %%
#
#
# S. pneumoniae 20k
sp_ppc001 = pd.read_csv("/nfs/research/jlees/anlei/analysis/s.pneumoniae/s.pneumoniae_20k_qc_fit/s.pneumoniae_20k_qc_fit0.001/s.pneumoniae_20k_qc_fit0.001_clusters.csv")
sp_ppc00325 = pd.read_csv("/nfs/research/jlees/anlei/analysis/s.pneumoniae/s.pneumoniae_20k_qc_fit/s.pneumoniae_20k_qc_fit0.0032500000000000003/s.pneumoniae_20k_qc_fit0.0032500000000000003_clusters.csv")
sp_ppc00525 = pd.read_csv("/nfs/research/jlees/anlei/analysis/s.pneumoniae/s.pneumoniae_20k_qc_fit/s.pneumoniae_20k_qc_fit0.00525/s.pneumoniae_20k_qc_fit0.00525_clusters.csv")
# %%
sp_stcc = pd.read_csv("/nfs/research/jlees/anlei/analysis/s.pneumoniae/s.pneumoniae_MLST/20k_MLST_qc_stcc.csv",na_values=["-","new_alleles","new_ST"])
print(sp_stcc.head)
print(len(sp_stcc))
# %%
sp_ppc001.rename(columns={"Taxon":"Genome","Cluster":"ppc_0.001"},inplace=True)
sp_ppc00325.rename(columns={"Taxon":"Genome","Cluster":"ppc_0.00325"},inplace=True)
sp_ppc00525.rename(columns={"Taxon":"Genome","Cluster":"ppc_0.00525"},inplace=True)
sp_ppc0 = pd.merge(sp_ppc00525,sp_ppc00325,"left")
sp_ppc = pd.merge(sp_ppc0,sp_ppc001,"left")
print(sp_ppc.head)
# %%
sp_merged = pd.merge(sp_ppc,sp_stcc,"left")
print(sp_merged)
# %%
sp_merged.to_csv("/nfs/research/jlees/anlei/analysis/s.pneumoniae/s.pneumoniae_20k_all_clusters.csv")
# %%
print("Not typed genomes",sp_merged["ST"].isna().sum())
sp_noNA = sp_merged[sp_merged["ST"].notna()]
print(len(sp_noNA))
# %%
pp001_clusters = sp_noNA["ppc_0.001"].values.tolist()
pp00325_clusters = sp_noNA["ppc_0.00325"].values.tolist()
pp00525_clusters = sp_noNA["ppc_0.00525"].values.tolist()
st_clusters = sp_noNA["ST"].values.tolist()
print("PopPUNK 0.001 against ST:\t",adjusted_rand_score(st_clusters,pp001_clusters))
print("PopPUNK 0.00325 against ST:\t",adjusted_rand_score(st_clusters,pp00325_clusters))
print("PopPUNK 0.00525 against ST:\t",adjusted_rand_score(st_clusters,pp00525_clusters))
# %%
#
#
# S. aureus
sa_ppc0035 = pd.read_csv("/nfs/research/jlees/anlei/analysis/s.aureus/s.aureus_20k_qc_fit/s.aureus_20k_qc_fit0.0035000/s.aureus_20k_qc_fit0.0035000_clusters.csv")
sa_ppc0075 = pd.read_csv("/nfs/research/jlees/anlei/analysis/s.aureus/s.aureus_20k_qc_fit/s.aureus_20k_qc_fit0.0075000/s.aureus_20k_qc_fit0.0075000_clusters.csv")
sa_st = pd.read_csv("/nfs/research/jlees/anlei/analysis/s.aureus/s.aureus_MLST/20k_MLST_ST.csv",na_values=["-","new_alleles","new_ST"])
# %%
sa_ppc0035.rename(columns={"Taxon":"Genome","Cluster":"ppc_0.0035"},inplace=True)
sa_ppc0075.rename(columns={"Taxon":"Genome","Cluster":"ppc_0.0075"},inplace=True)
sa_ppc = pd.merge(sa_ppc0075,sa_ppc0035,"left")
sa_merged = pd.merge(sa_ppc,sa_st,"left")
print(sa_merged)
sa_merged.to_csv("/nfs/research/jlees/anlei/analysis/s.aureus/s.aureus_20k_all_clusters.csv")

# %%
print("Not typed genomes",sa_merged["ST"].isna().sum())
sa_noNA = sa_merged[sa_merged["ST"].notna()]
print(len(sa_noNA))
# %%
pp0035_clusters = sa_noNA["ppc_0.0035"].values.tolist()
pp0075_clusters = sa_noNA["ppc_0.0075"].values.tolist()
st_clusters = sa_noNA["ST"].values.tolist()
print("PopPUNK 0.0035 against ST:\t",adjusted_rand_score(st_clusters,pp0035_clusters))
print("PopPUNK 0.0075 against ST:\t",adjusted_rand_score(st_clusters,pp0075_clusters))

# %%
#
#
# S. enterica
se_ppc00175 = pd.read_csv("/nfs/research/jlees/anlei/analysis/s.enterica/s.enterica_20k_qc_fit/s.enterica_20k_qc_fit0.00175/s.enterica_20k_qc_fit0.00175_clusters.csv")
se_ppc00225 = pd.read_csv("/nfs/research/jlees/anlei/analysis/s.enterica/s.enterica_20k_qc_fit/s.enterica_20k_qc_fit0.0022500000000000003/s.enterica_20k_qc_fit0.0022500000000000003_clusters.csv")
se_ppc00525 = pd.read_csv("/nfs/research/jlees/anlei/analysis/s.enterica/s.enterica_20k_qc_fit/s.enterica_20k_qc_fit0.0052499999999999995/s.enterica_20k_qc_fit0.0052499999999999995_clusters.csv")
se_stcc = pd.read_csv("/nfs/research/jlees/anlei/analysis/s.enterica/s.enterica_MLST/20k_MLST_stcc.csv",na_values=["-","new_alleles","new_ST"])
# %%
se_ppc00175.rename(columns={"Taxon":"Genome","Cluster":"ppc_0.00175"},inplace=True)
se_ppc00225.rename(columns={"Taxon":"Genome","Cluster":"ppc_0.00225"},inplace=True)
se_ppc00525.rename(columns={"Taxon":"Genome","Cluster":"ppc_0.00525"},inplace=True)
merge1 = pd.merge(se_ppc00525,se_ppc00225,"left")
se_ppc = pd.merge(merge1,se_ppc00175,"left")
se_merged = pd.merge(se_ppc,se_stcc,"left")
print(se_merged)
se_merged.to_csv("/nfs/research/jlees/anlei/analysis/s.enterica/s.enterica_20k_all_clusters.csv")
# %%
print("Not typed genomes",se_merged["ST"].isna().sum())
print("No clonal complex genomes",se_merged["clonal_complex"].isna().sum())
se_noNA = se_merged.dropna()
print(len(se_noNA))
# %%
pp00175_clusters = se_noNA["ppc_0.00175"].values.tolist()
pp00225_clusters = se_noNA["ppc_0.00225"].values.tolist()
pp00525_clusters = se_noNA["ppc_0.00525"].values.tolist()
st_clusters = se_noNA["ST"].values.tolist()
cc_clusters = se_noNA["clonal_complex"].values.tolist()
print("PopPUNK 0.00175 against ST:\t",adjusted_rand_score(st_clusters,pp00175_clusters))
print("PopPUNK 0.00175 against clonal complex:\t",adjusted_rand_score(cc_clusters,pp00175_clusters))
print("PopPUNK 0.00225 against ST:\t",adjusted_rand_score(st_clusters,pp00225_clusters))
print("PopPUNK 0.00225 against clonal complex:\t",adjusted_rand_score(cc_clusters,pp00225_clusters))
print("PopPUNK 0.00525 against ST:\t",adjusted_rand_score(st_clusters,pp00525_clusters))
print("PopPUNK 0.00525 against clonal complex:\t",adjusted_rand_score(cc_clusters,pp00525_clusters))

# %%
#
#
# M. tuberculosis
mt_ppc00035 = pd.read_csv("/nfs/research/jlees/anlei/analysis/m.tuberculosis/mycobacterium_tuberculosis_20000_qc_fit_20260107_092231/mycobacterium_tuberculosis_20000_qc_fit0.0003500/mycobacterium_tuberculosis_20000_qc_fit0.0003500_clusters.csv")
mt_st = pd.read_csv("/nfs/research/jlees/anlei/analysis/m.tuberculosis/m.tuberculosis_MLST/20k_MLST_ST.csv",na_values=["-","new_alleles","new_ST"])
# %%
mt_ppc00035.rename(columns={"Taxon":"Genome","Cluster":"ppc_0.00035"},inplace=True)
mt_merged = pd.merge(mt_ppc00035,mt_st,"left")
print(mt_merged)
mt_merged.to_csv("/nfs/research/jlees/anlei/analysis/m.tuberculosis/m.tuberculosis_20k_all_clusters.csv")
# %%
print("Not typed genomes",mt_merged["ST"].isna().sum())
mt_noNA = mt_merged.dropna()
pp00035_clusters = mt_noNA["ppc_0.00035"].values.tolist()
st_clusters = mt_noNA["ST"].values.tolist()
print("PopPUNK 0.00035 against ST:\t",adjusted_rand_score(st_clusters,pp00035_clusters))

# %%
#
#
# N. gonorrhoeae
ng_ppc001625 = pd.read_csv("/nfs/research/jlees/anlei/analysis/n.gonorrhoeae/neisseria_gonorrhoeae_20000_qc_fit/neisseria_gonorrhoeae_20000_qc_fit0.001625/neisseria_gonorrhoeae_20000_qc_fit0.001625_clusters.csv")
ng_ppc0023 = pd.read_csv("/nfs/research/jlees/anlei/analysis/n.gonorrhoeae/neisseria_gonorrhoeae_20000_qc_fit/neisseria_gonorrhoeae_20000_qc_fit0.0023000/neisseria_gonorrhoeae_20000_qc_fit0.0023000_clusters.csv")
ng_st = pd.read_csv("/nfs/research/jlees/anlei/analysis/n.gonorrhoeae/n.gonorrhoeae_MLST/20k_MLST_ST.csv",na_values=["-","new_alleles","new_ST"])
# %%
ng_ppc001625.rename(columns={"Taxon":"Genome","Cluster":"ppc_0.001625"},inplace=True)
ng_ppc0023.rename(columns={"Taxon":"Genome","Cluster":"ppc_0.0023"},inplace=True)
ng_ppc = pd.merge(ng_ppc0023,ng_ppc001625,"left")
ng_merged = pd.merge(ng_ppc,ng_st,"left")
print(ng_merged)
ng_merged.to_csv("/nfs/research/jlees/anlei/analysis/n.gonorrhoeae/n.gonorrhoeae_20k_all_clusters.csv")
# %%
print("Not typed genomes",ng_merged["ST"].isna().sum())
ng_noNA = ng_merged.dropna()
# %%
pp001625_clusters = ng_noNA["ppc_0.001625"].values.tolist()
pp0023_clusters = ng_noNA["ppc_0.0023"].values.tolist()
st_clusters = ng_noNA["ST"].values.tolist()
print("PopPUNK 0.001625 against ST:\t",adjusted_rand_score(st_clusters,pp001625_clusters))
print("PopPUNK 0.0023 against ST:\t",adjusted_rand_score(st_clusters,pp0023_clusters))

# %%
