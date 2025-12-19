#%%
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
#%%
dists = np.load("../analysis/e.coli/e.coli_5k_pi0.1_db/e.coli_5k_pi0.1_db.dists.npy")
print(dists)
print(len(dists))
#%%

pi = dists[:,0]
a = dists[:,1]

plt.hist(pi,bins=30)
plt.xlabel("π")
plt.ylabel("Frequency")
plt.show()

plt.hist(a,bins=30)
plt.xlabel("a")
plt.ylabel("Frequency")
plt.show()

#%%
import pypickle as ppk 
# %%
names = ppk.load("../analysis/e.coli/e.coli_5k_pi0.1_db/e.coli_5k_pi0.1_db.dists.pkl")
print(names[0])
print(names[1])
print(names[2])
#print(len(names[2]))
#print(len(names[1]))
# %%

# Zero distance checks
#
# Salmonella houtenae distances
jaccard_dists = np.load("../analysis/s.houtenae/s.houtanae_sketchlib_db/s.houtanae_sketchlib_db.jaccard.npy")
print(jaccard_dists)
# %%
print(len(jaccard_dists))
print(jaccard_dists.shape)
# %%
names = ["Query","Reference","Core","Accessory"]
full_dict = pd.read_csv("../analysis/s.houtenae/s.houtanae_sketchlib_db/s.houtanae_sketchlib_db.dists.out",sep="\t",names=names)
#zero_dict = pd.read_csv("../analysis/s.houtenae/s.houtanae_sketchlib_db/s.houtanae_sketchlib_db.dists.zeros",sep="\t",names=names)
print(zero_dict)
# %%
query_list = zero_dict["Query"].tolist()
ref_list = zero_dict["Reference"].tolist()
print(query_list.count("SAMN01933116"))
# %%
pp_zero_dict = pd.read_csv("../analysis/s.houtenae/salmonella_houtenae_db/salmonella_houtenae_db.dists.zeros",sep="\t",names=names)
print(pp_zero_dict)
# %%
pp_query_list = pp_zero_dict["Query"].tolist()
pp_ref_list = pp_zero_dict["Reference"].tolist()
print(pp_query_list.count("SAMN33020872")+pp_ref_list.count("SAMN33020872"))

# %%
# Plot Jaccard
k_list = [17,21,25,29,33,37,41]
plt.scatter(k_list,jaccard_dists[0])
plt.plot()
# %%
zero_index = []
#for row in full_dict:
# %%
print(full_dict.shape)
# %%
