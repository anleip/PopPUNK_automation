# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# %%
stats = pd.read_csv("../analysis/species_overview_for_analysis.tsv",sep="\t")
# %%
print("Total number of sequences in the database:")
print(sum(stats["count"]))
print("Total number of sepecies in the database:")
print(len(stats))
print("Number of sequences from species with 50 or less sequences")
print(sum(stats["count"][(stats["count"]<50)]))

print("Number of species with 50 sequences or more:")
print(len(stats[stats["count"]>50]))


# %%
stats_non_singular = stats[stats["count"]>1]
print("Number of species with more than one sequence")
print(len(stats_non_singular))
print(len(stats_non_singular)/len(stats)*100)
# %%
species_50_100 = stats[(stats["count"]>=50) & (stats["count"]<100)]
species_100_200 = stats[(stats["count"]>=100) & (stats["count"]<200)]
species_200_400 = stats[(stats["count"]>=200) & (stats["count"]<400)]
species_400_800 = stats[(stats["count"]>=400) & (stats["count"]<800)]
species_800_1600 = stats[(stats["count"]>=800) & (stats["count"]<1600)]
species_1600_and_more = stats[stats["count"]>1600]
species_2000_and_more = stats[stats["count"]>2000]
species_50_5000 = stats[(stats["count"]>=50) & (stats["count"]<=5000)]
species_more_than_5000 = stats[stats["count"]>5000]
# %%
dfs = [species_50_100,species_100_200,species_200_400,species_400_800,species_800_1600,species_1600_and_more]
frequency = np.array([])
for df in dfs:
    frequency = np.append(frequency,len(df))
percentage = frequency/len(stats)*100
percentage_in_non_singular = frequency/len(stats_non_singular)*100
print(percentage)
print(percentage_in_non_singular)
# %%
xlabel = ['50~100','100~200','200~400','400~800','800~1600','>1600']
plt.bar(xlabel,frequency)
plt.title("Species abundancy distribution")
plt.xlabel("Number of sequences")
plt.ylabel("Number of species")
plt.show()
# %%
print(species_50_100[80:])
# %%
print(species_100_200[40:])

# %%
print(species_200_400.iloc[30:])

# %%
print(species_400_800)

# %%
print(species_800_1600)

# %%
print(species_2000_and_more)
print(len(species_2000_and_more["species"]))
# %%
#species_50_5000["species"].to_csv("species_to_take_all.tsv",sep="\t",header=False,index=False)
# %%
#species_more_than_5000["species"].to_csv("species_to_sample_5000.tsv",sep="\t",header=False,index=False)

# %%
print(stats[stats["species"]=="Enterobacter hormaechei_A"])
# %%
