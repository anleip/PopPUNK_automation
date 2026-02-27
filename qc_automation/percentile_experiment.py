#%%
import numpy as np
import matplotlib.pyplot as plt

#%%

dists = np.load("../analysis/e.coli/e.coli_5k_db/e.coli_5k_db.dists.npy")
print(len(dists))
pi = dists[:,0]
a = dists[:,1]
# %%
pi_pcs = np.percentile(pi,np.linspace(0.1,100,1000))
a_pcs = np.percentile(a,np.linspace(0.1,100,1000))
# %%
def percentile_experiment(species="e.coli",data_size="5k",x=1.1,r=5):
    """Finds jumps in distances. Visualisation of cumulative plot. Find a strict maximum distance percentiles for QC.
    
    Args:
        species (str)
            Name of bacteria species in the form of g.species
        data_size (str)
            Number of genomes used in the raw sketch, in unit of k.
        x (float)
            Sensitivity. Multiplication factor threshold between consecutive 1 percentiles to be consider as a jump.
        r (float)
            Resolution of detection in terms of number of genomes.
    """
    dists = np.load(f"../analysis/{species}/{species}_{data_size}_db/{species}_{data_size}_db.dists.npy")
    #Look for the number of bins that gives a resolution of r genomes per step
    N = int(data_size[0])
    n = int(N*1000/r)   # n = number of bins
    y = int(n/100-1)    # y = (number of bins in one percentile) -1. n/100 needs to be an integer, so r must be a factor of 1000.
    percentiles = np.linspace(100/n,100,n)

    pi = dists[:,0]
    a = dists[:,1]
    pi_pcs = np.percentile(pi,percentiles)
    a_pcs = np.percentile(a,percentiles)
    
    #Jump detection
    pi_jumps = []; pi_jumps_pc=[]
    a_jumps = []; a_jumps_pc=[]
    j = 0
    for pcs in [pi_pcs,a_pcs]:    
        #Start searching for jumps from 80% upwards, stopping at i+1 = n.
        for i in range(int(len(pcs)*0.8),len(pcs)-1):
            if pcs[i-y]*x < pcs[i+1] and j == 0:
                pi_jumps.append(pcs[i])
                pi_jumps_pc.append(percentiles[i])
            if pcs[i-y]*x < pcs[i+1] and j == 1:
                a_jumps.append(pcs[i])
                a_jumps_pc.append(percentiles[i])
        j +=1
    
    #Quantile plot
    plt.plot(percentiles,pi_pcs)
    plt.vlines(pi_jumps_pc,ymin=0,ymax=pi_pcs[-1],colors="red",linestyles="dotted")
    plt.xlabel("Percentile (%)")
    plt.ylabel("π")
    plt.title(species + f" Core Distances, x={x}")
    plt.show()

    plt.plot(percentiles,a_pcs)
    plt.vlines(a_jumps_pc,ymin=0,ymax=a_pcs[-1],colors="red",linestyles="dotted")
    plt.xlabel("Percentile (%)")
    plt.ylabel("a")
    plt.title(species + f" Accessory Distances, x={x}")
    plt.show()

    print(f"Core distance jumps at \n{pi_jumps_pc} percentiles, with values \n{pi_jumps}")
    print("\n")
    print(f"Accessory distance jumps at \n{a_jumps_pc} percentiles, with values \n{a_jumps}")

    return pi_jumps_pc, a_jumps_pc


#percentile_experiment()

# %%
percentile_experiment(x=1.1)
# %%
percentile_experiment("s.pneumoniae","5k",x=1.2)
# %%
percentile_experiment("m.tuberculosis","2k",x=1.3)
# %%
percentile_experiment("p.aeruginosa","5k",x=1.1)
# %%
percentile_experiment("n.gonorrhoeae","2k", x=1.1)
# %%
percentile_experiment("n.meningitidis","5k")
# %%
percentile_experiment("s.enterica","5k",x=1.2)
# %%
percentile_experiment("l.monocytogenes","5k",x=1.2)
# %%
percentile_experiment("s.aureus","5k",x=1.1)
# %%
def percentile_experiment_universal(db_name,dist_name,species="e.coli",x=0.1,r=50):
    """Finds jumps in distances. Visualisation of cumulative plot. Find a strict maximum distance percentiles for QC.
    
    Args:
        db_name (str)
            Name of PopPUNK database. The folder that contains the .npy and .pkl files.
        dist_name (str)
            Name of the .npy and .pkl files. Prefix before .dist.xxx.
        species (str)
            Name of the species in the form of g.species
        step (int)
            Number of steps across two distances to compare between in junp detection. Smaller step, more sensitive to fluctuation in data.
        x (float)
            Sensitivity. Proportion increase threshold between distances of consecutive 1 percentiles to be consider as a jump.
        r (float)
            Resolution of detection in terms of number of distances.
    """
    import pypickle as ppk 
    dists = np.load(f"/nfs/research/jlees/anlei/analysis/{species}/{db_name}/{dist_name}.dists.npy")
    names = ppk.load(f"/nfs/research/jlees/anlei/analysis/{species}/{db_name}/{dist_name}.dists.pkl")
    L = len(dists)  # Total number of distance pairs
    N = len(names[0])   # Number of sampels
    #Look for the number of bins that gives a resolution of r genomes per step
    n = int(L/r)   # n = number of bins
    step = int(n//100)  # number of bins in 1%
    s = step-1
    y = 100*step*x/n+1  # y = modified x. Same sensitivity as x but across steps.
    percentiles = np.linspace(100/n,100,n)

    pi = dists[:,0]
    a = dists[:,1]
    pi_pcs = np.percentile(pi,percentiles)
    a_pcs = np.percentile(a,percentiles)
    
    #Jump detection
    pi_jumps = []; pi_jumps_pc=[]
    a_jumps = []; a_jumps_pc=[]
    j = 0
    for pcs in [pi_pcs,a_pcs]:    
        #Start searching for jumps from 75% upwards, stopping at i+1 = n.
        for i in range(int(len(pcs)*0.75),len(pcs)-1):
            if pcs[i-s]*y < pcs[i+1] and j == 0:
                pi_jumps.append(pcs[i])
                pi_jumps_pc.append(percentiles[i])
            if pcs[i-s]*y < pcs[i+1] and j == 1:
                a_jumps.append(pcs[i])
                a_jumps_pc.append(percentiles[i])
        j +=1
    
    #Quantile plot
    plt.plot(percentiles,pi_pcs)
    plt.vlines(pi_jumps_pc,ymin=0,ymax=pi_pcs[-1],colors="red",linestyles="dotted")
    plt.xlabel("Percentile (%)")
    plt.ylabel("π")
    plt.title(species + " Core Distances" + f" n={N}, x={x}")
    plt.show()

    plt.plot(percentiles,a_pcs)
    plt.vlines(a_jumps_pc,ymin=0,ymax=a_pcs[-1],colors="red",linestyles="dotted")
    plt.xlabel("Percentile (%)")
    plt.ylabel("a")
    plt.title(species + " Accessory Distances" + f" n={N}, x={x}")
    plt.show()

    

    if len(pi_jumps) != 0:
        max_pi = min(pi_jumps)
        print(f"Core distance jumps at \n{pi_jumps_pc} percentiles, with values \n{pi_jumps}")
    else:
        max_pi = max(pi)
        print("No trimming needed for core distance")

    if len(a_jumps) != 0:
        max_a = min(a_jumps)
        print(f"Accessory distance jumps at \n{a_jumps_pc} percentiles, with values \n{a_jumps}")
    else:
        max_a = max(a)
        print("No trimming needed for accessory distance")
    
    return max_pi, max_a


# %%

percentile_experiment_universal(db_name="salmonella_enterica_30k_reference_qc_pi0.06",dist_name="salmonella_enterica_30k_reference_qc_pi0.06",species="s.enterica",x=0.01)


# %%
# Comparing the coarse method and the fine method
percentile_experiment("l.monocytogenes","5k",x=1.2)
percentile_experiment_universal(db_name="l.monocytogenes_5k_db",dist_name="l.monocytogenes_5k_db",species="l.monocytogenes",x=0.2)
# %% s.aureus
# Repeated QC, see if more genomes get trimmed
# First extra QC
#percentile_experiment_universal(db_name="qc_5k_x0.1_r50",dist_name="qc_5k_x0.1_r50",species="l.monocytogenes")
# Second extra QC
#percentile_experiment_universal(db_name="qc_5k_x1.1_r5_2nd",dist_name="qc_5k_x1.1_r5_2nd",species="l.monocytogenes")
#%%
#percentile_experiment("s.aureus","5k")
#percentile_experiment_universal(db_name="s.aureus_5k_db",dist_name="s.aureus_5k_db",species="s.aureus",x=0.2)
# First extra QC
#percentile_experiment_universal(db_name="qc_5k_x0.2_r50",dist_name="qc_5k_x0.2_r50",species="s.aureus",x=0.2)
# Second extra QC
percentile_experiment_universal(db_name="qc_5k_x0.1_r50_2nd",dist_name="qc_5k_x0.1_r50_2nd",species="s.aureus",x=0.2)

# %% s.pneumoniae
percentile_experiment("s.pneumoniae","5k",x=1.2)
percentile_experiment_universal(db_name="s.pneumoniae_5k_db",dist_name="s.pneumoniae_5k_db",species="s.pneumoniae",x=0.2)
# %%
# First extra QC
#percentile_experiment_universal(db_name="qc_5k_x0.2_r50",dist_name="qc_5k_x0.2_r50",species="s.pneumoniae",x=0.2)
# Second extra QC detected a previously ignored jump.
percentile_experiment_universal(db_name="qc_5k_x0.2_r50_3rd",dist_name="qc_5k_x0.2_r50_3rd",species="s.pneumoniae",x=0.2)
# Third extra QC
percentile_experiment_universal(db_name="qc_5k_x0.1_r50_3rd",dist_name="qc_5k_x0.1_r50_3rd",species="s.pneumoniae",x=0.1)
# Fourth extra

# %%
percentile_experiment(x=1.1)
percentile_experiment_universal(db_name="qc_5k_x1.1_r5",dist_name="qc_5k_x1.1_r5",species="e.coli")
# %%
percentile_experiment("p.aeruginosa","5k",x=1.1)
percentile_experiment_universal(db_name="qc_5k_x1.1_r5",dist_name="qc_5k_x1.1_r5",species="p.aeruginosa")
# %%
percentile_experiment("m.tuberculosis","2k",x=1.1)
#percentile_experiment_universal(db_name="qc_2k_x1.1_r5",dist_name="qc_2k_x1.1_r5",species="m.tuberculosis")



# %%
# Experimenting with histogram plot
dists = np.load("../analysis/s.enterica/s.enterica_5k_db/s.enterica_5k_db.dists.npy")
pi = dists[:,0]
a = dists[:,1]

plt.hist(pi,1000)
plt.xlim(0,max(pi))
plt.show()
# %%
def hist_plot(db_name,species="e.coli",data_size="5k",bins=1000,accessory=False):
    """Plot histograms of distances from database after QC
    
    
    """
    
    dists = np.load(f"../analysis/{species}/{db_name}/{db_name}.dists.npy")
    pi = dists[:,0]
    a = dists[:,1]

    print(db_name)
    plt.hist(pi,bins,)
    plt.xlim(0,max(pi))
    plt.xlabel("π")
    plt.ylabel("Geonome count")
    plt.title(f"{species} {data_size} core distance histogram")
    plt.show()
    if accessory:
        plt.hist(a,bins,color="orange")
        plt.xlim(0,max(a))
        plt.xlabel("a")
        plt.ylabel("Geonome count")
        plt.title(f"{species} {data_size} accessory distance histogram")
        plt.show()
# %%
hist_plot("e.coli_5k_x1.1_r5_num0_db",accessory=True)
hist_plot("e.coli_5k_num1_x1.1_db")
hist_plot("e.coli_5k_num2_x1.1_db")
hist_plot("e.coli_5k_num3_x1.1_db")
# %%
hist_plot("s.pneumoniae_5k_num0_x1.1_db","s.pneumoniae",accessory=True)
# %%
hist_plot("m.tuberculosis_2k_num0_x1.1_db","m.tuberculosis","2k",accessory=True)
# %%
hist_plot("p.aeruginosa_5k_x1.1_db","p.aeruginosa",accessory=True)
# %%
hist_plot("n.gonorrhoeae_2k_num0_x1.1_db","n.gonorrhoeae","2k",accessory=True)
# %%
hist_plot("n.meningitidis_5k_x1.1_db","n.meningitidis",accessory=True)
# %%
hist_plot("s.enterica_5k_num0_x1.1_db","s.enterica",accessory=True)
# %%
hist_plot("l.monocytogenes_5k_x1.1_db","l.monocytogenes",accessory=True)
# %%
hist_plot("s.aureus_5k_x1.1_db","s.aureus",accessory=True)
# %%
