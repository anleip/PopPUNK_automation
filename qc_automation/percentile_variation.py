# %%
import numpy as np
import matplotlib.pyplot as plt

# %%
def percentile_variation(species="e.coli",data_size="5k",num=0,x=1.1,r=5):
    """Finds the cut off for distances. Find a strict maximum distance percentiles for QC.
    
    Args:
        species (str)
            Name of bacteria species in the form of g.species
        data_size (str)
            Number of genomes used in the raw sketch, in unit of k.
        num (int)
            Sample number
        x (float)
            Sensitivity. Multiplication factor threshold between consecutive 1 percentiles to be consider as a jump.
        r (float)
            Resolution of detection in terms of number of genomes.
    """
    dists = np.load(f"../analysis/{species}/{species}_{data_size}_num{num}_db/{species}_{data_size}_num{num}_db.dists.npy")
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

    return min(pi_jumps), min(a_jumps)

    
# %%
pi,a = percentile_variation("s.pneumoniae",num=3,x=1.2)
# %%
print(pi)
print(a)
# %%
