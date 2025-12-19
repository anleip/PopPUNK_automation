import numpy as np
import argparse


def percentile_jump_detection(db_name,species,x=0.2,r=50):
    """Finds jumps in distances. Find a strict maximum distance percentiles for QC.
    
    Args:
        species (str)
            Name of the species in the form of genus_species
        x (float)
            Sensitivity. Proportion increase threshold between distances of consecutive 1 percentiles to be consider as a jump.
        r (float)
            Resolution of detection in terms of number of distances.
    """ 
    dists = np.load(f"{db_name}/{db_name}.dists.npy")

    L = len(dists)  # Total number of distance pairs
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
    
    #print(f"Detecting maximum distance cut-offs for {species}")
    #print(f"Using x={x}, r={r}")

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

    
    if len(pi_jumps) != 0:
        max_pi = min(pi_jumps)
        #print(f"Core distance jumps at \n{pi_jumps_pc} percentiles, π = \n{pi_jumps}")
    else:
        max_pi = max(pi)
        #print("No trimming needed for core distance")

    if len(a_jumps) != 0:
        max_a = min(a_jumps)
        #print(f"Accessory distance jumps at \n{a_jumps_pc} percentiles, a = \n{a_jumps}")
    else:
        max_a = max(a)
        #print("No trimming needed for accessory distance")
    
    return max_pi, max_a


def main():
    parser = argparse.ArgumentParser(description="Find max core and max accessory for PopPUNK QC")
    parser.add_argument("db_name",type=str,help="Name of database to be QC")
    parser.add_argument("species",type=str,help="Species name (e.g. escherichia_coli)")
    parser.add_argument("--x", type=float,default=0.2,
    help="Sensitivity. Proportion increase threshold between distances of consecutive 1 percentiles to be consider as a jump.")
    parser.add_argument("--r",type=int,default=50,
    help="Resolution of detection in terms of number of distances.")
    
    args = parser.parse_args()      

    max_pi, max_a = percentile_jump_detection(args.db_name,args.species,args.x,args.r)
    print(f"max_pi={max_pi}")
    print(f"max_a={max_a}")
    #print(f"\nDone!")

if __name__ == "__main__":
    main()

