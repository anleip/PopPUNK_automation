import subprocess
from pathlib import Path
from datetime import datetime
import re
import argparse
import numpy as np
import math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def abbreviate_species(name: str) -> str:
    parts = name.strip().split()
    if len(parts) != 2:
        raise ValueError("Input must be in the form 'Genus species'")
    
    genus, species = parts
    return f"{genus[0].lower()}.{species.lower()}"

def decimals_from_step(step):
    s = f"{step:.10f}".rstrip("0")
    return len(s.split(".")[1])

def extractScores(thresholds, abb_species, N, db_name, decimals, run_id, suffix=None):
    """
    Run PopPUNK for each threshold, stream stderr to a log file,
    then extract 3 scores per threshold.
    """

    if suffix:
        log_file = Path(f"{abb_species}_{N}_scores_{suffix}_{run_id}.log")
        out_file = Path(f"{abb_species}_{N}_scores_{suffix}_{run_id}.tsv")
    else:
        log_file = Path(f"{abb_species}_{N}_scores_{run_id}.log")
        out_file = Path(f"{abb_species}_{N}_scores_{run_id}.tsv")

    # Run PopPUNK and stream stderr to log file
    with open(log_file, "w") as log:
        for threshold in thresholds:
            threshold_rounded = f"{threshold:.{decimals}f}"
            log.write(f"\n### Threshold {threshold_rounded} ###\n")
            log.flush()

            try:
                subprocess.run(
                    [
                        "bash",
                        "/nfs/research/jlees/anlei/bin/boundary_fitting/poppunk_fit_threshold.sh",
                        threshold_rounded,
                        str(abb_species),
                        str(db_name),
                        run_id,
                    ],
                    stdout=subprocess.DEVNULL,
                    stderr=log,      
                    text=True,
                    timeout=3600,
                    check=True,
                )
            except subprocess.TimeoutExpired as e:
                raise RuntimeError(f"Threshold {threshold_rounded} timed out") from e
            except subprocess.CalledProcessError as e:
                raise RuntimeError(f"Threshold {threshold_rounded} failed") from e

    # Parse scores from log, trusting every threshold would give exactly 3 scores.
    scores = []
    with open(log_file) as f:
        for line in f:
            if "Score" in line:
                m = re.search(r"[0-9]+\.[0-9]+", line)
                if m:
                    scores.append(float(m.group()))
    
    s = np.array(scores)
    s = s.reshape(len(thresholds),3)

    np.savetxt(out_file, s, delimiter="\t")

    return s


def plotNetworkScores(s,thresholds,species,abb_species,N,suffix=None):
    # Plot scores
    mean = np.mean(s,1)
    plt.scatter(x=thresholds,y=s[:,0],label="default score")
    plt.scatter(x=thresholds,y=s[:,1],label="score w/ betweenness")
    plt.scatter(x=thresholds,y=s[:,2],label="score w/ weighted betweenness")
    plt.plot(thresholds,mean,label="mean",color="red")
    plt.legend()
    plt.title(f"{species}, N = {N}")
    plt.xlabel("Core boundary")
    plt.ylabel("Network score")
    if suffix:
        plt.savefig(f"/nfs/research/jlees/anlei/bin/boundary_fitting/output/{abb_species}_{N}_scores_{suffix}.png")
    else:
        plt.savefig(f"/nfs/research/jlees/anlei/bin/boundary_fitting/output/{abb_species}_{N}_scores.png")


def main():
    threshold_list = np.linspace(0.0015,0.0055,17)
    parser = argparse.ArgumentParser(description='Calculate and plot network scores for a database')
    parser.add_argument("species",type=str)
    parser.add_argument("N",type=int)
    parser.add_argument("db_name",type=str)
    parser.add_argument("--suffix",type=str,default=None)
    parser.add_argument("--run-id", default=None)
    
    args = parser.parse_args()
    
    abb_species = abbreviate_species(args.species)
    step = threshold_list[1] - threshold_list[0]
    decimals = decimals_from_step(step)
    run_id = args.run_id or datetime.now().strftime("%Y%m%d_%H%M%S")
    s = extractScores(threshold_list, abb_species, args.N, args.db_name, decimals, run_id, args.suffix)
    plotNetworkScores(s,threshold_list, args.species, abb_species, args.N, args.suffix)

if __name__ == "__main__":
    main()