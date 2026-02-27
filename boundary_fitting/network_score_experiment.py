# %%
import subprocess
import re
import numpy as np
import matplotlib.pyplot as plt
# %%
# Use Salmonella enterica as an example
# A 20k-sample database

#dists20k = np.load("../analysis/s.enterica/s.enterica_20k_db/s.enterica_20k_db.dists.npy")
# %%
#pi = dists20k[:,0]
#print(max(pi))
# %%
# Run bash script and capture output
#result = subprocess.run(["sbatch", "./poppunk_fit_threshold.sh"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
# %%
#print(result.stderr)
# %%
#print(result.stdout)

# %%
def extractScores(thresholds):
    # Calculate scores using PopPUNK and extract output
    scores = []
    for threshold in thresholds:
        result = subprocess.run(["bash", "./poppunk_fit_threshold.sh", str(threshold)], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        print(result.stderr)
        score = []
        for line in result.stderr.splitlines():
            if "Score" in line:
                m = re.search(r"[0-9]+\.[0-9]+", line)
                if m:
                    score.append(float(m.group()))
        scores.append(score)
    s = np.array(scores)
    return s

def plotNetworkScores(s,thresholds,species,N):
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
    plt.show()

# %%
# Range 0 to 0.01, covering the clsuters closest to origin
threshold_list = np.linspace(0.0005,0.01,20)
s = extractScores(threshold_list)
plotNetworkScores(s,threshold_list,"Salmonella enterica",5000)

# %%
# Double resolution
threshold_list = np.linspace(0.0015,0.0055,17)
s = extractScores(threshold_list)
plotNetworkScores(s,threshold_list,"Salmonella enterica",5000)

# %%
print(np.linspace(0.0045,0.006,7))
# %%
s = np.loadtxt("/nfs/research/jlees/anlei/analysis/s.enterica/salmonella_enterica500_scores_zoomed0.001-0.012.tsv",delimiter="\t")
mean_s = np.mean(s,1)
print(mean_s)
# %%
plotNetworkScores(s,np.linspace(0.0010,0.0120,45),"Salmonella enterica",500)
# %%
# Extracting network scores from slurm output text file
with open("/nfs/research/jlees/anlei/analysis/l.monocytogenes/l.monocytogenes_20000_scores_20260108_101503.log") as f:
    result = f.read()
i = 0
score = []
for line in result.splitlines():
    if "Score" in line:
        m = re.search(r"[0-9]+\.[0-9]+", line)
        if m:
            score.append(float(m.group()))
            i+=1
print(i)
s = np.array(score)
scores = s.reshape(23,3)
print(scores)
#np.savetxt("/nfs/research/jlees/anlei/analysis/s.aureus/s.aureus_5000_scores_zoomed.tsv",s_new,delimiter="\t")
# %%
plotNetworkScores(scores,np.linspace(0.0005,0.0060,23),"Listeria monocytogenes",20000)
# %%
# Getting core boundary from maximum score
scores = np.loadtxt("/nfs/research/jlees/anlei/analysis/s.enterica/s.enterica_20000_scores_test_20260220_142524.tsv")
# %%
mean = np.mean(scores,1)
x = np.argmax(scores[:,0])
print(x)
y = np.argmax(scores[:,1])
print(y)
z = np.argmax(scores[:,2])
print(z)
m = np.argmax(mean)
print(m)
# %%
threshold_list = np.linspace(0.0015,0.0055,17)
xb = threshold_list[x]
yb = threshold_list[y]
zb = threshold_list[z]
mb = threshold_list[m]
print(f"{xb:.5f},{yb:.5f},{zb:.5f},{mb:.5f}")
print(scores[x],scores[y],scores[z],scores[m])

# %%
names = ["default","w/ betweenness","w/ weighted betweenness", "mean"]
xa = [scores[x,0],scores[y,0],scores[z,0],scores[m,0]]
ya = [scores[x,1],scores[y,1],scores[z,1],scores[m,1]]
za = [scores[x,2],scores[y,2],scores[z,2],scores[m,2]]
ma = [mean[x],mean[y],mean[z],mean[m]]
plt.scatter(names,xa,label="default score")
plt.scatter(names,ya,label="score w/ betweenness")
plt.scatter(names,za,label="score w/ weighted betweenness")
plt.plot(names,ma,color="red",label="mean")
plt.legend()
plt.show()
# %%
def decimals_from_step(step):
    s = f"{step:.10f}".rstrip("0")
    return len(s.split(".")[1])
# %%
threshold_list = np.linspace(0.0005,0.0065,25)
step = threshold_list[1] - threshold_list[0]
print(step)
decimals = decimals_from_step(step)
print(decimals)
# %%
a = 0.000250000000003
b = decimals_from_step(a)
print(b)
# %%
