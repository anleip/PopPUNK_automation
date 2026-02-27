#!/bin/bash

# Base parameters
TIME="0:20:00"
CPUS=16
MEM="5G"
#DB_PATH="/nfs/research/jlees/anlei/analysis/pubmlst"
DB_PATH="/hps/nobackup/rdf/metagenomics/research-team/johanna/projects/PhD_main/analysis/2025-04-09-TAC2-results/MLST/pubmlst"
SCHEME="neisseria"

mkdir -p result

# Loop from 001 to 020
for i in $(seq 20); do
  printf -v index "%03d" "$i"
  echo "Submitting job for batch_${index}..."

  sbatch --job-name=mlst_${index} \
         --time=${TIME} \
         --cpus-per-task=${CPUS} \
         --mem=${MEM} \
         --wrap="fastmlst batch_${index}/* --db_path ${DB_PATH} --scheme ${SCHEME} -t ${CPUS} -to result/batch_${index}.csv"

  echo "Job for batch_${index} submitted"
done

echo "All jobs submitted!"

#fastmlst batch_005/* --db_path /hps/nobackup/rdf/metagenomics/research-team/johanna/projects/PhD_main/analysis/2025-04-09-TAC2-results/MLST/pubmlst --scheme "ecoli#1" -t 16 -to result/batch_005.csv


