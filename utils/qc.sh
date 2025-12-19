#!/bin/bash

db_name=$1
species=$2

prefix=$(echo "${db_name}" | sed 's/_db//g')

# Find QC thresholds
printf "Looking for core and accessory distance cut-offs...\n"
eval $(python /nfs/research/jlees/anlei/bin/qc_automation/qc_automated.py "${db_name}" "${species}")
printf "\nMaximum core distance for QC: $max_pi"
printf "\nMaximum accessory distance for QC: $max_a"
printf "\nPerforming quality control using the maximum distances...\n"
# Perform QC
/nfs/research/jlees/anlei/bin/poppunk_edits/PopPUNK/poppunk-runner.py --qc-db \
--ref-db ${db_name} --max-pi-dist $max_pi --max-a-dist $max_a --max-zero-dist 10 --output ${prefix}_qc --threads 48
