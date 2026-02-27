#!/bin/bash

species=$1
N=$2

# Convert to genus_species
hyphened_species=$(echo "$species" | tr 'A-Z' 'a-z' | sed 's/ /_/g')

# Generate rfile
python /nfs/research/jlees/anlei/bin/PopPUNK_automation/utils/extract_species_subset.py "${species}" $N -o .
mv -f ${hyphened_species}.id ${hyphened_species}_${N}.id
mv -f ${hyphened_species}.rfile ${hyphened_species}_${N}.rfile

# Sketch
echo "Sketching $N genomes for ${species}..."
/nfs/research/jlees/anlei/bin/poppunk_edits/PopPUNK/poppunk-runner.py --create-db \
--output ${hyphened_species}_${N}_db --r-files ${hyphened_species}_${N}.rfile --sketch-size 100000 --min-k 17 --max-k 41 --plot-fit 3 --threads 48

# Find QC thresholds
printf "Looking for core and accessory distance cut-offs...\n"
eval $(python /nfs/research/jlees/anlei/bin/PopPUNK_automation/qc_automation/qc_automated.py "${hyphened_species}_${N}_db" "${hyphened_species}")
printf "\nMaximum core distance for QC: $max_pi"
printf "\nMaximum accessory distance for QC: $max_a"
printf "\nPerforming quality control using the maximum distances...\n"
# Perform QC
/nfs/research/jlees/anlei/bin/poppunk_edits/PopPUNK/poppunk-runner.py --qc-db \
--ref-db ${hyphened_species}_${N}_db --max-pi-dist $max_pi --max-a-dist $max_a --max-zero-dist 10 --output ${hyphened_species}_${N}_qc --threads 48
