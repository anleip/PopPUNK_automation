#!/bin/bash

# Current directory is the one contains rfiles folder
# Make an array of "genus_species" (strip .rfile from filenames)
# For cases where there is only one database for each species
species_array=($(ls *.rfile 2>/dev/null | xargs -n1 basename | sed 's/\.rfile$//'))

#echo "${species_array[@]}"

for species in ${species_array[@]}; do 
    N=$(wc -l < "${species}.rfile")
    echo "Sketching $N genomes for ${species}..."
    /nfs/research/jlees/anlei/bin/poppunk_edits/PopPUNK/poppunk-runner.py --create-db \
     --output ${species}_db --r-files ${species}.rfile --sketch-size 100000 --min-k 17 --max-k 41 --plot-fit 3 --threads 48
    
    printf "Looking for core and accessory distance cut-offs...\n"
    eval $(python /nfs/research/jlees/anlei/bin/qc_automation/qc_automated.py "${species}_db" "${species}")
    printf "\nMaximum core distance for QC: $max_pi"
    printf "\nMaximum accessory distance for QC: $max_a"
    printf "\nPerforming quality control using the maximum distances...\n"
    
    /nfs/research/jlees/anlei/bin/poppunk_edits/PopPUNK/poppunk-runner.py \
    --qc-db --ref-db ${species}_db --max-pi-dist $max_pi --max-a-dist $max_a --max-zero-dist 10 --output ${species}_qc --threads 48
done