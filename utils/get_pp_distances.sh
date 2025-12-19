

# Current directory is the one contains rfiles folder
# Make an array of "genus_species" (strip .rfile from filenames)
species_array=($(ls *.rfile 2>/dev/null | xargs -n1 basename | sed 's/\.rfile$//'))
# Convert .npy and .pkl files to text
for species in ${species_array[@]}; do 
    /nfs/research/jlees/anlei/bin/poppunk_edits/PopPUNK/scripts/poppunk_extract_distances.py --distances "${species}_qc/${species}_qc.dists" --output "../pp_distances/${species}.dists.out"

done