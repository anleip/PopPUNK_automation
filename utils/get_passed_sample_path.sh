# Current directory is the one contains rfiles folder
# Make an array of "genus_species" (strip .rfile from filenames)
species_array=($(ls *.rfile 2>/dev/null | xargs -n1 basename | sed 's/\.rfile$//'))
# Extract failed smaple names from qcreport.txt, remove lines from rfile, and retain only the second column (file paths).
for species in ${species_array[@]}; do 
    grep -v -F -f <(cut -f1 "${species}_qc/${species}_qc_qcreport.txt") "${species}.rfile" | cut -f2 > ../qc_passed_paths/"${species}_passed_paths.txt"
done