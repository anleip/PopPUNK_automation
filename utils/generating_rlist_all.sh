#mapfile -t species_list < /nfs/research/jlees/anlei/bin/species_to_take_all.tsv
#echo "${species_list[@]}"
while IFS= read -r species; do 
    /nfs/research/jlees/anlei/bin/extract_species_all.sh "${species}" /nfs/research/jlees/anlei/analysis/pp-mod 
done < /nfs/research/jlees/anlei/bin/species_to_take_all.tsv


#for species in ${species_list[@]}; do
#    /nfs/research/jlees/anlei/bin/extract_species_all.sh "${species}" /nfs/research/jlees/anlei/analysis/pp-mod/.
#done