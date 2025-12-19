
while IFS= read -r species; do   
    python /nfs/research/jlees/anlei/bin/extract_species_subset.py "${species}" 5000 -o /nfs/research/jlees/anlei/analysis/pp-mod-5000-2nd  
done < /nfs/research/jlees/anlei/bin/species_to_sample_5000.tsv