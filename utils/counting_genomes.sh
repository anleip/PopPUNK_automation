for f in /nfs/research/jlees/anlei/analysis/e.coli/e.coli_MLST/result/*; do
    lines=$(wc -l < "$f")
    echo $lines
done

# then run 
# ../bin/counting_genomes.sh | awk 'BEGIN{sum=0}{sum+=$1}END{print sum}'