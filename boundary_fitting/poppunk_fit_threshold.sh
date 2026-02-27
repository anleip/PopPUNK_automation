
#!/bin/bash

# Check if the user provided a threshold
if [ "$#" -ne 4 ]; then
  echo "Usage: $0 <threshold> <abb_species> <db_name> <run_id>"
  exit 1
fi

# Set the user-defined threshold
THRESHOLD=$1
abb_species=$2
db_name=$3
run_id=$4

OUTDIR="/nfs/research/jlees/anlei/analysis/${abb_species}/${db_name}_fit_${run_id}"
OUTNAME="${db_name}_fit${THRESHOLD}"

mkdir -p "$OUTDIR"

# Prevent accidental overwrite / hang
if [ -e "${OUTDIR}/${OUTNAME}" ]; then
  echo "ERROR: Output ${OUTDIR}/${OUTNAME} already exists" >&2
  exit 2
fi

# Run poppunk with the user-defined threshold
poppunk --fit-model threshold \
    --ref-db /nfs/research/jlees/anlei/analysis/"${abb_species}"/"${db_name}" \
    --threshold ${THRESHOLD} \
    --output "${OUTDIR}/${OUTNAME}" \
    --threads 32
