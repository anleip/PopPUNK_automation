#!/usr/bin/env bash

# Usage:
#   ./extract_species_samples.sh "Enterococcus faecalis"
#   ./extract_species_samples.sh "Enterococcus faecalis" /path/to/output_dir

# === CONFIGURATION ===
ANNOT_DIR="/nfs/research/jlees/anlei/analysis"
ANNOT_PATTERN="${ANNOT_DIR}/all_annotations/2kk_hq_annotations*.csv"
RFILE="${ANNOT_DIR}/2kk.rfile"

# === INPUT CHECKS ===
if [[ $# -lt 1 || $# -gt 2 ]]; then
    echo "Usage: $0 \"<species name>\" [output_directory]"
    echo "Example: $0 \"Enterococcus faecalis\" ./results"
    exit 1
fi

SPECIES="$1"
OUTDIR="${2:-$(pwd)}"   # Default to current directory if not provided

# Convert species name to safe lowercase filename, e.g. "Enterococcus faecalis" → "enterococcus_faecalis"
OUTBASE=$(echo "$SPECIES" | tr '[:upper:]' '[:lower:]' | tr ' ' '_')
OUTFILE="${OUTDIR}/${OUTBASE}.rfile"

# === MAIN ===
echo "Extracting sample IDs for species: \"$SPECIES\""
echo "Output will be saved to: $OUTFILE"

# Ensure output directory exists
mkdir -p "$OUTDIR"

# Empty the output file before starting
> "$OUTFILE"

# Temporary file for sample IDs — created in current working directory
TMP_IDS="${OUTDIR}/${OUTBASE}_tmp_ids.txt"
> "$TMP_IDS"

# Step 1: Extract sample IDs matching the species name
for file in $ANNOT_PATTERN; do
    if [[ -f "$file" ]]; then
        grep -F "$SPECIES" "$file" | cut -d',' -f1 >> "$TMP_IDS"
    fi
done

# Step 2: Deduplicate sample IDs
sort -u "$TMP_IDS" -o "$TMP_IDS"

# Step 3: Extract matching lines from 2kk.rfile using exact (whole word) match
echo "Matching all sample IDs..."
grep -wFf "$TMP_IDS" "$RFILE" >> "$OUTFILE"

# Cleanup
rm -f "$TMP_IDS"

echo "Done!"
