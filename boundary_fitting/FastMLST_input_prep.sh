#!/bin/bash

# Had a file of 20,000 genomes. Need to put them into folders of 1000 genomes.

set -euo pipefail

LIST="neisseria_gonorrhoeae_20k_paths.tsv"
DEST="n.gonorrhoeae_MLST"
BATCH_SIZE=1000

mkdir -p "$DEST"

batch=1

# Split the list into chunks of 1000 lines and process each chunk
split -l "$BATCH_SIZE" -d -a 3 "$LIST" /tmp/file_batch_

for chunk in /tmp/file_batch_*; do
    batch_dir=$(printf "%s/batch_%03d" "$DEST" "$batch")
    mkdir -p "$batch_dir"

    xargs -d '\n' -a "$chunk" cp -t "$batch_dir"

    ((batch++))
done

# Cleanup
rm /tmp/file_batch_*
