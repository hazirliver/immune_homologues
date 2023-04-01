#!/bin/bash

N=$1
FILES=$(find ./Temporary_files/string_db/ -type f -name "*_*.tsv")

# Split files into chunks of size N
CHUNKS=$(echo "$FILES" | xargs -n $N | sed 's/ /,/g')

# Loop through chunks and run N parallel instances of the script
for chunk in $CHUNKS; do
  PIDS=()
  for file in $(echo $chunk | tr ',' ' '); do
    python string_db_parser_double.py $file &
    PIDS+=("$!")
  done

  # Wait for all the PIDs in the PIDS array to finish
  for pid in "${PIDS[@]}"; do
    wait $pid
  done
done