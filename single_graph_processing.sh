#!/bin/bash

N=$1
FILES=$(find ./Temporary_files/string_db/ -type f)

# Split files into chunks of size N
CHUNKS=$(echo "$FILES" | xargs -n $N | sed 's/ /,/g')

# Loop through chunks and run N parallel instances of the script
for chunk in $CHUNKS; do
  for file in $(echo $chunk | tr ',' ' '); do
    python string_db_parser.py $file &
  done
  wait
done