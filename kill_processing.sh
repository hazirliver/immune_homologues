#!/bin/bash

# Get the list of process IDs for all processes containing "string_db"
pid_list=$(pgrep -f "string_db")

# If the list is empty, exit
if [[ -z "$pid_list" ]]; then
  echo "No processes containing 'string_db' found."
  exit 1
fi

# Kill all processes in the list
kill $pid_list

echo "All processes containing 'string_db' have been terminated."
exit 0
