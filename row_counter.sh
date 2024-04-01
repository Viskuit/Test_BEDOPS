#!/bin/bash

# Check if a directory path is provided as an argument
if [ $# -eq 0 ]; then
    echo "Usage: $0 <path_to_directory>"
    exit 1
fi

# Assign the first argument as the directory path
DIRECTORY_PATH=$1

# List CSV files in the specified directory, sort them numerically based on the numbers in filenames, and iterate over each file
for file in $(ls $DIRECTORY_PATH/*.csv | sort -V); do
    # Count the number of rows in the current CSV file
    num_rows=$(wc -l < "$file")
    
    # Print the filename and the number of rows
    echo "$file: $num_rows"
done
