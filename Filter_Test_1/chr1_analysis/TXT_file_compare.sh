#!/bin/bash

# Path to the directory containing the .txt files
directory_path="../chr1_analysis/results"

# Find the first .txt file to serve as a reference
reference_file=$(find "$directory_path" -type f -name "*.txt" | head -1)

# Flag to indicate if all files are the same
all_same=true

# Loop through all .txt files in the directory
find "$directory_path" -type f -name "*.txt" | while read -r file; do
    # Skip the reference file itself
    if [[ "$file" == "$reference_file" ]]; then
        continue
    fi

    # Compare the current file with the reference file
    if ! diff -q "$reference_file" "$file" > /dev/null; then
        all_same=false
        echo "File $file is different."
        # Exit the loop early since we found a difference
        break
    fi
done

# Final check
if $all_same; then
    echo "All .txt files are exactly the same."
else
    echo "Not all .txt files are the same."
fi
