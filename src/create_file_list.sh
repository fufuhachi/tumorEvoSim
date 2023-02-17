#!/bin/bash

# Input file containing list of folders
folder_list=$1

# Output file to store list of files
output_file="all_files.txt"

# Remove output file if it already exists
rm -f $output_file

# Loop through each folder in the input file
while read folder; do
    # Strip any trailing newline character
    folder=$(echo $folder | tr -d '\n')

    # Check if folder exists
    if [ -d "$folder" ]; then
        # Append list of files in folder to output file
        ls $folder >> $output_file
    else
        echo "Folder $folder not found. Skipping."
    fi
done < $folder_list
