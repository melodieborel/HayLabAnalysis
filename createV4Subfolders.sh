#!/bin/bash

# Description: This script finds and lists all directories containing .avi files and organizes them into subfolders if needed.

# Starting directory (default is the current directory)

start_directory="/crnldata/forgetting/Aurelie/CheeseboardExperiment/DAQ_data/"

echo "Searching for folders containing .avi files in '$start_directory'..."

# Use a for loop to process directories containing .avi files
avi_files=$(find "$start_directory" -type f -name "*.avi")

declare -A folder_count

for file in $avi_files; do
    dir=$(dirname "$file")
    # Track the number of .avi files in each directory
    folder_count["$dir"]=$((folder_count["$dir"] + 1))
done


# Process each directory to organize .avi files
for dir in "${!folder_count[@]}"; do
    count=${folder_count["$dir"]}

    # Check if the parent directory name is "V4miniscope"
    if [[ "$dir" == *"/My_V4_Miniscope"* && "$count" -gt 15 ]]; then
        echo "Directory '$dir' contains $count .avi files. Creating subfolders..."

        # Determine the name of the grandparent directory
        grandparent_dir=$(basename $(dirname "$dir"))

        # Create subfolders and distribute .avi files
        avi_files_in_dir=( $(ls "$dir"/*.avi | sort -V) )  # Sort files in numeric order
        subfolder_index=1
        file_count=0

        for avi_file in "${avi_files_in_dir[@]}"; do
            subfolder="$dir/${grandparent_dir}_$subfolder_index"

            # Create the subfolder if it doesn't exist
            mkdir -p "$subfolder"

            # Move the file into the subfolder
            mv "$avi_file" "$subfolder"

            file_count=$((file_count + 1))

            # Switch to a new subfolder after 15 files
            if [ "$file_count" -ge 15 ]; then
                subfolder_index=$((subfolder_index + 1))
                file_count=0
            fi
        done
    fi

done
