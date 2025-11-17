#!/bin/bash

# Define the starting directory
SRC="/crnldata/forgetting/Aurelie/MiniscopeOE_data/L2_3_mice/"

echo "Searching for folders in '$SRC'..." 

# Adjust Internal Field Separator to handle spaces in folder names
IFS=$'\n' 

DEST="/mnt/data/AurelieB_other/"

cd "$SRC"

find . \
    -type f \
    -path "*/Exploration_task/*" \
    -path "*/Cheeseboard/*" \
    ! -path "*/OpenEphys/*" \
    ! -name "*.avi" \
    -print0 | while IFS= read -r -d '' file; do
        echo "Copying: $file"
        # Create directory structure in DEST
        mkdir -p "$DEST/$(dirname "$file")"
        # Copy file
        cp "$file" "$DEST/$file"
    done


echo "Transfer complete."


# Reset IFS to default
IFS=$'\n'