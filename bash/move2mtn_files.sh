#!/bin/bash

# cd bash
# ./move2mtn_files.sh

# Define the starting directory
SRC="/crnldata/forgetting/Aurelie/MiniscopeOE_data/L1NDNF_mice/"

echo "Searching for folders in '$SRC'..." 

# Adjust Internal Field Separator to handle spaces in folder names
IFS=$'\n' 

DEST="/mnt/data/TheaM_other/"

cd "$SRC"

find . \
    -type f \
    ! -name "*.avi" \
    ! -path "*/BW_notused/*" \
    ! -path "*/NO_nodata/*" \
    ! -path "*/NR_notdata/*" \
    ! -path "*/Habituation/*" \
    ! -name "*.daq" \
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