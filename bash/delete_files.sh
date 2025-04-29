#!/bin/bash

# Define the starting directory
BASE_DIR="//10.69.168.1/crnldata/waking/audrey_hay/L1imaging/Analysed2025_AB/"

echo "Searching for files to delete in '$BASE_DIR'..." 

# Adjust Internal Field Separator to handle spaces in folder names
IFS=$'\n' 

# Find all 'OpenEphys' directories
find "$BASE_DIR" -type d -name "OpenEphys" | while read -r dir; do
  echo "Processing directory: $dir"
  
  # Delete all .xlsx files in the directory
  find "$dir" -type f -name "*.xlsx" -exec rm -v {} \;

  # Delete all files containing '_AB' in their name
  find "$dir" -type f -name "*AH*" -exec rm -v {} \;  

  # Delete all files containing '_AB' in their name
  find "$dir" -type f -name "*AB*" -exec rm -v {} \;  
  
  # Delete all files containing 'initial' in their name
  find "$dir" -type f -name "*initial*" -exec rm -v {} \;

# Delete all files containing 'intial' in their name
  find "$dir" -type f -name "*intial*" -exec rm -v {} \;
done

# Reset IFS to default
IFS=$'\n'
