#!/bin/bash

# Define the starting directory (default: current directory)
#START_DIR="."
START_DIR="/crnldata/forgetting/Aurelie/CheeseboardExperiment/"
START_DIR="/crnldata/forgetting/Aurelie/CheeseboardExperiment/DAQ_data/AB/Training/Lou/Cheeseboard/2024_12_02/12_28_00/"


# Find all directories named "My_V4_Miniscope"
echo "Searching for folders named 'My_V4_Miniscope' in $START_DIR..."
# Adjust Internal Field Separator to handle spaces in folder names
IFS=$'\n' 

# Loop through the results
for pathtofolder in $(find "$START_DIR" -type d -name "My_V4_Miniscope"); do
    echo "Found: $pathtofolder"

    rm -rf /mnt/data/minianAB/*
    cp -r "${pathtofolder}" /mnt/data/minianAB/
    srun --mem=80G --cpus-per-task=10 python /home/aurelie.brecier/HayLabAnalysis/python/MinianCluster.py
    # Check the exit status of srun
    if [ $? -ne 0 ]; then
        echo "Error: srun failed for $pathtofolder. Skipping..."
        continue
    fi
    cp -r /mnt/data/minianAB/ "${pathtofolder}"
    rm -rf /mnt/data/minianAB/*

done

# Reset IFS to default
IFS=$'\n'