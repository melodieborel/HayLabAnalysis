#!/bin/bash

# Define the starting directory (default: current directory)
#START_DIR="."
START_DIR="/crnldata/forgetting/Aurelie/CheeseboardExperiment/"
START_DIR="/crnldata/forgetting/Aurelie/CheeseboardExperiment/DAQ_data/AB/Habituation/Blue/SleepBefore/2024_11_28/13_33_25/"

# Find all directories named "My_V4_Miniscope"

echo "Searching for folders containing .avi files in '$START_DIR'..." # echo "Searching for folders named 'My_V4_Miniscope' in $START_DIR..."

# Adjust Internal Field Separator to handle spaces in folder names
IFS=$'\n' 

# Loop through the results

for pathtofolder in $(find "$START_DIR" -type f -name "*.avi" -exec dirname {} \; | sort -u); do

    if [[ "$pathtofolder" == *"My_V4_Miniscope"* ]]; then
        echo "Found folder: $pathtofolder"

        rm -rf /mnt/data/minianAB/* #empty mnt data
        #cp -r "${pathtofolder}" /mnt/data/minianAB/ #copy crnldata to mnt data 
        cp -r "${pathtofolder}/"* /mnt/data/minianAB/ #copy crnldata to mnt data 
        srun --mem=80G --cpus-per-task=10 python /home/aurelie.brecier/HayLabAnalysis/python/MinianCluster.py
        # Check the exit status of srun
        if [ $? -ne 0 ]; then
            echo "Error: srun failed for $pathtofolder. Skipping..."
            continue
        fi
        #cp -r /mnt/data/minianAB/ "${pathtofolder}" #copy mnt data to crnldata 
        cp -r /mnt/data/minianAB/ "${pathtofolder}/"*  #copy crnldata to mnt data 
        rm -rf /mnt/data/minianAB/* #empty mnt data

    fi
done

# Reset IFS to default
IFS=$'\n'
