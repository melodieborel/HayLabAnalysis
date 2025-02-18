#!/bin/bash

# Define the starting directory
START_DIR="/crnldata/forgetting/Aurelie/CheeseboardExperiment/"
START_DIR="/crnldata/forgetting/Aurelie/CheeseboardExperiment/DAQ_data/AB/Habituation/Green/"

echo "Searching for folders containing .avi files in '$START_DIR'..." 

# Adjust Internal Field Separator to handle spaces in folder names
IFS=$'\n' 

# Loop through all the folders containing .avi files 
for pathtofolder in $(find "$START_DIR" -type f -name "*.avi" -exec dirname {} \; | sort -u); do
    
    if [[ "$pathtofolder" == *"My_V4_Miniscope"* && ! -d "$pathtofolder/minian" ]]; then #only process miniscope movies that were not already processed
    # if [[ "$pathtofolder" == *"My_V4_Miniscope"* ]]; then  #only process miniscope movies 

        echo "Found folder: $pathtofolder"

        rm -rf /mnt/data/minianAB/* #empty mnt data
        cp -r "${pathtofolder}/"* /mnt/data/minianAB/ #copy crnldata to mnt data 
        
        srun --mem=250G --cpus-per-task=40 python /home/aurelie.brecier/HayLabAnalysis/python/pipelineMinian.py
        
        # Check the exit status of srun
        if [ $? -ne 0 ]; then
            echo "Error: srun failed for $pathtofolder. Skipping..."
            continue #break
        fi

        cp -r /mnt/data/minianAB/minian "${pathtofolder}/" #copy minian folder of mnt data to crnldata 
        rm -rf /mnt/data/minianAB/* #empty mnt data
    fi
done

# Reset IFS to default
IFS=$'\n'
