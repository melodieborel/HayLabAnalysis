#!/bin/bash

# Activate "minian" environement before
# cd in /HayLabAnalysis/bash
#(minian) aurelie.brecier@node14:~/HayLabAnalysis/bash$ ./autoProcessing_minian.sh

# Define the starting directory

#START_DIR="/crnldata/forgetting/Clementine/CheeseboardExperiment/DAQ_data/ClementineR/Training/GL/Cheesboard/2025_03_31/12_22_47/"
START_DIR="/crnldata/forgetting/Clementine/CheeseboardExperiment/DAQ_data/ClementineR/Training/YL/SleepAfter/2025_03_24/13_58_54/"

echo "Searching for folders containing .avi files in '$START_DIR'..." 

# Adjust Internal Field Separator to handle spaces in folder names
IFS=$'\n' 

# Loop through all the folders containing .avi files 
for pathtofolder in $(find "$START_DIR" -type f -name "*.avi" -exec dirname {} \; | sort -u); do
    
    if [[ "$pathtofolder" == *"V4_Miniscope"* && ! -d "$pathtofolder/minian" ]]; then #only process miniscope movies that were not already processed
    # if [[ "$pathtofolder" == *"V4_Miniscope"* ]]; then  # process all miniscope movies 

        echo "Found folder: $pathtofolder"
        
        rm -rf /mnt/data/ClemR_minian/* #empty mnt data
 
        # Use '/' as a delimiter to split the path and count the parts
        pathtofolder_len=$(echo "$pathtofolder" | awk -F'/' '{print NF}')

        mouse_name=$(echo "$pathtofolder" | awk -F'/' '{print $9}')
        session_name=$(echo "$pathtofolder" | awk -F'/' '{print $12}')

        echo "$mouse_name"
        echo "$session_name"

        mkdir -p "/mnt/data/ClemR_minian/$mouse_name/$session_name/"

        cp -r "${pathtofolder}/"* /mnt/data/ClemR_minian/$mouse_name/$session_name/ #copy crnldata to mnt data 
        
        #srun --mem=250G --cpus-per-task=40 python /home/aurelie.brecier/HayLabAnalysis/python/MINI_1_detect_units_pipeline.py
        srun --mem=90G --cpus-per-task=40 python /home/clementine.robein/HayLabAnalysis/python/MINI_1_detect_units_pipeline.py #&>/dev/null

        # Check the exit status of srun
        if [ $? -ne 0 ]; then
            echo "Error: srun failed for $pathtofolder. Skipping..."
            continue #break
        fi

        cp -r /mnt/data/ClemR_minian/$mouse_name/$session_name/minian "${pathtofolder}/" #copy minian folder of mnt data to crnldata 
        rm -rf /mnt/data/ClemR_minian/* #empty mnt data
    fi
done

# Reset IFS to default
IFS=$'\n'