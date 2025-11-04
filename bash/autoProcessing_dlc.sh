#!/bin/bash

# Activate "dlc" environement before
# cd in /HayLabAnalysis/bash
#(dlc) aurelie.brecier@node14:~/HayLabAnalysis/bash$ ./autoProcessing_dlc.sh

# Define the starting directory
START_DIR="/crnldata/forgetting/Carla/Cheeseboard/"
echo "Searching for folders containing .mp4 files in '$START_DIR'..."

# Adjust Internal Field Separator to handle spaces in folder names
IFS=$'\n'

# Loop through all the folders containing .mp4 files
for pathtofolder in $(find "$START_DIR" -type f -name "*.mp4" -exec dirname {} \; | sort -u); do
    
    # Check if folder contains "Cheeseboard"
    if [[ "$pathtofolder" == *"Cheeseboard"* ]]; then
        
        # Check if any .h5 file already exists in this folder
        h5_exists=false
        for mp4_file in "$pathtofolder"/*.mp4; do
            if [[ -f "$mp4_file" ]]; then
                # Get base name without extension
                video_basename=$(basename "$mp4_file" .mp4)
                h5_file="${pathtofolder}/${video_basename}DLC_Resnet50_CheeseboardFeb6shuffle1_snapshot_010.h5"
                
                if [[ -f "$h5_file" ]]; then
                    h5_exists=true
                    break
                fi
            fi
        done
        
        # Only process if no .h5 file exists
        if [[ "$h5_exists" == false ]]; then
            echo "Found folder to process: $pathtofolder"
            rm -rf /mnt/data/DLC_Carla/* #empty mnt data
            cp -r "${pathtofolder}/"* /mnt/data/DLC_Carla/ #copy crnldata to mnt data
           
            srun --mem=90G --cpus-per-task=20 python /home/carla.burnet-merlin/HayLabAnalysis/python/DLC_2_analyze_videos.py
            
            # Check the exit status of srun
            if [ $? -ne 0 ]; then
                echo "Error: srun failed for $pathtofolder. Skipping..."
                continue
            fi
            
            cp -r /mnt/data/DLC_Carla/* "${pathtofolder}/" #copy dlc folder of mnt data to crnldata
            rm -rf /mnt/data/DLC_Carla/* #empty mnt data
        else
            echo "Skipping $pathtofolder - .h5 files already exist"
        fi
    fi
done

# Reset IFS to default
IFS=$' \t\n'