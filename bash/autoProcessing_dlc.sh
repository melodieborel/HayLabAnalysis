#!/bin/bash

# Activate "dlc" environement before
# cd in /HayLabAnalysis/bash
#(dlc) aurelie.brecier@node14:~/HayLabAnalysis/bash$ ./autoProcessing_dlc.sh

# Define the starting directory
START_DIR="/crnldata/forgetting/Aurelie/MiniscopeOE_data/"

echo "Searching for folders containing .avi files in '$START_DIR'..." 

# Adjust Internal Field Separator to handle spaces in folder names
IFS=$'\n' 

# Loop through all the folders containing .avi files 
for pathtofolder in $(find "$START_DIR" -type f -name "*.avi" -exec dirname {} \; | sort -u); do
    
    #if [[ "$pathtofolder" == *"My_First_WebCam"* && "$pathtofolder" == *"/Cheeseboard/"* && ! -d "$pathtofolder/plot-poses" ]]; then #only process cheeseboard movies that were not already processed
    #if [[ "$pathtofolder" == *"My_First_WebCam"* && "$pathtofolder" == *"/MemoryTest/"* && "$pathtofolder" == *"/Cheeseboard/"* && ! -d "$pathtofolder/plot-poses" ]]; then #only process cheeseboard movies that were not already processed
    if [[ "$pathtofolder" == *"My_First_WebCam"* && "$pathtofolder" == *"/Cheeseboard/"* ]]; then  #only process cheeseboard movies 

        echo "Found folder: $pathtofolder"

        rm -rf /mnt/data/AurelieB_dlc/* #empty mnt data
        cp -r "${pathtofolder}/"* /mnt/data/AurelieB_dlc/ #copy crnldata to mnt data 
        
        srun --mem=90G --cpus-per-task=10 python /home/aurelie.brecier/HayLabAnalysis/python/DLC_2_analyze_videos.py
        #srun --partition=GPU --mem=20G --cpus-per-task=4 --gres=gpu:1g.20gb:1 python /home/aurelie.brecier/HayLabAnalysis/python/DLC_2_analyze_videos.py
        
        # Check the exit status of srun
        if [ $? -ne 0 ]; then
            echo "Error: srun failed for $pathtofolder. Skipping..."
            continue #break
        fi

        cp -r /mnt/data/AurelieB_dlc/* "${pathtofolder}/" #copy dlc folder of mnt data to crnldata 
        rm -rf /mnt/data/AurelieB_dlc/* #empty mnt data
    fi
done

# Reset IFS to default
IFS=$'\n'