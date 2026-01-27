#!/bin/bash

# Activate "minian" environement before
# cd in /HayLabAnalysis/bash
#(minian) aurelie.brecier@node14:~/HayLabAnalysis/bash$ ./autoProcessing_minian.sh

# Define the starting directory
START_DIR="/crnldata/forgetting/Aurelie/MiniscopeOE_data/L1NDNF_mice/OW/Allocentric_task/Training/"

echo "Searching for folders containing .avi files in '$START_DIR'..." 

# Adjust Internal Field Separator to handle spaces in folder names
IFS=$'\n' 

# Loop through all the folders containing .avi files 
for pathtofolder in $(find "$START_DIR" -type f -name "*.avi" -exec dirname {} \; | sort -u); do
    
     if [[ "$pathtofolder" == *"V4_Miniscope"* && ! -d "$pathtofolder/minian" ]]; then #only process miniscope movies that were not already processed
    #if [[ "$pathtofolder" == *"V4_Miniscope"* ]]; then  # process all miniscope movies 

        echo "Found folder: $pathtofolder"
        
        rm -rf /mnt/data/TheaM_minian/* #empty mnt data
 
        # Use '/' as a delimiter to split the path and count the parts
        pathtofolder_len=$(echo "$pathtofolder" | awk -F'/' '{print NF}')

        mouse_name=$(echo "$pathtofolder" | awk -F'/' '{print $7}')
        session_name=$(echo "$pathtofolder" | awk -F'/' '{print $9}')

        #if [ "$pathtofolder_len" -eq 12 ]; then          
        #    session_name=$(echo "$pathtofolder" | awk -F'/' '{print $11}')
        #elif [ "$pathtofolder_len" -eq 11 ]; then
        #    text=$(echo "$pathtofolder" | awk -F'/' '{print $10}')
        #    session_name="${text##*_}"      
        #fi        

        echo "$mouse_name"
        echo "$session_name"

        mkdir -p "/mnt/data/TheaM_minian/$mouse_name/$session_name/"

        cp -r "${pathtofolder}/"* /mnt/data/TheaM_minian/$mouse_name/$session_name/ #copy crnldata to mnt data 
        
        #srun --mem=250G --cpus-per-task=40 python /home/thea.michel/HayLabAnalysis/python/MINI_1_detect_units_remote.py #&>/dev/null
        #srun --mem=80G --cpus-per-task=10 python /home/thea.michel/HayLabAnalysis/python/MINI_1_detect_units_remote.py 
        srun --partition=GPU --mem=20G --cpus-per-task=4 --gres=gpu:1g.20gb:1 python /home/thea.michel/HayLabAnalysis/python/MINI_1_detect_units_remote.py


        # Check the exit status of srun
        if [ $? -ne 0 ]; then
            echo "Error: srun failed for $pathtofolder. Skipping..."
            continue #break
        fi

        cp -r /mnt/data/TheaM_minian/$mouse_name/$session_name/minian "${pathtofolder}/" #copy minian folder of mnt data to crnldata 
        rm -rf /mnt/data/TheaM_minian/* #empty mnt data
    fi
done

# Reset IFS to default
IFS=$'\n'