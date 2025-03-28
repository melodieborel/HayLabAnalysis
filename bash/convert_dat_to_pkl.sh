#!/bin/bash


# Activate "minian" environement before
# cd in /HayLabAnalysis/bash
#(minian) aurelie.brecier@node14:~/HayLabAnalysis/bash$ ./concert_dat_to_pkl.sh


# Define the starting directory
START_DIR="/crnldata/forgetting/Aurelie/CheeseboardExperiment/"
START_DIR="/crnldata/waking/audrey_hay/L1imaging/Analysed2025_AB/"


echo "Searching for folders containing .dat files in '$START_DIR'..." 

# Adjust Internal Field Separator to handle spaces in folder names
IFS=$'\n' 

# Loop through all the folders containing .dat files 
for pathtofolder in $(find "$START_DIR" -type d -name "OpenEphys" | sort -u); do 
 
    if [[ ! -f "$pathtofolder/DataFrame_rawdataDS.pkl" ]]; then #only process if DataFrame_rawdataDS.pkl does not exist already

        echo "Found folder: $pathtofolder"

        rm -rf /mnt/data/AurelieB_other/* #empty mnt data
        cp -r "${pathtofolder}/"* /mnt/data/AurelieB_other/ #copy crnldata to mnt data 
        
        #srun --mem=250G --cpus-per-task=40 python /home/aurelie.brecier/HayLabAnalysis/python/OE_2_downsample_rec.py
        srun --mem=250G --cpus-per-task=40 python /home/aurelie.brecier/HayLabAnalysis/python/OE_2_downsample_rec.py

        # Check the exit status of srun
        if [ $? -ne 0 ]; then
            echo "Error: srun failed for $pathtofolder. Skipping..."
            continue #break
        fi

        cp -r /mnt/data/AurelieB_other/* "${pathtofolder}/" #copy mnt data to crnldata 
        rm -rf /mnt/data/AurelieB_other/* #empty mnt data
    fi
done

# Reset IFS to default
IFS=$'\n'
