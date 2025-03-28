#!/bin/bash

# cd in /HayLabAnalysis/bash
#Activate dlc env (ffmpeg installed)
#(dlc) aurelie.brecier@node14:~/HayLabAnalysis/bash$ 

# Define the starting directory
START_DIR="/crnldata/forgetting/Aurelie/CheeseboardExperiment/"
START_DIR="/crnldata/forgetting/Aurelie/CheeseboardExperiment/DAQ_data/AB/Training/"


echo "Searching for folders containing .avi files in '$START_DIR'..." 

# Adjust Internal Field Separator to handle spaces in folder names
IFS=$'\n' 

# Loop through all the folders containing .avi files 
for pathtofolder in $(find "$START_DIR" -type f -name "*.avi" -exec dirname {} \; | sort -u); do
    
    #if [[ "$pathtofolder" == *"My_First_WebCam"* && "$pathtofolder" == *"/Cheeseboard/"* && ! -d "$pathtofolder/plot-poses" ]]; then #only process cheeseboard movies that were not already processed
    if [[ "$pathtofolder" == *"My_First_WebCam"* && "$pathtofolder" == *"/Cheeseboard/"* ]]; then  #only process cheeseboard movies 

        echo "Found folder: $pathtofolder"

        rm -rf /mnt/data/AurelieB_dlc/* #empty mnt data
        cp -r "${pathtofolder}/"* /mnt/data/AurelieB_dlc/ #copy crnldata to mnt data 


        for file in /mnt/data/AurelieB_dlc/*.avi; do
        
            # Define output file name
            output_file="${file%.avi}_compressed.avi"

            echo "Compressing: $file → $output_file"

            # Compress the .avi file using ffmpeg
            
            ffmpeg -i "$file" -vcodec libx264 -crf 23 -preset medium -acodec aac -b:a 128k "$output_file" -loglevel quiet -y > /dev/null 2>&1

            if [[ -s "$output_file" ]]; then # Check if the compressed file is empty = failed compresssion
                rm -f $file
                echo "Done: $output_file"
            else
                rm -f $output_file
                echo "/!\ Compression failed /!\ "
                continue
            fi          
        done
        
        rm -rf ${pathtofolder}/* #empty folder    
        cp -r /mnt/data/AurelieB_dlc/* "${pathtofolder}/" #copy dlc folder of mnt data to crnldata 
        rm -rf /mnt/data/AurelieB_dlc/* #empty mnt data
    fi
done

# Reset IFS to default
IFS=$'\n'



