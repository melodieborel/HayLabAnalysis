#!/bin/bash

# Define the starting directory
START_DIR="/crnldata/forgetting/Aurelie/CheeseboardExperiment/"
START_DIR="/crnldata/forgetting/Aurelie/CheeseboardExperiment/DAQ_data/AB/Training/"

echo "Searching for folders containing .avi files in '$START_DIR'..." 

# Adjust Internal Field Separator to handle spaces in folder names
IFS=$'\n' 

# Loop through all the folders containing .avi files 
for pathtofolder in $(find "$START_DIR" -type f -name "*.avi" -exec dirname {} \; | sort -u); do
    
    if [[ "$pathtofolder" == *"My_First_WebCam"* && "$pathtofolder" == *"/Cheeseboard/"* && ! -d "$pathtofolder/plot-poses" ]]; then #only process cheeseboard movies that were not already processed
    #if [[ "$pathtofolder" == *"My_First_WebCam"* && "$pathtofolder" == *"Cheeseboard"* ]]; then  #only process cheeseboard movies 

        echo "Found folder: $pathtofolder"

        rm -rf /mnt/data/AurelieB_other/* #empty mnt data
        cp -r "${pathtofolder}/"* /mnt/data/AurelieB_other/ #copy crnldata to mnt data 

        # Define folder containing AVI files
        INPUT_FOLDER="/mnt/data/AurelieB_other/"
        OUTPUT_FILE="$INPUT_FOLDER/output.avi"
        FILE_LIST="file_list.txt"

        # Check for ffmpeg
        if ! command -v ffmpeg &> /dev/null; then
            echo "ffmpeg is not installed. Please install it and try again."
            exit 1
        fi

        # Create a file list
        rm -f "$FILE_LIST"
        for file in $(ls "$INPUT_FOLDER"/*.avi 2>/dev/null | sort -V); do
            echo "file '$file'" >> "$FILE_LIST"
        done

        # Check if the list is not empty
        if [ ! -s "$FILE_LIST" ]; then
            echo "No .avi files found in $INPUT_FOLDER."
            rm -f "$FILE_LIST"
            exit 1
        fi

        # Remove DLC video analysis
        #find "${INPUT_FOLDER}" -type f \( -name "*.h5" -o -name "*.pickle" \) -exec rm -f {} \;
        #find "${INPUT_FOLDER}" -type d -name "plot-poses" -exec rm -rf {} \;

        # Merge using ffmpeg, suppress all output
        ffmpeg -f concat -safe 0 -i "$FILE_LIST" -c copy "$OUTPUT_FILE" > /dev/null 2>&1

        # Cleanup: Remove the .avi files from the input folder after merging, but keep output.avi
        find "$INPUT_FOLDER" -type f -name "*.avi" ! -name "output.avi" -exec rm -f {} \;

        # Remove the temporary file list
        rm -f "$FILE_LIST"

        # Rename avi file 
        mv "$INPUT_FOLDER/output.avi" "$INPUT_FOLDER/0_compressed.avi"
        
        echo "Merging complete. Output file: 0.avi"

        rm -rf ${pathtofolder}/* #empty folder    
        cp -r /mnt/data/AurelieB_other/* "${pathtofolder}/" #copy dlc folder of mnt data to crnldata 
        rm -rf /mnt/data/AurelieB_other/* #empty mnt data
    fi
done

# Reset IFS to default
IFS=$'\n'