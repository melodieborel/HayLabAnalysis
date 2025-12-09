#!/bin/bash

# Define the starting directory
START_DIR="//10.69.168.1/crnldata/forgetting/Aurelie/MiniscopeOE_data/L2_3_mice/PC/Exploration_task/"

echo "Searching for folders containing .avi files in '$START_DIR'..." 

# Adjust Internal Field Separator to handle spaces in folder names
IFS=$'\n' 


# Loop through all the folders containing .avi files 
for pathtofolder in $(find "$START_DIR" -type f -name "*.avi" -exec dirname {} \; | sort -u); do
    
    #if [[ "$pathtofolder" == *"My_First_WebCam"* && "$pathtofolder" == *"/Cheeseboard/"* && ! -d "$pathtofolder/plot-poses" ]]; then #only process cheeseboard movies that were not already processed
    #if [[ "$pathtofolder" == *"My_First_WebCam"* && "$pathtofolder" == *"Cheeseboard"* ]]; then  #only process cheeseboard movies 
    if [[ "$pathtofolder" == *"My_First_WebCam"* ]]; then  #only process cheeseboard movies 
        count=$(find "$pathtofolder" -maxdepth 1 -type f -name "*.avi" | wc -l)
        if [ "$count" -gt 1 ]; then
            echo "Folder $pathtofolder has $count .avi files"

            echo "Found folder: $pathtofolder"

            # Define folder containing AVI files
            INPUT_FOLDER="$pathtofolder"
            OUTPUT_FILE="$pathtofolder/0_compressed.avi"
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

            # Merge using ffmpeg, suppress all prompt output
            ffmpeg -f concat -safe 0 -i "$FILE_LIST" -c copy "$OUTPUT_FILE" > /dev/null 2>&1

            # Cleanup: Remove the .avi files from the input folder after merging, but keep 0_compressed.avi
            find "$INPUT_FOLDER" -type f -name "*.avi" ! -name "0_compressed.avi" -exec rm -f {} \;

            # Remove the temporary file list
            rm -f "$FILE_LIST"

            echo "Merging complete. Output file: 0_compressed.avi"

        fi
    fi
done

# Reset IFS to default
IFS=$'\n'