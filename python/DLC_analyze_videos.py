# Activate "dlc" environement before
# (dlc) aurelie.brecier@node14:~$ srun --mem=90G --cpus-per-task=40 python /home/aurelie.brecier/HayLabAnalysis/python/DLC_analyze_videos.py

import os
import shutil
import subprocess
from pathlib import Path
import deeplabcut
import time
totst = time.time()

# Define start directory
START_DIR = "/crnldata/forgetting/Aurelie/CheeseboardExperiment/DAQ_data/AB/Training/testtest/"

print(f"Searching for folders containing .avi files in '{START_DIR}'...")

# Find unique folders containing .avi files
folders_with_avi = set()
for root, _, files in os.walk(START_DIR):
    if any(file.endswith(".avi") for file in files):
        folders_with_avi.add(root)

# Process only folders that match the criteria
for folder in sorted(folders_with_avi):
    if "My_First_WebCam" in folder and "/Cheeseboard/" in folder:
        print(f"Found folder: {folder}")

        mnt_data_path = "/mnt/data/AurelieB_dlc"

        # Empty /mnt/data/AurelieB_dlc
        shutil.rmtree(mnt_data_path, ignore_errors=True)
        os.makedirs(mnt_data_path, exist_ok=True)

        # Copy contents from found folder to /mnt/data/AurelieB_dlc
        for item in os.listdir(folder):
            src_path = os.path.join(folder, item)
            dest_path = os.path.join(mnt_data_path, item)
            if os.path.isdir(src_path):
                shutil.copytree(src_path, dest_path, dirs_exist_ok=True)
            else:
                shutil.copy2(src_path, dest_path)

        # Process data
        st = time.time()

        #path_config_file='/home/aurelie.brecier/Cheeseboard-Aurelie-2025-02-19/config.yaml'
        path_config_file='/home/aurelie.brecier/CheeseboardMiniscope-AurelieB-2025-02-21/config.yaml'

        #videofile_path = ["/mnt/data/DLC_Aurelie/0.avi"]
        videofile_path = ["/mnt/data/AurelieB_dlc/"] #will analyse all .avi files in this directory

        deeplabcut.analyze_videos(path_config_file,videofile_path, videotype='.avi')#, save_as_csv=True)
        deeplabcut.filterpredictions(path_config_file, videofile_path, shuffle=1, p_bound=0.01, ARdegree=3, MAdegree=1)

        #deeplabcut.create_labeled_video(path_config_file,videofile_path)
        #deeplabcut.plot_trajectories(path_config_file,videofile_path, imagetype=".svg")

        elapsed_time = time.time() - st
        print('Execution time:', time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))


        try:
            result = subprocess.run(["srun", "your_command_here"], check=True)
        except subprocess.CalledProcessError:
            print(f"Error: srun failed for {folder}. Skipping...")
            continue

        # Copy processed data back to original folder
        for item in os.listdir(mnt_data_path):
            src_path = os.path.join(mnt_data_path, item)
            dest_path = os.path.join(folder, item)
            if os.path.isdir(src_path):
                shutil.copytree(src_path, dest_path, dirs_exist_ok=True)
            else:
                shutil.copy2(src_path, dest_path)

        # Empty /mnt/data/AurelieB_dlc again
        shutil.rmtree(mnt_data_path, ignore_errors=True)

totelapsed_time = time.time() - totst
print('Total execution time:', time.strftime("%H:%M:%S", time.gmtime(totelapsed_time)))



