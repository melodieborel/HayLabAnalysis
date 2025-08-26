##### Called by autoProcessing_dlc.sh #####

import deeplabcut
import time

st = time.time()

path_config_file='/home/aurelie.brecier/CheeseboardMiniscope-AurelieB-2025-02-21/config.yaml'
videofile_path = ["/mnt/data/AurelieB_dlc/"] #will analyse all .avi files in this directory
#videofile_path = ["//10.69.168.1/crnldata/forgetting/Aurelie/MiniscopeOE_data/"] #will analyse all .avi files in this directory

deeplabcut.analyze_videos(path_config_file,videofile_path, videotype='.avi')#, save_as_csv=True)
#deeplabcut.filterpredictions(path_config_file, videofile_path, shuffle=1, p_bound=0.01, ARdegree=3, MAdegree=1)
#deeplabcut.create_labeled_video(path_config_file,videofile_path)
#deeplabcut.plot_trajectories(path_config_file,videofile_path, imagetype=".svg")

elapsed_time = time.time() - st
print('Execution time:', time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))