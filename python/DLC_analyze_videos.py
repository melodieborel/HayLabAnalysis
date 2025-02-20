import deeplabcut
import time

st = time.time()

path_config_file='/home/aurelie.brecier/DeepLabCut/Cheeseboard-Aurelie-2025-02-19/config.yaml'

#videofile_path = ["/mnt/data/DLC_Aurelie/0.avi"]
videofile_path = ["/mnt/data/AurelieB_dlc/"] #will analyse all .avi files in this directory

deeplabcut.analyze_videos(path_config_file,videofile_path, videotype='.avi', save_as_csv=True)
deeplabcut.create_labeled_video(path_config_file,videofile_path)
deeplabcut.plot_trajectories(path_config_file,videofile_path)

elapsed_time = time.time() - st
print('Execution time:', time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))