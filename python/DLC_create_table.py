
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import random
import math
import os
from pathlib import Path
import re
from scipy import stats

folder_path = "//10.69.168.1/crnldata/forgetting/Aurelie/CheeseboardExperiment/DAQ_data/AB" #will analyse all .avi files in this directory
counter=0
day=1
trial=1
previousmice=0
previous_session_time=0
Summary_table=pd.DataFrame()

for root, dirs, files in os.walk(folder_path):
    h5_files = [f for f in files if f.endswith(('.h5'))]
    
    for h5_file in h5_files:
        filename = os.path.join(root, h5_file)
        filename = os.path.normpath(filename)

        # Load HDF5 file
        df = pd.read_hdf(filename)
        directory = os.path.dirname(filename)
        timestamps_path = Path(directory,'timeStamps.csv')
        if timestamps_path.exists():
            timestamps = pd.read_csv(timestamps_path)
            frame_rate = round(1/(np.mean(np.diff(timestamps.iloc[:,1]))/1000))  # fps
        else:
            frame_rate = 16  # fps /!\ CHANGE ACCORDING TO YOUR DATA
        with open("//10.69.168.1/crnldata/forgetting/Aurelie/CheeseboardExperiment/DAQ_data/AB/Reward_position.txt", "r") as file:
            text = file.read()  
        numbers = re.findall(r"[-+]?\d*\.\d+|\d+", text)
        reward_x, reward_y = map(float, numbers)

        session_type = os.path.basename(os.path.abspath(os.path.join(filename, "../../../../../..")))
        mice = os.path.basename(os.path.abspath(os.path.join(filename, "../../../../..")))
        sesstion_time = os.path.basename(os.path.abspath(os.path.join(filename, "../../..")))
        trial_time = os.path.basename(os.path.abspath(os.path.join(filename, "../..")))

        if mice == previousmice:
            if sesstion_time==previous_session_time:
                trial+=1
            else:
                day+=1
                trial=1
        else:
            day=1
            trial=1

        # Define parameters
        pixel_to_cm = 2.25  
        reward_zone = 8 * pixel_to_cm if "Training" in filename else 20 * pixel_to_cm  # 8 cm for training, 20 cm for test
        table_center_x, table_center_y = 313, 283  # Center of the cheeseboard table on the video
        table_radius = 270 / 2
        min_stay_at_reward = 5 * frame_rate  # 5 seconds

        # Define functions
        def calculate_relative_distance(x1, y1, x2, y2):
            return math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)

        def calculate_distance_run(x_coords, y_coords):
            distances = np.sqrt(np.diff(x_coords) ** 2 + np.diff(y_coords) ** 2)
            for i in range(1, len(distances) - 1):
                if np.isnan(distances[i]):
                    neighbors = [distances[i-1], distances[i+1]]
                    distances[i] = np.nanmean([x for x in neighbors if not np.isnan(x)])
            total_distance_cm = np.nansum(distances) / pixel_to_cm  # Convert to cm
            return total_distance_cm, distances

        def remove_abnormal_speed(x, y, max_speed):
            distances = np.sqrt(np.diff(x) ** 2 + np.diff(y) ** 2)  
            for i in range(1, len(distances) - 1):
                if np.isnan(distances[i]):
                    neighbors = [distances[i-1], distances[i+1]]
                    distances[i] = np.nanmean([x for x in neighbors if not np.isnan(x)])
            mask = np.insert(distances <= max_speed, 0, True)  # Keep first point
            return x[mask], y[mask]

        # Remove uncertain location predictions (likelihood < 0.9)
        df.iloc[:, 0] = df.apply(lambda row: row.iloc[0] if row.iloc[-1] > 0.9 else np.nan, axis=1)
        df.iloc[:, 1] = df.apply(lambda row: row.iloc[1] if row.iloc[-1] > 0.9 else np.nan, axis=1)

        # Separate the individual's positions into x and y coordinates
        X = df.iloc[:, 0]
        individual_x = np.array(X.values)
        Y = df.iloc[:, 1]
        individual_y = np.array(Y.values)

        # Define when the mouse is on the cheeseboard (start)
        for i, x in enumerate(individual_x):
            y = individual_y[i]
            if calculate_relative_distance(x, y, table_center_x, table_center_y) >= table_radius:
                individual_x[i] = np.nan
                individual_y[i] = np.nan

        x_start = individual_x[~np.isnan(individual_x)][0] if np.any(~np.isnan(individual_x)) else None # first non nan value
        y_start = individual_y[~np.isnan(individual_y)][0] if np.any(~np.isnan(individual_y)) else None # first non nan value

        start_frame = np.where(individual_x == x_start)[0].item()
        for i in range(len(individual_x)-1, 0, -1): # Find the last non-NaN value which is not isolated
            if not np.isnan(individual_x[i]) and not np.isnan(individual_x[i-1]):
                last_frame = i
                break
            
        if timestamps_path.exists():
            start_time = timestamps.iloc[start_frame,1].item() / 1000
            end_time = timestamps.iloc[-1,1].item() / 1000
            duration_trial = end_time - start_time
        else:
            duration_trial = (last_frame - start_frame) / frame_rate

        total_distance, distances = calculate_distance_run(individual_x[start_frame:last_frame], individual_y[start_frame:last_frame])

        # Define when the mouse eats the reward
        found_reward_frame = np.nan
        consecutive_count = 0
        for i, (x, y) in enumerate(zip(individual_x, individual_y)):
            if calculate_relative_distance(x, y, reward_x, reward_y) <= reward_zone:
                consecutive_count += 1
            else:
                consecutive_count = 0  
            if consecutive_count > min_stay_at_reward:
                found_reward_frame = (i - min_stay_at_reward)
                break
            
        if not np.isnan(found_reward_frame):
            if timestamps_path.exists():
                found_reward_time = timestamps.iloc[found_reward_frame,1].item() / 1000
                latency = (found_reward_time - start_time)
            else:
                latency = (found_reward_frame - start_frame) / frame_rate
            distance_to_reward, distances_to_reward = calculate_distance_run(individual_x[start_frame:found_reward_frame], individual_y[start_frame:found_reward_frame])
            enter_reward_zone = 0
            consecutive_count = 0
            for i, (x, y) in enumerate(zip(individual_x, individual_y)):
                if calculate_relative_distance(x, y, reward_x, reward_y) <= reward_zone:
                    if consecutive_count == 0:
                        enter_reward_zone += 1
                        consecutive_count = 1
                else:
                    consecutive_count = 0   
            crossings_per_m=enter_reward_zone / (round(distance_to_reward) / 100)
        
        else:
            latency=np.nan
            distance_to_reward=np.nan
            distances_to_reward=np.nan
            enter_reward_zone=np.nan
            crossings_per_m=np.nan

        # Define the time spent inside the reward zone
        time_spent_in_zone = 0

        for i, (x, y) in enumerate(zip(individual_x, individual_y)):
            if calculate_relative_distance(x, y, reward_x, reward_y) <= reward_zone:
                time_spent_in_zone += 1
        time_spent_in_zone = time_spent_in_zone / frame_rate

        Summary_table.loc[counter, 'session_type'] = session_type
        Summary_table.loc[counter, 'mice'] = mice
        Summary_table.loc[counter, 'session'] = int(day)
        Summary_table.loc[counter, 'session_time'] = sesstion_time
        Summary_table.loc[counter, 'trial'] = int(trial)
        Summary_table.loc[counter, 'trial_time'] = trial_time
        Summary_table.loc[counter, 'reward_location_pix'] = str([reward_x,reward_y])
        Summary_table.loc[counter, 'duration_trial_s'] = round(duration_trial,2)
        Summary_table.loc[counter, 'total_distance_cm'] = round(total_distance,2)
        Summary_table.loc[counter, 'latency_to_reward_s'] = round(latency,2)
        Summary_table.loc[counter, 'time_spent_at_reward_s'] = round(time_spent_in_zone,2)
        Summary_table.loc[counter, 'distance_to_reward_cm'] = round(distance_to_reward,2)
        Summary_table.loc[counter, 'crossings'] = enter_reward_zone
        Summary_table.loc[counter, 'crossings_per_m'] = round(crossings_per_m, 2)

        previousmice = mice
        previous_session_time = sesstion_time
        counter+=1

filenameOut = f'{folder_path}/Summary_table.xlsx'
with pd.ExcelWriter(filenameOut) as writer:
    Summary_table.to_excel(writer, index=False)