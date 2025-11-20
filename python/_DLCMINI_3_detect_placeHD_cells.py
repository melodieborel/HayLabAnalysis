local = False
all_expe_types =['Cheeseboard']

if local:
    dir = "//10.69.168.1/crnldata/forgetting/Aurelie/MiniscopeOE_data/L2_3_mice/"
else: 
    dir = "/mnt/data/AurelieB_other/"

import os
import re
import numpy as np
from scipy import signal
import math 
import json
from collections import Counter
from pathlib import Path
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, Cursor
import pickle
import sys 
from datetime import datetime
import shutil
from ast import literal_eval
from scipy.signal import find_peaks
from scipy.stats import zscore
from scipy import stats
from itertools import groupby
from IPython.display import display
from scipy.interpolate import griddata
import numpy as np
from scipy.signal import butter, filtfilt, hilbert
import matplotlib.pyplot as plt
from scipy import interpolate
import numpy.matlib
from sklearn.decomposition import PCA
from scipy import stats
import numpy as np
import pickle
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import resample
from scipy.signal import resample_poly
from math import gcd
import warnings
warnings.filterwarnings("ignore")
import sys
from scipy.interpolate import interp1d
from collections import defaultdict
import bisect
from scipy.ndimage import gaussian_filter
import random
from scipy.ndimage import gaussian_filter1d
from matplotlib.collections import LineCollection
from itertools import islice
from concurrent.futures import ProcessPoolExecutor
import multiprocessing


class Tee:
    def __init__(self, *files):
        self.files = files
    def write(self, obj):
        for f in self.files:
            f.write(obj)
            f.flush()
    def flush(self):
        for f in self.files:
            f.flush()

if local: 
    os.chdir("C:/Users/Manip2/SCRIPTS/minian/")
else: 
    sys.path.append("/home/aurelie.brecier/minian/")

    
minian_path = os.path.join(os.path.abspath('.'),'minian')
print("The folder used for minian procedures is : {}".format(minian_path))
sys.path.append(minian_path)

from minian.utilities import (
    TaskAnnotation,
    get_optimal_chk,
    load_videos,
    open_minian,
    save_minian,
)


pixel_to_cm = 2.25  
table_center_x, table_center_y = 313, 283  # Center of the cheeseboard table on the video
table_center_x, table_center_y = 300, 270  # Center of the cheeseboard table on the video
table_center_x, table_center_y = 315, 275  # Center of the cheeseboard table on the video

table_radius = 290 / 2
square_size = pixel_to_cm * 6

bin_size = 10
bin_edges = np.arange(0, 361, bin_size)
bin_centers = np.radians(bin_edges[:-1] + bin_size/2)


def find_closest_index_sorted(arr, target):
    idx = bisect.bisect_left(arr, target)  # Find the insertion point
    if idx == 0:
        return 0
    if idx == len(arr):
        return len(arr) - 1
    before = idx - 1
    after = idx
    return before if abs(arr[before] - target) <= abs(arr[after] - target) else after

def calculate_relative_distance(x1, y1, x2, y2):
    return math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)

def find_long_non_nan_sequences(arr, min_length=100):
    mask = ~np.isnan(arr)  # True for non-NaN values
    diff = np.diff(np.concatenate(([0], mask.astype(int), [0])))  # Add padding to detect edges
    starts = np.where(diff == 1)[0]  # Where a sequence starts
    ends = np.where(diff == -1)[0]   # Where a sequence ends
    sequences = [arr[start:end] for start, end in zip(starts, ends) if (end - start) > min_length]
    return sequences

def replace_high_speed_points_with_nan(x, y, speed_threshold):
    x = np.array(x, dtype='float')
    y = np.array(y, dtype='float')
    # Compute speed between consecutive points
    dx = np.diff(x)
    dy = np.diff(y)
    speeds = np.sqrt(dx**2 + dy**2)
    # Create mask for speed exceeding threshold
    high_speed_mask = speeds > speed_threshold
    # We mark i+1 as NaN if speed between them is too high
    x_out = x.copy()
    y_out = y.copy()
    for i in range(len(high_speed_mask)):
        if high_speed_mask[i]:
            # Only mark the faster of the two points
            if i > 0 and i < len(x) - 1:
                if speeds[i] > speeds[i - 1]:
                    x_out[i + 1] = np.nan
                    y_out[i + 1] = np.nan
                else:
                    x_out[i] = np.nan
                    y_out[i] = np.nan
    return x_out, y_out

def interpolate_2d_path(x, y, kind='linear', fill='extrapolate'):
    x = np.array(x, dtype='float')
    y = np.array(y, dtype='float')
    indices = np.arange(len(x))
    valid_mask = ~np.isnan(x) & ~np.isnan(y)
    if np.sum(valid_mask) < 2:
        raise ValueError("Not enough valid points to interpolate/extrapolate.")
    interp_x = interp1d(indices[valid_mask], x[valid_mask], kind=kind, fill_value=fill, bounds_error=False)
    interp_y = interp1d(indices[valid_mask], y[valid_mask], kind=kind, fill_value=fill, bounds_error=False)
    x_filled = x.copy()
    y_filled = y.copy()
    nan_mask = np.isnan(x) | np.isnan(y)
    x_filled[nan_mask] = interp_x(indices[nan_mask])
    y_filled[nan_mask] = interp_y(indices[nan_mask])
    return x_filled, y_filled

def limit_speed(x, y, max_speed):
    dx = np.diff(x.copy())
    dy = np.diff(y.copy())
    speeds = np.sqrt(dx**2 + dy**2)
    for i,t in enumerate(speeds):
        if t > max_speed:        
            x[i+1] = x[i] 
            y[i+1] = y[i] 
            x[i+2] = x[i] 
            y[i+2] = y[i] 
    return x, y

def remove_short_sequences(arr, max_len=10):
    arr = np.array(arr, dtype='float')
    result = arr.copy()
    is_value = ~np.isnan(arr)
    i = 0
    while i < len(arr):
        if is_value[i]:
            start = i
            while i < len(arr) and is_value[i]:
                i += 1
            end = i
            seq_len = end - start
            # Check if surrounded by NaNs and short enough
            if seq_len <= max_len:
                left_nan = (start == 0) or np.isnan(arr[start - 1])
                right_nan = (end == len(arr)) or np.isnan(arr[end])  # safe for edge
                if left_nan and right_nan:
                    result[start:end] = np.nan
        else:
            i += 1
    return result

# Conversion des quaternions en angles d'Euler
def quaternion_to_euler(qw, qx, qy, qz):
    sinr_cosp = 2 * (qw * qx + qy * qz)
    cosr_cosp = 1 - 2 * (qx * qx + qy * qy)
    roll = np.arctan2(sinr_cosp, cosr_cosp)

    sinp = 2 * (qw * qy - qz * qx)
    sinp = np.clip(sinp, -1.0, 1.0)
    pitch = np.arcsin(sinp)

    siny_cosp = 2 * (qw * qz + qx * qy)
    cosy_cosp = 1 - 2 * (qy * qy + qz * qz)
    yaw = np.arctan2(siny_cosp, cosy_cosp)

    return roll, pitch, yaw

def direction(x1, y1, x2, y2):
    dx = x2 - x1
    dy = y2 - y1
    if dx == 0 and dy == 0:
        return np.nan  # undefined direction
    angle_rad = math.atan2(dy, dx)
    angle_deg = math.degrees(angle_rad)
    return (angle_deg + 360) % 360


def runpermutations(args):
    worker_id, first_perm, last_perm = args
    print(worker_id, 'will run permuatations', first_perm, 'til', last_perm)

    bin_size = 10
    bin_edges = np.arange(0, 361, bin_size)
    bin_centers = np.radians(bin_edges[:-1] + bin_size/2)
    sigma = 1  # Standard deviation for Gaussian kernel
    max_row = int(table_radius*2/square_size)+1
    max_col = int(table_radius*2/square_size)+1

    ShuffledSpatial_info={}
    ShuffledHeadDirection_info={}
    ShuffledRoll_info={}
    ShuffledPitch_info={}
    ShuffledYaw_info={}
    ShuffledSpeed_info={}
    ShuffledAHV_info={}

    with open("/mnt/data/AurelieB_other/SpatialInfo.pkl", "rb") as f:
        SpatialInfo = pickle.load(f)

    for nr in SpatialInfo.keys():
        ShuffledSpatial_info[nr]={}
        ShuffledHeadDirection_info[nr]={}
        ShuffledRoll_info[nr]={}
        ShuffledPitch_info[nr]={}
        ShuffledYaw_info[nr]={}
        ShuffledSpeed_info[nr]={}
        ShuffledAHV_info[nr]={}

    for permutation_nb in range(first_perm, last_perm):

        AllCellDict_PC_Sh={}
        AllTimeDict_PC_Sh={}

        AllCellDict_HD_Sh={}
        AllTimeDict_HD_Sh={}

        AllCellDict_Roll_Sh={}
        AllTimeDict_Roll_Sh={}

        AllCellDict_Pitch_Sh={}
        AllTimeDict_Pitch_Sh={}

        AllCellDict_Yaw_Sh={}
        AllTimeDict_Yaw_Sh={}

        AllCellDict_Speed_Sh={}
        AllTimeDict_Speed_Sh={}

        AllCellDict_AHV_Sh={}
        AllTimeDict_AHV_Sh={}

        for dpath in Path(dir).glob('**/Exploration_task/mappingsAB.pkl'):

            mappfile = open(dpath.parents[0]/ f'mappingsAB.pkl', 'rb')
            mapping = pickle.load(mappfile)
            mapping_sess = mapping['session']   
            mice = dpath.parents[1].parts[-1]
            
            sess=0

            minian_folders = [f for f in dpath.parents[0].rglob('minian') if f.is_dir()]

            for minianpath in minian_folders: # for each minian folders found in this mouse 

                if any(p in all_expe_types for p in minianpath.parts): # have to be to the expe_types

                    try: 
                        session_path=minianpath.parents[1]        
                        tsmini= pd.read_csv(list(session_path.glob('*V4_Miniscope/timeStamps.csv'))[0])['Time Stamp (ms)']
                        V4subfolder=False
                        session_type=minianpath.parents[2].name
                        session_date=minianpath.parents[3].name
                        session_time=minianpath.parents[1].name
                        session=session_time
                
                    except:
                        session_path=minianpath.parents[2]      
                        tsmini= pd.read_csv(list(session_path.glob('*V4_Miniscope/timeStamps.csv'))[0])['Time Stamp (ms)']                      
                        V4subfolder=True
                        session_type=minianpath.parents[4].name
                        session_date=minianpath.parents[5].name
                        session_time=minianpath.parents[0].name
                        session=session_time

                    # Minian data 
                                
                    minian_ds = open_minian(minianpath)
                    Co = minian_ds['C']  # calcium traces
                    minian_freq=round(1/np.mean(np.diff(np.array(tsmini)/1000)))
                    try: 
                        TodropFile = minianpath / f'TodropFileAB.json'
                        with open(TodropFile, 'r') as f:
                            unit_to_drop = json.load(f)
                    except:
                        TodropFile = minianpath.parent / f'TodropFileAB.json'
                        with open(TodropFile, 'r') as f:
                            unit_to_drop = json.load(f)

                    # Extraction des quaternions

                    vestibular_df= pd.read_csv(list(session_path.glob('*V4_Miniscope/headOrientation.csv'))[0])
                    tsvest = vestibular_df['Time Stamp (ms)']
                    qw = vestibular_df['qw'].to_numpy()
                    qx = vestibular_df['qx'].to_numpy()
                    qy = vestibular_df['qy'].to_numpy()
                    qz = vestibular_df['qz'].to_numpy()
                    roll, pitch, yaw = quaternion_to_euler(qw, qx, qy, qz)
                    roll_deg = np.degrees(roll) # nose goes clockwise or anticlockwise
                    pitch_deg = np.degrees(pitch) # nose up or down
                    yaw_deg = np.degrees(yaw) # nose moves from side to side   
                    last_valid = tsvest.last_valid_index() + 1
                    roll_deg = roll_deg[:last_valid]
                    pitch_deg = pitch_deg[:last_valid]
                    yaw_deg = yaw_deg[:last_valid]

                    if V4subfolder:
                        V4subfolder_id = int(minianpath.parent.name[-1]) - 1 
                        ts_start = V4subfolder_id*15*1000 # cause 15 videos per subfolders and 1000 frames per videos
                        ts_stop = np.shape(Co)[1] + ts_start 
                        tsmini_sub=tsmini[ts_start:ts_stop]
                        tsmini_sub=tsmini_sub.reset_index(drop=True)                  
                        tsvest_sub = tsvest[ts_start:ts_stop]
                        tsvest_sub = tsvest_sub.reset_index(drop=True)    
                        roll_deg= roll_deg[ts_start:ts_stop]
                        pitch_deg= pitch_deg[ts_start:ts_stop]
                        yaw_deg= yaw_deg[ts_start:ts_stop]            
                    else:
                        tsmini_sub=tsmini  
                        tsvest_sub=tsvest  
                    

                    # DeepLabCut data

                    dlcfile=''
                    dlcpath=Path(f'{Path(session_path)}/My_First_WebCam/')
                    for file in os.listdir(dlcpath):
                        if file.endswith(('.h5')):
                            dlcfile=file
                            break
                    dlc_path = os.path.join(dlcpath, dlcfile)
                    df = pd.read_hdf(dlc_path)
                    directory = os.path.dirname(dlc_path)
                    timestamps_path = Path(directory,'timeStamps.csv')
                    if timestamps_path.exists():
                        timestamps = pd.read_csv(timestamps_path)
                        tswebcam = np.array(timestamps['Time Stamp (ms)'])
                        frame_rate = round(1/(np.mean(np.diff(timestamps.iloc[:,1]))/1000))  # fps
                    else:
                        frame_rate = 16  # fps /!\ CHANGE ACCORDING TO YOUR DATA
                    
                    # Nose 

                    df_nose= df.copy()
                    df_nose.iloc[:, 0] = df_nose.apply(lambda row: row.iloc[0] if row.iloc[2] > 0.5 else np.nan, axis=1)
                    df_nose.iloc[:, 1] = df_nose.apply(lambda row: row.iloc[1] if row.iloc[2] > 0.5 else np.nan, axis=1)
                    X = df_nose.iloc[:, 0]
                    Y = df_nose.iloc[:, 1]        
                    nose_xO= np.array(X.values) 
                    nose_yO = np.array(Y.values)
                    for i, x in enumerate(nose_xO):# Define when the mouse is on the cheeseboard (start)
                        y = nose_yO[i]
                        if calculate_relative_distance(x, y, table_center_x, table_center_y) >= table_radius:
                            nose_xO[i] = np.nan
                            nose_yO[i] = np.nan
                    nose_xOO = remove_short_sequences(nose_xO, max_len=3)
                    nose_yOO = remove_short_sequences(nose_yO, max_len=3)
                    x_start = find_long_non_nan_sequences(nose_xOO)[0][0] # first value of the first long non nan sequence
                    y_start = find_long_non_nan_sequences(nose_yOO)[0][0] # first value of the first long non nan sequence
                    start_frame = np.where(nose_xOO == x_start)[0][0].item()
                    
                    nose_xOO[:start_frame]=np.nan # remove any path before the real start
                    nose_yOO[:start_frame]=np.nan # remove any path before the real start
                    nose_x1, nose_y1 = replace_high_speed_points_with_nan(nose_xOO, nose_yOO, speed_threshold=10)
                    last_frame = len(nose_x1)
                    nose_x2, nose_y2 = interpolate_2d_path(nose_x1[start_frame:last_frame], nose_y1[start_frame:last_frame], kind='nearest')
                    max_speed = 20
                    nose_x3, nose_y3 = limit_speed(nose_x2, nose_y2, max_speed=max_speed)
                    nose_x = np.concatenate((nose_x1[:start_frame], nose_x3))
                    nose_y = np.concatenate((nose_y1[:start_frame], nose_y3))
                    
                    # Neck 

                    df_neck= df.copy()
                    df_neck.iloc[:, 3] = df_neck.apply(lambda row: row.iloc[3] if row.iloc[5] > 0.5 else np.nan, axis=1)
                    df_neck.iloc[:, 4] = df_neck.apply(lambda row: row.iloc[4] if row.iloc[5] > 0.5 else np.nan, axis=1)
                    Xt = df_neck.iloc[:, 3]
                    Yt = df_neck.iloc[:, 4]        
                    neck_xO= np.array(Xt.values) 
                    neck_yO = np.array(Yt.values)
                    for i, x in enumerate(neck_xO):# Define when the mouse is on the cheeseboard (start)
                        y = neck_yO[i]
                        if calculate_relative_distance(x, y, table_center_x, table_center_y) >= table_radius:
                            neck_xO[i] = np.nan
                            neck_yO[i] = np.nan
                    neck_xOO = remove_short_sequences(neck_xO, max_len=3)
                    neck_yOO = remove_short_sequences(neck_yO, max_len=3)
                    neck_xOO[:start_frame]=np.nan # remove any path before the real start
                    neck_yOO[:start_frame]=np.nan # remove any path before the real start
                    neck_x1, neck_y1 = replace_high_speed_points_with_nan(neck_xOO, neck_yOO, speed_threshold=10)
                    neck_x2, neck_y2 = interpolate_2d_path(neck_x1[start_frame:last_frame], neck_y1[start_frame:last_frame], kind='nearest')
                    neck_x3, neck_y3 = limit_speed(neck_x2, neck_y2, max_speed=max_speed)
                    neck_x = np.concatenate((neck_x1[:start_frame], neck_x3))
                    neck_y = np.concatenate((neck_y1[:start_frame], neck_y3))
                                

                    # Keep only crossregistered cells

                    C_sel=Co.drop_sel(unit_id=unit_to_drop)
                    indexMappList=mapping_sess[session]
                    kept_uniq_unit_List=[]
                    for unit in C_sel['unit_id'].values:
                        indexMapp = np.where(indexMappList == unit)[0]
                        kept_uniq_unit_List.append(str(indexMapp))
                    nb_unit=len(C_sel['unit_id'].values)
                    if nb_unit == 0:
                        continue


                    # If the miniscope recording started after the webcam recording, cut the webcam data
                    if tsmini_sub.iloc[0] > tswebcam[start_frame]:
                        Newstart_frame = np.where(tswebcam >= tsmini_sub.iloc[0].item())[0][1].item()
                        start_frame = Newstart_frame 
                    # If the miniscope recording is shorter than the webcam recording, cut the webcam data
                    if tsmini_sub.iloc[-1] < tswebcam[last_frame-1]:
                        Newlast_frame = np.where(tswebcam <= tsmini_sub.iloc[-1].item())[0][-1].item()
                        last_frame = Newlast_frame

                    # Align data

                    nose_x_ = nose_x[start_frame:last_frame]
                    nose_y_ = nose_y[start_frame:last_frame]
                    neck_x_ = neck_x[start_frame:last_frame]
                    neck_y_ = neck_y[start_frame:last_frame]            
                    angles_deg = np.array([direction(x1, y1, x2, y2) for x1, y1, x2, y2 in zip(neck_x_, neck_y_, nose_x_, nose_y_)])

                    ahv_deg= np.diff(angles_deg)
                    ahv_deg_= np.concatenate(([0], ahv_deg))

                    dx= np.diff(nose_x_)
                    dy= np.diff(nose_y_)
                    speed = np.sqrt(dx**2 + dy**2)
                    speed_= np.concatenate(([0], speed))
                    speed_= np.where(speed_>max_speed, 0, speed_)

                    # Keep only crossregistered cells

                    Carray = C_sel.to_numpy().T
                    
                    # Randomly shift the calcium trace to break temporal correlation with behavior
                    offset = random.randint(20*minian_freq, len(Carray)) # Random shift between 20s and the length of the session
                    Carray = np.roll(Carray, shift=offset, axis=0) 

                    Calcium = pd.DataFrame(Carray.T, index=C_sel['unit_id'].values.tolist()) 

                    n = int(np.floor(2 * table_radius / square_size))

                    for unit in range(nb_unit): 
                        indexMapp = np.where(mapping_sess[session] == Calcium.index[unit])[0]
                        if len(indexMapp)>0 : # The neuron needs to be in the cross-registration
                            unique_unit = f'{mice}{str(indexMapp[0])}'
                            sess_unit = f'{mice}{str(indexMapp[0])}_s{sess}_{unit}'
                            Carray_unit = Carray[:,unit]

                            nr_tot_act_biaised_PC = defaultdict(int) # Neuron activity per square
                            counts_PC = defaultdict(int) # Count visits per square

                            nr_tot_act_biaised_HD = np.zeros(len(bin_centers))
                            counts_HD = np.zeros(len(bin_centers)) 

                            nr_tot_act_biaised_Roll = np.zeros(len(bin_centers))
                            counts_Roll = np.zeros(len(bin_centers))
                            
                            nr_tot_act_biaised_Pitch = np.zeros(len(bin_centers))
                            counts_Pitch = np.zeros(len(bin_centers))           

                            nr_tot_act_biaised_Yaw = np.zeros(len(bin_centers))
                            counts_Yaw = np.zeros(len(bin_centers))
                            
                            nr_tot_act_biaised_Speed = np.zeros(max_speed)
                            counts_Speed = np.zeros(max_speed)

                            nr_tot_act_biaised_AHV = np.zeros(len(bin_centers)*2)
                            counts_AHV = np.zeros(len(bin_centers)*2)

                            for idx, (px, py, angle, _speed, _ahv_deg) in enumerate(zip(nose_x_, nose_y_, angles_deg, speed_, ahv_deg_)):
                                if np.sqrt((px - table_center_x)**2 + (py - table_center_y)**2) > table_radius:
                                    continue  # skip points outside circle
                                closest_point = find_closest_index_sorted(tsmini_sub, tswebcam[start_frame+idx])

                                ix = int(np.floor((px - (table_center_x - n/2 * square_size)) / square_size))
                                iy = int(np.floor((py - (table_center_y - n/2 * square_size)) / square_size))
                                nr_tot_act_biaised_PC[(ix, iy)] += Carray_unit[closest_point]
                                counts_PC[(ix, iy)] += 1/frame_rate 

                                bin_idx = int(angle // bin_size)
                                nr_tot_act_biaised_HD[bin_idx] += Carray_unit[closest_point]
                                counts_HD[bin_idx] += 1/frame_rate
                                                            
                                AHVbin_idx = int(_ahv_deg // bin_size*2)
                                nr_tot_act_biaised_AHV[AHVbin_idx] += Carray_unit[closest_point]
                                counts_AHV[AHVbin_idx] += 1/frame_rate

                                speedbin_idx = min(int(_speed), max_speed-1) # between 0 and 20 
                                nr_tot_act_biaised_Speed[speedbin_idx] += Carray_unit[closest_point]
                                counts_Speed[speedbin_idx] += 1/frame_rate  

                            if len(roll_deg) != len(Carray_unit):
                                print('ISSSUE with length Calcium array & vestibular infos')

                            for idx, (_roll, _pitch, _yaw) in enumerate(zip(roll_deg, pitch_deg, yaw_deg)):

                                rollbin_idx = int(_roll // bin_size)
                                nr_tot_act_biaised_Roll[rollbin_idx] += Carray_unit[idx]
                                counts_Roll[rollbin_idx] += 1/minian_freq      

                                pitchbin_idx = int(_pitch // bin_size)
                                nr_tot_act_biaised_Pitch[pitchbin_idx] += Carray_unit[idx]
                                counts_Pitch[pitchbin_idx] += 1/minian_freq

                                yawbin_idx = int(_yaw // bin_size)
                                nr_tot_act_biaised_Yaw[yawbin_idx] += Carray_unit[idx]
                                counts_Yaw[yawbin_idx] += 1/minian_freq

                            if unique_unit in AllCellDict_PC_Sh:
                                AllCellDict_PC_Sh[unique_unit][session] = nr_tot_act_biaised_PC
                                AllTimeDict_PC_Sh[unique_unit][session] = counts_PC

                                AllCellDict_HD_Sh[unique_unit][session] = nr_tot_act_biaised_HD 
                                AllTimeDict_HD_Sh[unique_unit][session] = counts_HD

                                AllCellDict_Roll_Sh[unique_unit][session] = nr_tot_act_biaised_Roll 
                                AllTimeDict_Roll_Sh[unique_unit][session] = counts_Roll
                                
                                AllCellDict_Pitch_Sh[unique_unit][session] = nr_tot_act_biaised_Pitch 
                                AllTimeDict_Pitch_Sh[unique_unit][session] = counts_Pitch

                                AllCellDict_Yaw_Sh[unique_unit][session] = nr_tot_act_biaised_Yaw 
                                AllTimeDict_Yaw_Sh[unique_unit][session] = counts_Yaw

                                AllCellDict_Speed_Sh[unique_unit][session] = nr_tot_act_biaised_Speed 
                                AllTimeDict_Speed_Sh[unique_unit][session] = counts_Speed

                                AllCellDict_AHV_Sh[unique_unit][session] = nr_tot_act_biaised_AHV 
                                AllTimeDict_AHV_Sh[unique_unit][session] = counts_AHV

                            else :
                                AllCellDict_PC_Sh[unique_unit] = {}
                                AllCellDict_PC_Sh[unique_unit][session] = nr_tot_act_biaised_PC 
                                AllTimeDict_PC_Sh[unique_unit] = {}
                                AllTimeDict_PC_Sh[unique_unit][session] = counts_PC
                                
                                AllCellDict_HD_Sh[unique_unit] = {}
                                AllCellDict_HD_Sh[unique_unit][session] = nr_tot_act_biaised_HD 
                                AllTimeDict_HD_Sh[unique_unit] = {}
                                AllTimeDict_HD_Sh[unique_unit][session] = counts_HD

                                AllCellDict_Roll_Sh[unique_unit] = {}
                                AllCellDict_Roll_Sh[unique_unit][session] = nr_tot_act_biaised_Roll 
                                AllTimeDict_Roll_Sh[unique_unit] = {}
                                AllTimeDict_Roll_Sh[unique_unit][session] = counts_Roll

                                AllCellDict_Pitch_Sh[unique_unit] = {}
                                AllCellDict_Pitch_Sh[unique_unit][session] = nr_tot_act_biaised_Pitch 
                                AllTimeDict_Pitch_Sh[unique_unit] = {}
                                AllTimeDict_Pitch_Sh[unique_unit][session] = counts_Pitch

                                AllCellDict_Yaw_Sh[unique_unit] = {}
                                AllCellDict_Yaw_Sh[unique_unit][session] = nr_tot_act_biaised_Yaw 
                                AllTimeDict_Yaw_Sh[unique_unit] = {}
                                AllTimeDict_Yaw_Sh[unique_unit][session] = counts_Yaw

                                AllCellDict_Speed_Sh[unique_unit] = {}
                                AllCellDict_Speed_Sh[unique_unit][session] = nr_tot_act_biaised_Speed 
                                AllTimeDict_Speed_Sh[unique_unit] = {}
                                AllTimeDict_Speed_Sh[unique_unit][session] = counts_Speed

                                AllCellDict_AHV_Sh[unique_unit] = {}
                                AllCellDict_AHV_Sh[unique_unit][session] = nr_tot_act_biaised_AHV
                                AllTimeDict_AHV_Sh[unique_unit] = {}
                                AllTimeDict_AHV_Sh[unique_unit][session] = counts_AHV
                
                    sess+=1

        #########################################################################################
                        # Average activity map across sessions for each cell #
        ##########################################################################################
        
        # Compute Place cells activity
        for uniquecell in AllCellDict_PC_Sh.keys():
            sums = {}
            for subdict in AllCellDict_PC_Sh[uniquecell].values():
                for key, value in subdict.items():
                    sums[key] = sums.get(key, 0) + value
            nr_tot_act =  {key: sums[key] for key in sums}
            nr_tot_act_array = [[None for _ in range(max_col)] for _ in range(max_row)] # Initialize 2D array
            for (row, col), value in nr_tot_act.items():# Fill array
                nr_tot_act_array[row][col] = value    
            nr_tot_act_array = [[np.nan if v is None else v for v in row] for row in nr_tot_act_array] # Replace None with np.nan
            nr_tot_act_array = np.array(nr_tot_act_array, dtype=float)
            nr_tot_act_array=nr_tot_act_array.T

            array_for_filter = np.nan_to_num(nr_tot_act_array, nan=0.0) # Replace nan with 0 for Gaussian filtering
            act_smoothed_array = gaussian_filter(array_for_filter, sigma=sigma)    # Apply 2D Gaussian filter

            sums = {}
            for subdict in AllTimeDict_PC_Sh[uniquecell].values():
                for key, value in subdict.items():
                    sums[key] = sums.get(key, 0) + value
            nr_tot_time = {key: sums[key] for key in sums}
            nr_tot_time_array = [[None for _ in range(max_col)] for _ in range(max_row)] # Initialize 2D array
            for (row, col), value in nr_tot_time.items():# Fill array
                nr_tot_time_array[row][col] = value
            nr_tot_time_array = [[np.nan if v is None else v for v in row] for row in nr_tot_time_array] # Replace None with np.nan
            nr_tot_time_array = np.array(nr_tot_time_array, dtype=float)
            nr_tot_time_array=nr_tot_time_array.T

            array_for_filter = np.nan_to_num(nr_tot_time_array, nan=0.0) # Replace nan with 0 for Gaussian filtering
            time_smoothed_array = gaussian_filter(array_for_filter, sigma=sigma)    # Apply 2D Gaussian filter
            time_smoothed_array[time_smoothed_array < .5] = np.nan  # do not consider if spent less than 0.1 second in the square

            rows, cols = max_row, max_col
            center_row, center_col = rows // 2, cols // 2
            radius = min(rows, cols) // 2   
            Y, X = np.ogrid[:rows, :cols]# Create a circular mask
            dist_from_center = np.sqrt((X - center_col)**2 + (Y - center_row)**2)
            mask = dist_from_center < radius    
            act_smoothed_array_mask = np.where(mask, act_smoothed_array, np.nan)# Apply mask: set values outside circle to NaN
            time_smoothed_array_mask = np.where(mask, time_smoothed_array, np.nan)# Apply mask: set values outside circle to NaN
            
            meanact_smoothed_array_mask = np.divide(act_smoothed_array_mask, time_smoothed_array_mask)

            mean_act=np.nanmean(meanact_smoothed_array_mask)
            sum_time=np.nansum(time_smoothed_array_mask)
            ShuffledSpatial_info[uniquecell][permutation_nb] = 0
            for rows in np.arange(np.shape(meanact_smoothed_array_mask)[0]) : 
                for cols in np.arange(np.shape(meanact_smoothed_array_mask)[1]) :             
                    bin_act=meanact_smoothed_array_mask[rows,cols]
                    p_bin=time_smoothed_array_mask[rows,cols]/sum_time 
                    if ~ np.isnan(bin_act) and ~ np.isnan(p_bin): 
                        if bin_act!=0:             
                            spatial_info_bin=p_bin*(bin_act/mean_act)*math.log2(bin_act/mean_act)
                        else:
                            spatial_info_bin=0    
                        ShuffledSpatial_info[uniquecell][permutation_nb] += spatial_info_bin

        # Compute HD activity
        for uniquecell, inner_dict in AllCellDict_HD_Sh.items():
            first_array = next(iter(inner_dict.values()))
            sumact_array = np.zeros_like(first_array)
            for arr in inner_dict.values():
                sumact_array += arr
            inner_dict= AllTimeDict_HD_Sh[uniquecell]
            timesum_array = np.zeros_like(first_array)
            for arr in inner_dict.values():
                timesum_array += arr
            timesum_array = np.where(timesum_array < 0.5, np.nan, timesum_array)
            
            ShuffledHeadDirection_info[uniquecell][permutation_nb] = np.divide(sumact_array,timesum_array)    

        # Compute Roll activity
        for uniquecell, inner_dict in AllCellDict_Roll_Sh.items():
            first_array = next(iter(inner_dict.values()))
            sumact_array = np.zeros_like(first_array)
            for arr in inner_dict.values():
                sumact_array += arr
            inner_dict= AllTimeDict_Roll_Sh[uniquecell]
            timesum_array = np.zeros_like(first_array)
            for arr in inner_dict.values():
                timesum_array += arr
            timesum_array = np.where(timesum_array < 0.5, np.nan, timesum_array)
            
            ShuffledRoll_info[uniquecell][permutation_nb] = np.divide(sumact_array,timesum_array)    

        # Compute Pitch activity
        for uniquecell, inner_dict in AllCellDict_Pitch_Sh.items():
            first_array = next(iter(inner_dict.values()))
            sumact_array = np.zeros_like(first_array)
            for arr in inner_dict.values():
                sumact_array += arr
            inner_dict= AllTimeDict_Pitch_Sh[uniquecell]
            timesum_array = np.zeros_like(first_array)
            for arr in inner_dict.values():
                timesum_array += arr
            timesum_array = np.where(timesum_array < 0.5, np.nan, timesum_array)
            
            ShuffledPitch_info[uniquecell][permutation_nb] = np.divide(sumact_array,timesum_array)    

        # Compute Yaw activity
        for uniquecell, inner_dict in AllCellDict_Yaw_Sh.items():
            first_array = next(iter(inner_dict.values()))
            sumact_array = np.zeros_like(first_array)
            for arr in inner_dict.values():
                sumact_array += arr
            inner_dict= AllTimeDict_Yaw_Sh[uniquecell]
            timesum_array = np.zeros_like(first_array)
            for arr in inner_dict.values():
                timesum_array += arr
            timesum_array = np.where(timesum_array < 0.5, np.nan, timesum_array)
            
            ShuffledYaw_info[uniquecell][permutation_nb] = np.divide(sumact_array,timesum_array)    

        # Compute AHV activity
        for uniquecell, inner_dict in AllCellDict_AHV_Sh.items():
            first_array = next(iter(inner_dict.values()))
            sumact_array = np.zeros_like(first_array)
            for arr in inner_dict.values():
                sumact_array += arr
            inner_dict= AllTimeDict_AHV_Sh[uniquecell]
            timesum_array = np.zeros_like(first_array)
            for arr in inner_dict.values():
                timesum_array += arr
            timesum_array = np.where(timesum_array < 0.5, np.nan, timesum_array)
        
            ShuffledAHV_info[uniquecell][permutation_nb] = np.divide(sumact_array,timesum_array)    

        # Compute Speed activity
        for uniquecell, inner_dict in AllCellDict_Speed_Sh.items():
            first_array = next(iter(inner_dict.values()))
            sumact_array = np.zeros_like(first_array)
            for arr in inner_dict.values():
                sumact_array += arr
            inner_dict= AllTimeDict_Speed_Sh[uniquecell]
            timesum_array = np.zeros_like(first_array)
            for arr in inner_dict.values():
                timesum_array += arr
            timesum_array = np.where(timesum_array < 0.5, np.nan, timesum_array)

            ShuffledSpeed_info[uniquecell][permutation_nb] = np.divide(sumact_array,timesum_array)    

    with open(os.path.join(folder_to_save, f"ShuffledSpatialInfo_{worker_id}.pkl"), "wb") as f:
        pickle.dump(ShuffledSpatial_info, f)
    with open(os.path.join(folder_to_save, f"ShuffledHeadDirection_{worker_id}.pkl"), "wb") as f:
        pickle.dump(ShuffledHeadDirection_info, f)
    with open(os.path.join(folder_to_save, f"ShuffledRoll_{worker_id}.pkl"), "wb") as f:
        pickle.dump(ShuffledRoll_info, f)
    with open(os.path.join(folder_to_save, f"ShuffledPitch_{worker_id}.pkl"), "wb") as f:
        pickle.dump(ShuffledPitch_info, f)
    with open(os.path.join(folder_to_save, f"ShuffledYaw_{worker_id}.pkl"), "wb") as f:
        pickle.dump(ShuffledYaw_info, f)
    with open(os.path.join(folder_to_save, f"ShuffledSpeed_{worker_id}.pkl"), "wb") as f:
        pickle.dump(ShuffledSpeed_info, f)
    with open(os.path.join(folder_to_save, f"ShuffledAHV_{worker_id}.pkl"), "wb") as f:
        pickle.dump(ShuffledAHV_info, f)


if __name__ == "__main__":
    
    nb_tot_permutation=1000

    FolderNameSave=str(datetime.now())[:19] # Get the current date and time
    FolderNameSave = FolderNameSave.replace(" ", "_").replace(".", "_").replace(":", "_")
    destination_folder= f"/mnt/data/AurelieB_other/0_NeuronIdentity_{nb_tot_permutation}Permutations_perCPUs/" # _{FolderNameSave}" 
    os.makedirs(destination_folder, exist_ok=True)
    folder_to_save=Path(destination_folder)
    print(folder_to_save)


    cpus = multiprocessing.cpu_count()
    print(cpus, 'CPUs found, dividing the', nb_tot_permutation, 'permutations to each worker.')
        
    base = nb_tot_permutation // cpus        # base tasks per worker
    remainder = nb_tot_permutation % cpus    # extra tasks to distribute
    tasks = []
    start = 0
    for i in range(cpus):
        end = start + base + (1 if i < remainder else 0)
        tasks.append((i, start, end))
        start = end

    with ProcessPoolExecutor(cpus) as ex:
        list(ex.map(runpermutations, tasks))
    
    destination_folder2= f"/mnt/data/AurelieB_other/0_NeuronIdentity_AllPermutations/" # _{FolderNameSave}" 
    os.makedirs(destination_folder2, exist_ok=True)
    folder_to_save2=Path(destination_folder2)
    print(folder_to_save2)
    
    Vars= ['ShuffledAHV', 'ShuffledHeadDirection',  'ShuffledPitch', 'ShuffledRoll', 'ShuffledSpeed', 'ShuffledYaw', 'ShuffledSpatial']
    for Var in Vars:
        pkl_files = list(folder_to_save.glob(f"{Var}*.pkl"))
        pooled = defaultdict(lambda: defaultdict(list))
        for file in pkl_files:
            with open(file, "rb") as f:
                data = pickle.load(f)  # assume dict of dicts
                for key, subdict in data.items():
                    for subkey, value in subdict.items():
                        pooled[key][subkey]=value
        pooled = {k: dict(v) for k, v in pooled.items()}

        with open(os.path.join(folder_to_save2, f"{Var}Info.pkl"), "wb") as f:
            pickle.dump(pooled, f)




