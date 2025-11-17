
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
table_radius = 290 / 2
square_size = pixel_to_cm * 6

bin_size = 10
bin_edges = np.arange(0, 361, bin_size)
bin_centers = np.radians(bin_edges[:-1] + bin_size/2)



def remove_outliers_avg_filter(data):
    data = np.array(data, dtype=float)  # Ensure NumPy array with float type
    filtered_data = np.copy(data)  # Copy to avoid modifying original data
    for i in range(len(data)):
        if not np.isnan(data[i]):  # Skip valid values
            continue
        # Find the closest previous non-NaN value
        prev_idx = i - 1
        while prev_idx >= 0 and np.isnan(data[prev_idx]):
            prev_idx -= 1        
        # Find the closest next non-NaN value
        next_idx = i + 1
        while next_idx < len(data) and np.isnan(data[next_idx]):
            next_idx += 1
        # Compute average if both values exist
        if prev_idx >= 0 and next_idx < len(data):
            filtered_data[i] = (data[prev_idx] + data[next_idx]) / 2
        # If neither exists, NaN remains
    return filtered_data

def find_closest_index_sorted(arr, target):
    idx = bisect.bisect_left(arr, target)  # Find the insertion point
    if idx == 0:
        return 0
    if idx == len(arr):
        return len(arr) - 1
    before = idx - 1
    after = idx
    return before if abs(arr[before] - target) <= abs(arr[after] - target) else after

# Sample callback function
def update_my_folder(chooser):
    global dpath
    dpath = chooser.selected
    %store dpath
    return 

def detect_longest_lowest_sequence(arr, margin=0):
    min_val = np.nanmin(arr)  # Find minimum value
    threshold = min_val + margin  # Define threshold based on margin    
    # Get indices where values are within the threshold
    min_indices = np.where(arr <= threshold)[0]
    # Identify consecutive sequences
    longest_sequence = None
    if len(min_indices) > 0:
        start = min_indices[0]
        max_duration = 0  # Track longest duration        
        for i in range(1, len(min_indices)):
            if min_indices[i] != min_indices[i - 1] + 1:  # Not consecutive
                duration = min_indices[i - 1] - start + 1
                if duration > max_duration:
                    max_duration = duration
                    longest_sequence = (start, min_indices[i - 1], duration)
                start = min_indices[i]  # Reset start index        
        # Check last detected sequence
        duration = min_indices[-1] - start + 1
        if duration > max_duration:
            longest_sequence = (start, min_indices[-1], duration)
    return min_val, threshold, longest_sequence


def resample_matrix(data, orig_rate, target_rate, axis=0):
    if orig_rate == target_rate:
        return data.copy()
    # Compute integer up/down factors using GCD
    up = int(target_rate)
    down = int(orig_rate)
    factor = gcd(up, down)
    up //= factor
    down //= factor
    return resample_poly(data, up=up, down=down, axis=axis)

def calculate_relative_distance(x1, y1, x2, y2):
    return math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)

def calculate_distance_run(x_coords, y_coords):
    distances = np.sqrt(np.diff(x_coords) ** 2 + np.diff(y_coords) ** 2)
    for i in range(1, len(distances) - 1):
        if np.isnan(distances[i]):
            neighbors = [distances[i-1], distances[i+1]]
            distances[i] = np.mean([x for x in neighbors if not np.isnan(x)])
    total_distance_cm = np.nansum(distances) / pixel_to_cm  # Convert to cm
    return total_distance_cm, distances

def find_long_non_nan_sequences(arr, min_length=100):
    mask = ~np.isnan(arr)  # True for non-NaN values
    diff = np.diff(np.concatenate(([0], mask.astype(int), [0])))  # Add padding to detect edges
    starts = np.where(diff == 1)[0]  # Where a sequence starts
    ends = np.where(diff == -1)[0]   # Where a sequence ends
    sequences = [arr[start:end] for start, end in zip(starts, ends) if (end - start) > min_length]
    return sequences

def remove_outliers_median_filter(data, window=1):
    data = np.array(data, dtype=float)  # Ensure NumPy array with float type
    filtered_data = np.copy(data)  # Copy to avoid modifying original data
    half_window = window // 2
    for i in range(len(data)):
        # Define window range, ensuring it doesn't exceed bounds
        start = max(0, i - half_window)
        end = min(len(data), i + half_window + 1)
        # Extract local values in window
        local_values = data[start:end]
        # Check if the window contains at least one non-NaN value
        if np.all(np.isnan(local_values)):
            median_value = np.nan  # Keep NaN if no valid numbers
        else:
            median_value = np.nanmedian(local_values)  # Compute median ignoring NaNs
        # Replace only if the current value is not NaN
        if not np.isnan(data[i]):
            filtered_data[i] = median_value
    return filtered_data

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
    return x, y, speeds

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

def Convert(string):
            li = list(string.split(", "))
            li2 = len(li)
            return li2

def marcenkopastur(significance):
    nbins = significance.nbins
    nneurons = significance.nneurons
    tracywidom = significance.tracywidom
    q = float(nbins)/float(nneurons)
    lambdaMax = pow((1+np.sqrt(1/q)),2)
    lambdaMax += tracywidom*pow(nneurons,-2./3)
    return lambdaMax

def getlambdacontrol(zactmat_):
    significance_ = PCA()
    significance_.fit(zactmat_.T)
    lambdamax_ = np.max(significance_.explained_variance_)
    return lambdamax_

def binshuffling(zactmat,significance):
    np.random.seed()
    lambdamax_ = np.zeros(significance.nshu)
    for shui in range(significance.nshu):
        zactmat_ = np.copy(zactmat)
        for (neuroni,activity) in enumerate(zactmat_):
            randomorder = np.argsort(np.random.rand(significance.nbins))
            zactmat_[neuroni,:] = activity[randomorder]
        lambdamax_[shui] = getlambdacontrol(zactmat_)
    lambdaMax = np.percentile(lambdamax_,significance.percentile)
    return lambdaMax

def circshuffling(zactmat,significance):
    np.random.seed()
    lambdamax_ = np.zeros(significance.nshu)
    for shui in range(significance.nshu):
        zactmat_ = np.copy(zactmat)
        for (neuroni,activity) in enumerate(zactmat_):
            cut = int(np.random.randint(significance.nbins*2))
            zactmat_[neuroni,:] = np.roll(activity,cut)
        lambdamax_[shui] = getlambdacontrol(zactmat_)
    lambdaMax = np.percentile(lambdamax_,significance.percentile)
    return lambdaMax

def runSignificance(zactmat,significance):
    if significance.nullhyp == 'mp':
        lambdaMax = marcenkopastur(significance)
    elif significance.nullhyp == 'bin':
        lambdaMax = binshuffling(zactmat,significance)
    elif significance.nullhyp == 'circ':
        lambdaMax = circshuffling(zactmat,significance)
    else:
        print('ERROR !')
        print('    nyll hypothesis method '+str(nullhyp)+' not understood')
        significance.nassemblies = np.nan
    nassemblies = np.sum(significance.explained_variance_>lambdaMax)
    significance.nassemblies = nassemblies
    return significance

def extractPatterns(actmat,significance,method):
    nassemblies = significance.nassemblies
    if method == 'pca':
        idxs = np.argsort(-significance.explained_variance_)[0:nassemblies]
        patterns = significance.components_[idxs,:]
    elif method == 'ica':
        from sklearn.decomposition import FastICA
        ica = FastICA(n_components=nassemblies, max_iter=1000)
        ica.fit(actmat.T)
        patterns = ica.components_
    else:
        print('ERROR !')
        print('    assembly extraction method '+str(method)+' not understood')
        patterns = np.nan
    if patterns is not np.nan:
        patterns = patterns.reshape(nassemblies,-1)
        norms = np.linalg.norm(patterns,axis=1)
        patterns /= np.matlib.repmat(norms,np.size(patterns,1),1).T
    return patterns

def runPatterns(actmat, method='ica', nullhyp = 'mp', nshu = 1000, percentile = 99, tracywidom = False):
    nneurons = np.size(actmat,0)
    nbins = np.size(actmat,1)
    silentneurons = np.var(actmat,axis=1)==0
    actmat_ = actmat[~silentneurons,:]
    zactmat_ = stats.zscore(actmat_,axis=1)
    significance = PCA()
    significance.fit(zactmat_.T)
    significance.nneurons = nneurons
    significance.nbins = nbins
    significance.nshu = nshu
    significance.percentile = percentile
    significance.tracywidom = tracywidom
    significance.nullhyp = nullhyp
    significance = runSignificance(zactmat_,significance)
    if np.isnan(significance.nassemblies):
        patterns = []
        zactmat = []
        significance = []
        #return
    if significance.nassemblies<1:
        print('WARNING 1!')
        print('    no assembly detecded!')
        patterns = []
        zactmat = []
        significance = []
    else:
        patterns_ = extractPatterns(zactmat_,significance,method)
        if patterns_ is np.nan:
            patterns = []
            zactmat = []
            significance = []
            #return
        patterns = np.zeros((np.size(patterns_,0),nneurons))
        patterns[:,~silentneurons] = patterns_
        zactmat = np.copy(actmat)
        zactmat[~silentneurons,:] = zactmat_
    return patterns,significance,zactmat

def computeAssemblyActivity(patterns,zactmat,zerodiag = True):
    if len(patterns) == 0:
        print('WARNING 2!')
        print('    no assembly detecded!')
        assemblyAct = []
    else:
        nassemblies = len(patterns)
        nbins = np.size(zactmat,1)
        assemblyAct = np.zeros((nassemblies,nbins))
        for (assemblyi,pattern) in enumerate(patterns):
            projMat = np.outer(pattern,pattern)
            projMat -= zerodiag*np.diag(np.diag(projMat))
            for bini in range(nbins):
                assemblyAct[assemblyi,bini] = np.dot(np.dot(zactmat[:,bini],projMat),zactmat[:,bini])
    return assemblyAct

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




AllCellDict_PC={}
AllTimeDict_PC={}
AllCellDictSess_PC={}
AllTimeDictSess_PC={}

AllCellDict_HD={}
AllTimeDict_HD={}
AllCellDictSess_HD={}
AllTimeDictSess_HD={}

AllCellDict_Roll={}
AllTimeDict_Roll={}
AllCellDictSess_Roll={}
AllTimeDictSess_Roll={}

AllCellDict_Pitch={}
AllTimeDict_Pitch={}
AllCellDictSess_Pitch={}
AllTimeDictSess_Pitch={}

AllCellDict_Yaw={}
AllTimeDict_Yaw={}
AllCellDictSess_Yaw={}
AllTimeDictSess_Yaw={}

AllCellDict_Speed={}
AllTimeDict_Speed={}
AllCellDictSess_Speed={}
AllTimeDictSess_Speed={}

for dpath in Path(dir).glob('**/Exploration_task/mappingsAB.pkl'):

    mappfile = open(dpath.parents[0]/ f'mappingsAB.pkl', 'rb')
    mapping = pickle.load(mappfile)
    mapping_sess = mapping['session']   
        
    centfile = open(dpath.parents[0]/ f'centsAB.pkl', 'rb')
    centroids = pickle.load(centfile) 

    mice = dpath.parents[1].parts[-1]
    
    print(f"Processing mouse {mice} ")

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
                session_type=minianpath.parents[3].name
                session_date=minianpath.parents[4].name
                session_time=minianpath.parents[0].name
                session=session_time
            print(f"Processing {session_type} session: {session} on the {session_date} ")

            # Minian data 
                        
            minian_ds = open_minian(minianpath)
            Co = minian_ds['C']  # calcium traces
            minian_freq=round(1/np.mean(np.diff(np.array(tsmini)/1000)))
            print('... miniscope sample rate =', minian_freq, 'Hz')            
            try: 
                TodropFile = minianpath / f'TodropFileAB.json'
                with open(TodropFile, 'r') as f:
                    unit_to_drop = json.load(f)
            except:
                TodropFile = minianpath.parent / f'TodropFileAB.json'
                with open(TodropFile, 'r') as f:
                    unit_to_drop = json.load(f)

            if V4subfolder:
                V4subfolder_id = int(minianpath.parent.name[-1]) - 1 
                ts_start = V4subfolder_id*15*1000 # cause 15 videos per subfolders and 1000 frames per videos
                ts_stop = np.shape(Co)[1] + ts_start
                tsmini_sub=tsmini[ts_start:ts_stop]
                tsmini_sub=tsmini_sub.reset_index(drop=True)    
            else:
                tsmini_sub=tsmini


            # Extraction des quaternions

            vestibular_df= pd.read_csv(list(session_path.glob('*V4_Miniscope/headOrientation.csv'))[0])
            qw = vestibular_df['qw'].to_numpy()
            qx = vestibular_df['qx'].to_numpy()
            qy = vestibular_df['qy'].to_numpy()
            qz = vestibular_df['qz'].to_numpy()
            roll, pitch, yaw = quaternion_to_euler(qw, qx, qy, qz)
            roll_deg = np.degrees(roll) # nose goes clockwise or anticlockwise
            pitch_deg = np.degrees(pitch) # nose up or down
            yaw_deg = np.degrees(yaw) # nose moves from side to side                
            

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
            nose_x3, nose_y3, speed = limit_speed(nose_x2, nose_y2, max_speed=max_speed)
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
            neck_x3, neck_y3, speed = limit_speed(neck_x2, neck_y2, max_speed=max_speed)
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
                print(f'... no cells kept in the session: {session}')
                continue
            print(f'... {nb_unit} kept cells in the session')


            # If the miniscope recording started after the webcam recording, cut the webcam data
            if tsmini_sub.iloc[0] > tswebcam[start_frame]:
                Newstart_frame = np.where(tswebcam >= tsmini_sub.iloc[0].item())[0][1].item()
                print(f'... webcam data cut to match miniscope length, new start at frame {Newstart_frame} (instead of {start_frame})')
                start_frame = Newstart_frame 
            # If the miniscope recording is shorter than the webcam recording, cut the webcam data
            if tsmini_sub.iloc[-1] < tswebcam[last_frame-1]:
                Newlast_frame = np.where(tswebcam <= tsmini_sub.iloc[-1].item())[0][-1].item()
                print(f'... webcam data cut to match miniscope length, new end at frame {Newlast_frame} (instead of {last_frame})')
                last_frame = Newlast_frame
            

            # Align data

            nose_x_ = nose_x[start_frame:last_frame]
            nose_y_ = nose_y[start_frame:last_frame]
            neck_x_ = neck_x[start_frame:last_frame]
            neck_y_ = neck_y[start_frame:last_frame]            
            angles_deg = np.array([direction(x1, y1, x2, y2) for x1, y1, x2, y2 in zip(neck_x_, neck_y_, nose_x_, nose_y_)])

            speed= np.concatenate(([0], speed))
            speed_= speed[start_frame:last_frame]

            closest_start= find_closest_index_sorted(tsmini_sub, tswebcam[start_frame])
            closest_end= find_closest_index_sorted(tsmini_sub, tswebcam[last_frame-1])

            roll_deg_ = roll_deg[closest_start:closest_end]
            pitch_deg_ = pitch_deg[closest_start:closest_end]
            yaw_deg_ = yaw_deg[closest_start:closest_end]

            # Keep only crossregistered cells

            Carray = C_sel.to_numpy().T
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

                    for idx, (px, py, angle, _roll, _pitch, _yaw, _speed) in enumerate(zip(nose_x_, nose_y_, angles_deg, roll_deg_, pitch_deg_, yaw_deg_, speed_)):
                        if np.sqrt((px - table_center_x)**2 + (py - table_center_y)**2) > table_radius:
                            continue  # skip points outside circle
                        closest_point = find_closest_index_sorted(tsmini_sub, tswebcam[start_frame+idx])

                        ix = int(np.floor((px - (table_center_x - n/2 * square_size)) / square_size))
                        iy = int(np.floor((py - (table_center_y - n/2 * square_size)) / square_size))
                        nr_tot_act_biaised_PC[(ix, iy)] += Carray_unit[closest_point]
                        counts_PC[(ix, iy)] += 1/frame_rate 

                        bin_idx = int(angle // bin_size)
                        nr_tot_act_biaised_HD[bin_idx] += Carray_unit[closest_point]
                        counts_HD[bin_idx] += 1    

                        rollbin_idx = int(_roll // bin_size)
                        nr_tot_act_biaised_Roll[rollbin_idx] += Carray_unit[closest_point]
                        counts_Roll[rollbin_idx] += 1      

                        pitchbin_idx = int(_pitch // bin_size)
                        nr_tot_act_biaised_Pitch[pitchbin_idx] += Carray_unit[closest_point]
                        counts_Pitch[pitchbin_idx] += 1

                        yawbin_idx = int(_yaw // bin_size)
                        nr_tot_act_biaised_Yaw[yawbin_idx] += Carray_unit[closest_point]
                        counts_Yaw[yawbin_idx] += 1  
                        
                        speedbin_idx = min(int(_speed), 19) # between 0 and 20 
                        nr_tot_act_biaised_Speed[speedbin_idx] += Carray_unit[closest_point]
                        counts_Speed[speedbin_idx] += 1  


                    if unique_unit in AllCellDict_PC:
                        AllCellDict_PC[unique_unit][session] = nr_tot_act_biaised_PC
                        AllTimeDict_PC[unique_unit][session] = counts_PC

                        AllCellDict_HD[unique_unit][session] = nr_tot_act_biaised_HD 
                        AllTimeDict_HD[unique_unit][session] = counts_HD

                        AllCellDict_Roll[unique_unit][session] = nr_tot_act_biaised_Roll 
                        AllTimeDict_Roll[unique_unit][session] = counts_Roll
                        
                        AllCellDict_Pitch[unique_unit][session] = nr_tot_act_biaised_Pitch 
                        AllTimeDict_Pitch[unique_unit][session] = counts_Pitch

                        AllCellDict_Yaw[unique_unit][session] = nr_tot_act_biaised_Yaw 
                        AllTimeDict_Yaw[unique_unit][session] = counts_Yaw

                        AllCellDict_Speed[unique_unit][session] = nr_tot_act_biaised_Speed 
                        AllTimeDict_Speed[unique_unit][session] = counts_Speed

                    else :
                        AllCellDict_PC[unique_unit] = {}
                        AllCellDict_PC[unique_unit][session] = nr_tot_act_biaised_PC 
                        AllTimeDict_PC[unique_unit] = {}
                        AllTimeDict_PC[unique_unit][session] = counts_PC
                        
                        AllCellDict_HD[unique_unit] = {}
                        AllCellDict_HD[unique_unit][session] = nr_tot_act_biaised_HD 
                        AllTimeDict_HD[unique_unit] = {}
                        AllTimeDict_HD[unique_unit][session] = counts_HD

                        AllCellDict_Roll[unique_unit] = {}
                        AllCellDict_Roll[unique_unit][session] = nr_tot_act_biaised_Roll 
                        AllTimeDict_Roll[unique_unit] = {}
                        AllTimeDict_Roll[unique_unit][session] = counts_Roll

                        AllCellDict_Pitch[unique_unit] = {}
                        AllCellDict_Pitch[unique_unit][session] = nr_tot_act_biaised_Pitch 
                        AllTimeDict_Pitch[unique_unit] = {}
                        AllTimeDict_Pitch[unique_unit][session] = counts_Pitch

                        AllCellDict_Yaw[unique_unit] = {}
                        AllCellDict_Yaw[unique_unit][session] = nr_tot_act_biaised_Yaw 
                        AllTimeDict_Yaw[unique_unit] = {}
                        AllTimeDict_Yaw[unique_unit][session] = counts_Yaw

                        AllCellDict_Speed[unique_unit] = {}
                        AllCellDict_Speed[unique_unit][session] = nr_tot_act_biaised_Speed 
                        AllTimeDict_Speed[unique_unit] = {}
                        AllTimeDict_Speed[unique_unit][session] = counts_Speed

                    AllCellDictSess_PC[sess_unit] = nr_tot_act_biaised_PC 
                    AllTimeDictSess_PC[sess_unit] = counts_PC

                    AllCellDictSess_HD[sess_unit] = nr_tot_act_biaised_HD 
                    AllTimeDictSess_HD[sess_unit] = counts_HD

                    AllCellDictSess_Roll[sess_unit] = nr_tot_act_biaised_Roll 
                    AllTimeDictSess_Roll[sess_unit] = counts_Roll

                    AllCellDictSess_Pitch[sess_unit] = nr_tot_act_biaised_Pitch 
                    AllTimeDictSess_Pitch[sess_unit] = counts_Pitch

                    AllCellDictSess_Yaw[sess_unit] = nr_tot_act_biaised_Yaw 
                    AllTimeDictSess_Yaw[sess_unit] = counts_Yaw

                    AllCellDictSess_Speed[sess_unit] = nr_tot_act_biaised_Speed 
                    AllTimeDictSess_Speed[sess_unit] = counts_Speed


            sess+=1

    print(f'{len(AllCellDict_PC.keys())} unique cells found so far')  
    
#########################################################################################
                # Average activity map across sessions for each cell #
##########################################################################################

SumActDict_PC={}
for uniquecell in AllCellDict_PC.keys():
    sums = {}
    for subdict in AllCellDict_PC[uniquecell].values():
        for key, value in subdict.items():
            sums[key] = sums.get(key, 0) + value
    SumActDict_PC[uniquecell] = {key: sums[key] for key in sums}
SumTimeDict_PC={}
for uniquecell in AllTimeDict_PC.keys():
    sums = {}
    for subdict in AllTimeDict_PC[uniquecell].values():
        for key, value in subdict.items():
            sums[key] = sums.get(key, 0) + value
    SumTimeDict_PC[uniquecell] = {key: sums[key] for key in sums}

SumActDict_HD={}
for outer_key, inner_dict in AllCellDict_HD.items():
    first_array = next(iter(inner_dict.values()))
    sum_array = np.zeros_like(first_array)
    for arr in inner_dict.values():
        sum_array += arr
    SumActDict_HD[outer_key] = sum_array
SumTimeDict_HD={}
for outer_key, inner_dict in AllTimeDict_HD.items():
    first_array = next(iter(inner_dict.values()))
    sum_array = np.zeros_like(first_array)
    for arr in inner_dict.values():
        sum_array += arr
    SumTimeDict_HD[outer_key] = sum_array

SumActDict_Roll={}
for outer_key, inner_dict in AllCellDict_Roll.items():
    first_array = next(iter(inner_dict.values()))
    sum_array = np.zeros_like(first_array)
    for arr in inner_dict.values():
        sum_array += arr
    SumActDict_Roll[outer_key] = sum_array
SumTimeDict_Roll={}
for outer_key, inner_dict in AllTimeDict_Roll.items():
    first_array = next(iter(inner_dict.values()))
    sum_array = np.zeros_like(first_array)
    for arr in inner_dict.values():
        sum_array += arr
    SumTimeDict_Roll[outer_key] = sum_array

SumActDict_Pitch={}
for outer_key, inner_dict in AllCellDict_Pitch.items():
    first_array = next(iter(inner_dict.values()))
    sum_array = np.zeros_like(first_array)
    for arr in inner_dict.values():
        sum_array += arr
    SumActDict_Pitch[outer_key] = sum_array
SumTimeDict_Pitch={}
for outer_key, inner_dict in AllTimeDict_Pitch.items():
    first_array = next(iter(inner_dict.values()))
    sum_array = np.zeros_like(first_array)
    for arr in inner_dict.values():
        sum_array += arr
    SumTimeDict_Pitch[outer_key] = sum_array

SumActDict_Yaw={}
for outer_key, inner_dict in AllCellDict_Yaw.items():
    first_array = next(iter(inner_dict.values()))
    sum_array = np.zeros_like(first_array)
    for arr in inner_dict.values():
        sum_array += arr
    SumActDict_Yaw[outer_key] = sum_array
SumTimeDict_Yaw={}
for outer_key, inner_dict in AllTimeDict_Yaw.items():
    first_array = next(iter(inner_dict.values()))
    sum_array = np.zeros_like(first_array)
    for arr in inner_dict.values():
        sum_array += arr
    SumTimeDict_Yaw[outer_key] = sum_array

SumActDict_Speed={}
for outer_key, inner_dict in AllCellDict_Speed.items():
    first_array = next(iter(inner_dict.values()))
    sum_array = np.zeros_like(first_array)
    for arr in inner_dict.values():
        sum_array += arr
    SumActDict_Speed[outer_key] = sum_array
SumTimeDict_Speed={}
for outer_key, inner_dict in AllTimeDict_Speed.items():
    first_array = next(iter(inner_dict.values()))
    sum_array = np.zeros_like(first_array)
    for arr in inner_dict.values():
        sum_array += arr
    SumTimeDict_Speed[outer_key] = sum_array


print(f'Total unique cell = {len(SumActDict_PC.keys())}')  



FolderNameSave=str(datetime.now())[:19] # Get the current date and time
FolderNameSave = FolderNameSave.replace(" ", "_").replace(".", "_").replace(":", "_")
if local:
    destination_folder= f"//10.69.168.1/crnldata/forgetting/Aurelie/MiniscopeOE_analysis/Exploration_task/0_NeuronIdentity_{FolderNameSave}" 
else: 
    destination_folder= f"/crnldata/forgetting/Aurelie/MiniscopeOE_analysis/Exploration_task/0_NeuronIdentity_{FolderNameSave}" 
os.makedirs(destination_folder)
folder_to_save=Path(destination_folder)

plot = False


sorted_keys = sorted(AllCellDictSess_PC.keys())
AllCellDictSess_PC = {key: AllCellDictSess_PC[key]for rank, key in enumerate(sorted_keys)}

if plot: 
    nb_subplot=len(AllCellDictSess_PC.keys())
    rows = int(np.ceil(np.sqrt(nb_subplot)))  # Rows: ceil(sqrt(X))
    cols = int(np.ceil(np.sqrt(nb_subplot)))  # Columns: floor(sqrt(X))

    # Create the figure with a 2x2 grid
    fig, axs = plt.subplots(rows, cols, figsize=(int(nb_subplot/10), int(nb_subplot/10)))
    axs = axs.flatten()
    plt.tight_layout()

sigma = 1  # Standard deviation for Gaussian kernel
max_row = int(table_radius*2/square_size)+1
max_col = int(table_radius*2/square_size)+1

Smoothed_ActSess={}
Smoothed_TimeSess={}
Spatial_infoSess = {}

for nsubplot, nr in enumerate(AllCellDictSess_PC.keys()):

    nr_tot_act = AllCellDictSess_PC[nr]    
    nr_tot_act_array = [[None for _ in range(max_col)] for _ in range(max_row)] # Initialize 2D array
    for (row, col), value in nr_tot_act.items(): # Fill array
        nr_tot_act_array[row][col] = value    
    nr_tot_act_array = [[np.nan if v is None else v for v in row] for row in nr_tot_act_array] # Replace None with np.nan
    nr_tot_act_array = np.array(nr_tot_act_array, dtype=float)
    nr_tot_act_array=nr_tot_act_array.T

    array_for_filter = np.nan_to_num(nr_tot_act_array, nan=0.0) # Replace nan with 0 for Gaussian filtering
    act_smoothed_array = gaussian_filter(array_for_filter, sigma=sigma)    # Apply 2D Gaussian filter


    nr_tot_time = AllTimeDictSess_PC[nr]
    nr_tot_time_array = [[None for _ in range(max_col)] for _ in range(max_row)] # Initialize 2D array
    for (row, col), value in nr_tot_time.items():# Fill array
        nr_tot_time_array[row][col] = value
    nr_tot_time_array = [[np.nan if v is None else v for v in row] for row in nr_tot_time_array] # Replace None with np.nan
    nr_tot_time_array = np.array(nr_tot_time_array, dtype=float)
    nr_tot_time_array=nr_tot_time_array.T

    array_for_filter = np.nan_to_num(nr_tot_time_array, nan=0.0) # Replace nan with 0 for Gaussian filtering
    time_smoothed_array = gaussian_filter(array_for_filter, sigma=sigma)    # Apply 2D Gaussian filter
    time_smoothed_array[time_smoothed_array < .5] = np.nan  # do not consider if spent less than 0.1 second in the square

    # Ratio
    act_smoothed_array = np.divide(act_smoothed_array, time_smoothed_array)

    
    # --- Mask outside circle ---
    rows, cols = act_smoothed_array.shape
    center_row, center_col = rows // 2, cols // 2
    radius = min(rows, cols) // 2   
    Y, X = np.ogrid[:rows, :cols]# Create a circular mask
    dist_from_center = np.sqrt((X - center_col)**2 + (Y - center_row)**2)
    mask = dist_from_center < radius    
    act_smoothed_array_mask = np.where(mask, act_smoothed_array, np.nan)# Apply mask: set values outside circle to NaN
    Smoothed_ActSess[nr]= act_smoothed_array_mask

    time_smoothed_array_mask = np.where(mask, time_smoothed_array, np.nan)# Apply mask: set values outside circle to NaN
    Smoothed_TimeSess[nr]= time_smoothed_array_mask


    # --- Spatial information ---
    mean_act=np.nanmean(Smoothed_ActSess[nr])
    sum_time=np.nansum(Smoothed_TimeSess[nr])
    Spatial_infoSess[nr] = 0
    for rows in np.arange(np.shape(Smoothed_ActSess[nr])[0]) : 
        for cols in np.arange(np.shape(Smoothed_ActSess[nr])[1]) :             
            bin_act=Smoothed_ActSess[nr][rows,cols]
            p_bin=Smoothed_TimeSess[nr][rows,cols]/sum_time 
            if ~ np.isnan(bin_act): 
                if bin_act!=0:             
                    spatial_info_bin=p_bin*(bin_act/mean_act)*math.log2(bin_act/mean_act)
                else:
                    spatial_info_bin=0    
                Spatial_infoSess[nr] += spatial_info_bin


    if plot: 
        # Mask outside circle 
        rows, cols = np.full((max_row, max_col),1).shape
        center_row, center_col = rows // 2, cols // 2
        radius = min(rows, cols) // 2   
        Y, X = np.ogrid[:rows, :cols]# Create a circular mask
        dist_from_center = np.sqrt((X - center_col)**2 + (Y - center_row)**2)
        mask = dist_from_center < radius    
        masked_array = np.where(mask, np.full((max_row, max_col),1), np.nan)# Apply mask: set values outside circle to NaN
        
        # Draw map
        axs[nsubplot].imshow(masked_array, cmap='Greys_r', origin='upper')
        axs[nsubplot].imshow(act_smoothed_array_mask, cmap='jet', origin='upper')
        axs[nsubplot].set_aspect('equal')
        axs[nsubplot].set_title(f'{nr} ({np.round(Spatial_infoSess[nr], 1)})', pad=0, loc='left', fontsize=15)
        axs[nsubplot].set_title(f"{nr.split('_')[0]}_{nr.split('_')[1]}", pad=0, loc='left', fontsize=max(7,int(nb_subplot/20)))

            
        # Remove box (spines)
        for spine in axs[nsubplot].spines.values():
            spine.set_visible(False)
        axs[nsubplot].set_xticks([])  # No x-axis ticks
        axs[nsubplot].set_yticks([])  # No y-axis ticks

if plot: 
    # Adjust layout to avoid clipping
    fig.subplots_adjust(wspace=.1, hspace=.2)
    plt.show()

sorted_keys = sorted(SumActDict_PC.keys())
SumActDict_PC = {key: SumActDict_PC[key]for rank, key in enumerate(sorted_keys)}

if plot:
    nb_subplot=len(SumActDict_PC.keys())
    rows = int(np.ceil(np.sqrt(nb_subplot)))  # Rows: ceil(sqrt(X))
    cols = int(np.ceil(np.sqrt(nb_subplot)))  # Columns: floor(sqrt(X))

    # Create the figure with a 2x2 grid
    fig, axs = plt.subplots(rows, cols, figsize=(int(nb_subplot/10), int(nb_subplot/10)))
    axs = axs.flatten()
    plt.tight_layout()

sigma = 1  # Standard deviation for Gaussian kernel
max_row = int(table_radius*2/square_size)+1
max_col = int(table_radius*2/square_size)+1

Smoothed_Act={}
Smoothed_Time={}
Spatial_info = {}

for nsubplot, nr in enumerate(SumActDict_PC.keys()):

    nr_tot_act = SumActDict_PC[nr]    
    nr_tot_act_array = [[None for _ in range(max_col)] for _ in range(max_row)] # Initialize 2D array
    for (row, col), value in nr_tot_act.items():# Fill array
        nr_tot_act_array[row][col] = value    
    nr_tot_act_array = [[np.nan if v is None else v for v in row] for row in nr_tot_act_array] # Replace None with np.nan
    nr_tot_act_array = np.array(nr_tot_act_array, dtype=float)
    nr_tot_act_array=nr_tot_act_array.T

    array_for_filter = np.nan_to_num(nr_tot_act_array, nan=0.0) # Replace nan with 0 for Gaussian filtering
    act_smoothed_array = gaussian_filter(array_for_filter, sigma=sigma)    # Apply 2D Gaussian filter


    nr_tot_time = SumTimeDict_PC[nr]
    nr_tot_time_array = [[None for _ in range(max_col)] for _ in range(max_row)] # Initialize 2D array
    for (row, col), value in nr_tot_time.items():# Fill array
        nr_tot_time_array[row][col] = value
    nr_tot_time_array = [[np.nan if v is None else v for v in row] for row in nr_tot_time_array] # Replace None with np.nan
    nr_tot_time_array = np.array(nr_tot_time_array, dtype=float)
    nr_tot_time_array=nr_tot_time_array.T

    array_for_filter = np.nan_to_num(nr_tot_time_array, nan=0.0) # Replace nan with 0 for Gaussian filtering
    time_smoothed_array = gaussian_filter(array_for_filter, sigma=sigma)    # Apply 2D Gaussian filter
    time_smoothed_array[time_smoothed_array < .5] = np.nan  # do not consider if spent less than 0.1 second in the square

    # Ratio
    act_smoothed_array = np.divide(act_smoothed_array, time_smoothed_array)

    # --- Mask outside circle ---
    rows, cols = act_smoothed_array.shape
    center_row, center_col = rows // 2, cols // 2
    radius = min(rows, cols) // 2   
    Y, X = np.ogrid[:rows, :cols]# Create a circular mask
    dist_from_center = np.sqrt((X - center_col)**2 + (Y - center_row)**2)
    mask = dist_from_center < radius    
    act_smoothed_array_mask = np.where(mask, act_smoothed_array, np.nan)# Apply mask: set values outside circle to NaN
    Smoothed_Act[nr]= act_smoothed_array_mask

    time_smoothed_array_mask = np.where(mask, time_smoothed_array, np.nan)# Apply mask: set values outside circle to NaN
    Smoothed_Time[nr]= time_smoothed_array_mask


    # --- Spatial information ---
    mean_act=np.nanmean(Smoothed_Act[nr])
    sum_time=np.nansum(Smoothed_Time[nr])
    Spatial_info[nr] = 0
    for rows in np.arange(np.shape(Smoothed_Act[nr])[0]) : 
        for cols in np.arange(np.shape(Smoothed_Act[nr])[1]) :             
            bin_act=Smoothed_Act[nr][rows,cols]
            p_bin=Smoothed_Time[nr][rows,cols]/sum_time 
            if ~ np.isnan(bin_act): 
                if bin_act!=0:             
                    spatial_info_bin=p_bin*(bin_act/mean_act)*math.log2(bin_act/mean_act)
                else:
                    spatial_info_bin=0    
                Spatial_info[nr] += spatial_info_bin

    if plot: 
        # Mask outside circle 
        rows, cols = np.full((max_row, max_col),1).shape
        center_row, center_col = rows // 2, cols // 2
        radius = min(rows, cols) // 2   
        Y, X = np.ogrid[:rows, :cols]# Create a circular mask
        dist_from_center = np.sqrt((X - center_col)**2 + (Y - center_row)**2)
        mask = dist_from_center < radius    
        masked_array = np.where(mask, np.full((max_row, max_col),1), np.nan)# Apply mask: set values outside circle to NaN
        
        # Draw map
        axs[nsubplot].imshow(masked_array, cmap='Greys_r', origin='upper')
        axs[nsubplot].imshow(act_smoothed_array_mask, cmap='jet', origin='upper')
        axs[nsubplot].set_aspect('equal')
        axs[nsubplot].set_title(f'{nr} ({np.round(Spatial_info[nr], 1)})', pad=0, loc='left', fontsize=max(7,int(nb_subplot/10)))

            
        # Remove box (spines)
        for spine in axs[nsubplot].spines.values():
            spine.set_visible(False)
        axs[nsubplot].set_xticks([])  # No x-axis ticks
        axs[nsubplot].set_yticks([])  # No y-axis ticks


if plot: 
    # Adjust layout to avoid clipping
    fig.subplots_adjust(wspace=.1, hspace=.2)
    plt.show()

with open(os.path.join(folder_to_save, f"SpatialInfo.pkl"), "wb") as f:
    pickle.dump(Spatial_info, f)


sorted_keys = sorted(AllCellDictSess_HD.keys())
AllCellDictSess_HD = {key: AllCellDictSess_HD[key] for rank, key in enumerate(sorted_keys)}

if plot: 
    nb_subplot=len(AllCellDictSess_HD.keys())
    rows = int(np.ceil(np.sqrt(nb_subplot)))  # Rows: ceil(sqrt(X))
    cols = int(np.ceil(np.sqrt(nb_subplot)))  # Columns: floor(sqrt(X))

    # Create the figure with a 2x2 grid
    fig, axs = plt.subplots(rows, cols, figsize=(int(nb_subplot/10), int(nb_subplot/10)),  subplot_kw={'projection': 'polar'})
    axs = axs.flatten()
    plt.tight_layout()

bin_size = 10
bin_edges = np.arange(0, 361, bin_size)
bin_centers = np.radians(bin_edges[:-1] + bin_size/2)

HeadDirection_infoSess={}
for nsubplot, nr in enumerate(AllCellDictSess_HD.keys()):

    HeadDirection_infoSess[nr] = np.divide(AllCellDictSess_HD[nr], AllTimeDictSess_HD[nr])
    
    if plot: 
        ax = axs[nsubplot]
        ax.set_theta_zero_location("E")  # 0° at East
        ax.set_theta_direction(-1)       # Clockwise
        ax.set_title(f"{nr.split('_')[0]}_{nr.split('_')[1]}", pad=0, loc='left', fontsize=max(7,int(nb_subplot/10)))
        
        #ax.bar(bin_centers, HeadDirection_infoSess[nr], width=np.radians(bin_size), alpha=0.6, color= 'k')   
                
        theta = bin_centers
        theta = np.append(theta, theta[0])   
        r = HeadDirection_infoSess[nr]
        r_smooth = gaussian_filter1d(r, sigma=1, mode='wrap')  # wrap for circular smoothing
        r_smooth = np.append(r_smooth, r_smooth[0])
        pref_idx = np.argmax(r_smooth)    # Find preferred angle
        pref_angle = theta[pref_idx]    
        color = plt.cm.hsv(pref_angle / (2*np.pi))# pref_angle ranges 0 to 2*pi -> normalize to 0-1 for colormap
        ax.plot(theta, r_smooth, color=color)
        
        #r = np.append(r, r[0])
        #ax.plot(theta, r, linestyle='-', alpha=0.7, color= 'r')

        ax.set_xticks([])
        ax.set_yticks([])

if plot: 
    # Adjust layout to avoid clipping
    fig.subplots_adjust(wspace=.1, hspace=.2)
    plt.show()


sorted_keys = sorted(SumActDict_HD.keys())
SumActDict_HD = {key: SumActDict_HD[key]for rank, key in enumerate(sorted_keys)}

if plot: 
    nb_subplot=len(SumActDict_HD.keys())
    rows = int(np.ceil(np.sqrt(nb_subplot)))  # Rows: ceil(sqrt(X))
    cols = int(np.ceil(np.sqrt(nb_subplot)))  # Columns: floor(sqrt(X))

    # Create the figure with a 2x2 grid
    fig, axs = plt.subplots(rows, cols, figsize=(int(nb_subplot/10), int(nb_subplot/10)),  subplot_kw={'projection': 'polar'})
    axs = axs.flatten()
    plt.tight_layout()

bin_size = 10
bin_edges = np.arange(0, 361, bin_size)
bin_centers = np.radians(bin_edges[:-1] + bin_size/2)

HeadDirection_info={}
for nsubplot, nr in enumerate(SumActDict_HD.keys()):

    HeadDirection_info[nr] = np.divide(SumActDict_HD[nr], SumTimeDict_HD[nr])
    
    if plot: 
        ax = axs[nsubplot]
        ax.set_theta_zero_location("E")  # 0° at East
        ax.set_theta_direction(-1)        # Clockwise
        ax.set_title(f'{nr}', pad=0, loc='left', fontsize=max(7,int(nb_subplot/10)))

        #ax.bar(bin_centers, HeadDirection_info[nr], width=np.radians(bin_size), alpha=0.6, color= 'k')   
        
        theta = bin_centers
        theta = np.append(theta, theta[0])  
        r = HeadDirection_info[nr]
        r_smooth = gaussian_filter1d(r, sigma=1, mode='wrap')  # wrap for circular smoothing
        r_smooth = np.append(r_smooth, r_smooth[0])
        pref_idx = np.argmax(r_smooth)    # Find preferred angle
        pref_angle = theta[pref_idx]    
        color = plt.cm.hsv(pref_angle / (2*np.pi))# pref_angle ranges 0 to 2*pi -> normalize to 0-1 for colormap
        ax.plot(theta, r_smooth, color=color)
        
        #r = np.append(r, r[0])
        #ax.plot(theta, r, linestyle='-', alpha=0.7, color= 'r')
        
        ax.set_xticks([])
        ax.set_yticks([])

if plot: 
    fig.subplots_adjust(wspace=.1, hspace=.2)
    plt.show()


with open(os.path.join(folder_to_save, f"HeadDirectionInfo.pkl"), "wb") as f:
    pickle.dump(HeadDirection_info, f)


sorted_keys = sorted(AllCellDictSess_Roll.keys())
AllCellDictSess_Roll = {key: AllCellDictSess_Roll[key] for rank, key in enumerate(sorted_keys)}

if plot: 
    nb_subplot=len(AllCellDictSess_Roll.keys())
    rows = int(np.ceil(np.sqrt(nb_subplot)))  # Rows: ceil(sqrt(X))
    cols = int(np.ceil(np.sqrt(nb_subplot)))  # Columns: floor(sqrt(X))

    # Create the figure with a 2x2 grid
    fig, axs = plt.subplots(rows, cols, figsize=(int(nb_subplot/10), int(nb_subplot/10)),  subplot_kw={'projection': 'polar'})
    axs = axs.flatten()
    plt.tight_layout()

bin_size = 10
bin_edges = np.arange(0, 361, bin_size)
bin_centers = np.radians(bin_edges[:-1] + bin_size/2)

Roll_infoSess={}
for nsubplot, nr in enumerate(AllCellDictSess_Roll.keys()):

    Roll_infoSess[nr] = np.divide(AllCellDictSess_Roll[nr], AllTimeDictSess_Roll[nr])
    
    if plot: 
        ax = axs[nsubplot]
        ax.set_theta_zero_location("E")  # 0° at East
        ax.set_theta_direction(-1)       # Clockwise
        ax.set_title(f"{nr.split('_')[0]}_{nr.split('_')[1]}", pad=0, loc='left', fontsize=max(7,int(nb_subplot/10)))
        
        #ax.bar(bin_centers, HeadDirection_infoSess[nr], width=np.radians(bin_size), alpha=0.6, color= 'k')   
                
        theta = bin_centers
        theta = np.append(theta, theta[0])   
        r = Roll_infoSess[nr]
        r_smooth = gaussian_filter1d(r, sigma=1, mode='wrap')  # wrap for circular smoothing
        r_smooth = np.append(r_smooth, r_smooth[0])
        pref_idx = np.argmax(r_smooth)    # Find preferred angle
        pref_angle = theta[pref_idx]    
        color = plt.cm.hsv(pref_angle / (2*np.pi))# pref_angle ranges 0 to 2*pi -> normalize to 0-1 for colormap
        ax.plot(theta, r_smooth, color=color)
        
        #r = np.append(r, r[0])
        #ax.plot(theta, r, linestyle='-', alpha=0.7, color= 'r')

        ax.set_xticks([])
        ax.set_yticks([])

if plot: 
    # Adjust layout to avoid clipping
    fig.subplots_adjust(wspace=.1, hspace=.2)
    plt.show()

sorted_keys = sorted(SumActDict_Roll.keys())
SumActDict_Roll = {key: SumActDict_Roll[key]for rank, key in enumerate(sorted_keys)}

if plot: 
    nb_subplot=len(SumActDict_Roll.keys())
    rows = int(np.ceil(np.sqrt(nb_subplot)))  # Rows: ceil(sqrt(X))
    cols = int(np.ceil(np.sqrt(nb_subplot)))  # Columns: floor(sqrt(X))

    # Create the figure with a 2x2 grid
    fig, axs = plt.subplots(rows, cols, figsize=(int(nb_subplot/10), int(nb_subplot/10)),  subplot_kw={'projection': 'polar'})
    axs = axs.flatten()
    plt.tight_layout()

bin_size = 10
bin_edges = np.arange(0, 361, bin_size)
bin_centers = np.radians(bin_edges[:-1] + bin_size/2)

Roll_info={}
for nsubplot, nr in enumerate(SumActDict_Roll.keys()):

    Roll_info[nr] = np.divide(SumActDict_Roll[nr], SumTimeDict_Roll[nr])
    
    if plot: 
        ax = axs[nsubplot]
        ax.set_theta_zero_location("E")  # 0° at East
        ax.set_theta_direction(-1)        # Clockwise
        ax.set_title(f'{nr}', pad=0, loc='left', fontsize=max(7,int(nb_subplot/10)))

        #ax.bar(bin_centers, HRoll_info[nr], width=np.radians(bin_size), alpha=0.6, color= 'k')   
        
        theta = bin_centers
        theta = np.append(theta, theta[0])  
        r = Roll_info[nr]
        r_smooth = gaussian_filter1d(r, sigma=1, mode='wrap')  # wrap for circular smoothing
        r_smooth = np.append(r_smooth, r_smooth[0])
        pref_idx = np.argmax(r_smooth)    # Find preferred angle
        pref_angle = theta[pref_idx]    
        color = plt.cm.hsv(pref_angle / (2*np.pi))# pref_angle ranges 0 to 2*pi -> normalize to 0-1 for colormap
        ax.plot(theta, r_smooth, color=color)
        
        #r = np.append(r, r[0])
        #ax.plot(theta, r, linestyle='-', alpha=0.7, color= 'r')
        
        ax.set_xticks([])
        ax.set_yticks([])

if plot: 
    fig.subplots_adjust(wspace=.1, hspace=.2)
    plt.show()

with open(os.path.join(folder_to_save, f"RollInfo.pkl"), "wb") as f:
    pickle.dump(Roll_info, f)


sorted_keys = sorted(AllCellDictSess_Pitch.keys())
AllCellDictSess_Pitch = {key: AllCellDictSess_Pitch[key] for rank, key in enumerate(sorted_keys)}

if plot: 
    nb_subplot=len(AllCellDictSess_Pitch.keys())
    rows = int(np.ceil(np.sqrt(nb_subplot)))  # Rows: ceil(sqrt(X))
    cols = int(np.ceil(np.sqrt(nb_subplot)))  # Columns: floor(sqrt(X))

    # Create the figure with a 2x2 grid
    fig, axs = plt.subplots(rows, cols, figsize=(int(nb_subplot/10), int(nb_subplot/10)),  subplot_kw={'projection': 'polar'})
    axs = axs.flatten()
    plt.tight_layout()

bin_size = 10
bin_edges = np.arange(0, 361, bin_size)
bin_centers = np.radians(bin_edges[:-1] + bin_size/2)

Pitch_infoSess={}
for nsubplot, nr in enumerate(AllCellDictSess_Pitch.keys()):

    Pitch_infoSess[nr] = np.divide(AllCellDictSess_Pitch[nr], AllTimeDictSess_Pitch[nr])
    
    if plot: 
        ax = axs[nsubplot]
        ax.set_theta_zero_location("E")  # 0° at East
        ax.set_theta_direction(-1)       # Clockwise
        ax.set_title(f"{nr.split('_')[0]}_{nr.split('_')[1]}", pad=0, loc='left', fontsize=max(7,int(nb_subplot/10)))
        
        #ax.bar(bin_centers, HeadDirection_infoSess[nr], width=np.radians(bin_size), alpha=0.6, color= 'k')   
                
        theta = bin_centers
        theta = np.append(theta, theta[0])   
        r = Pitch_infoSess[nr]
        r_smooth = gaussian_filter1d(r, sigma=1, mode='wrap')  # wrap for circular smoothing
        r_smooth = np.append(r_smooth, r_smooth[0])
        pref_idx = np.argmax(r_smooth)    # Find preferred angle
        pref_angle = theta[pref_idx]    
        color = plt.cm.hsv(pref_angle / (2*np.pi))# pref_angle ranges 0 to 2*pi -> normalize to 0-1 for colormap
        ax.plot(theta, r_smooth, color=color)
        
        #r = np.append(r, r[0])
        #ax.plot(theta, r, linestyle='-', alpha=0.7, color= 'r')

        ax.set_xticks([])
        ax.set_yticks([])

if plot: 
    # Adjust layout to avoid clipping
    fig.subplots_adjust(wspace=.1, hspace=.2)
    plt.show()


sorted_keys = sorted(SumActDict_Pitch.keys())
SumActDict_Pitch = {key: SumActDict_Pitch[key]for rank, key in enumerate(sorted_keys)}

if plot: 
    nb_subplot=len(SumActDict_Pitch.keys())
    rows = int(np.ceil(np.sqrt(nb_subplot)))  # Rows: ceil(sqrt(X))
    cols = int(np.ceil(np.sqrt(nb_subplot)))  # Columns: floor(sqrt(X))

    # Create the figure with a 2x2 grid
    fig, axs = plt.subplots(rows, cols, figsize=(int(nb_subplot/10), int(nb_subplot/10)),  subplot_kw={'projection': 'polar'})
    axs = axs.flatten()
    plt.tight_layout()

bin_size = 10
bin_edges = np.arange(0, 361, bin_size)
bin_centers = np.radians(bin_edges[:-1] + bin_size/2)

Pitch_info={}
for nsubplot, nr in enumerate(SumActDict_Pitch.keys()):

    Pitch_info[nr] = np.divide(SumActDict_Pitch[nr], SumTimeDict_Pitch[nr])
    
    if plot: 
        ax = axs[nsubplot]
        ax.set_theta_zero_location("E")  # 0° at East
        ax.set_theta_direction(-1)        # Clockwise
        ax.set_title(f'{nr}', pad=0, loc='left', fontsize=max(7,int(nb_subplot/10)))

        #ax.bar(bin_centers, Pitch_info[nr], width=np.radians(bin_size), alpha=0.6, color= 'k')   
        
        theta = bin_centers
        theta = np.append(theta, theta[0])  
        r = Pitch_info[nr]
        r_smooth = gaussian_filter1d(r, sigma=1, mode='wrap')  # wrap for circular smoothing
        r_smooth = np.append(r_smooth, r_smooth[0])
        pref_idx = np.argmax(r_smooth)    # Find preferred angle
        pref_angle = theta[pref_idx]    
        color = plt.cm.hsv(pref_angle / (2*np.pi))# pref_angle ranges 0 to 2*pi -> normalize to 0-1 for colormap
        ax.plot(theta, r_smooth, color=color)
        
        #r = np.append(r, r[0])
        #ax.plot(theta, r, linestyle='-', alpha=0.7, color= 'r')
        
        ax.set_xticks([])
        ax.set_yticks([])

if plot: 
    fig.subplots_adjust(wspace=.1, hspace=.2)
    plt.show()

with open(os.path.join(folder_to_save, f"PitchInfo.pkl"), "wb") as f:
    pickle.dump(Pitch_info, f)


sorted_keys = sorted(AllCellDictSess_Yaw.keys())
AllCellDictSess_Yaw = {key: AllCellDictSess_Yaw[key] for rank, key in enumerate(sorted_keys)}

if plot: 
    nb_subplot=len(AllCellDictSess_Yaw.keys())
    rows = int(np.ceil(np.sqrt(nb_subplot)))  # Rows: ceil(sqrt(X))
    cols = int(np.ceil(np.sqrt(nb_subplot)))  # Columns: floor(sqrt(X))

    # Create the figure with a 2x2 grid
    fig, axs = plt.subplots(rows, cols, figsize=(int(nb_subplot/10), int(nb_subplot/10)),  subplot_kw={'projection': 'polar'})
    axs = axs.flatten()
    plt.tight_layout()

bin_size = 10
bin_edges = np.arange(0, 361, bin_size)
bin_centers = np.radians(bin_edges[:-1] + bin_size/2)

Yaw_infoSess={}
for nsubplot, nr in enumerate(AllCellDictSess_Yaw.keys()):

    Yaw_infoSess[nr] = np.divide(AllCellDictSess_Yaw[nr], AllTimeDictSess_Yaw[nr])
    
    if plot: 
        ax = axs[nsubplot]
        ax.set_theta_zero_location("E")  # 0° at East
        ax.set_theta_direction(-1)       # Clockwise
        ax.set_title(f"{nr.split('_')[0]}_{nr.split('_')[1]}", pad=0, loc='left', fontsize=max(7,int(nb_subplot/10)))
        
        #ax.bar(bin_centers, HeadDirection_infoSess[nr], width=np.radians(bin_size), alpha=0.6, color= 'k')   
                
        theta = bin_centers
        theta = np.append(theta, theta[0])   
        r = Yaw_infoSess[nr]
        r_smooth = gaussian_filter1d(r, sigma=1, mode='wrap')  # wrap for circular smoothing
        r_smooth = np.append(r_smooth, r_smooth[0])
        pref_idx = np.argmax(r_smooth)    # Find preferred angle
        pref_angle = theta[pref_idx]    
        color = plt.cm.hsv(pref_angle / (2*np.pi))# pref_angle ranges 0 to 2*pi -> normalize to 0-1 for colormap
        ax.plot(theta, r_smooth, color=color)
        
        #r = np.append(r, r[0])
        #ax.plot(theta, r, linestyle='-', alpha=0.7, color= 'r')

        ax.set_xticks([])
        ax.set_yticks([])

if plot: 
    # Adjust layout to avoid clipping
    fig.subplots_adjust(wspace=.1, hspace=.2)
    plt.show()


sorted_keys = sorted(SumActDict_Yaw.keys())
SumActDict_Yaw = {key: SumActDict_Yaw[key]for rank, key in enumerate(sorted_keys)}

if plot: 
    nb_subplot=len(SumActDict_Yaw.keys())
    rows = int(np.ceil(np.sqrt(nb_subplot)))  # Rows: ceil(sqrt(X))
    cols = int(np.ceil(np.sqrt(nb_subplot)))  # Columns: floor(sqrt(X))

    # Create the figure with a 2x2 grid
    fig, axs = plt.subplots(rows, cols, figsize=(int(nb_subplot/10), int(nb_subplot/10)),  subplot_kw={'projection': 'polar'})
    axs = axs.flatten()
    plt.tight_layout()

bin_size = 10
bin_edges = np.arange(0, 361, bin_size)
bin_centers = np.radians(bin_edges[:-1] + bin_size/2)

Yaw_info={}
for nsubplot, nr in enumerate(SumActDict_Yaw.keys()):

    Yaw_info[nr] = np.divide(SumActDict_Yaw[nr], SumTimeDict_Yaw[nr])
    
    if plot: 
        ax = axs[nsubplot]
        ax.set_theta_zero_location("E")  # 0° at East
        ax.set_theta_direction(-1)        # Clockwise
        ax.set_title(f'{nr}', pad=0, loc='left', fontsize=max(7,int(nb_subplot/10)))

        #ax.bar(bin_centers, Yaw_info[nr], width=np.radians(bin_size), alpha=0.6, color= 'k')   
        
        theta = bin_centers
        theta = np.append(theta, theta[0])  
        r = Yaw_info[nr]
        r_smooth = gaussian_filter1d(r, sigma=1, mode='wrap')  # wrap for circular smoothing
        r_smooth = np.append(r_smooth, r_smooth[0])
        pref_idx = np.argmax(r_smooth)    # Find preferred angle
        pref_angle = theta[pref_idx]    
        color = plt.cm.hsv(pref_angle / (2*np.pi))# pref_angle ranges 0 to 2*pi -> normalize to 0-1 for colormap
        ax.plot(theta, r_smooth, color=color)
        
        #r = np.append(r, r[0])
        #ax.plot(theta, r, linestyle='-', alpha=0.7, color= 'r')
        
        ax.set_xticks([])
        ax.set_yticks([])

if plot: 
    fig.subplots_adjust(wspace=.1, hspace=.2)
    plt.show()

with open(os.path.join(folder_to_save, f"YawInfo.pkl"), "wb") as f:
    pickle.dump(Yaw_info, f)


sorted_keys = sorted(AllCellDictSess_Speed.keys())
AllCellDictSess_Speed = {key: AllCellDictSess_Speed[key] for rank, key in enumerate(sorted_keys)}

if plot: 
    nb_subplot=len(AllCellDictSess_Speed.keys())
    rows = int(np.ceil(np.sqrt(nb_subplot)))  # Rows: ceil(sqrt(X))
    cols = int(np.ceil(np.sqrt(nb_subplot)))  # Columns: floor(sqrt(X))

    # Create the figure with a 2x2 grid
    fig, axs = plt.subplots(rows, cols, figsize=(int(nb_subplot/10), int(nb_subplot/10)),  subplot_kw={'projection': 'polar'})
    axs = axs.flatten()
    plt.tight_layout()

bin_size = 1
bin_edges = np.arange(0, 20, bin_size)
bin_centers = bin_edges[:-1] + bin_size/2


Speed_infoSess={}
for nsubplot, nr in enumerate(AllCellDictSess_Speed.keys()):

    Speed_infoSess[nr] = np.divide(AllCellDictSess_Speed[nr], AllTimeDictSess_Speed[nr])
    
    if plot: 
        ax = axs[nsubplot]        
        ax.set_title(f"{nr.split('_')[0]}_{nr.split('_')[1]}", pad=0, loc='left', fontsize=max(7,int(nb_subplot/10)))        
        ax.bar(bin_centers, HeadDirection_infoSess[nr], width=bin_size, alpha=0.6, color= 'k')   
        
        ax.set_xticks([])
        ax.set_yticks([])

if plot: 
    # Adjust layout to avoid clipping
    fig.subplots_adjust(wspace=.1, hspace=.2)
    plt.show()


sorted_keys = sorted(SumActDict_Speed.keys())
SumActDict_Speed = {key: SumActDict_Speed[key]for rank, key in enumerate(sorted_keys)}

if plot: 
    nb_subplot=len(SumActDict_Speed.keys())
    rows = int(np.ceil(np.sqrt(nb_subplot)))  # Rows: ceil(sqrt(X))
    cols = int(np.ceil(np.sqrt(nb_subplot)))  # Columns: floor(sqrt(X))

    # Create the figure with a 2x2 grid
    fig, axs = plt.subplots(rows, cols, figsize=(int(nb_subplot/10), int(nb_subplot/10)),  subplot_kw={'projection': 'polar'})
    axs = axs.flatten()
    plt.tight_layout()

bin_size = 1
bin_edges = np.arange(0, 20, bin_size)
bin_centers = bin_edges[:-1] + bin_size/2

Speed_info={}
for nsubplot, nr in enumerate(SumActDict_Speed.keys()):

    Speed_info[nr] = np.divide(SumActDict_Speed[nr], SumTimeDict_Speed[nr])
    
    if plot: 
        ax = axs[nsubplot]
        ax.set_title(f'{nr}', pad=0, loc='left', fontsize=max(7,int(nb_subplot/10)))
        ax.bar(bin_centers, Speed_info[nr], width=bin_size, alpha=0.6, color= 'k')   
        ax.set_xticks([])
        ax.set_yticks([])

if plot: 
    fig.subplots_adjust(wspace=.1, hspace=.2)
    plt.show()


with open(os.path.join(folder_to_save, f"SpeedInfo.pkl"), "wb") as f:
    pickle.dump(Speed_info, f)


ShuffledSpatial_info={}
ShuffledHeadDirection_info={}

ShuffledRoll_info={}
ShuffledPitch_info={}
ShuffledYaw_info={}
ShuffledSpeed_info={}

for nr in SumActDict_PC.keys():
    ShuffledSpatial_info[nr]={}
    ShuffledHeadDirection_info[nr]={}
    ShuffledRoll_info[nr]={}
    ShuffledPitch_info[nr]={}
    ShuffledYaw_info[nr]={}
    ShuffledSpeed_info[nr]={}

nb_tot_permutation=1000

for permutation_nb in range(nb_tot_permutation):

    AllCellDict_PC_Sh={}
    AllTimeDict_PC_Sh={}
    AllCellDictSess_PC_Sh={}
    AllTimeDictSess_PC_Sh={}

    AllCellDict_HD_Sh={}
    AllTimeDict_HD_Sh={}
    AllCellDictSess_HD_Sh={}
    AllTimeDictSess_HD_Sh={}

    AllCellDict_Roll_Sh={}
    AllTimeDict_Roll_Sh={}
    AllCellDictSess_Roll_Sh={}
    AllTimeDictSess_Roll_Sh={}

    AllCellDict_Pitch_Sh={}
    AllTimeDict_Pitch_Sh={}
    AllCellDictSess_Pitch_Sh={}
    AllTimeDictSess_Pitch_Sh={}

    AllCellDict_Yaw_Sh={}
    AllTimeDict_Yaw_Sh={}
    AllCellDictSess_Yaw_Sh={}
    AllTimeDictSess_Yaw_Sh={}

    AllCellDict_Speed_Sh={}
    AllTimeDict_Speed_Sh={}
    AllCellDictSess_Speed_Sh={}
    AllTimeDictSess_Speed_Sh={}

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
                    session_type=minianpath.parents[3].name
                    session_date=minianpath.parents[4].name
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

                if V4subfolder:
                    V4subfolder_id = int(minianpath.parent.name[-1]) - 1 
                    ts_start = V4subfolder_id*15*1000 # cause 15 videos per subfolders and 1000 frames per videos
                    ts_stop = np.shape(Co)[1] + ts_start
                    tsmini_sub=tsmini[ts_start:ts_stop]
                    tsmini_sub=tsmini_sub.reset_index(drop=True)    
                else:
                    tsmini_sub=tsmini


                # Extraction des quaternions

                vestibular_df= pd.read_csv(list(session_path.glob('*V4_Miniscope/headOrientation.csv'))[0])
                qw = vestibular_df['qw'].to_numpy()
                qx = vestibular_df['qx'].to_numpy()
                qy = vestibular_df['qy'].to_numpy()
                qz = vestibular_df['qz'].to_numpy()
                roll, pitch, yaw = quaternion_to_euler(qw, qx, qy, qz)
                roll_deg = np.degrees(roll) # nose goes clockwise or anticlockwise
                pitch_deg = np.degrees(pitch) # nose up or down
                yaw_deg = np.degrees(yaw) # nose moves from side to side                
                

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
                nose_x3, nose_y3, speed = limit_speed(nose_x2, nose_y2, max_speed=max_speed)
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
                neck_x3, neck_y3, speed = limit_speed(neck_x2, neck_y2, max_speed=max_speed)
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

                speed= np.concatenate(([0], speed))
                speed_= speed[start_frame:last_frame]

                closest_start= find_closest_index_sorted(tsmini_sub, tswebcam[start_frame])
                closest_end= find_closest_index_sorted(tsmini_sub, tswebcam[last_frame-1])

                roll_deg_ = roll_deg[closest_start:closest_end]
                pitch_deg_ = pitch_deg[closest_start:closest_end]
                yaw_deg_ = yaw_deg[closest_start:closest_end]


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

                        for idx, (px, py, angle, _roll, _pitch, _yaw, _speed) in enumerate(zip(nose_x_, nose_y_, angles_deg, roll_deg_, pitch_deg_, yaw_deg_, speed_)):
                            if np.sqrt((px - table_center_x)**2 + (py - table_center_y)**2) > table_radius:
                                continue  # skip points outside circle
                            closest_point = find_closest_index_sorted(tsmini_sub, tswebcam[start_frame+idx])

                            ix = int(np.floor((px - (table_center_x - n/2 * square_size)) / square_size))
                            iy = int(np.floor((py - (table_center_y - n/2 * square_size)) / square_size))
                            nr_tot_act_biaised_PC[(ix, iy)] += Carray_unit[closest_point]
                            counts_PC[(ix, iy)] += 1/frame_rate 

                            bin_idx = int(angle // bin_size)
                            nr_tot_act_biaised_HD[bin_idx] += Carray_unit[closest_point]
                            counts_HD[bin_idx] += 1    

                            rollbin_idx = int(_roll // bin_size)
                            nr_tot_act_biaised_Roll[rollbin_idx] += Carray_unit[closest_point]
                            counts_Roll[rollbin_idx] += 1      

                            pitchbin_idx = int(_pitch // bin_size)
                            nr_tot_act_biaised_Pitch[pitchbin_idx] += Carray_unit[closest_point]
                            counts_Pitch[pitchbin_idx] += 1

                            yawbin_idx = int(_yaw // bin_size)
                            nr_tot_act_biaised_Yaw[yawbin_idx] += Carray_unit[closest_point]
                            counts_Yaw[yawbin_idx] += 1  
                            
                            speedbin_idx = min(int(_speed), 19) # between 0 and 20 
                            nr_tot_act_biaised_Speed[speedbin_idx] += Carray_unit[closest_point]
                            counts_Speed[speedbin_idx] += 1  


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

                        AllCellDictSess_PC_Sh[sess_unit] = nr_tot_act_biaised_PC 
                        AllTimeDictSess_PC_Sh[sess_unit] = counts_PC

                        AllCellDictSess_HD_Sh[sess_unit] = nr_tot_act_biaised_HD 
                        AllTimeDictSess_HD_Sh[sess_unit] = counts_HD

                        AllCellDictSess_Roll_Sh[sess_unit] = nr_tot_act_biaised_Roll 
                        AllTimeDictSess_Roll_Sh[sess_unit] = counts_Roll

                        AllCellDictSess_Pitch_Sh[sess_unit] = nr_tot_act_biaised_Pitch 
                        AllTimeDictSess_Pitch_Sh[sess_unit] = counts_Pitch

                        AllCellDictSess_Yaw_Sh[sess_unit] = nr_tot_act_biaised_Yaw 
                        AllTimeDictSess_Yaw_Sh[sess_unit] = counts_Yaw

                        AllCellDictSess_Speed_Sh[sess_unit] = nr_tot_act_biaised_Speed 
                        AllTimeDictSess_Speed_Sh[sess_unit] = counts_Speed
            
                sess+=1

    #########################################################################################
                    # Average activity map across sessions for each cell #
    ##########################################################################################

    SumActDict_PC_Sh={}
    for uniquecell in AllCellDict_PC_Sh.keys():
        sums = {}
        for subdict in AllCellDict_PC_Sh[uniquecell].values():
            for key, value in subdict.items():
                sums[key] = sums.get(key, 0) + value
        SumActDict_PC_Sh[uniquecell] = {key: sums[key] for key in sums}
    SumTimeDict_PC_Sh={}
    for uniquecell in AllTimeDict_PC_Sh.keys():
        sums = {}
        for subdict in AllTimeDict_PC_Sh[uniquecell].values():
            for key, value in subdict.items():
                sums[key] = sums.get(key, 0) + value
        SumTimeDict_PC_Sh[uniquecell] = {key: sums[key] for key in sums}
    
    SumActDict_HD_Sh={}
    for outer_key, inner_dict in AllCellDict_HD_Sh.items():
        first_array = next(iter(inner_dict.values()))
        sum_array = np.zeros_like(first_array)
        for arr in inner_dict.values():
            sum_array += arr
        SumActDict_HD_Sh[outer_key] = sum_array
    SumTimeDict_HD_Sh={}
    for outer_key, inner_dict in AllTimeDict_HD_Sh.items():
        first_array = next(iter(inner_dict.values()))
        sum_array = np.zeros_like(first_array)
        for arr in inner_dict.values():
            sum_array += arr
        SumTimeDict_HD_Sh[outer_key] = sum_array


    SumActDict_Roll_Sh={}
    for outer_key, inner_dict in AllCellDict_Roll_Sh.items():
        first_array = next(iter(inner_dict.values()))
        sum_array = np.zeros_like(first_array)
        for arr in inner_dict.values():
            sum_array += arr
        SumActDict_Roll_Sh[outer_key] = sum_array
    SumTimeDict_Roll_Sh={}
    for outer_key, inner_dict in AllTimeDict_Roll_Sh.items():
        first_array = next(iter(inner_dict.values()))
        sum_array = np.zeros_like(first_array)
        for arr in inner_dict.values():
            sum_array += arr
        SumTimeDict_Roll_Sh[outer_key] = sum_array

    SumActDict_Pitch_Sh={}
    for outer_key, inner_dict in AllCellDict_Pitch_Sh.items():
        first_array = next(iter(inner_dict.values()))
        sum_array = np.zeros_like(first_array)
        for arr in inner_dict.values():
            sum_array += arr
        SumActDict_Pitch_Sh[outer_key] = sum_array
    SumTimeDict_Pitch_Sh={}
    for outer_key, inner_dict in AllTimeDict_Pitch_Sh.items():
        first_array = next(iter(inner_dict.values()))
        sum_array = np.zeros_like(first_array)
        for arr in inner_dict.values():
            sum_array += arr
        SumTimeDict_Pitch_Sh[outer_key] = sum_array

    SumActDict_Yaw_Sh={}
    for outer_key, inner_dict in AllCellDict_Yaw_Sh.items():
        first_array = next(iter(inner_dict.values()))
        sum_array = np.zeros_like(first_array)
        for arr in inner_dict.values():
            sum_array += arr
        SumActDict_Yaw_Sh[outer_key] = sum_array
    SumTimeDict_Yaw_Sh={}
    for outer_key, inner_dict in AllTimeDict_Yaw_Sh.items():
        first_array = next(iter(inner_dict.values()))
        sum_array = np.zeros_like(first_array)
        for arr in inner_dict.values():
            sum_array += arr
        SumTimeDict_Yaw_Sh[outer_key] = sum_array

    SumActDict_Speed_Sh={}
    for outer_key, inner_dict in AllCellDict_Speed_Sh.items():
        first_array = next(iter(inner_dict.values()))
        sum_array = np.zeros_like(first_array)
        for arr in inner_dict.values():
            sum_array += arr
        SumActDict_Speed_Sh[outer_key] = sum_array
    SumTimeDict_Speed_Sh={}
    for outer_key, inner_dict in AllTimeDict_Speed_Sh.items():
        first_array = next(iter(inner_dict.values()))
        sum_array = np.zeros_like(first_array)
        for arr in inner_dict.values():
            sum_array += arr
        SumTimeDict_Speed_Sh[outer_key] = sum_array



    sorted_keys = sorted(SumActDict_PC_Sh.keys())
    SumActDict_PC_Sh = {key: SumActDict_PC_Sh[key]for rank, key in enumerate(sorted_keys)}

    sigma = 1  # Standard deviation for Gaussian kernel
    max_row = int(table_radius*2/square_size)+1
    max_col = int(table_radius*2/square_size)+1

    Smoothed_Act_PC_Sh={}
    Smoothed_Time_PC_Sh={}
    
    for nsubplot, nr in enumerate(SumActDict_PC_Sh.keys()):

        nr_tot_act = SumActDict_PC_Sh[nr]    
        nr_tot_act_array = [[None for _ in range(max_col)] for _ in range(max_row)] # Initialize 2D array
        for (row, col), value in nr_tot_act.items():# Fill array
            nr_tot_act_array[row][col] = value    
        nr_tot_act_array = [[np.nan if v is None else v for v in row] for row in nr_tot_act_array] # Replace None with np.nan
        nr_tot_act_array = np.array(nr_tot_act_array, dtype=float)
        nr_tot_act_array=nr_tot_act_array.T

        array_for_filter = np.nan_to_num(nr_tot_act_array, nan=0.0) # Replace nan with 0 for Gaussian filtering
        act_smoothed_array = gaussian_filter(array_for_filter, sigma=sigma)    # Apply 2D Gaussian filter


        nr_tot_time = SumTimeDict_PC_Sh[nr]
        nr_tot_time_array = [[None for _ in range(max_col)] for _ in range(max_row)] # Initialize 2D array
        for (row, col), value in nr_tot_time.items():# Fill array
            nr_tot_time_array[row][col] = value
        nr_tot_time_array = [[np.nan if v is None else v for v in row] for row in nr_tot_time_array] # Replace None with np.nan
        nr_tot_time_array = np.array(nr_tot_time_array, dtype=float)
        nr_tot_time_array=nr_tot_time_array.T

        array_for_filter = np.nan_to_num(nr_tot_time_array, nan=0.0) # Replace nan with 0 for Gaussian filtering
        time_smoothed_array = gaussian_filter(array_for_filter, sigma=sigma)    # Apply 2D Gaussian filter
        time_smoothed_array[time_smoothed_array < .5] = np.nan  # do not consider if spent less than 0.1 second in the square

        # Ratio
        act_smoothed_array = np.divide(act_smoothed_array, time_smoothed_array)

        # --- Mask outside circle ---
        rows, cols = act_smoothed_array.shape
        center_row, center_col = rows // 2, cols // 2
        radius = min(rows, cols) // 2   
        Y, X = np.ogrid[:rows, :cols]# Create a circular mask
        dist_from_center = np.sqrt((X - center_col)**2 + (Y - center_row)**2)
        mask = dist_from_center < radius    
        act_smoothed_array_mask = np.where(mask, act_smoothed_array, np.nan)# Apply mask: set values outside circle to NaN
        Smoothed_Act_PC_Sh[nr]= act_smoothed_array_mask

        time_smoothed_array_mask = np.where(mask, time_smoothed_array, np.nan)# Apply mask: set values outside circle to NaN
        Smoothed_Time_PC_Sh[nr]= time_smoothed_array_mask


        # --- Spatial information ---
        mean_act=np.nanmean(Smoothed_Act_PC_Sh[nr])
        sum_time=np.nansum(Smoothed_Time_PC_Sh[nr])

        ShuffledSpatial_info[nr][permutation_nb] = 0

        for rows in np.arange(np.shape(Smoothed_Act_PC_Sh[nr])[0]) : 
            for cols in np.arange(np.shape(Smoothed_Act_PC_Sh[nr])[1]) :             
                bin_act=Smoothed_Act_PC_Sh[nr][rows,cols]
                p_bin=Smoothed_Time_PC_Sh[nr][rows,cols]/sum_time 
                if ~ np.isnan(bin_act): 
                    if bin_act!=0:             
                        spatial_info_bin=p_bin*(bin_act/mean_act)*math.log2(bin_act/mean_act)
                    else:
                        spatial_info_bin=0 

                    ShuffledSpatial_info[nr][permutation_nb] += spatial_info_bin
        
        ShuffledHeadDirection_info[nr][permutation_nb] = np.divide(SumActDict_HD_Sh[nr], SumTimeDict_HD_Sh[nr])    
        ShuffledRoll_info[nr][permutation_nb] = np.divide(SumActDict_Roll_Sh[nr], SumTimeDict_Roll_Sh[nr])    
        ShuffledPitch_info[nr][permutation_nb] = np.divide(SumActDict_Pitch_Sh[nr], SumTimeDict_Pitch_Sh[nr])    
        ShuffledYaw_info[nr][permutation_nb] = np.divide(SumActDict_Yaw_Sh[nr], SumTimeDict_Yaw_Sh[nr])    
        ShuffledSpeed_info[nr][permutation_nb] = np.divide(SumActDict_Speed_Sh[nr], SumTimeDict_Speed_Sh[nr])    

    print(f'Permutation n° {permutation_nb} done')  


    if permutation_nb % 100 == 0:
        with open(os.path.join(folder_to_save, f"ShuffledSpatialInfo_{permutation_nb}.pkl"), "wb") as f:
            pickle.dump(ShuffledSpatial_info, f)
        with open(os.path.join(folder_to_save, f"ShuffledHeadDirection_{permutation_nb}.pkl"), "wb") as f:
            pickle.dump(ShuffledHeadDirection_info, f)
        with open(os.path.join(folder_to_save, f"ShuffledRoll_{permutation_nb}.pkl"), "wb") as f:
            pickle.dump(ShuffledRoll_info, f)
        with open(os.path.join(folder_to_save, f"ShuffledPitch_{permutation_nb}.pkl"), "wb") as f:
            pickle.dump(ShuffledPitch_info, f)
        with open(os.path.join(folder_to_save, f"ShuffledYaw_{permutation_nb}.pkl"), "wb") as f:
            pickle.dump(ShuffledYaw_info, f)
        with open(os.path.join(folder_to_save, f"ShuffledSpeed_{permutation_nb}.pkl"), "wb") as f:
            pickle.dump(ShuffledSpeed_info, f)

with open(os.path.join(folder_to_save, f"ShuffledSpatialInfo.pkl"), "wb") as f:
    pickle.dump(ShuffledSpatial_info, f)
with open(os.path.join(folder_to_save, f"ShuffledHeadDirection.pkl"), "wb") as f:
    pickle.dump(ShuffledHeadDirection_info, f)
with open(os.path.join(folder_to_save, f"ShuffledRoll.pkl"), "wb") as f:
    pickle.dump(ShuffledRoll_info, f)
with open(os.path.join(folder_to_save, f"ShuffledPitch.pkl"), "wb") as f:
    pickle.dump(ShuffledPitch_info, f)
with open(os.path.join(folder_to_save, f"ShuffledYaw.pkl"), "wb") as f:
    pickle.dump(ShuffledYaw_info, f)
with open(os.path.join(folder_to_save, f"ShuffledSpeed.pkl"), "wb") as f:
    pickle.dump(ShuffledSpeed_info, f)


