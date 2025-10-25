# # Associate Ca2+ signal with sleep stages for each session & subsessions using crossregistration

#######################################################################################
                            # Define Experiment type #
#######################################################################################

DrugExperiment=0 # =1 if CGP Experiment // DrugExperiment=0 if Baseline Experiment

AnalysisID='_' 

Cell_Assembly_folder= "CellAssemblies_2025-10-25_18_51_09"

local = True
if local:
    dir= "//10.69.168.1/crnldata/forgetting/Aurelie/MiniscopeOE_data/L2_3_mice/"
else: 
    dir= "/crnldata/forgetting/Aurelie/MiniscopeOE_data/L2_3_mice/"

mapp = {1: 'AW',  2: 'QW', 3: 'NREM',  4: 'IS', 5: 'REM',  6: 'undefined'}
drugs=['baseline']

#######################################################################################
                                # Load packages #
#######################################################################################

import os
import numpy as np
from scipy import signal
import quantities as pq
import math 
import neo
import json
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
import numpy as np
from scipy.stats import zscore

warnings.filterwarnings("ignore")
import sys
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


minian_path = os.path.join(os.path.abspath('..'),'minian')
print("The folder used for minian procedures is : {}".format(minian_path))
sys.path.append(minian_path)

from minian.utilities import (
    TaskAnnotation,
    get_optimal_chk,
    load_videos,
    open_minian,
    save_minian,
)
#######################################################################################
                                # Define functions #
#######################################################################################

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

def bin_sum_fractional(x, Fs1, Fs2):
    N = len(x)
    bin_width = Fs1 / Fs2
    y = []
    start = 0
    while start < N:
        end = start + bin_width
        i_start = int(np.floor(start))
        i_end = int(np.floor(end))
        if i_end >= N:
            i_end = N - 1
        total = 0.0
        total += x[i_start] * (1 - (start - i_start))
        for i in range(i_start + 1, i_end):
            total += x[i]
        if i_end > i_start:
            total += x[i_end] * (end - i_end)
        y.append(total)
        start = end
    return np.array(y)

#######################################################################################
                # Load sleep score and Ca2+ time series numpy arrays #
#######################################################################################

all_expe_types=['baseline','preCGP', 'postCGP'] if DrugExperiment else ['baseline', 'preCGP']

# Get the current date and time
FolderNameSave=str(datetime.now())[:19]
FolderNameSave = FolderNameSave.replace(" ", "_").replace(".", "_").replace(":", "_")

destination_folder= f"//10.69.168.1/crnldata/forgetting/Aurelie/MiniscopeOE_analysis/PlaceCells_experiment/Reactivation_{FolderNameSave}{AnalysisID}" if local else f"/crnldata/forgetting/Aurelie/MiniscopeOE_analysis/PlaceCells_experiment/Reactivation_{FolderNameSave}{AnalysisID}"
os.makedirs(destination_folder)
folder_to_save=Path(destination_folder)

logfile = open(f"{destination_folder}/output_log.txt", 'w')
sys.stdout = Tee(sys.stdout, logfile)  # print goes to both


# Copy the script file to the destination folder
source_script = "C:/Users/Manip2/SCRIPTS/HayLabAnalysis/python/_MINI&OE_3_estimate_reactivation_strenght.py" if local else "/home/aurelie.brecier/HayLabAnalysis/python/_MINI&OE_3_estimate_reactivation_strenght.py"
destination_file_path = f"{destination_folder}/_MINI&OE_3_estimate_reactivation_strenght.txt"
shutil.copy(source_script, destination_file_path)

data = {}
counter=0
CellAssembly_GlobalResults= pd.DataFrame(data, columns=['Mice','NeuronType','Session', 'Session_Date', 'Session_Time', 'Assembly_ID', 'Assembly_size', 
                                                        'Cells_in_Assembly','ExpeType', 'Drug', 'Substate', 'SubstateDuration', 'Session_ID', 
                                                        'Avg_ReactStr', 'EventFreq', 'EventTime' ])

filenameOut =  f'//10.69.168.1/crnldata/forgetting/Aurelie/MiniscopeOE_analysis/PlaceCells_experiment/{Cell_Assembly_folder}/CellAssembly_dict.pkl'
with open(filenameOut, 'rb') as pickle_file:
    CellAssembly_dict = pickle.load(pickle_file)

filenameOut =  f'//10.69.168.1/crnldata/forgetting/Aurelie/MiniscopeOE_analysis/PlaceCells_experiment/{Cell_Assembly_folder}/CellAssembly_members.pkl'
with open(filenameOut, 'rb') as pickle_file:
    CellAssembly_members = pickle.load(pickle_file)

for dpath in Path(dir).glob('**/PlaceCells_experiment/mappingsAB.pkl'):

    mappfile = open(dpath.parents[0]/ f'mappingsAB.pkl', 'rb')
    mapping = pickle.load(mappfile)
    mapping_sess = mapping['session']   

    centroidfile = open(dpath.parents[0]/ f'centsAB.pkl', 'rb')
    centroids = pickle.load(centroidfile)
    centroids_sess = centroids['session']   

    mice = dpath.parents[1].parts[-1]
    NeuronType = dpath.parents[2].parts[-1]
    
    print(f"####################################################################################")
    print(f"################################### {mice} ####################################")
    print(f"####################################################################################")

    nb_minian_total=0
    dict_Path={}
    dict_Calcium = {}
    dict_Deconv = {}
    dict_Scoring = {}
    dict_Stamps = {}
    dict_TodropFile = {}
    dict_StampsMiniscope = {}

    previousEndTime=0
    InitialStartTime=0

    minian_folders = [f for f in dpath.parents[0].rglob('minian') if f.is_dir()]

    for minianpath in minian_folders: # for each minian folders found in this mouse 
        drug= 'baseline' 

        if 1==1: # any(p in all_expe_types for p in minianpath.parts): # have to be to the expe_types

            if any(minianpath.parents[1].glob('*V4_Miniscope/timeStamps.csv')): 
                session_path=minianpath.parents[1]        
                tsmini= pd.read_csv(list(session_path.glob('*V4_Miniscope/timeStamps.csv'))[0])['Time Stamp (ms)']
                V4subfolder=False               
                if any(session_path.parent.rglob('OpenEphys/EphyViewer_ManualScor5s_6Stages.csv')):
                    session_type=minianpath.parents[3].name
                    session_date=minianpath.parents[4].name
                    session_time=minianpath.parents[1].name
                else:
                    session_type=minianpath.parents[2].name
                    session_date=minianpath.parents[3].name
                    session_time=minianpath.parents[1].name
                session=session_time        
            elif any(minianpath.parents[2].glob('*V4_Miniscope/timeStamps.csv')):
                session_path=minianpath.parents[2]      
                tsmini= pd.read_csv(list(session_path.glob('*V4_Miniscope/timeStamps.csv'))[0])['Time Stamp (ms)']                      
                V4subfolder=True
                session_type=minianpath.parents[4].name
                session_date=minianpath.parents[5].name
                session_time=minianpath.parents[0].name
                session=session_time
            print(f"Processing {session_type} session: {session} on the {session_date}, subfolders = {V4subfolder}")

            minian_ds = open_minian(minianpath)
            dict_Calcium[session] = minian_ds['C'] # calcium traces 
            dict_Deconv[session] = minian_ds['S'] # estimated spikes deconvolved activity
            
            try: 
                try: 
                    with open(minianpath.parent / f'TodropFileAB.json', 'r') as f:
                        unit_to_drop = json.load(f)
                        dict_TodropFile[session]  = unit_to_drop
                except: 
                    with open(minianpath/ f'TodropFileAB.json', 'r') as f:
                        unit_to_drop = json.load(f)
                        dict_TodropFile[session]  = unit_to_drop
                nb_minian_total+=1
            except:
                print(' !!!! Minian session not validated !!!!')
                continue
            print(session_path)

            
            # Start time & freq miniscope
            
            dict_StampsMiniscope[session]  = pd.read_csv(list(session_path.glob('*V4_Miniscope/timeStamps.csv'))[0])
            tsmini=dict_StampsMiniscope[session]['Time Stamp (ms)']
            minian_freq=round(1/np.mean(np.diff(np.array(tsmini)/1000)))  

            freqLFP=1000    
            numbdropfr= 0         
            
            nb_minian_total+=1

            #######################################################################################
            # Distribute Ca2+ intensity & spikes to vigilance states for each sessions #
            #######################################################################################
                
            # Adjust the StartTime if subsessions

            dict_Stamps[session]  = pd.read_excel(session_path/ f'SynchroFile.xlsx')     
            StartTime = list(dict_Stamps[session][0])[0] 
            if math.isnan(StartTime):   
                StartTime = 0 # No OpenEphys file found, so Miniscope start relative to OpenEphys start equal 0

            if InitialStartTime==0:
                InitialStartTime=StartTime    
                firstframe=0
                StartTimeMiniscope=0 # start time of miniscope rec of that subsesssions relative to the start of the mniscope recording
            else:
                if StartTime == InitialStartTime: # just a subsession
                    StartTime = previousEndTime + 1/minian_freq #  +1 frame in seconds
                    StartTimeMiniscope= StartTime-InitialStartTime
                else:  
                    InitialStartTime=StartTime # this is a new session
                    firstframe=0
                    StartTimeMiniscope=0 
                    

            # Remove bad units from recordings

            C=dict_Calcium[session]
            D=dict_Deconv[session] 
            unit_to_drop=dict_TodropFile[session]    

            C_sel=C.drop_sel(unit_id=unit_to_drop)
            D_sel=D.drop_sel(unit_id=unit_to_drop)

            CalciumSub = pd.DataFrame(C_sel, index=C_sel['unit_id'])
            DeconvSub = pd.DataFrame(D_sel, index=D_sel['unit_id'])

            indexMappList=mapping_sess[session]
            kept_uniq_unit_List=[]
            for unit in CalciumSub.index:
                indexMapp = np.where(indexMappList == unit)[0]
                kept_uniq_unit_List.append(f"{mice}{str(indexMapp).replace('[','').replace(']','')}")

            nb_unit=len(CalciumSub)
            if nb_unit==0:
                print(f'no cells kept in the session: {session}')
                continue  # next iteration
            

            # Realigned the traces to the recorded timestamps 

            timestamps =  np.array(tsmini[firstframe:firstframe+len(CalciumSub.T)])/freqLFP
            new_timestamps= np.arange(timestamps[0], timestamps[-1], 1/minian_freq)
            
            Calcium = pd.DataFrame(index=CalciumSub.index, columns=new_timestamps)
            for feature in CalciumSub.index:
                interpolator = interpolate.interp1d(timestamps, CalciumSub.loc[feature], kind='linear')
                Calcium.loc[feature] = interpolator(new_timestamps)
            Carray=Calcium.values.T.astype(float) # Calcium activity

            Deconv = pd.DataFrame(index=DeconvSub.index, columns=new_timestamps)
            for feature in DeconvSub.index:
                interpolator = interpolate.interp1d(timestamps, DeconvSub.loc[feature], kind='linear')
                Deconv.loc[feature] = interpolator(new_timestamps)
            Darray=Deconv.values.T.astype(float) # Deconvolved activity

            Sarray= np.zeros((np.shape(Darray)[0], np.shape(Darray)[1])) # Spike activity
            for i in np.arange(np.shape(Darray)[1]):
                Darray_unit =Darray[:,i]
                peaks, _ = find_peaks(Darray_unit)
                Sarray_unit=np.zeros(len(Darray_unit))
                Sarray_unit[peaks]=1   
                Sarray[:,i]=Sarray_unit
            Spike= pd.DataFrame(Sarray.T, index=CalciumSub.index, columns=new_timestamps)

            firstframe+=len(CalciumSub.T)
            rec_dur=len(CalciumSub.T)
            rec_dur_sec= timestamps[-1]- timestamps[0] 
            EndTime = StartTime + rec_dur_sec # (upd_rec_dur/minian_freq) # in seconds
            previousEndTime=EndTime 

            if any(session_path.parent.rglob('OpenEphys/EphyViewer_ManualScor5s_6Stages.csv')):
                dict_Scoring[session]  = pd.read_csv(session_path.parent/ f'OpenEphys/EphyViewer_ManualScor5s_6Stages.csv')
                               
                # Upscale scoring to miniscope frequency
                SleepScored=dict_Scoring[session]
                Bool= np.sum(SleepScored['duration'])/len(SleepScored) == 5

                if not Bool:
                    print('/!\ Scoring freq=', np.sum(SleepScored['duration'])/len(SleepScored), 'Hz !!!!')
                    continue 

                SleepScored['label']= SleepScored['label'].str.extract(r'(\d+)', expand=False)
                SleepScoredTS=np.array(SleepScored['label'])

                scale_factor=minian_freq/0.2  #cause scoring was done in 5 seconds bin, ie 0.2 Hz              
                SleepScoredTS_upscaled = np.repeat(SleepScoredTS, scale_factor, axis=0)
                StartTime_frame=round(StartTime*minian_freq)
                SleepScoredTS_upscaled_ministart=SleepScoredTS_upscaled[StartTime_frame:StartTime_frame+rec_dur]
                #SleepScoredTS_upscaled_ministart[SleepScoredTS_upscaled_ministart == '2'] = '1'

                # Determine each substate identity and duration
                SleepScoredTS_upscaled_ministart=SleepScoredTS_upscaled_ministart.astype(int)
                substates_duration = [len(list(group)) for key, group in groupby(SleepScoredTS_upscaled_ministart)]
                substates_identity = [key for key, _ in groupby(SleepScoredTS_upscaled_ministart)]
                substates_end = np.array(substates_duration).cumsum()        
                substates_start =np.append([0],substates_end[:-1]) #substates_start =np.append([1],substates_end+1) create 1 second gap
                substates_identity = [mapp[num] for num in substates_identity]
                substates = pd.DataFrame(list(zip(substates_identity, substates_duration, substates_start, substates_end)), columns=['Identity', 'Duration', 'Start','End'])
    
            else:
                if session_type == 'Cheeseboard':
                    print(f"No scoring file found, assuming it's all Wake cause it's a Cheeseboard session")
                    substates = pd.DataFrame(data={'Identity': ['AW'], 'Duration': [Carray.shape[0]], 'Start': [StartTime], 'End': [Carray.shape[0]]})
                    SleepScoredTS_upscaled_ministart = np.ones(Carray.shape[0]).astype(int) # all awake
                    SleepScoredTS = np.ones(int(Carray.shape[0]*1.5)).astype(int) # all awake
                else: 
                    print(f'!!!! No scoring file found in a {session_type} session !!!!')
                    continue

            print(session, ': starts at', round(StartTime,1), 's & ends at', round(EndTime,1), 's (', round(rec_dur_sec,1), 's duration, ', numbdropfr, 'dropped frames, minian frequency =', minian_freq, 'Hz, experiment type = ', session_type, ')...') 
            sentence1= f"... kept values = {kept_uniq_unit_List}"
            print(sentence1)

            # Define cell assemblies
            target_rate = 20 #Hz == 50ms bins
            #Array_bin = resample_matrix(Carray, orig_rate=minian_freq, target_rate=target_rate)
            Array_bin = bin_sum_fractional(Sarray, minian_freq, target_rate)            
            
            # Add missing cells from indexMapp
            dN_=  zscore(Array_bin)
            df_array=pd.DataFrame(dN_.T, index=kept_uniq_unit_List)
            dN = df_array.reindex(CellAssembly_dict[mice].index, fill_value=0).to_numpy().T

            patterns = CellAssembly_dict[mice].to_numpy().T
            all_cells = CellAssembly_dict[mice].index.tolist()
            assemblies_ID = CellAssembly_dict[mice].columns.tolist()

            if len(patterns)>0:
                patterns_th = patterns.copy()
                for ass in np.arange(np.shape(patterns)[0]):
                    
                    assembly_ID = assemblies_ID[ass]
                    cells_in_assembly = CellAssembly_members[mice][assembly_ID]

                    template = np.add.outer(abs(patterns[ass]),abs(patterns[ass]))      
                    template = np.nan_to_num(template, nan=0)
                    template = template - np.diag(np.diag(template))
                    tmp = dN @ template 
                    Reactivation_Strength= np.nansum(tmp * dN, axis=1) 
                    Reactivation_Strength = zscore(Reactivation_Strength)

                    scale_factor=target_rate/0.2  #cause scoring was done in 5 seconds bin, ie 0.2 Hz   
                    SleepScoredTS_binned = np.repeat(SleepScoredTS, scale_factor, axis=0)
                    StartTime_binned=round(StartTime*target_rate)
                    SleepScoredTS_binned=SleepScoredTS_binned[StartTime_binned:StartTime_binned+(rec_dur_sec*target_rate).astype(int)] 
                    SleepScoredTS_binned=SleepScoredTS_binned.astype(int)           

                    for m in mapp:  
                        Bool = (SleepScoredTS_binned == m)
                        Reactivation_Strength_VigSpe = Reactivation_Strength.copy()
                        Reactivation_Strength_VigSpe = Reactivation_Strength_VigSpe[0:np.shape(SleepScoredTS_binned)[0]] # if Calcium imaging longer than LFP rec
                        Reactivation_Strength_VigSpe[~Bool] = np.nan
                        sizeVigSt=len(Reactivation_Strength_VigSpe[Bool])
                        mean_act_ass = np.nanmean(Reactivation_Strength_VigSpe)
                        
                        peaks, properties = find_peaks(Reactivation_Strength_VigSpe, prominence=2) # 2 standard deviations away from the mean

                        CellAssembly_GlobalResults.loc[counter, 'Mice'] = mice
                        CellAssembly_GlobalResults.loc[counter, 'NeuronType'] = NeuronType
                        
                        CellAssembly_GlobalResults.loc[counter, 'Session'] = session
                        CellAssembly_GlobalResults.loc[counter, 'Session_Date'] = session_date 
                        CellAssembly_GlobalResults.loc[counter, 'Session_Time'] = session_time                    

                        CellAssembly_GlobalResults.loc[counter, 'Assembly_ID'] = assembly_ID
                        CellAssembly_GlobalResults.loc[counter, 'Assembly_size'] = len(cells_in_assembly)
                        CellAssembly_GlobalResults.loc[counter, 'Cells_in_Assembly'] = cells_in_assembly
                        
                        CellAssembly_GlobalResults.loc[counter, 'ExpeType'] =  session_type
                        CellAssembly_GlobalResults.loc[counter, 'Drug'] =  drug

                        CellAssembly_GlobalResults.loc[counter, 'Substate'] = mapp[m]
                        CellAssembly_GlobalResults.loc[counter, 'SubstateDuration'] = sizeVigSt/minian_freq
                        CellAssembly_GlobalResults.loc[counter, 'Session_ID'] = session_date + '_' + session_time

                        CellAssembly_GlobalResults.loc[counter, 'Avg_ReactStr'] = mean_act_ass
                        CellAssembly_GlobalResults.loc[counter, 'EventFreq'] = len(peaks)/(sizeVigSt/minian_freq) if sizeVigSt>0 else 0
                        CellAssembly_GlobalResults.loc[counter, 'EventTime'] = str(np.round(peaks/minian_freq+StartTime, 2))

                        counter+=1
            
    filenameOut = folder_to_save / f'CellAssembly_Global_{mice}.pkl'
    with open(filenameOut, 'wb') as pickle_file:
        pickle.dump(CellAssembly_GlobalResults, pickle_file)

filenameOut = folder_to_save / f'CellAssembly_Global.pkl'
with open(filenameOut, 'wb') as pickle_file:
    pickle.dump(CellAssembly_GlobalResults, pickle_file)

sys.stdout = sys.__stdout__
logfile.close()