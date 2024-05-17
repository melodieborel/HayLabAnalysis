# # Associate Ca2+ signal with sleep stages for each session & subsessions using crossregistration


#######################################################################################
                                # Load packages #
#######################################################################################

import os
import numpy as np
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

from itertools import groupby
from IPython.display import display

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
#######################################################################################
                                # Define functions #
#######################################################################################

def Convert(string):
            li = list(string.split(", "))
            li2 = len(li)
            return li2

#######################################################################################
                # Load sleep score and Ca2+ time series numpy arrays #
#######################################################################################

MiceList=['BlackLinesOK', 'BlueLinesOK', 'GreenDotsOK', 'GreenLinesOK', 'Purple', 'RedLinesOK','ThreeColDotsOK', 'ThreeBlueCrossesOK']
MiceList=['GreenDotsOK']

# Get the current date and time
FolderNameSave=str(datetime.now())
FolderNameSave = FolderNameSave.replace(" ", "_").replace(".", "_").replace(":", "_")
destination_folder= f"//10.69.168.1/crnldata/waking/audrey_hay/L1imaging/AB_Analysis/Analysis_VigStates_{FolderNameSave}"
os.makedirs(destination_folder)
folder_to_save=Path(destination_folder)

# Copy the script file to the destination folder
source_script = "C:/Users/Manip2/SCRIPTS/Code python audrey/code python aurelie/interfaceJupyter/python/12_13AB_AssociateMinianSleepScore_FullAuto.py"
destination_file_path = f"{destination_folder}/12_13AB_AssociateMinianSleepScore_FullAuto.txt"
shutil.copy(source_script, destination_file_path)

for micename in MiceList:
    # Load sleep score and Ca2+ time series numpy arrays
    dpath0 = "//10.69.168.1/crnldata/waking/audrey_hay/L1imaging/AnalysedMarch2023/Gaelle/Baseline_recording_ABmodified/"
    dpath=dpath0 + micename
    print(f"####################################################################################")
    print(f"################################### {micename} ####################################")
    print(f"####################################################################################")
    print(f"Path to the folder : {dpath}")
    folder_base = Path(dpath)

    nb_sessions = sum(1 for p in folder_base.iterdir() if p.is_dir() and p.name.startswith("session"))
    try:
        mfile = open(folder_base / f'mappingsAB.pkl', 'rb')
        mapping = pickle.load(mfile)
    except:
        mfile = open(folder_base / f'mappings.pkl', 'rb')
        mapping = pickle.load(mfile)

    sessions = []
    subsessions = []
    nb_minian_total=0
    dict_Calcium = {}
    dict_Spike = {}
    dict_Scoring = {}
    dict_Stamps = {}
    dict_TodropFile = {}
    dict_StampsMiniscope = {}

    sessions = [folder.name for folder in folder_base.iterdir() if folder.is_dir() and "session" in folder.name]

    for session in sessions: #range(1, nb_sessions+1):
        #session= 'session' + str(y)
        print(session)
        folder_mini = folder_base / session / f'V4_Miniscope'
        nb_subsessions = sum(1 for p in folder_mini.iterdir() if p.is_dir() and p.name.startswith("session"))
        ScoringFile = folder_base / session/ f'OpenEphys/ScoredSleep.npy'
        StampsFile = folder_base / session/ f'SynchroFile.xlsx'
        StampsMiniscopeFile = folder_mini / f'timeStamps.csv'

        if nb_subsessions!=0:
            for x in range(1, nb_subsessions+1):            
                subsession= session + str(x)
                subsessions.append(subsession)    
                minian_ds = open_minian(folder_mini / subsession / f'minian')      # OR minianAB
                dict_Calcium[subsession] = minian_ds['C'] # calcium traces 
                dict_Spike[subsession] = minian_ds['S'] # estimated spikes
                dict_Scoring[subsession]  = np.load(ScoringFile)
                dict_Stamps[subsession]  = pd.read_excel(StampsFile)
                dict_StampsMiniscope[subsession]  = pd.read_csv(StampsMiniscopeFile)
                try:
                    TodropFile = folder_mini / subsession / f'minian/TodropFileAB.json'
                    with open(TodropFile, 'r') as f:
                        unit_to_drop = json.load(f)
                        dict_TodropFile[subsession]  = unit_to_drop
                except:
                    TodropFile = folder_mini / subsession / f'minian/TodropFile.json'
                    with open(TodropFile, 'r') as f:
                        unit_to_drop = json.load(f)
                        dict_TodropFile[subsession]  = unit_to_drop
                nb_minian_total+=1
        else:
            minian_ds = open_minian(folder_mini / f'minian')            # OR minianAB
            dict_Calcium[session] = minian_ds['C'] # calcium traces 
            dict_Spike[session] = minian_ds['S'] # estimated spikes
            dict_Scoring[session]  = np.load(ScoringFile) 
            dict_Stamps[session]  = pd.read_excel(StampsFile)
            dict_StampsMiniscope[session]  = pd.read_csv(StampsMiniscopeFile)
            try:
                TodropFile = folder_mini / f'minian/TodropFileAB.json'
                with open(TodropFile, 'r') as f:
                    unit_to_drop = json.load(f)
                    dict_TodropFile[session]  = unit_to_drop
            except:
                TodropFile = folder_mini / f'minian/TodropFile.json'
                with open(TodropFile, 'r') as f:
                    unit_to_drop = json.load(f)
                    dict_TodropFile[session]  = unit_to_drop
            nb_minian_total+=1  

    #######################################################################################
                            # Cross registration results #
    #######################################################################################
    
    B = mapping['session']    
    if os.path.basename(folder_base) == 'Purple':
        index = B.columns
        B.columns = index.str.replace('part', 'session2')

    #######################################################################################
      # Distribute Ca2+ intensity & spikes to vigilance states for each sessions/subsessions #
    #######################################################################################
    
    data = {}
    counter=0
    VigilanceState_GlobalResults= pd.DataFrame(data, columns=['Mice','Session', 'Session_Time', 'Unique_Unit','UnitNumber','UnitValue', 'Substate','SubstateNumber','DurationSubstate', 'CalciumActivity', 'AUC_calcium', 'Avg_CalciumActivity','Avg_AUC_calcium','SpikeActivity', 'AUC_spike','Avg_SpikeActivity','Avg_AUC_spike'])

    previousEndTime=0
    InitialStartTime=0

    for session in list(dict_Stamps.keys()):

        # Start time & freq miniscope
        StartTime = list(dict_Stamps[session][0])[0]
        minian_freq=list(dict_Stamps[session][0])[2]
       
        C=dict_Calcium[session]
        Cupd = C.loc[:, :]
        rec_dur = Cupd.shape[1]
        S=dict_Spike[session] 
        Supd = S.loc[:, :] 

        # Adjust the StartTime if subsessions

        if InitialStartTime==0:
            InitialStartTime=StartTime    
            StartTimeMiniscope=0 # start time of miniscope rec of that subsesssions relative to the start of the mniscope recording
        else:
            if StartTime == InitialStartTime:
                StartTime = previousEndTime + 1/minian_freq #  +1 frame in seconds
                StartTimeMiniscope= StartTime-InitialStartTime
            else:  
                InitialStartTime=StartTime
                StartTimeMiniscope=0   

        # Deal with dropped frames (failure to acquire miniscope images)

        list_droppedframes = literal_eval(dict_Stamps[session][0][3])    

        numbdropfr= 0   
        upd_rec_dur=rec_dur
        droppedframes_inrec=[]
        for item in list_droppedframes: 
            if item < (int(StartTimeMiniscope*minian_freq) + upd_rec_dur) and item > int(StartTimeMiniscope*minian_freq):
                droppedframes_inrec.append(item-int(StartTimeMiniscope*minian_freq))
                upd_rec_dur+=1 #add the dropped frame to the recording length
                numbdropfr+=1                        

        EndTime = StartTime + (upd_rec_dur/minian_freq) # in seconds
        previousEndTime=EndTime     

        print(session, ': starts at', round(StartTime,1), 's & ends at', round(EndTime,1), 's (', round(upd_rec_dur/minian_freq,1), 's duration, ', numbdropfr, 'dropped frames, minian frequency =', minian_freq, ')...') 

        # Upscale scoring to miniscope frequency

        scale_factor=minian_freq/0.2  #cause scoring was done in 5 seconds bin, ie 0.2 Hz
        SleepScoredTS=dict_Scoring[session]
        SleepScoredTS_upscaled = np.repeat(SleepScoredTS, scale_factor, axis=0)
        StartTime_frame=int(StartTime*minian_freq)
        SleepScoredTS_upscaled_ministart=SleepScoredTS_upscaled[int(StartTime_frame):int(StartTime_frame)+rec_dur]

        # Identify start time and upscale scoring to miniscope acquisition frequency

        unit_to_drop=dict_TodropFile[session]
        D = C['unit_id']
        copyD = list(D.copy())
        for r in range(len(unit_to_drop)):
            elem = unit_to_drop[r]
            copyD.remove(elem)
        unit_to_keep = copyD

        C_upd = Cupd.loc[unit_to_keep,:]
        S_upd = Supd.loc[unit_to_keep,:]
        nb_unit = C_upd.shape[0]
        print(len(C_upd.unit_id), 'selected units in', session)
            
        # Determine each substate identity and duration
        array=SleepScoredTS_upscaled_ministart
        substates_duration = [len(list(group)) for key, group in groupby(array)]
        substates_identity = [key for key, _ in groupby(array)]
        substates_end = np.array(substates_duration).cumsum()        
        substates_start =np.append([0],substates_end[:-1]) #substates_start =np.append([1],substates_end+1) create 1 second gap
        mapp = {0: 'NREM', 0.5: 'N2', 1: 'REM', 1.5: 'Wake'}
        substates_identity = [mapp[num] for num in substates_identity]
        substates = pd.DataFrame(list(zip(substates_identity, substates_duration, substates_start, substates_end)), columns=['Identity', 'Duration', 'Start','End'])

        C_upd_unit_id = C_upd['unit_id'].values 
        Carray = Cupd.values.T
        Sarray = Supd.values.T
                   
        C_upd_unit_id = Cupd['unit_id'].values
        kept_uniq_unit_List=[]
        for unit in range(nb_unit):
            indexMapp = np.where(B[session] == C_upd_unit_id[unit])[0]
            kept_uniq_unit_List.append(str(indexMapp))

        sentence1= f"... kept values = {kept_uniq_unit_List}"
        print(sentence1) 

         # Replace dropped frame in calcium and spike traces with the previous value

        TimeStamps_miniscope=list(dict_StampsMiniscope[session]["Time Stamp (ms)"]) # + (StartTime*1000))
        diff_TimeStamps = np.concatenate(([0], np.diff(TimeStamps_miniscope)))
        rounded_diff_TimeStamps = np.round(diff_TimeStamps / minian_freq)

        print(diff_TimeStamps)
        print(rounded_diff_TimeStamps)

        for droppedframe in droppedframes_inrec: 
            repeat_count = rounded_diff_TimeStamps[droppedframe]
            print('frame=', droppedframe, 'repeated',  repeat_count, 'times')
            row_to_repeat = np.tile(Carray[droppedframe], (repeat_count, 1))
            Carray = np.vstack((Carray[:droppedframe], row_to_repeat, Carray[droppedframe:]))
            row_to_repeat = np.tile(Sarray[droppedframe], (repeat_count, 1))
            Sarray = np.vstack((Sarray[:droppedframe], row_to_repeat, Sarray[droppedframe:]))

        for unit in range(nb_unit): 

            Carray_unit =Carray[:,unit]
            Sarray_unit =Sarray[:,unit]

            for index in range(len(substates)):
                ca_input_sub=Carray_unit[substates.Start[index]:substates.End[index]]
                sp_input_sub=Sarray_unit[substates.Start[index]:substates.End[index]]

                peaks, _ = find_peaks(sp_input_sub, height=np.std(sp_input_sub))
                SpTrace=np.zeros(len(sp_input_sub))
                SpTrace[peaks]=1

                VigilanceState_GlobalResults.loc[counter, 'Mice'] = os.path.basename(folder_base)
                VigilanceState_GlobalResults.loc[counter, 'Session'] = session
                VigilanceState_GlobalResults.loc[counter, 'Session_Time'] = None 
                
                indexMapp = np.where(B[session] == C_upd_unit_id[unit])[0]
                VigilanceState_GlobalResults.loc[counter, 'Unique_Unit'] = indexMapp if len(indexMapp)>0 else None
                VigilanceState_GlobalResults.loc[counter, 'UnitNumber'] = unit
                VigilanceState_GlobalResults.loc[counter, 'UnitValue'] = C_upd_unit_id[unit]
                VigilanceState_GlobalResults.loc[counter, 'Substate'] = substates.Identity[index]
                VigilanceState_GlobalResults.loc[counter, 'SubstateNumber'] = substates.index[index]
                VigilanceState_GlobalResults.loc[counter, 'DurationSubstate'] = substates.Duration[index]
                VigilanceState_GlobalResults.loc[counter, 'CalciumActivity'] = ca_input_sub.mean()
                VigilanceState_GlobalResults.loc[counter, 'AUC_calcium'] = np.trapz(ca_input_sub,np.arange(0,len(ca_input_sub),1))
                VigilanceState_GlobalResults.loc[counter, 'Avg_CalciumActivity'] = Carray_unit.mean()
                VigilanceState_GlobalResults.loc[counter, 'Avg_AUC_calcium'] = np.trapz(Carray_unit,np.arange(0,len(Carray_unit),1))
                VigilanceState_GlobalResults.loc[counter, 'SpikeActivityHz'] = SpTrace.sum()/substates.Duration[index]
                VigilanceState_GlobalResults.loc[counter, 'SpikeActivity'] = SpTrace.sum()
                VigilanceState_GlobalResults.loc[counter, 'AUC_spike'] = np.trapz(sp_input_sub,np.arange(0,len(sp_input_sub),1))
                VigilanceState_GlobalResults.loc[counter, 'Avg_SpikeActivity'] = Sarray_unit.mean()
                VigilanceState_GlobalResults.loc[counter, 'Avg_AUC_spike'] = np.trapz(Sarray_unit,np.arange(0,len(Sarray_unit),1))
                counter+=1

    mice=os.path.basename(folder_base) 
    filenameOut = folder_to_save / f'VigilanceState_GlobalResultsAB_{mice}.xlsx'
    writer = pd.ExcelWriter(filenameOut)
    VigilanceState_GlobalResults.to_excel(writer)
    writer.close()

