#!/usr/bin/env python
# coding: utf-8

# 
# # Associate Ca2+ signal with sleep stages for each session & subsessions using crossregistration
# 
# Load packages

# In[1]:


import os
os.chdir("C:/Users/Manip2/SCRIPTS/Code python audrey/code python aurelie/interfaceJupyter/minian")

# In[2]:


import quantities as pq
import numpy as np
import math 
import neo
import json
from pathlib import Path
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, Cursor
import pickle
get_ipython().run_line_magic('matplotlib', 'widget')

from itertools import groupby
from ephyviewer import mkQApp, MainViewer, TraceViewer
from IPython.display import display
from ipyfilechooser import FileChooser


from minian.utilities import (
    TaskAnnotation,
    get_optimal_chk,
    load_videos,
    open_minian,
    save_minian,
)


# Load sleep score and Ca2+ time series numpy arrays

# In[3]:


dpath = "//10.69.168.1/crnldata/waking/audrey_hay/L1imaging/AnalysedMarch2023/Gaelle/Baseline_recording"
try:
    get_ipython().run_line_magic('store', '-r dpath')
except:
    print("data not in strore")
    #dpath = "/Users/mb/Documents/Syntuitio/AudreyHay/PlanB/ExampleRedLines/2022_08_06/13_30_01/My_V4_Miniscope/"
    dpath ="//10.69.168.1/crnldata/waking/audrey_hay/L1imaging/AnalysedMarch2023/Gaelle/Baseline_recording"

# Set up Initial Basic Parameters#
minian_path = "."

fc1 = FileChooser(dpath,select_default=True, show_only_dirs = True, title = "<b>Folder with videos</b>")
display(fc1)

# Sample callback function
def update_my_folder(chooser):
    global dpath
    dpath = chooser.selected
    get_ipython().run_line_magic('store', 'dpath')
    return 

# Register callback function
fc1.register_callback(update_my_folder)


# In[15]:


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

for y in range(1, nb_sessions+1):
    session= 'session' + str(y)
    print(session)
    sessions.append(session)
    folder_mini = folder_base / f'session{y}/V4_Miniscope'
    nb_subsessions = sum(1 for p in folder_mini.iterdir() if p.is_dir() and p.name.startswith("session"))
    ScoringFile = folder_base / f'session{y}/OpenEphys/ScoredSleep.npy'
    StampsFile = folder_base / f'session{y}/SynchroFile.xlsx'

    if nb_subsessions!=0:
        for x in range(1, nb_subsessions+1):            
            subsession= "session"  + str(y) + str(x)
            subsessions.append(subsession)    
            minian_ds = open_minian(folder_mini / subsession / f'minian')      # OR minianAB
            dict_Calcium[subsession] = minian_ds['C'] # calcium traces 
            dict_Spike[subsession] = minian_ds['S'] # estimated spikes
            dict_Scoring[subsession]  = np.load(ScoringFile)
            dict_Stamps[subsession]  = pd.read_excel(StampsFile)
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
            print(nb_minian_total)
    else:
        minian_ds = open_minian(folder_mini / f'minian')            # OR minianAB
        dict_Calcium[session] = minian_ds['C'] # calcium traces 
        dict_Spike[session] = minian_ds['S'] # estimated spikes
        dict_Scoring[session]  = np.load(ScoringFile) 
        dict_Stamps[session]  = pd.read_excel(StampsFile)
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
        print(nb_minian_total)


# Cross registration results

# In[16]:


B = mapping['session']
for c in range(len(B)):
    print('unit n°', c)
    for sess in B.keys(): #list(dict_Stamps.keys()):
        print('= unit', int(B[sess][c]), 'in', sess) if math.isnan (float(B[sess][c])) == False else None


# Distribute Ca2+ intensity & spikes to vigilance state for each sessions/subsessions

# In[35]:


data = {}
counter=0
VigilanceState_GlobalResults= pd.DataFrame(data, columns=['Mice','Session', 'Session_Time', 'Unique_Unit','UnitNumber','UnitValue', 'Substate','SubstateNumber','DurationSubstate', 'CalciumActivity', 'AUC_calcium', 'Avg_CalciumActivity','Avg_AUC_calcium','SpikeActivity', 'AUC_spike','Avg_SpikeActivity','Avg_AUC_spike'])

for i in list(dict_Stamps.keys()):

    # Start time & freq miniscope
    StartTime = list(dict_Stamps[i][0])[0]
    minian_freq=list(dict_Stamps[i][0])[2]

    First_frame = StartTime*minian_freq
    C=dict_Calcium[i]
    Cupd = C.loc[:, :] #C.loc[:, First_frame:]
    nb_unit = Cupd.shape[0]
    rec_dur = Cupd.shape[1]
    S=dict_Spike[i] 
    Supd = S.loc[:, :] #S.loc[:, First_frame:]

    # Upscale scoring to miniscope frequency
    scale_factor=minian_freq/0.2  #cause scoring was done in 5 seconds bin, ie 0.2 Hz
    SleepScoredTS=dict_Scoring[i]
    SleepScoredTS_upscaled = np.repeat(SleepScoredTS, scale_factor, axis=0)
    StartTime_frame=int(StartTime*minian_freq)
    SleepScoredTS_upscaled_ministart=SleepScoredTS_upscaled[int(StartTime_frame):int(StartTime_frame)+rec_dur]

    # Identify start time and upscale scoring to miniscope acquisition frequency
    unit_to_drop=dict_TodropFile[i]
    D = C['unit_id']
    copyD = list(D.copy())
    for r in range(len(unit_to_drop)):
        elem = unit_to_drop[r]
        copyD.remove(elem)
    unit_to_keep = copyD

    C_upd = Cupd.loc[unit_to_keep,:]
    S_upd = Supd.loc[unit_to_keep,:]
    nb_unit = C_upd.shape[0]
    print(len(C_upd.unit_id), 'selected units in', i)

        
    # Determine each substate identity and duration
    array=SleepScoredTS_upscaled_ministart
    substates_duration = [len(list(group)) for key, group in groupby(array)]
    substates_identity = [key for key, _ in groupby(array)]
    substates_end = np.array(substates_duration).cumsum()
    substates_start =np.append([1],substates_end+1)
    mapp = {0: 'NREM', 0.5: 'N2', 1: 'REM', 1.5: 'Wake'}
    substates_identity = [mapp[num] for num in substates_identity]
    substates = pd.DataFrame(list(zip(substates_identity, substates_duration, substates_start, substates_end)), columns=['Identity', 'Duration', 'Start','End'])

    C_upd_unit_id = C_upd['unit_id'].values #added by AB
    calciumtraces_all = C_upd.to_series()
    spikestraces_all = S_upd.to_series()

    for unit in range(nb_unit): # nb_unit):
        ca_input = np.array(calciumtraces_all)[(unit)*rec_dur:(unit+1)*rec_dur]
        sp_input = np.array(spikestraces_all)[(unit)*rec_dur:(unit+1)*rec_dur]
        for index in range(len(substates)):
            ca_input_sub=ca_input[substates.Start[index]:substates.End[index]]
            sp_input_sub=sp_input[substates.Start[index]:substates.End[index]]
            VigilanceState_GlobalResults.loc[counter, 'Mice'] = os.path.basename(folder_base)
            VigilanceState_GlobalResults.loc[counter, 'Session'] = i 
            VigilanceState_GlobalResults.loc[counter, 'Session_Time'] = None 
            
            keys_list = list(dict_Stamps.keys())
            key_to_find = i
            position = keys_list.index(key_to_find)
            keys_list2 = list(B.keys())
            key_at_position = keys_list2[position]
            indexMapp = np.where(B[key_at_position] == C_upd_unit_id[unit])[0]

            VigilanceState_GlobalResults.loc[counter, 'Unique_Unit'] = indexMapp
            VigilanceState_GlobalResults.loc[counter, 'UnitNumber'] = unit
            VigilanceState_GlobalResults.loc[counter, 'UnitValue'] = C_upd_unit_id[unit]
            VigilanceState_GlobalResults.loc[counter, 'Substate'] = substates.Identity[index]
            VigilanceState_GlobalResults.loc[counter, 'SubstateNumber'] = substates.index[index]
            VigilanceState_GlobalResults.loc[counter, 'DurationSubstate'] = substates.Duration[index]
            VigilanceState_GlobalResults.loc[counter, 'CalciumActivity'] = ca_input_sub.mean()
            VigilanceState_GlobalResults.loc[counter, 'AUC_calcium'] = np.trapz(ca_input_sub,np.arange(0,len(ca_input_sub),1))
            VigilanceState_GlobalResults.loc[counter, 'Avg_CalciumActivity'] = ca_input.mean()
            VigilanceState_GlobalResults.loc[counter, 'Avg_AUC_calcium'] = np.trapz(ca_input,np.arange(0,len(ca_input),1))
            VigilanceState_GlobalResults.loc[counter, 'SpikeActivity'] = sp_input_sub.mean()
            VigilanceState_GlobalResults.loc[counter, 'AUC_spike'] = np.trapz(sp_input_sub,np.arange(0,len(sp_input_sub),1))
            VigilanceState_GlobalResults.loc[counter, 'Avg_SpikeActivity'] = sp_input.mean()
            VigilanceState_GlobalResults.loc[counter, 'Avg_AUC_spike'] = np.trapz(sp_input,np.arange(0,len(sp_input),1))
            counter+=1


# In[36]:


VigilanceState_GlobalResults


# Save to Excel

# In[37]:


mice=os.path.basename(folder_base) 
filenameOut = folder_base / f'VigilanceState_GlobalResultsAB_{mice}.xlsx'
writer = pd.ExcelWriter(filenameOut)
VigilanceState_GlobalResults.to_excel(writer)
writer.close()

