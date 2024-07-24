# # Associate Ca2+ signal with spindles for each session & subsessions using crossregistration

#######################################################################################
                            # Define Experiment type #
#######################################################################################

#DrugExperiment=0 #if Baseline Experiment =1#if CGP Experiment
DrugExperiment=1

Method=0 # 1=AB 0=AH
AnalysisID='_FINAL' 

suffix='_AB' if Method else '_AH'

#######################################################################################
                                # Load packages #
#######################################################################################

import os
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
from scipy.interpolate import interp2d
from scipy.signal import find_peaks
from scipy.stats import zscore
import pickle
import os
from scipy.interpolate import griddata
import logging
import sys 
import shutil
from bisect import bisect_left
from ast import literal_eval

from itertools import groupby
from ephyviewer import mkQApp, MainViewer, TraceViewer
from IPython.display import display
from ipyfilechooser import FileChooser
from datetime import datetime

import warnings
warnings.filterwarnings("ignore")

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

def take_closest(myList, myNumber):        
    #Assumes myList is sorted. Returns closest value to myNumber.
    #If two numbers are equally close, return the smallest number.        
    pos = bisect_left(myList, myNumber)
    if pos == 0:
        return 0
    if pos == len(myList):
        return len(myList)
    before = myList[pos - 1]
    after = myList[pos]
    if after - myNumber < myNumber - before:
        return after
    else:
        return before

def take_closest2(myList, myNumber):
    value2 = 10000000
    for ind in range(len(myList)):
        value = abs(myList[ind]-myNumber)
        if value < value2:
            value2 = value
            index = myList[ind]
    return index

def take_closest3(myList, myNumber):
    pos = bisect_left(myList, myNumber)
    if pos == 0:
        return myList[0]
    if pos == len(myList):
        return myList[-1]
    before = myList[pos - 1]
    after = myList[pos]
    if after - myNumber < myNumber - before:
        dummy = myList.index(after)
        return dummy
    else:
        dummy = myList.index(before)
        return dummy

def is_between(myList, starttime, endtime):
    IsTrue='False'
    for ind in range(len(myList)):
        if starttime <= myList[ind] <= endtime:
            IsTrue='True'
    return IsTrue

def is_overlapping(starttime, endtime, starttimeList, endtimeList):
    IsTrue='False'
    for ind in range(len(starttimeList)):
        if starttime<=starttimeList[ind] and starttimeList[ind]<=endtime: # event n°2 begins after the start n°1               
            if (endtime-starttimeList[ind])>=int(0.5*(endtime-starttime)): # overlapp > to 50% of the duration of the event n°1
                IsTrue='True'
        elif starttime<=endtimeList[ind] and endtimeList[ind]<=endtime: # event n°2 ends before the end n°1 
            if (endtimeList[ind]-starttime)>=int(0.5*(endtime-starttime)): # overlapp > to 50% of the duration of the event n°1
                IsTrue='True'
    return IsTrue

def Convert(string):
    li = list(string.split(", "))
    li2 = len(li)
    return li2

def find_session_folders(root_path):
    sessions = []
    sessions_path=[]
    # Iterate through items in the root_path
    for item in os.listdir(root_path):
        item_path = os.path.join(root_path, item)
        if os.path.isdir(item_path):
            # Check if the directory name contains "session"
            if "session" in item:
                sessions.append(item)
                sessions_path.append(item_path)
            else:
                # Check the subdirectories of the current directory
                for sub_item in os.listdir(item_path):
                    sub_item_path = os.path.join(item_path, sub_item)
                    if os.path.isdir(sub_item_path) and "session" in sub_item:
                        sessions.append(sub_item)
                        sessions_path.append(sub_item_path)
                        
    return sessions, sessions_path

    
#######################################################################################
                # Load sleep score and Ca2+ time series numpy arrays #
#######################################################################################

MiceList=['BlackLinesOK', 'BlueLinesOK', 'GreenDotsOK','Purple' ,'ThreeColDotsOK'] if DrugExperiment else ['BlackLinesOK', 'BlueLinesOK', 'GreenDotsOK', 'GreenLinesOK', 'Purple', 'RedLinesOK','ThreeColDotsOK', 'ThreeBlueCrossesOK']

# Get the current date and time
FolderNameSave=str(datetime.now())[:19]
FolderNameSave = FolderNameSave.replace(" ", "_").replace(".", "_").replace(":", "_").replace("-", "_")
destination_folder= f"//10.69.168.1/crnldata/waking/audrey_hay/L1imaging/AnalysedMarch2023/Gaelle/CGP/AB_Analysis/Osc_{FolderNameSave}{suffix}{AnalysisID}" if DrugExperiment else f"//10.69.168.1/crnldata/waking/audrey_hay/L1imaging/AnalysedMarch2023/Gaelle/Baseline_recording_ABmodified/AB_Analysis/Osc_{FolderNameSave}{suffix}{AnalysisID}"

os.makedirs(destination_folder)
folder_to_save=Path(destination_folder)

# Copy the script file to the destination folder
source_script = "C:/Users/Manip2/SCRIPTS/CodePythonAudrey/CodePythonAurelie/HayLabAnalysis/python/14_16_AssociationCa2DownStatesSpindSWR_FullAuto.py"
destination_file_path = f"{destination_folder}/14_16_AssociationCa2DownStatesSpindSWR_FullAuto.txt"
shutil.copy(source_script, destination_file_path)

for mice in MiceList:

    dpath0 = "//10.69.168.1/crnldata/waking/audrey_hay/L1imaging/AnalysedMarch2023/Gaelle/CGP/" if DrugExperiment else "//10.69.168.1/crnldata/waking/audrey_hay/L1imaging/AnalysedMarch2023/Gaelle/Baseline_recording_ABmodified/"
    dpath=dpath0 + mice
    print(f"####################################################################################")
    print(f"################################### {mice} ####################################")
    print(f"####################################################################################")
    print(f"Path to the folder : {dpath}")
    folder_base = Path(dpath)

    try:
        mfile = open(folder_base / f'mappingsAB.pkl', 'rb')
        mapping = pickle.load(mfile)
        print('mappingsAB.pkl opened')
    except:
        mfile = open(folder_base / f'mappings.pkl', 'rb')
        mapping = pickle.load(mfile)
        print('mappings.pkl opened')

    subsessions = []
    dict_Calcium = {}
    dict_Spike = {}
    dict_SWRprop = {}
    dict_Spindleprop_PFC = {}
    dict_Spindleprop_S1 = {}
    dict_Stamps = {}
    dict_StampsMiniscope = {}
    dict_TodropFile = {}
    dict_DSprop_S1={}
    dict_DSprop_PFC={}
    dict_Path={}

    sessions, sessions_path = find_session_folders(folder_base)
    nb_sessions=len(sessions)

    for sess,session in enumerate(sessions):  
        session_path=Path(sessions_path[sess])
        folder_mini = session_path / f'V4_Miniscope'
        nb_subsessions = sum(1 for p in folder_mini.iterdir() if p.is_dir() and p.name.startswith("session"))
        SWRproperties = session_path / f'OpenEphys/SWRproperties_8sd_AB.xlsx' if Method else session_path / f'OpenEphys/SWRproperties.csv'
        Spindleproperties_PFC = session_path / f'OpenEphys/Spindlesproperties_PFC_7sd_AB.xlsx' if Method else session_path / f'OpenEphys/Spindleproperties_PFC.csv'
        Spindleproperties_S1 = session_path / f'OpenEphys/Spindlesproperties_S1_7sd_AB.xlsx' if Method else session_path / f'OpenEphys/Spindleproperties_S1.csv'
        DownStatesproperties_S1 = session_path / f'OpenEphys/DownStatesproperties_S1_2Pro3Height_AB.xlsx' 
        DownStatesproperties_PFC = session_path / f'OpenEphys/DownStatesproperties_PFC_2Pro3Height_AB.xlsx' 
        StampsFile = session_path / f'SynchroFile.xlsx'
        StampsMiniscopeFile = folder_mini / f'timeStamps.csv'
        if nb_subsessions!=0:
            for x in range(1, nb_subsessions+1):            
                subsession= session + str(x)
                subsessions.append(subsession)    
                minian_ds = open_minian(folder_mini / subsession / f'minian')      # OR minianAB
                dict_SWRprop[subsession]  = pd.read_excel(SWRproperties) if Method else pd.read_csv(SWRproperties)
                dict_Spindleprop_PFC[subsession]  = pd.read_excel(Spindleproperties_PFC) if Method else pd.read_csv(Spindleproperties_PFC)
                dict_Spindleprop_S1[subsession]  = pd.read_excel(Spindleproperties_S1) if Method else pd.read_csv(Spindleproperties_S1)           
                dict_DSprop_S1[subsession]  = pd.read_excel(DownStatesproperties_S1)           
                dict_DSprop_PFC[subsession]  = pd.read_excel(DownStatesproperties_PFC)       
                dict_Path[subsession] = session_path
                dict_Calcium[subsession] = minian_ds['C'] # calcium traces 
                dict_Spike[subsession] = minian_ds['S'] # estimated spikes
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
        else:
            minian_ds = open_minian(folder_mini / f'minian')            # OR minianAB
            dict_Path[session] = session_path
            dict_Calcium[session] = minian_ds['C'] # calcium traces 
            dict_Spike[session] = minian_ds['S'] # estimated spikes
            dict_SWRprop[session]  = pd.read_excel(SWRproperties) if Method else pd.read_csv(SWRproperties)
            dict_Spindleprop_PFC[session]  = pd.read_excel(Spindleproperties_PFC) if Method else pd.read_csv(Spindleproperties_PFC)
            dict_Spindleprop_S1[session]  = pd.read_excel(Spindleproperties_S1) if Method else pd.read_csv(Spindleproperties_S1)
            dict_DSprop_S1[session]  = pd.read_excel(DownStatesproperties_S1)     
            dict_DSprop_PFC[session]  = pd.read_excel(DownStatesproperties_PFC)        
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
    
    #######################################################################################
                             # Detect Global Spindles #
    #######################################################################################

    for sess,session in enumerate(sessions):  
        session_path=Path(sessions_path[sess])
        folder_mini = session_path / f'V4_Miniscope'
        nb_subsessions = sum(1 for p in folder_mini.iterdir() if p.is_dir() and p.name.startswith("session"))
        if nb_subsessions!=0:
            for x in range(1, nb_subsessions+1):            
                subsession= session + str(x)
                listSpdlPFC= dict_Spindleprop_PFC[subsession]
                listSpdlS1= dict_Spindleprop_S1[subsession]
                listSpdlPFCstarts=listSpdlPFC["start time"]
                listSpdlPFCends=listSpdlPFC["end time"]
                listSpdlS1starts=listSpdlS1["start time"]
                listSpdlS1ends=listSpdlS1["end time"]
                listSpdlS1['GlobalSpindle']=None
                listSpdlPFC['GlobalSpindle']=None
                for ss in range(len(listSpdlPFCstarts)): # for PFC 
                    startPFC=listSpdlPFCstarts[ss]
                    endPFC=listSpdlPFCends[ss]
                    Istrue=is_overlapping(startPFC, endPFC, listSpdlS1starts, listSpdlS1ends)
                    listSpdlPFC.loc[ss, 'GlobalSpindle'] =Istrue
                dict_Spindleprop_PFC[subsession]=listSpdlPFC
                if x==1 : #no need to overwritte for each subsession
                    filenameOut = session_path / f'OpenEphys/Spindlesproperties_PFC_7sd_AB.xlsx'  if Method else session_path / f'OpenEphys/Spindleproperties_PFC.csv'
                    print(filenameOut)
                    #dict_Spindleprop_PFC[subsession].to_excel(filenameOut) if Method else dict_Spindleprop_PFC[subsession].to_csv(filenameOut)
                for ss in range(len(listSpdlS1starts)): # for S1 
                    startS1=listSpdlS1starts[ss]
                    endS1=listSpdlS1ends[ss]
                    Istrue=is_overlapping(startS1, endS1, listSpdlPFCstarts, listSpdlPFCends)
                    listSpdlS1.loc[ss, 'GlobalSpindle'] =Istrue       
                dict_Spindleprop_S1[subsession]=listSpdlS1
                if x==1 : #no need to overwritte for each subsession
                    filenameOut = session_path / f'OpenEphys/Spindlesproperties_S1_7sd_AB.xlsx'  if Method else session_path / f'OpenEphys/Spindleproperties_S1.csv'
                    print(filenameOut)
                    #dict_Spindleprop_S1[subsession].to_excel(filenameOut) if Method else dict_Spindleprop_S1[subsession].to_csv(filenameOut) 
        else: 
            listSpdlPFC= dict_Spindleprop_PFC[session]
            listSpdlS1= dict_Spindleprop_S1[session]
            listSpdlPFCstarts=listSpdlPFC["start time"]
            listSpdlPFCends=listSpdlPFC["end time"]
            listSpdlS1starts=listSpdlS1["start time"]
            listSpdlS1ends=listSpdlS1["end time"]
            listSpdlS1['GlobalSpindle']=None
            listSpdlPFC['GlobalSpindle']=None
            for ss in range(len(listSpdlPFCstarts)): # for PFC 
                startPFC=listSpdlPFCstarts[ss]
                endPFC=listSpdlPFCends[ss]
                Istrue=is_overlapping(startPFC, endPFC, listSpdlS1starts, listSpdlS1ends)
                listSpdlPFC.loc[ss, 'GlobalSpindle'] =Istrue
            dict_Spindleprop_PFC[session]=listSpdlPFC
            filenameOut = session_path / f'OpenEphys/Spindlesproperties_PFC_7sd_AB.xlsx' if Method else session_path / f'OpenEphys/Spindleproperties_PFC.csv'
            print(filenameOut)
            #dict_Spindleprop_PFC[session].to_excel(filenameOut) if Method else dict_Spindleprop_PFC[session].to_csv(filenameOut)
            for ss in range(len(listSpdlS1starts)): # for S1 
                startS1=listSpdlS1starts[ss]
                endS1=listSpdlS1ends[ss]
                Istrue=is_overlapping(startS1, endS1, listSpdlPFCstarts, listSpdlPFCends)
                listSpdlS1.loc[ss, 'GlobalSpindle'] =Istrue       
            dict_Spindleprop_S1[session]=listSpdlS1
            filenameOut = session_path / f'OpenEphys/Spindlesproperties_S1_7sd_AB.xlsx' if Method else session_path / f'OpenEphys/Spindleproperties_S1.csv'
            print(filenameOut)
            #dict_Spindleprop_S1[session].to_excel(filenameOut) if Method else dict_Spindleprop_S1[session].to_csv(filenameOut)
    
    #######################################################################################
                            # Cross registration results #
    #######################################################################################

    B = mapping['session']
    if mice == 'Purple' and DrugExperiment==0:
        index = B.columns
        B.columns = index.str.replace('part', 'session2')

    #######################################################################################
      # Distribute Ca2+ intensity & spikes to oscillations for each sessions/subsessions #
    #######################################################################################

    CortexList= ['PFC', 'S1']

    for Cortex in CortexList:

        data = {}        
        before = 1000 # Max distance in ms between a SWR and a spindle to be considered as Precoupled
        after = 1000 # Max distance in ms between a spindle and a SWR to be considered as Postcoupled
        durationSpdl = 2 # number of sec before and after the Spdl onset taken into acount
        durationSWR = 0.5 # number of sec before and after the SWR onset taken into acount
        durationDS = 1 # number of sec before and after the DownStates onset taken into acount
        counter=0
        counter2=0
        counter3=0

        Drugs=['Baseline', 'CGP'] if DrugExperiment else ['Baseline']

        norm_freq=20 # final miniscope frequency used for all recordings

        Spindles_GlobalResults= pd.DataFrame(data, columns=['Mice', 'Session','Session_Time','Unique_Unit','UnitNumber','UnitValue','Drug', 'SpdlStatut','SpdlNumber','SpdlDuration','SWR_inside_Spdl','GlobalSpindle','CalciumActivityPreference', 'CalciumActivityBefore','CalciumActivityDuring','CalciumActivityAfter','AUC_calciumBefore','AUC_calciumDuring','AUC_calciumAfter','SpikeActivityPreference','SpikeActivityBefore','SpikeActivityDuring','SpikeActivityAfter'])
        SWR_GlobalResults= pd.DataFrame(data, columns=['Mice', 'Session','Session_Time','Unique_Unit','UnitNumber','UnitValue','Drug','SWRStatut','SWRNumber','SWRDuration','SWR_inside_Spdl','CalciumActivityPreference', 'CalciumActivityBefore','CalciumActivityDuring','CalciumActivityAfter','AUC_calciumBefore','AUC_calciumDuring','AUC_calciumAfter','SpikeActivityPreference','SpikeActivityBefore','SpikeActivityDuring','SpikeActivityAfter'])
        DS_GlobalResults= pd.DataFrame(data, columns=['Mice', 'Session','Session_Time','Unique_Unit','UnitNumber','UnitValue','Drug','DSStatut','DSNumber','DSDuration','DS_inside_Spdl','CalciumActivityPreference', 'CalciumActivityBefore','CalciumActivityDuring','CalciumActivityAfter','AUC_calciumBefore','AUC_calciumDuring','AUC_calciumAfter','SpikeActivityPreference','SpikeActivityBefore','SpikeActivityDuring','SpikeActivityAfter'])

        dict_All_ActivityCa_ds_Baseline={}
        dict_All_ActivityCa_ds_Precoupled_Baseline={}
        dict_All_ActivityCa_ds_Postcoupled_Baseline={}
        dict_All_ActivityCa_ds_Uncoupled_Baseline={}
                
        dict_All_ActivitySp_ds_Baseline={}
        dict_All_ActivitySp_ds_Precoupled_Baseline={}
        dict_All_ActivitySp_ds_Postcoupled_Baseline={}
        dict_All_ActivitySp_ds_Uncoupled_Baseline={}
        
        dict_All_ActivityCa_ds_CGP={}
        dict_All_ActivityCa_ds_Precoupled_CGP={}
        dict_All_ActivityCa_ds_Postcoupled_CGP={}
        dict_All_ActivityCa_ds_Uncoupled_CGP={}
                
        dict_All_ActivitySp_ds_CGP={}
        dict_All_ActivitySp_ds_Precoupled_CGP={}
        dict_All_ActivitySp_ds_Postcoupled_CGP={}
        dict_All_ActivitySp_ds_Uncoupled_CGP={}
        
        dict_All_ActivityCa_spin_Baseline={}
        dict_All_ActivityCa_spin_Precoupled_Baseline={}
        dict_All_ActivityCa_spin_Postcoupled_Baseline={}
        dict_All_ActivityCa_spin_Uncoupled_Baseline={}
        dict_All_ActivityCa_GlobalSpdl_Baseline={}
        dict_All_ActivityCa_LocalSpdl_Baseline={}        
        
        dict_All_ActivitySp_spin_Baseline={}
        dict_All_ActivitySp_spin_Precoupled_Baseline={}
        dict_All_ActivitySp_spin_Postcoupled_Baseline={}
        dict_All_ActivitySp_spin_Uncoupled_Baseline={}
        dict_All_ActivitySp_GlobalSpdl_Baseline={}
        dict_All_ActivitySp_LocalSpdl_Baseline={}

        dict_All_ActivityCa_spin_CGP={}
        dict_All_ActivityCa_spin_Precoupled_CGP={}
        dict_All_ActivityCa_spin_Postcoupled_CGP={}
        dict_All_ActivityCa_spin_Uncoupled_CGP={}
        dict_All_ActivityCa_GlobalSpdl_CGP={}
        dict_All_ActivityCa_LocalSpdl_CGP={}

        dict_All_ActivitySp_spin_CGP={}
        dict_All_ActivitySp_spin_Precoupled_CGP={}
        dict_All_ActivitySp_spin_Postcoupled_CGP={}
        dict_All_ActivitySp_spin_Uncoupled_CGP={}
        dict_All_ActivitySp_GlobalSpdl_CGP={}
        dict_All_ActivitySp_LocalSpdl_CGP={}
        
        dict_All_ActivityCa_swr_Baseline={}
        dict_All_ActivityCa_swr_Precoupled_Baseline={}
        dict_All_ActivityCa_swr_Postcoupled_Baseline={}
        dict_All_ActivityCa_swr_Uncoupled_Baseline={}
        
        dict_All_ActivitySp_swr_Baseline={}
        dict_All_ActivitySp_swr_Precoupled_Baseline={}
        dict_All_ActivitySp_swr_Postcoupled_Baseline={}
        dict_All_ActivitySp_swr_Uncoupled_Baseline={}

        dict_All_ActivityCa_swr_CGP={}
        dict_All_ActivityCa_swr_Precoupled_CGP={}
        dict_All_ActivityCa_swr_Postcoupled_CGP={}
        dict_All_ActivityCa_swr_Uncoupled_CGP={}
        
        dict_All_ActivitySp_swr_CGP={}
        dict_All_ActivitySp_swr_Precoupled_CGP={}
        dict_All_ActivitySp_swr_Postcoupled_CGP={}
        dict_All_ActivitySp_swr_Uncoupled_CGP={}

        previousEndTime=0
        InitialStartTime=0

        for session in list(dict_Stamps.keys()):    
            cPreCoupled=0
            cPostCoupled=0
            cUnCoupled=0
            cGlobal=0
            cLocal=0

            cPreCoupledSWR=0
            cPostCoupledSWR=0
            cUnCoupledSWR=0                        
                        
            cPreCoupledDS=0
            cPostCoupledDS=0
            cUnCoupledDS=0
            
            drug=os.path.basename(os.path.dirname(dict_Path[session])) if DrugExperiment else 'Baseline'

            # Start time & freq miniscope

            StartTime = list(dict_Stamps[session][0])[0] # in seconds
            minian_freq=list(dict_Stamps[session][0])[2] # in Hz

            if minian_freq>=20: # should only remove 1 session                

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

                print(session, ': starts at', round(StartTime,1), 's & ends at', round(EndTime,1), 's (', round(upd_rec_dur/minian_freq,1), 's duration, ', numbdropfr, 'dropped frames, minian frequency =', minian_freq, 'Hz, drug = ', drug, ')...') 
                
                # Remove bad units from recordings

                AA = C['unit_id']
                copyAA = list(AA.copy())
                unit_to_drop=dict_TodropFile[session]    
                for u in unit_to_drop: # ugly way to do it, need to be improved to only use unit_to_drop
                    copyAA.remove(u)
                unit_to_keep = copyAA
                Cupd = Cupd.loc[unit_to_keep,:]
                Carray = Cupd.values.T
                Supd = Supd.loc[unit_to_keep,:]
                Sarray = Supd.values.T
                nb_unit = Cupd.shape[0]
                units = range(nb_unit)
                              
                C_upd_unit_id = Cupd['unit_id'].values
                kept_uniq_unit_List=[]
                for unit in units:
                    indexMapp = np.where(B[session] == C_upd_unit_id[unit])[0]
                    kept_uniq_unit_List.append(str(indexMapp))

                sentence1= f"... kept values = {kept_uniq_unit_List}"
                print(sentence1)     
                
                if nb_unit==0:
                    continue  # next iteration
                
                # Zscore traces
                #Carray=zscore(Carray, axis=0)
                #Sarray=zscore(Sarray, axis=0)

                # Replace dropped frame in calcium and spike traces with the previous value

                for droppedframe in droppedframes_inrec: 
                    row_to_repeat = Carray[droppedframe]  
                    Carray = np.vstack((Carray[:droppedframe], row_to_repeat, Carray[droppedframe:]))
                    row_to_repeat = Sarray[droppedframe]  
                    Sarray = np.vstack((Sarray[:droppedframe], row_to_repeat, Sarray[droppedframe:]))

                # Align Oscillations to miniscope start 

                dictnameSpdl=f'dict_Spindleprop_{Cortex}' #change dict according to the cortex 
                dictnameDS=f'dict_DSprop_{Cortex}' #change dict according to the cortex 
                dictSpdl = globals()[dictnameSpdl]
                dictDS = globals()[dictnameDS]

                SpipropO=dictSpdl[session]
                SpipropM=SpipropO.copy()
                SWRpropO=dict_SWRprop[session]
                SWRpropM=SWRpropO.copy()
                DSpropO=dictDS[session]
                DSpropM=DSpropO.copy()
                SpipropM[["peak time", "start time", "end time"]] = SpipropM[["peak time", "start time", "end time"]]-(StartTime*1000)
                SWRpropM[["peak time", "start time", "end time"]] = SWRpropM[["peak time", "start time", "end time"]]-(StartTime*1000)        
                DSpropM[["peak time", "start time", "end time"]] = DSpropM[["peak time", "start time", "end time"]]-(StartTime*1000)        

                timeSpdl = range(int(durationSpdl*2*minian_freq))
                HalfSpdl = int(durationSpdl*minian_freq)
                
                timeSWR = range(int(durationSWR*2*minian_freq))
                HalfSWR = int(durationSWR*minian_freq)

                timeDS = range(int(durationDS*2*minian_freq))
                HalfDS = int(durationDS*minian_freq)

                TimeStamps_miniscope=list(dict_StampsMiniscope[session]["Time Stamp (ms)"]) # + (StartTime*1000))

                SpipropTrunc = SpipropM[SpipropM["start time"]>0]
                SpipropTrunc = SpipropTrunc[SpipropTrunc["start time"]< (EndTime-StartTime)*1000]
                SWRpropTrunc = SWRpropM[SWRpropM["start time"]>0]
                SWRpropTrunc = SWRpropTrunc[SWRpropTrunc["start time"] < (EndTime-StartTime)*1000]
                DSpropTrunc = DSpropM[DSpropM["start time"]>0]
                DSpropTrunc = DSpropTrunc[DSpropTrunc["start time"] < (EndTime-StartTime)*1000]

                nb_spindle = SpipropTrunc.shape[0]
                nb_swr = SWRpropTrunc.shape[0]
                nb_ds = DSpropTrunc.shape[0]

                for unit in units: # for each kept units (cause Cseries/Sseries only have kept units)

                    Carray_unit =Carray[:,unit]
                    Sarray_unit =Sarray[:,unit]
                    #peaks, _ = find_peaks(Sarray_unit)#, height=np.std(SpTrace))
                    #Sarray_unit=np.zeros(len(Sarray_unit))
                    #Sarray_unit[peaks]=1

                    #######################################################################################
                                                        # for SPDLs #
                    #######################################################################################

                    ActivityCa_Spin = [] #For each unit  
                    ActivityCa_Spin_Precoupled= [] #For each unit 
                    ActivityCa_Spin_Postcoupled= [] #For each unit 
                    ActivityCa_Spin_Uncoupled= [] #For each unit 
                    ActivitySp_Spin = [] #For each unit  
                    ActivitySp_Spin_Precoupled= [] #For each unit 
                    ActivitySp_Spin_Postcoupled= [] #For each unit 
                    ActivitySp_Spin_Uncoupled= [] #For each unit 
                    ActivityCa_GlobalSpdl= [] #For each unit 
                    ActivitySp_GlobalSpdl= [] #For each unit 
                    ActivityCa_LocalSpdl= [] #For each unit 
                    ActivitySp_LocalSpdl= [] #For each unit 

                    startSpiList = list(pd.Series(SpipropTrunc["start time"]))
                    endSpiList = list(pd.Series(SpipropTrunc["end time"]))
                    GlobalSpdlList = list(pd.Series(SpipropTrunc["GlobalSpindle"]))

                    for Pspin in range(nb_spindle): 
                        
                        # Get the calcium and spike trace associated with the spdl

                        startSpi=startSpiList[Pspin]
                        endSpi=endSpiList[Pspin]                        

                        TooEarlySpdl=startSpi/1000<durationSpdl # too close to the begining of the recording
                        TooLateSpdl=startSpi/1000+durationSpdl>round((upd_rec_dur)/minian_freq,1) # too close to the end of the recording
                        if TooEarlySpdl or TooLateSpdl:
                            print("/!\ Spindle too close to the begining/end of the recording,", session, ", Spdl n°", Pspin, ", Start Spdl =", round(startSpi/1000,1), "s") if unit==0 else None            
                        else:

                            Frame_Spindle_start = int(startSpi/1000*minian_freq)                            
                            CaTrace = list(Carray_unit[Frame_Spindle_start-HalfSpdl:Frame_Spindle_start+HalfSpdl])
                            SpTrace = list(Sarray_unit[Frame_Spindle_start-HalfSpdl:Frame_Spindle_start+HalfSpdl]) 
                            
                            ActivityCa_Spin.append(CaTrace)
                            ActivitySp_Spin.append(SpTrace)

                            # Define if that spindle is coupled with a SWR or not

                            Spdl_statut=[]
                            startSWRList = list(pd.Series(SWRpropTrunc["start time"]))
                            if len(startSWRList)>0:
                                startClosest_SWR = take_closest2(startSWRList, startSpi)
                                distance = startClosest_SWR - startSpi
                                if (distance > (- before)) and (distance <  0):
                                    Spdl_statut = 'PreCoupled'
                                    cPreCoupled+=1 if unit==0 else 0
                                    ActivityCa_Spin_Precoupled.append(CaTrace)
                                    ActivitySp_Spin_Precoupled.append(SpTrace)
                                elif (distance > (0)) and (distance <  after):
                                    Spdl_statut = 'PostCoupled'
                                    cPostCoupled+=1 if unit==0 else 0
                                    ActivityCa_Spin_Postcoupled.append(CaTrace)
                                    ActivitySp_Spin_Postcoupled.append(SpTrace)
                                else:
                                    Spdl_statut= 'UnCoupled'
                                    cUnCoupled+=1 if unit==0 else 0
                                    ActivityCa_Spin_Uncoupled.append(CaTrace)
                                    ActivitySp_Spin_Uncoupled.append(SpTrace)
                            else:
                                Spdl_statut= 'UnCoupled'
                                cUnCoupled+=1 if unit==0 else 0
                                ActivityCa_Spin_Uncoupled.append(CaTrace)
                                ActivitySp_Spin_Uncoupled.append(SpTrace)

                            # Define if that Spindle is local or global

                            if GlobalSpdlList[Pspin]=='True':
                                ActivityCa_GlobalSpdl.append(CaTrace)
                                ActivitySp_GlobalSpdl.append(SpTrace)
                                cGlobal+=1 if unit==0 else 0
                            else:
                                ActivityCa_LocalSpdl.append(CaTrace)
                                ActivitySp_LocalSpdl.append(SpTrace)
                                cLocal+=1 if unit==0 else 0

                            # Fill the big summary table Spindles_GlobalResults

                            Spindles_GlobalResults.loc[counter, 'Mice'] = mice
                            Spindles_GlobalResults.loc[counter, 'Session'] = session
                            Spindles_GlobalResults.loc[counter, 'Session_Time'] = None 

                            indexMapp = np.where(B[session] == C_upd_unit_id[unit])[0]
                            Spindles_GlobalResults.loc[counter, 'Unique_Unit'] = indexMapp 
                            Spindles_GlobalResults.loc[counter, 'UnitNumber'] = unit 
                            Spindles_GlobalResults.loc[counter, 'UnitValue'] = C_upd_unit_id[unit] 
                            
                            Spindles_GlobalResults.loc[counter, 'Drug'] =  os.path.basename(os.path.dirname(dict_Path[session])) if DrugExperiment else 'Baseline'

                            Spindles_GlobalResults.loc[counter, 'SpdlStatut'] = Spdl_statut
                            Spindles_GlobalResults.loc[counter, 'SpdlNumber'] = Pspin
                            Spindles_GlobalResults.loc[counter, 'SpdlDuration'] = endSpi- startSpi
                            IsTrue=is_between(startSWRList,startSpi, endSpi)
                            Spindles_GlobalResults.loc[counter, 'SWR_inside_Spdl'] = IsTrue

                            Spindles_GlobalResults.loc[counter, 'GlobalSpindle'] = GlobalSpdlList[Pspin]
                            
                            # Activity before/ during/after oscillation

                            durOsc=int((endSpi- startSpi)/1000*minian_freq)
                            TooEarlySpdl=startSpi/1000<durOsc/minian_freq # too close to the begining of the recording
                            TooLateSpdl=startSpi/1000+(durOsc/minian_freq*2)>round((upd_rec_dur)/minian_freq,1) # too close to the end of the recording
                            if TooEarlySpdl or TooLateSpdl:
                                print("/!\ Spindle too close to the begining/end of the recording,", session, ", Spdl n°", Pspin, ", Start Spdl =", round(startSpi/1000,1), "s, recdur=", round((upd_rec_dur)/minian_freq,1)) if unit==0 else None            
                            else:                                
                                CaTrace = list(Carray_unit[Frame_Spindle_start-durOsc:Frame_Spindle_start+durOsc*2])
                                SpTrace = list(Sarray_unit[Frame_Spindle_start-durOsc:Frame_Spindle_start+durOsc*2]) 
                            
                                ActBefore=np.mean(CaTrace[:durOsc],0)
                                ActDuring=np.mean(CaTrace[durOsc:durOsc*2],0)
                                ActAfter=np.mean(CaTrace[durOsc*2:durOsc*3],0)
                                        
                                if ActBefore > ActDuring and ActBefore > ActAfter:
                                    pref='Before'
                                elif ActAfter > ActDuring and ActAfter > ActBefore:
                                    pref='After' 
                                elif ActDuring > ActAfter and ActDuring > ActBefore:
                                    pref='During' 
                                else:
                                    pref='None'
                                Spindles_GlobalResults.loc[counter, 'CalciumActivityPreference'] = pref
                                Spindles_GlobalResults.loc[counter, 'CalciumActivityBefore'] = ActBefore
                                Spindles_GlobalResults.loc[counter, 'CalciumActivityDuring'] = ActDuring
                                Spindles_GlobalResults.loc[counter, 'CalciumActivityAfter'] = ActAfter
                                Spindles_GlobalResults.loc[counter, 'AUC_calciumBefore'] = np.trapz(CaTrace[:durOsc],np.arange(0,len(CaTrace[:durOsc]),1))
                                Spindles_GlobalResults.loc[counter, 'AUC_calciumDuring'] = np.trapz(CaTrace[durOsc:durOsc*2],np.arange(0,len(CaTrace[durOsc:durOsc*2]),1))          
                                Spindles_GlobalResults.loc[counter, 'AUC_calciumAfter'] = np.trapz(CaTrace[durOsc*2:durOsc*3],np.arange(0,len(CaTrace[durOsc*2:durOsc*3]),1))          

                                ActBefore=np.mean(SpTrace[:durOsc],0)
                                ActDuring=np.mean(SpTrace[durOsc:durOsc*2],0)
                                ActAfter=np.mean(SpTrace[durOsc*2:durOsc*3],0)

                                if ActBefore > ActDuring and ActBefore > ActAfter:
                                    pref='Before'
                                elif ActAfter > ActDuring and ActAfter > ActBefore:
                                    pref='After' 
                                elif ActDuring > ActAfter and ActDuring > ActBefore:
                                    pref='During' 
                                else:
                                    pref='None'
                                Spindles_GlobalResults.loc[counter, 'SpikeActivityPreference'] = pref
                                Spindles_GlobalResults.loc[counter, 'SpikeActivityBefore'] = np.mean(SpTrace[:durOsc],0)
                                Spindles_GlobalResults.loc[counter, 'SpikeActivityDuring'] = np.mean(SpTrace[durOsc:durOsc*2],0)
                                Spindles_GlobalResults.loc[counter, 'SpikeActivityAfter'] = np.mean(SpTrace[durOsc*2:durOsc*3],0)                         
                            counter+=1     

                    ## Peristimulus Time Histogram 
                    # All Ca traces for each spindles per Unique unit (according to cross-registration)
                    
                    list_ActivityCa= [f'ActivityCa_Spin', f'ActivityCa_Spin_Precoupled', f'ActivityCa_Spin_Postcoupled', f'ActivityCa_Spin_Uncoupled', f'ActivityCa_GlobalSpdl', f'ActivityCa_LocalSpdl']
                    list_dict_All_ActivityCa= [f'dict_All_ActivityCa_spin_{drug}', f'dict_All_ActivityCa_spin_Precoupled_{drug}', f'dict_All_ActivityCa_spin_Postcoupled_{drug}', f'dict_All_ActivityCa_spin_Uncoupled_{drug}', f'dict_All_ActivityCa_GlobalSpdl_{drug}', f'dict_All_ActivityCa_LocalSpdl_{drug}']
                    for it, ActivityCaNames in enumerate(list_ActivityCa): # for each Spdl types
                        if len(indexMapp) > 0: #not empty --> cause some units are not in the cross registration..! Need to know why 
                            ActivityCa = locals()[ActivityCaNames]
                            dict_All_ActivityCa = locals()[list_dict_All_ActivityCa[it]]       
                            if len(ActivityCa)>0 :                                
                                if np.shape(np.array(ActivityCa))[1] == int(norm_freq*durationSpdl*2):  #normalize traces to the same frequency rate         
                                    ActivityCa= np.reshape(np.array(ActivityCa), (-1, len(np.array(ActivityCa)))) if np.ndim(ActivityCa) == 1 else np.array(ActivityCa)    
                                    key=mice + str(indexMapp).replace('[','').replace(']','')
                                    dict_All_ActivityCa[key] = np.append(dict_All_ActivityCa[key], np.array(ActivityCa), axis=0) if key in dict_All_ActivityCa else np.array(ActivityCa)
                                else:
                                    dataO = np.array(ActivityCa)
                                    data= np.repeat(dataO, 2, axis=0) if dataO.shape[0] == 1 else dataO
                                    x_mesh, y_mesh = np.meshgrid(np.arange(data.shape[1]), np.arange(data.shape[0]))
                                    x_new_mesh, y_new_mesh = np.meshgrid(np.linspace(0, data.shape[1] - 1, int(norm_freq*durationSpdl*2)), np.linspace(0, data.shape[0] - 1, np.shape(data)[0]))
                                    resampled_dataO = griddata((x_mesh.flatten(), y_mesh.flatten()), data.flatten(), (x_new_mesh, y_new_mesh), method='linear')
                                    resampled_data= resampled_dataO[0,:] if dataO.shape[0] == 1 else resampled_dataO
                                    resampled_data= np.reshape(resampled_data, (-1, len(resampled_data))) if np.ndim(resampled_data) == 1 else resampled_data
                                    key=mice + str(indexMapp).replace('[','').replace(']','')
                                    dict_All_ActivityCa[key] = np.append(dict_All_ActivityCa[key], np.array(resampled_data), axis=0) if key in dict_All_ActivityCa else np.array(resampled_data)
                        else: 
                            print(f"/!\ Cell idx {unit} not in the cross registration") if it==1 else None
                
                    # All Sp traces for each spindles per Unique unit (according to cross-registration)

                    list_ActivitySp= [f'ActivitySp_Spin', f'ActivitySp_Spin_Precoupled', f'ActivitySp_Spin_Postcoupled', f'ActivitySp_Spin_Uncoupled', f'ActivitySp_GlobalSpdl', f'ActivitySp_LocalSpdl']
                    list_dict_All_ActivitySp= [f'dict_All_ActivitySp_spin_{drug}', f'dict_All_ActivitySp_spin_Precoupled_{drug}', f'dict_All_ActivitySp_spin_Postcoupled_{drug}', f'dict_All_ActivitySp_spin_Uncoupled_{drug}', f'dict_All_ActivitySp_GlobalSpdl_{drug}', f'dict_All_ActivitySp_LocalSpdl_{drug}']
                    for it, ActivitySpNames in enumerate(list_ActivitySp): # for each Spdl types
                        if len(indexMapp) > 0: #not empty --> cause some units are not in the cross registration..! Need to know why 
                            ActivitySp = locals()[ActivitySpNames]
                            dict_All_ActivitySp = locals()[list_dict_All_ActivitySp[it]]       
                            if len(ActivitySp)>0 :    
                                if np.shape(np.array(ActivitySp))[1] == int(norm_freq*durationSpdl*2):  #normalize traces to the same frequency rate         
                                    ActivitySp= np.reshape(np.array(ActivitySp), (-1, len(np.array(ActivitySp)))) if np.ndim(ActivitySp) == 1 else np.array(ActivitySp)    
                                    key=mice + str(indexMapp).replace('[','').replace(']','')
                                    dict_All_ActivitySp[key] = np.append(dict_All_ActivitySp[key], np.array(ActivitySp), axis=0) if key in dict_All_ActivitySp else np.array(ActivitySp)
                                else:
                                    dataO = np.array(ActivitySp)
                                    data= np.repeat(dataO, 2, axis=0) if dataO.shape[0] == 1 else dataO
                                    x_mesh, y_mesh = np.meshgrid(np.arange(data.shape[1]), np.arange(data.shape[0]))
                                    x_new_mesh, y_new_mesh = np.meshgrid(np.linspace(0, data.shape[1] - 1, int(norm_freq*durationSpdl*2)), np.linspace(0, data.shape[0] - 1, np.shape(data)[0]))
                                    resampled_dataO = griddata((x_mesh.flatten(), y_mesh.flatten()), data.flatten(), (x_new_mesh, y_new_mesh), method='nearest')
                                    resampled_data= resampled_dataO[0,:] if dataO.shape[0] == 1 else resampled_dataO
                                    resampled_data= np.reshape(resampled_data, (-1, len(resampled_data))) if np.ndim(resampled_data) == 1 else resampled_data
                                    key=mice + str(indexMapp).replace('[','').replace(']','')
                                    dict_All_ActivitySp[key] = np.append(dict_All_ActivitySp[key], np.array(resampled_data), axis=0) if key in dict_All_ActivitySp else np.array(resampled_data)
                                
                    #######################################################################################
                                                        # for SWRs #
                    #######################################################################################

                    ActivityCa_swr = [] #For each unit  
                    ActivityCa_swr_Precoupled= [] #For each unit 
                    ActivityCa_swr_Postcoupled= [] #For each unit 
                    ActivityCa_swr_Uncoupled= [] #For each unit 

                    ActivitySp_swr = [] #For each unit  
                    ActivitySp_swr_Precoupled= [] #For each unit 
                    ActivitySp_swr_Postcoupled= [] #For each unit 
                    ActivitySp_swr_Uncoupled= [] #For each unit 

                    startSwrList = list(pd.Series(SWRpropTrunc["start time"]))
                    endSwrList = list(pd.Series(SWRpropTrunc["end time"]))

                    for Pswr in range(nb_swr): 

                        # Get the calcium and spike trace associated with the SWR

                        startSwr=startSwrList[Pswr]
                        endSwr=endSwrList[Pswr]
                        
                        TooEarlySWR=startSwr/1000<durationSWR # too close to the begining of the recording
                        TooLateSWR=startSwr/1000+durationSWR>round((upd_rec_dur)/minian_freq,1) # too close to the end of the recording
                        if TooEarlySWR or TooLateSWR:
                            print("/!\ SWR too close to the begining/end of the recording,", session, ", SWR n°", Pswr, ", Start SWR =",  round(startSwr/1000,1), "s") if unit==0 else None 
                        else:

                            Frame_SWR_start = int(startSwr/1000*minian_freq)
                            CaTrace = list(Carray_unit[Frame_SWR_start-HalfSWR:Frame_SWR_start+HalfSWR])
                            SpTrace = list(Sarray_unit[Frame_SWR_start-HalfSWR:Frame_SWR_start+HalfSWR]) 

                            ActivityCa_swr.append(CaTrace) 
                            ActivitySp_swr.append(SpTrace) 

                            # Define if that SWR is coupled with a SPDL or not

                            SWR_statut=[]
                            startSpiList = list(pd.Series(SpipropTrunc["start time"]))
                            endSpiList = list(pd.Series(SpipropTrunc["end time"]))
                            if len(startSpiList)>0:
                                startClosest_Spi = take_closest2(startSpiList, startSwr)# + StartTimeIndexSpi])
                                indexSpi = startSpiList.index(startClosest_Spi)
                                endClosest_Spi=endSpiList[indexSpi]
                                distance = startClosest_Spi - startSwr #  + StartTimeIndexSpi]  
                                IsTrue = 'False'             
                                if (distance > (- before)) and (distance <  0):
                                    SWR_statut = 'Postcoupled'
                                    cPostCoupledSWR+=1 if unit==0 else 0
                                    ActivityCa_swr_Postcoupled.append(CaTrace)
                                    ActivitySp_swr_Postcoupled.append(SpTrace)
                                    if startSwr<endClosest_Spi:
                                        IsTrue = 'True' #SWR inside the Spindle
                                elif (distance > (0)) and (distance <  after):
                                    SWR_statut = 'Precoupled'
                                    cPreCoupledSWR+=1 if unit==0 else 0
                                    ActivityCa_swr_Precoupled.append(CaTrace)
                                    ActivitySp_swr_Precoupled.append(SpTrace)
                                else:
                                    SWR_statut= 'UnCoupled'
                                    cUnCoupledSWR+=1 if unit==0 else 0
                                    ActivityCa_swr_Uncoupled.append(CaTrace)
                                    ActivitySp_swr_Uncoupled.append(SpTrace)
                            else: 
                                SWR_statut= 'UnCoupled'
                                cUnCoupledSWR+=1 if unit==0 else 0
                                ActivityCa_swr_Uncoupled.append(CaTrace)
                                ActivitySp_swr_Uncoupled.append(SpTrace)

                            # Fill the big summary table SWR_GlobalResults

                            SWR_GlobalResults.loc[counter2, 'Mice'] = mice
                            SWR_GlobalResults.loc[counter2, 'Session'] = session
                            SWR_GlobalResults.loc[counter2, 'Session_Time'] = None 
                            indexMapp = np.where(B[session] == C_upd_unit_id[unit])[0]
                            SWR_GlobalResults.loc[counter2, 'Unique_Unit'] = indexMapp 
                            SWR_GlobalResults.loc[counter2, 'UnitNumber'] = unit 
                            SWR_GlobalResults.loc[counter2, 'UnitValue'] = C_upd_unit_id[unit] 
                            
                            SWR_GlobalResults.loc[counter2, 'Drug'] = os.path.basename(os.path.dirname(dict_Path[session])) if DrugExperiment else 'Baseline'

                            SWR_GlobalResults.loc[counter2, 'SWRStatut'] = SWR_statut
                            SWR_GlobalResults.loc[counter2, 'SWRNumber'] = Pswr
                            SWR_GlobalResults.loc[counter2, 'SWRDuration'] = endSwr- startSwr
                            SWR_GlobalResults.loc[counter2, 'SWR_inside_Spdl'] = IsTrue

                            # Activity before/ during/after oscillation

                            durOsc=int((endSwr- startSwr)/1000*minian_freq)
                            TooEarlySWR=startSwr/1000<durOsc/minian_freq # too close to the begining of the recording
                            TooLateSWR=startSwr/1000+(durOsc/minian_freq*2)>round((upd_rec_dur)/minian_freq,1) # too close to the end of the recording
                            if TooEarlySWR or TooLateSWR:
                                print("/!\ SWR too close to the begining/end of the recording,", session, ", SWR n°", Pswr, ", Start SWR =", round(startSwr/1000,1), "s, recdur=", round((upd_rec_dur)/minian_freq,1)) if unit==0 else None            
                            else:                                
                                CaTrace = list(Carray_unit[Frame_SWR_start-durOsc:Frame_SWR_start+durOsc*2])
                                SpTrace = list(Sarray_unit[Frame_SWR_start-durOsc:Frame_SWR_start+durOsc*2]) 

                                ActBefore=np.mean(CaTrace[:durOsc],0)
                                ActDuring=np.mean(CaTrace[durOsc:durOsc*2],0)
                                ActAfter=np.mean(CaTrace[durOsc*2:durOsc*3],0)
                                        
                                if ActBefore > ActDuring and ActBefore > ActAfter:
                                    pref='Before'
                                elif ActAfter > ActDuring and ActAfter > ActBefore:
                                    pref='After' 
                                elif ActDuring > ActAfter and ActDuring > ActBefore:
                                    pref='During' 
                                else:
                                    pref='None'
                                SWR_GlobalResults.loc[counter2, 'CalciumActivityPreference'] = pref
                                SWR_GlobalResults.loc[counter2, 'CalciumActivityBefore'] = ActBefore
                                SWR_GlobalResults.loc[counter2, 'CalciumActivityDuring'] = ActDuring
                                SWR_GlobalResults.loc[counter2, 'CalciumActivityAfter'] = ActAfter
                                SWR_GlobalResults.loc[counter2, 'AUC_calciumBefore'] = np.trapz(CaTrace[:durOsc],np.arange(0,len(CaTrace[:durOsc]),1))
                                SWR_GlobalResults.loc[counter2, 'AUC_calciumDuring'] = np.trapz(CaTrace[durOsc:durOsc*2],np.arange(0,len(CaTrace[durOsc:durOsc*2]),1))          
                                SWR_GlobalResults.loc[counter2, 'AUC_calciumAfter'] = np.trapz(CaTrace[durOsc*2:durOsc*3],np.arange(0,len(CaTrace[durOsc*2:durOsc*3]),1))          
                            
                                ActBefore=np.mean(SpTrace[:durOsc],0)
                                ActDuring=np.mean(SpTrace[durOsc:durOsc*2],0)
                                ActAfter=np.mean(SpTrace[durOsc*2:durOsc*3],0)

                                if ActBefore > ActDuring and ActBefore > ActAfter:
                                    pref='Before'
                                elif ActAfter > ActDuring and ActAfter > ActBefore:
                                    pref='After' 
                                elif ActDuring > ActAfter and ActDuring > ActBefore:
                                    pref='During' 
                                else:
                                    pref='None'
                                SWR_GlobalResults.loc[counter2, 'SpikeActivityPreference'] = pref
                                SWR_GlobalResults.loc[counter2, 'SpikeActivityBefore'] = np.mean(SpTrace[:durOsc],0)
                                SWR_GlobalResults.loc[counter2, 'SpikeActivityDuring'] = np.mean(SpTrace[durOsc:durOsc*2],0)
                                SWR_GlobalResults.loc[counter2, 'SpikeActivityAfter'] = np.mean(SpTrace[durOsc*2:durOsc*3],0)
                            counter2+=1    

                    ## Peristimulus Time Histogram 
                    # All Ca traces for each SWR per Unique unit (according to cross-registration) 
                    
                    list_ActivityCa= [f'ActivityCa_swr', f'ActivityCa_swr_Precoupled', f'ActivityCa_swr_Postcoupled', f'ActivityCa_swr_Uncoupled']
                    list_dict_All_ActivityCa= [f'dict_All_ActivityCa_swr_{drug}', f'dict_All_ActivityCa_swr_Precoupled_{drug}', f'dict_All_ActivityCa_swr_Postcoupled_{drug}', f'dict_All_ActivityCa_swr_Uncoupled_{drug}']
                    for it, ActivityCaNames in enumerate(list_ActivityCa): 
                        if len(indexMapp) > 0: #not empty --> cause some units are not in the cross registration..! Need to know why 
                            ActivityCa = locals()[ActivityCaNames]
                            dict_All_ActivityCa = locals()[list_dict_All_ActivityCa[it]]                
                            if len(ActivityCa)>0 :                                  
                                if np.shape(np.array(ActivityCa))[1] == int(norm_freq*durationSWR*2):   #normalize traces to the same frequency rate    
                                    ActivityCa= np.reshape(np.array(ActivityCa), (-1, len(np.array(ActivityCa)))) if np.ndim(ActivityCa) == 1 else np.array(ActivityCa)    
                                    key=mice + str(indexMapp).replace('[','').replace(']','')
                                    dict_All_ActivityCa[key] = np.append(dict_All_ActivityCa[key], np.array(ActivityCa), axis=0) if key in dict_All_ActivityCa else np.array(ActivityCa)
                                else:
                                    dataO = np.array(ActivityCa)
                                    data= np.repeat(dataO, 2, axis=0) if dataO.shape[0] == 1 else dataO
                                    x_mesh, y_mesh = np.meshgrid(np.arange(data.shape[1]), np.arange(data.shape[0]))
                                    x_new_mesh, y_new_mesh = np.meshgrid(np.linspace(0, data.shape[1] - 1, int(norm_freq*durationSWR*2)), np.linspace(0, data.shape[0] - 1, np.shape(data)[0]))
                                    resampled_dataO = griddata((x_mesh.flatten(), y_mesh.flatten()), data.flatten(), (x_new_mesh, y_new_mesh), method='linear')
                                    resampled_data= resampled_dataO[0,:] if dataO.shape[0] == 1 else resampled_dataO
                                    resampled_data= np.reshape(resampled_data, (-1, len(resampled_data))) if np.ndim(resampled_data) == 1 else resampled_data
                                    key=mice + str(indexMapp).replace('[','').replace(']','')
                                    dict_All_ActivityCa[key] = np.append(dict_All_ActivityCa[key], np.array(resampled_data), axis=0) if key in dict_All_ActivityCa else np.array(resampled_data)
                    
                    # All Sp traces for each SWR per Unique unit (according to cross-registration)

                    list_ActivitySp= [f'ActivitySp_swr', f'ActivitySp_swr_Precoupled', f'ActivitySp_swr_Postcoupled', f'ActivitySp_swr_Uncoupled']
                    list_dict_All_ActivitySp= [f'dict_All_ActivitySp_swr_{drug}', f'dict_All_ActivitySp_swr_Precoupled_{drug}', f'dict_All_ActivitySp_swr_Postcoupled_{drug}', f'dict_All_ActivitySp_swr_Uncoupled_{drug}']
                    for it, ActivitySpNames in enumerate(list_ActivitySp): 
                        if len(indexMapp) > 0: #not empty --> cause some units are not in the cross registration..! Need to know why 

                            ActivitySp = locals()[ActivitySpNames]
                            dict_All_ActivitySp = locals()[list_dict_All_ActivitySp[it]]                
                            if len(ActivitySp)>0 :  
                                if np.shape(np.array(ActivitySp))[1] == int(norm_freq*durationSWR*2):   
                                    ActivitySp= np.reshape(np.array(ActivitySp), (-1, len(np.array(ActivitySp)))) if np.ndim(ActivitySp) == 1 else np.array(ActivitySp)    
                                    key=mice + str(indexMapp).replace('[','').replace(']','')
                                    dict_All_ActivitySp[key] = np.append(dict_All_ActivitySp[key], np.array(ActivitySp), axis=0) if key in dict_All_ActivitySp else np.array(ActivitySp)
                                else: #normalize traces to the same frequency rate    
                                    dataO = np.array(ActivitySp)
                                    data= np.repeat(dataO, 2, axis=0) if dataO.shape[0] == 1 else dataO
                                    x_mesh, y_mesh = np.meshgrid(np.arange(data.shape[1]), np.arange(data.shape[0]))
                                    x_new_mesh, y_new_mesh = np.meshgrid(np.linspace(0, data.shape[1] - 1, int(norm_freq*durationSWR*2)), np.linspace(0, data.shape[0] - 1, np.shape(data)[0]))
                                    resampled_dataO = griddata((x_mesh.flatten(), y_mesh.flatten()), data.flatten(), (x_new_mesh, y_new_mesh), method='nearest')
                                    resampled_data= resampled_dataO[0,:] if dataO.shape[0] == 1 else resampled_dataO
                                    resampled_data= np.reshape(resampled_data, (-1, len(resampled_data))) if np.ndim(resampled_data) == 1 else resampled_data
                                    key=mice + str(indexMapp).replace('[','').replace(']','')
                                    dict_All_ActivitySp[key] = np.append(dict_All_ActivitySp[key], np.array(resampled_data), axis=0) if key in dict_All_ActivitySp else np.array(resampled_data)
                    
                    #######################################################################################
                                                    # for Down States #
                    #######################################################################################

                    ActivityCa_ds = [] #For each unit  
                    ActivityCa_ds_Precoupled= [] #For each unit 
                    ActivityCa_ds_Postcoupled= [] #For each unit 
                    ActivityCa_ds_Uncoupled= [] #For each unit 

                    ActivitySp_ds = [] #For each unit  
                    ActivitySp_ds_Precoupled= [] #For each unit 
                    ActivitySp_ds_Postcoupled= [] #For each unit 
                    ActivitySp_ds_Uncoupled= [] #For each unit 

                    startdsList = list(pd.Series(DSpropTrunc["start time"]))
                    enddsList = list(pd.Series(DSpropTrunc["end time"]))

                    for Pds in range(nb_ds): 

                        # Get the calcium and spike trace associated with the DS

                        startds=startdsList[Pds]
                        endds=enddsList[Pds]
                        
                        TooEarlyDS=startds/1000<durationDS # too close to the begining of the recording                        
                        TooLateDS=startds/1000+durationDS>round((upd_rec_dur)/minian_freq,1) # too close to the end of the recording
                        if TooEarlyDS or TooLateDS:
                            print("/!\ DS too close to the begining/end of the recording,", session, ", DS n°", Pds, ", Start DS =",  round(startds/1000,1), "s") if unit==0 else None 
                        else:

                            Frame_DS_start = int(startds/1000*minian_freq)
                            CaTrace = list(Carray_unit[Frame_DS_start-HalfDS:Frame_DS_start+HalfDS])
                            SpTrace = list(Sarray_unit[Frame_DS_start-HalfDS:Frame_DS_start+HalfDS]) 

                            ActivityCa_ds.append(CaTrace) 
                            ActivitySp_ds.append(SpTrace) 

                            # Define if that DS is coupled with a SPDL or not

                            DS_statut=[]
                            startSpiList = list(pd.Series(SpipropTrunc["start time"]))
                            endSpiList = list(pd.Series(SpipropTrunc["end time"]))
                            if len(startSpiList)>0:
                                startClosest_Spi = take_closest2(startSpiList, startds)# + StartTimeIndexSpi])
                                indexSpi = startSpiList.index(startClosest_Spi)
                                endClosest_Spi=endSpiList[indexSpi]
                                distance = startClosest_Spi - startds #  + StartTimeIndexSpi]  
                                IsTrue = 'False'             
                                if (distance > (- before)) and (distance <  0):
                                    DS_statut = 'Postcoupled'
                                    cPostCoupledDS+=1 if unit==0 else 0
                                    ActivityCa_ds_Postcoupled.append(CaTrace)
                                    ActivitySp_ds_Postcoupled.append(SpTrace)
                                    if startds<endClosest_Spi:
                                        IsTrue = 'True' #DS inside the Spindle
                                elif (distance > (0)) and (distance <  after):
                                    DS_statut = 'Precoupled'
                                    cPreCoupledDS+=1 if unit==0 else 0
                                    ActivityCa_ds_Precoupled.append(CaTrace)
                                    ActivitySp_ds_Precoupled.append(SpTrace)
                                else:
                                    DS_statut= 'UnCoupled'
                                    cUnCoupledDS+=1 if unit==0 else 0
                                    ActivityCa_ds_Uncoupled.append(CaTrace)
                                    ActivitySp_ds_Uncoupled.append(SpTrace)
                            else: 
                                DS_statut= 'UnCoupled'
                                cUnCoupledDS+=1 if unit==0 else 0
                                ActivityCa_ds_Uncoupled.append(CaTrace)
                                ActivitySp_ds_Uncoupled.append(SpTrace)

                            # Fill the big summary table DS_GlobalResults

                            DS_GlobalResults.loc[counter3, 'Mice'] = mice
                            DS_GlobalResults.loc[counter3, 'Session'] = session
                            DS_GlobalResults.loc[counter3, 'Session_Time'] = None 
                            indexMapp = np.where(B[session] == C_upd_unit_id[unit])[0]
                            DS_GlobalResults.loc[counter3, 'Unique_Unit'] = indexMapp 
                            DS_GlobalResults.loc[counter3, 'UnitNumber'] = unit 
                            DS_GlobalResults.loc[counter3, 'UnitValue'] = C_upd_unit_id[unit] 

                            DS_GlobalResults.loc[counter3, 'Drug'] =  os.path.basename(os.path.dirname(dict_Path[session])) if DrugExperiment else 'Baseline'

                            DS_GlobalResults.loc[counter3, 'DSStatut'] = DS_statut
                            DS_GlobalResults.loc[counter3, 'DSNumber'] = Pds
                            DS_GlobalResults.loc[counter3, 'DSDuration'] = endds- startds
                            DS_GlobalResults.loc[counter3, 'DS_inside_Spdl'] = IsTrue
                            
                            # Activity before/ during/after oscillation

                            durOsc=int((endds- startds)/1000*minian_freq)
                            TooEarlyDS=startds/1000<durOsc/minian_freq # too close to the begining of the recording
                            TooLateDS=startds/1000+(durOsc/minian_freq*2)>round((upd_rec_dur)/minian_freq,1) # too close to the end of the recording
                            if TooEarlyDS or TooLateDS:
                                print("/!\ DS too close to the begining/end of the recording,", session, ", DS n°", Pds, ", Start DS =", round(startds/1000,1), "s, recdur=", round((upd_rec_dur)/minian_freq,1))  if unit==0 else None            
                            else:                                
                                CaTrace = list(Carray_unit[Frame_DS_start-durOsc:Frame_DS_start+durOsc*2])
                                SpTrace = list(Sarray_unit[Frame_DS_start-durOsc:Frame_DS_start+durOsc*2]) 

                                ActBefore=np.mean(CaTrace[:durOsc],0)
                                ActDuring=np.mean(CaTrace[durOsc:durOsc*2],0)
                                ActAfter=np.mean(CaTrace[durOsc*2:durOsc*3],0)
                                        
                                if ActBefore > ActDuring and ActBefore > ActAfter:
                                    pref='Before'
                                elif ActAfter > ActDuring and ActAfter > ActBefore:
                                    pref='After' 
                                elif ActDuring > ActAfter and ActDuring > ActBefore:
                                    pref='During' 
                                else:
                                    pref='None'
                                DS_GlobalResults.loc[counter3, 'CalciumActivityPreference'] = pref
                                DS_GlobalResults.loc[counter3, 'CalciumActivityBefore'] = ActBefore
                                DS_GlobalResults.loc[counter3, 'CalciumActivityDuring'] = ActDuring
                                DS_GlobalResults.loc[counter3, 'CalciumActivityAfter'] = ActAfter
                                DS_GlobalResults.loc[counter3, 'AUC_calciumBefore'] = np.trapz(CaTrace[:durOsc],np.arange(0,len(CaTrace[:durOsc]),1))
                                DS_GlobalResults.loc[counter3, 'AUC_calciumDuring'] = np.trapz(CaTrace[durOsc:durOsc*2],np.arange(0,len(CaTrace[durOsc:durOsc*2]),1))          
                                DS_GlobalResults.loc[counter3, 'AUC_calciumAfter'] = np.trapz(CaTrace[durOsc*2:durOsc*3],np.arange(0,len(CaTrace[durOsc*2:durOsc*3]),1))          
                            
                                ActBefore=np.mean(SpTrace[:durOsc],0)
                                ActDuring=np.mean(SpTrace[durOsc:durOsc*2],0)
                                ActAfter=np.mean(SpTrace[durOsc*2:durOsc*3],0)

                                if ActBefore > ActDuring and ActBefore > ActAfter:
                                    pref='Before'
                                elif ActAfter > ActDuring and ActAfter > ActBefore:
                                    pref='After' 
                                elif ActDuring > ActAfter and ActDuring > ActBefore:
                                    pref='During' 
                                else:
                                    pref='None'
                                DS_GlobalResults.loc[counter3, 'SpikeActivityPreference'] = pref
                                DS_GlobalResults.loc[counter3, 'SpikeActivityBefore'] = np.mean(SpTrace[:durOsc],0)
                                DS_GlobalResults.loc[counter3, 'SpikeActivityDuring'] = np.mean(SpTrace[durOsc:durOsc*2],0)
                                DS_GlobalResults.loc[counter3, 'SpikeActivityAfter'] = np.mean(SpTrace[durOsc*2:durOsc*3],0)
                            counter3+=1     

                    ## Peristimulus Time Histogram 
                    # All Ca traces for each DS per Unique unit (according to cross-registration) 

                    list_ActivityCa= [f'ActivityCa_ds', f'ActivityCa_ds_Precoupled', f'ActivityCa_ds_Postcoupled', f'ActivityCa_ds_Uncoupled']
                    list_dict_All_ActivityCa= [f'dict_All_ActivityCa_ds_{drug}', f'dict_All_ActivityCa_ds_Precoupled_{drug}', f'dict_All_ActivityCa_ds_Postcoupled_{drug}', f'dict_All_ActivityCa_ds_Uncoupled_{drug}']
                    for it, ActivityCaNames in enumerate(list_ActivityCa): 
                        if len(indexMapp) > 0: #not empty --> cause some units are not in the cross registration..! Need to know why 
                            ActivityCa = locals()[ActivityCaNames]
                            dict_All_ActivityCa = locals()[list_dict_All_ActivityCa[it]]                
                            if len(ActivityCa)>0 :                                  
                                if np.shape(np.array(ActivityCa))[1] == int(norm_freq*durationDS*2):   #normalize traces to the same frequency rate    
                                    ActivityCa= np.reshape(np.array(ActivityCa), (-1, len(np.array(ActivityCa)))) if np.ndim(ActivityCa) == 1 else np.array(ActivityCa)    
                                    key=mice + str(indexMapp).replace('[','').replace(']','')
                                    dict_All_ActivityCa[key] = np.append(dict_All_ActivityCa[key], np.array(ActivityCa), axis=0) if key in dict_All_ActivityCa else np.array(ActivityCa)
                                else:
                                    dataO = np.array(ActivityCa)
                                    data= np.repeat(dataO, 2, axis=0) if dataO.shape[0] == 1 else dataO
                                    x_mesh, y_mesh = np.meshgrid(np.arange(data.shape[1]), np.arange(data.shape[0]))
                                    x_new_mesh, y_new_mesh = np.meshgrid(np.linspace(0, data.shape[1] - 1, int(norm_freq*durationDS*2)), np.linspace(0, data.shape[0] - 1, np.shape(data)[0]))
                                    resampled_dataO = griddata((x_mesh.flatten(), y_mesh.flatten()), data.flatten(), (x_new_mesh, y_new_mesh), method='linear')
                                    resampled_data= resampled_dataO[0,:] if dataO.shape[0] == 1 else resampled_dataO
                                    resampled_data= np.reshape(resampled_data, (-1, len(resampled_data))) if np.ndim(resampled_data) == 1 else resampled_data
                                    key=mice + str(indexMapp).replace('[','').replace(']','')
                                    dict_All_ActivityCa[key] = np.append(dict_All_ActivityCa[key], np.array(resampled_data), axis=0) if key in dict_All_ActivityCa else np.array(resampled_data)
                    
                    # All Sp traces for each DS per Unique unit (according to cross-registration)

                    list_ActivitySp= [f'ActivitySp_ds', f'ActivitySp_ds_Precoupled', f'ActivitySp_ds_Postcoupled', f'ActivitySp_ds_Uncoupled']
                    list_dict_All_ActivitySp= [f'dict_All_ActivitySp_ds_{drug}', f'dict_All_ActivitySp_ds_Precoupled_{drug}', f'dict_All_ActivitySp_ds_Postcoupled_{drug}', f'dict_All_ActivitySp_ds_Uncoupled_{drug}']
                    for it, ActivitySpNames in enumerate(list_ActivitySp): 
                        if len(indexMapp) > 0: #not empty --> cause some units are not in the cross registration..! Need to know why 

                            ActivitySp = locals()[ActivitySpNames]
                            dict_All_ActivitySp = locals()[list_dict_All_ActivitySp[it]]                
                            if len(ActivitySp)>0 :  
                                if np.shape(np.array(ActivitySp))[1] == int(norm_freq*durationDS*2):   
                                    ActivitySp= np.reshape(np.array(ActivitySp), (-1, len(np.array(ActivitySp)))) if np.ndim(ActivitySp) == 1 else np.array(ActivitySp)    
                                    key=mice + str(indexMapp).replace('[','').replace(']','')
                                    dict_All_ActivitySp[key] = np.append(dict_All_ActivitySp[key], np.array(ActivitySp), axis=0) if key in dict_All_ActivitySp else np.array(ActivitySp)
                                else: #normalize traces to the same frequency rate    
                                    dataO = np.array(ActivitySp)
                                    data= np.repeat(dataO, 2, axis=0) if dataO.shape[0] == 1 else dataO
                                    x_mesh, y_mesh = np.meshgrid(np.arange(data.shape[1]), np.arange(data.shape[0]))
                                    x_new_mesh, y_new_mesh = np.meshgrid(np.linspace(0, data.shape[1] - 1, int(norm_freq*durationDS*2)), np.linspace(0, data.shape[0] - 1, np.shape(data)[0]))
                                    resampled_dataO = griddata((x_mesh.flatten(), y_mesh.flatten()), data.flatten(), (x_new_mesh, y_new_mesh), method='nearest')
                                    resampled_data= resampled_dataO[0,:] if dataO.shape[0] == 1 else resampled_dataO
                                    resampled_data= np.reshape(resampled_data, (-1, len(resampled_data))) if np.ndim(resampled_data) == 1 else resampled_data
                                    key=mice + str(indexMapp).replace('[','').replace(']','')
                                    dict_All_ActivitySp[key] = np.append(dict_All_ActivitySp[key], np.array(resampled_data), axis=0) if key in dict_All_ActivitySp else np.array(resampled_data)

            else:
                print(f'/!\ {session} not taken into account cause minian frequency = {minian_freq}')
            sentence2=f"... in {Cortex}: {nb_ds} down states ({cPreCoupledDS} Pre, {cPostCoupledDS} Post & {cUnCoupledDS} Uncoupled DS), {nb_spindle} spindles ({cPreCoupled} Pre, {cPostCoupled} Post & {cUnCoupled} Uncoupled Spdl // {cGlobal} Global & {cLocal} Local) and {nb_swr} SWR detected ({cPreCoupledSWR} Pre, {cPostCoupledSWR} Post & {cUnCoupledSWR} Uncoupled SWR)"
            print(sentence2)       

        #######################################################################################
                                # Save Spindles analysis #
        #######################################################################################

        # Save the big summary table Spindles_GlobalResults

        filenameOut = folder_to_save / f'Spindles_{Cortex}_Global_{mice}.xlsx'
        writer = pd.ExcelWriter(filenameOut)
        Spindles_GlobalResults.to_excel(writer)
        writer.close()

        # Do average Calcium results for Spindles Peristimulus Time Histogram 
        
        filenameOut = folder_to_save / f'Spindles_{Cortex}_CalciumAvg_{mice}.xlsx'
        excel_writer = pd.ExcelWriter(filenameOut)

        for Drug in Drugs:

            dict_All_ActivityCa_spin=locals()[f'dict_All_ActivityCa_spin_{Drug}']
            dict_All_ActivityCa_spin_Uncoupled=locals()[f'dict_All_ActivityCa_spin_Uncoupled_{Drug}']
            dict_All_ActivityCa_spin_Precoupled=locals()[f'dict_All_ActivityCa_spin_Precoupled_{Drug}']
            dict_All_ActivityCa_spin_Postcoupled=locals()[f'dict_All_ActivityCa_spin_Postcoupled_{Drug}']
            dict_All_ActivityCa_GlobalSpdl=locals()[f'dict_All_ActivityCa_GlobalSpdl_{Drug}']
            dict_All_ActivityCa_LocalSpdl=locals()[f'dict_All_ActivityCa_LocalSpdl_{Drug}']

            AVG_dict_All_ActivityCa_spin = {key: np.mean(matrix,0) for key, matrix in dict_All_ActivityCa_spin.items()}
            AVG_dict_All_ActivityCa_spin_Uncoupled = {key: np.mean(matrix,0) for key, matrix in dict_All_ActivityCa_spin_Uncoupled.items()}
            AVG_dict_All_ActivityCa_spin_Precoupled = {key: np.mean(matrix,0) for key, matrix in dict_All_ActivityCa_spin_Precoupled.items()}
            AVG_dict_All_ActivityCa_spin_Postcoupled = {key: np.mean(matrix,0) for key, matrix in dict_All_ActivityCa_spin_Postcoupled.items()}
            AVG_dict_All_ActivityCa_spin_GlobalSpdl = {key: np.mean(matrix,0) for key, matrix in dict_All_ActivityCa_GlobalSpdl.items()}
            AVG_dict_All_ActivityCa_spin_LocalSpdl = {key: np.mean(matrix,0) for key, matrix in dict_All_ActivityCa_LocalSpdl.items()}

            Array=pd.DataFrame(AVG_dict_All_ActivityCa_spin).T
            ArrayUn=pd.DataFrame(AVG_dict_All_ActivityCa_spin_Uncoupled).T
            ArrayPre=pd.DataFrame(AVG_dict_All_ActivityCa_spin_Precoupled).T
            ArrayPost=pd.DataFrame(AVG_dict_All_ActivityCa_spin_Postcoupled).T
            ArrayGlobalSpdl=pd.DataFrame(AVG_dict_All_ActivityCa_spin_GlobalSpdl).T
            ArrayLocalSpdl=pd.DataFrame(AVG_dict_All_ActivityCa_spin_LocalSpdl).T

            Array.to_excel(excel_writer, sheet_name=f'{Drug}_All_Spindles', index=True, header=False)
            ArrayUn.to_excel(excel_writer, sheet_name=f'{Drug}_Uncoupled_Spindles', index=True, header=False)
            ArrayPre.to_excel(excel_writer, sheet_name=f'{Drug}_Precoupled_Spindles', index=True, header=False)
            ArrayPost.to_excel(excel_writer, sheet_name=f'{Drug}_Postcoupled_Spindles', index=True, header=False)
            ArrayGlobalSpdl.to_excel(excel_writer, sheet_name=f'{Drug}_Global_Spindles', index=True, header=False)
            ArrayLocalSpdl.to_excel(excel_writer, sheet_name=f'{Drug}_Local_Spindles', index=True, header=False)

        excel_writer.close()

        # Do average Spike results for Spindles Peristimulus Time Histogram 

        filenameOut = folder_to_save / f'Spindles_{Cortex}_SpikeAvg_{mice}.xlsx'
        excel_writer = pd.ExcelWriter(filenameOut)
        
        for Drug in Drugs:

            dict_All_ActivitySp_spin=locals()[f'dict_All_ActivitySp_spin_{Drug}']
            dict_All_ActivitySp_spin_Uncoupled=locals()[f'dict_All_ActivitySp_spin_Uncoupled_{Drug}']
            dict_All_ActivitySp_spin_Precoupled=locals()[f'dict_All_ActivitySp_spin_Precoupled_{Drug}']
            dict_All_ActivitySp_spin_Postcoupled=locals()[f'dict_All_ActivitySp_spin_Postcoupled_{Drug}']
            dict_All_ActivitySp_GlobalSpdl=locals()[f'dict_All_ActivitySp_GlobalSpdl_{Drug}']
            dict_All_ActivitySp_LocalSpdl=locals()[f'dict_All_ActivitySp_LocalSpdl_{Drug}']

            AVG_dict_All_ActivitySp_spin = {key: np.mean(matrix,0) for key, matrix in dict_All_ActivitySp_spin.items()}
            AVG_dict_All_ActivitySp_spin_Uncoupled = {key: np.mean(matrix,0) for key, matrix in dict_All_ActivitySp_spin_Uncoupled.items()}
            AVG_dict_All_ActivitySp_spin_Precoupled = {key: np.mean(matrix,0) for key, matrix in dict_All_ActivitySp_spin_Precoupled.items()}
            AVG_dict_All_ActivitySp_spin_Postcoupled = {key: np.mean(matrix,0) for key, matrix in dict_All_ActivitySp_spin_Postcoupled.items()}
            AVG_dict_All_ActivitySp_spin_GlobalSpdl = {key: np.mean(matrix,0) for key, matrix in dict_All_ActivitySp_GlobalSpdl.items()}
            AVG_dict_All_ActivitySp_spin_LocalSpdl = {key: np.mean(matrix,0) for key, matrix in dict_All_ActivitySp_LocalSpdl.items()}

            Array=pd.DataFrame(AVG_dict_All_ActivitySp_spin).T
            ArrayUn=pd.DataFrame(AVG_dict_All_ActivityCa_spin_Uncoupled).T
            ArrayPre=pd.DataFrame(AVG_dict_All_ActivitySp_spin_Precoupled).T
            ArrayPost=pd.DataFrame(AVG_dict_All_ActivitySp_spin_Postcoupled).T
            ArrayGlobalSpdl=pd.DataFrame(AVG_dict_All_ActivitySp_spin_GlobalSpdl).T
            ArrayLocalSpdl=pd.DataFrame(AVG_dict_All_ActivitySp_spin_LocalSpdl).T

            Array.to_excel(excel_writer, sheet_name=f'{Drug}_All_Spindles', index=True, header=False)
            ArrayUn.to_excel(excel_writer, sheet_name=f'{Drug}_Uncoupled_Spindles', index=True, header=False)
            ArrayPre.to_excel(excel_writer, sheet_name=f'{Drug}_Precoupled_Spindles', index=True, header=False)
            ArrayPost.to_excel(excel_writer, sheet_name=f'{Drug}_Postcoupled_Spindles', index=True, header=False)
            ArrayGlobalSpdl.to_excel(excel_writer, sheet_name=f'{Drug}_Global_Spindles', index=True, header=False)
            ArrayLocalSpdl.to_excel(excel_writer, sheet_name=f'{Drug}_Local_Spindles', index=True, header=False)

        excel_writer.close()

        #######################################################################################
                                        # Save SWR analysis #
        #######################################################################################

        # Save the big summary table SWR_GlobalResults

        filenameOut = folder_to_save / f'SWR_{Cortex}_Global_{mice}.xlsx'
        writer = pd.ExcelWriter(filenameOut)
        SWR_GlobalResults.to_excel(writer)
        writer.close()

        # Do average Calcium results for SWR Peristimulus Time Histogram 

        filenameOut = folder_to_save / f'SWR_{Cortex}_CalciumAvg_{mice}.xlsx'
        excel_writer = pd.ExcelWriter(filenameOut)

        for Drug in Drugs:

            dict_All_ActivityCa_swr=locals()[f'dict_All_ActivityCa_swr_{Drug}']
            dict_All_ActivityCa_swr_Uncoupled=locals()[f'dict_All_ActivityCa_swr_Uncoupled_{Drug}']
            dict_All_ActivityCa_swr_Precoupled=locals()[f'dict_All_ActivityCa_swr_Precoupled_{Drug}']
            dict_All_ActivityCa_swr_Postcoupled=locals()[f'dict_All_ActivityCa_swr_Postcoupled_{Drug}']

            AVG_dict_All_ActivityCa_swr = {key: np.mean(matrix,0) for key, matrix in dict_All_ActivityCa_swr.items()}
            AVG_dict_All_ActivityCa_swr_Uncoupled = {key: np.mean(matrix,0) for key, matrix in dict_All_ActivityCa_swr_Uncoupled.items()}
            AVG_dict_All_ActivityCa_swr_Precoupled = {key: np.mean(matrix,0) for key, matrix in dict_All_ActivityCa_swr_Precoupled.items()}
            AVG_dict_All_ActivityCa_swr_Postcoupled = {key: np.mean(matrix,0) for key, matrix in dict_All_ActivityCa_swr_Postcoupled.items()}

            Array=pd.DataFrame(AVG_dict_All_ActivityCa_swr).T
            ArrayUn=pd.DataFrame(AVG_dict_All_ActivityCa_swr_Uncoupled).T
            ArrayPre=pd.DataFrame(AVG_dict_All_ActivityCa_swr_Precoupled).T
            ArrayPost=pd.DataFrame(AVG_dict_All_ActivityCa_swr_Postcoupled).T

            Array.to_excel(excel_writer, sheet_name=f'{Drug}_All_SWR', index=True, header=False)
            ArrayUn.to_excel(excel_writer, sheet_name=f'{Drug}_Uncoupled_SWR', index=True, header=False)
            ArrayPre.to_excel(excel_writer, sheet_name=f'{Drug}_Precoupled_SWR', index=True, header=False)
            ArrayPost.to_excel(excel_writer, sheet_name=f'{Drug}_Postcoupled_SWR', index=True, header=False)

        excel_writer.close()

        # Do average Spike results for SWR Peristimulus Time Histogram 

        filenameOut = folder_to_save / f'SWR_{Cortex}_SpikeAvg_{mice}.xlsx'
        excel_writer = pd.ExcelWriter(filenameOut)

        for Drug in Drugs:

            dict_All_ActivitySp_swr=locals()[f'dict_All_ActivitySp_swr_{Drug}']
            dict_All_ActivitySp_swr_Uncoupled=locals()[f'dict_All_ActivitySp_swr_Uncoupled_{Drug}']
            dict_All_ActivitySp_swr_Precoupled=locals()[f'dict_All_ActivitySp_swr_Precoupled_{Drug}']
            dict_All_ActivitySp_swr_Postcoupled=locals()[f'dict_All_ActivitySp_swr_Postcoupled_{Drug}']

            AVG_dict_All_ActivitySp_swr = {key: np.mean(matrix,0) for key, matrix in dict_All_ActivitySp_swr.items()}
            AVG_dict_All_ActivitySp_swr_Uncoupled = {key: np.mean(matrix,0) for key, matrix in dict_All_ActivitySp_swr_Uncoupled.items()}
            AVG_dict_All_ActivitySp_swr_Precoupled = {key: np.mean(matrix,0) for key, matrix in dict_All_ActivitySp_swr_Precoupled.items()}
            AVG_dict_All_ActivitySp_swr_Postcoupled = {key: np.mean(matrix,0) for key, matrix in dict_All_ActivitySp_swr_Postcoupled.items()}

            Array=pd.DataFrame(AVG_dict_All_ActivitySp_swr).T
            ArrayUn=pd.DataFrame(AVG_dict_All_ActivitySp_swr_Uncoupled).T
            ArrayPre=pd.DataFrame(AVG_dict_All_ActivitySp_swr_Precoupled).T
            ArrayPost=pd.DataFrame(AVG_dict_All_ActivitySp_swr_Postcoupled).T

            Array.to_excel(excel_writer, sheet_name=f'{Drug}_All_SWR', index=True, header=False)
            ArrayUn.to_excel(excel_writer, sheet_name=f'{Drug}_Uncoupled_SWR', index=True, header=False)
            ArrayPre.to_excel(excel_writer, sheet_name=f'{Drug}_Precoupled_SWR', index=True, header=False)
            ArrayPost.to_excel(excel_writer, sheet_name=f'{Drug}_Postcoupled_SWR', index=True, header=False)

        excel_writer.close() 

        #######################################################################################
                                        # Save DownStates analysis #
        #######################################################################################

        # Save the big summary table DS_GlobalResults

        filenameOut = folder_to_save / f'DS_{Cortex}_Global_{mice}.xlsx'
        writer = pd.ExcelWriter(filenameOut)
        DS_GlobalResults.to_excel(writer)
        writer.close()

        # Do average Calcium results for DS Peristimulus Time Histogram 
        
        filenameOut = folder_to_save / f'DS_{Cortex}_CalciumAvg_{mice}.xlsx'
        excel_writer = pd.ExcelWriter(filenameOut)

        for Drug in Drugs:

            dict_All_ActivityCa_ds=locals()[f'dict_All_ActivityCa_ds_{Drug}']
            dict_All_ActivityCa_ds_Uncoupled=locals()[f'dict_All_ActivityCa_ds_Uncoupled_{Drug}']
            dict_All_ActivityCa_ds_Precoupled=locals()[f'dict_All_ActivityCa_ds_Precoupled_{Drug}']
            dict_All_ActivityCa_ds_Postcoupled=locals()[f'dict_All_ActivityCa_ds_Postcoupled_{Drug}']

            AVG_dict_All_ActivityCa_ds = {key: np.mean(matrix,0) for key, matrix in dict_All_ActivityCa_ds.items()}
            AVG_dict_All_ActivityCa_ds_Uncoupled = {key: np.mean(matrix,0) for key, matrix in dict_All_ActivityCa_ds_Uncoupled.items()}
            AVG_dict_All_ActivityCa_ds_Precoupled = {key: np.mean(matrix,0) for key, matrix in dict_All_ActivityCa_ds_Precoupled.items()}
            AVG_dict_All_ActivityCa_ds_Postcoupled = {key: np.mean(matrix,0) for key, matrix in dict_All_ActivityCa_ds_Postcoupled.items()}

            Array=pd.DataFrame(AVG_dict_All_ActivityCa_ds).T
            ArrayUn=pd.DataFrame(AVG_dict_All_ActivityCa_ds_Uncoupled).T
            ArrayPre=pd.DataFrame(AVG_dict_All_ActivityCa_ds_Precoupled).T
            ArrayPost=pd.DataFrame(AVG_dict_All_ActivityCa_ds_Postcoupled).T

            Array.to_excel(excel_writer, sheet_name=f'{Drug}_All_DS', index=True, header=False)
            ArrayUn.to_excel(excel_writer, sheet_name=f'{Drug}_Uncoupled_DS', index=True, header=False)
            ArrayPre.to_excel(excel_writer, sheet_name=f'{Drug}_Precoupled_DS', index=True, header=False)
            ArrayPost.to_excel(excel_writer, sheet_name=f'{Drug}_Postcoupled_DS', index=True, header=False)

        excel_writer.close()

        # Do average Spike results for DS Peristimulus Time Histogram 
        filenameOut = folder_to_save / f'DS_{Cortex}_SpikeAvg_{mice}.xlsx'
        excel_writer = pd.ExcelWriter(filenameOut)

        for Drug in Drugs:

            dict_All_ActivitySp_ds=locals()[f'dict_All_ActivitySp_ds_{Drug}']
            dict_All_ActivitySp_ds_Uncoupled=locals()[f'dict_All_ActivitySp_ds_Uncoupled_{Drug}']
            dict_All_ActivitySp_ds_Precoupled=locals()[f'dict_All_ActivitySp_ds_Precoupled_{Drug}']
            dict_All_ActivitySp_ds_Postcoupled=locals()[f'dict_All_ActivitySp_ds_Postcoupled_{Drug}']

            AVG_dict_All_ActivitySp_ds = {key: np.mean(matrix,0) for key, matrix in dict_All_ActivitySp_ds.items()}
            AVG_dict_All_ActivitySp_ds_Uncoupled = {key: np.mean(matrix,0) for key, matrix in dict_All_ActivitySp_ds_Uncoupled.items()}
            AVG_dict_All_ActivitySp_ds_Precoupled = {key: np.mean(matrix,0) for key, matrix in dict_All_ActivitySp_ds_Precoupled.items()}
            AVG_dict_All_ActivitySp_ds_Postcoupled = {key: np.mean(matrix,0) for key, matrix in dict_All_ActivitySp_ds_Postcoupled.items()}

            Array=pd.DataFrame(AVG_dict_All_ActivitySp_ds).T
            ArrayUn=pd.DataFrame(AVG_dict_All_ActivitySp_ds_Uncoupled).T
            ArrayPre=pd.DataFrame(AVG_dict_All_ActivitySp_ds_Precoupled).T
            ArrayPost=pd.DataFrame(AVG_dict_All_ActivitySp_ds_Postcoupled).T

            Array.to_excel(excel_writer, sheet_name=f'{Drug}_All_DS', index=True, header=False)
            ArrayUn.to_excel(excel_writer, sheet_name=f'{Drug}_Uncoupled_DS', index=True, header=False)
            ArrayPre.to_excel(excel_writer, sheet_name=f'{Drug}_Precoupled_DS', index=True, header=False)
            ArrayPost.to_excel(excel_writer, sheet_name=f'{Drug}_Postcoupled_DS', index=True, header=False)

        excel_writer.close()  

    sentence3=f"Nb of unique units for {mice} = {len(dict_All_ActivityCa_spin)}"
    print(sentence3)    