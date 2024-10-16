# # Associate Ca2+ signal with spindles for each session & subsessions using crossregistration

#######################################################################################
                            # Define Experiment type #
#######################################################################################

#DrugExperiment=0 #if Baseline Experiment =1#if CGP Experiment
DrugExperimentList=[0,1]

saveexcel=0

Method=0 # 1=AB 0=AH
AnalysisID='_correctScor' 

suffix='_AB' if Method else '_AH'

CTX=['S1', 'PFC', 'S1PFC']
Coupling=['', 'UnCoupled', 'Coupled']

before = 500 # Max distance in ms between a SWR and a spindle to be considered as Precoupled
after = 1000 # Max distance in ms between a spindle and a SWR to be considered as Postcoupled
durationSpdl = 5 # number of sec before and after the Spdl onset taken into acount
durationSWR = 3 # number of sec before and after the SWR onset taken into acount
durationDS = 3 # number of sec before and after the DS onset taken into acount

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
from scipy import interpolate

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

def is_between(myList, starttime, endtime):
    IsTrue=False
    for ind in range(len(myList)):
        if starttime <= myList[ind] <= endtime:
            IsTrue=True
    return IsTrue

def is_overlapping(starttime, endtime, starttimeList, endtimeList):
    IsTrue='False'
    for ind in starttimeList.index: #range(len(starttimeList)):
        if starttime<=starttimeList[ind] and starttimeList[ind]<=endtime: # event n°2 begins after the start n°1               
            if (endtime-starttimeList[ind])>=int(0.5*(endtime-starttime)): # overlapp > to 50% of the duration of the event n°1
                IsTrue='True'
                break                
        elif starttime<=endtimeList[ind] and endtimeList[ind]<=endtime: # event n°2 ends before the end n°1 
            if (endtimeList[ind]-starttime)>=int(0.5*(endtime-starttime)): # overlapp > to 50% of the duration of the event n°1
                IsTrue='True'
                break
    return IsTrue, ind

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

for DrugExperiment in DrugExperimentList: 

    Drugs=['Baseline', 'CGP'] if DrugExperiment else ['Baseline']

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

        mfile = open(folder_base / f'mappingsAB_ALL.pkl', 'rb')
        mapping = pickle.load(mfile)
        print('mappingsAB_ALL.pkl opened')

        subsessions = []
        dict_Calcium = {}
        dict_Spike = {}
        dict_SWRprop = {}
        dict_DSprop={}
        dict_Spindleprop = {}
        dict_Stamps = {}
        dict_StampsMiniscope = {}
        dict_TodropFile = {}
        dict_Path={}

        sessions, sessions_path = find_session_folders(folder_base)
        nb_sessions=len(sessions)

        for sess,session in enumerate(sessions):  
            session_path=Path(sessions_path[sess])
            folder_mini = session_path / f'V4_Miniscope'
            nb_subsessions = sum(1 for p in folder_mini.iterdir() if p.is_dir() and p.name.startswith("session"))
            SWRproperties = session_path / f'OpenEphys/SWRproperties_8sd_AB.xlsx' if Method else session_path / f'OpenEphys/SWRproperties.csv'
            Spindleproperties = session_path / f'OpenEphys/Spindlesproperties_S1&PFC_7sd_AB.xlsx' if Method else session_path / f'OpenEphys/Spindleproperties_S1&PFC.csv'
            DownStatesproperties = session_path / f'OpenEphys/DownStatesproperties_S1&PFC.csv' 
            StampsFile = session_path / f'SynchroFileCorrect.xlsx'
            StampsMiniscopeFile = folder_mini / f'timeStamps.csv'
            if nb_subsessions!=0:
                for x in range(1, nb_subsessions+1):            
                    subsession= session + str(x)
                    subsessions.append(subsession)    
                    minian_ds = open_minian(folder_mini / subsession / f'minian')      # OR minianAB
                    SWRlist= pd.read_excel(SWRproperties) if Method else pd.read_csv(SWRproperties)
                    SWRlist['toKeep'] = SWRlist['toKeep'].astype(str)  if DrugExperiment else 'True'
                    dict_SWRprop[subsession]  =SWRlist[SWRlist['toKeep'].isin(['VRAI', 'True'])]
                    Spdllist = pd.read_excel(Spindleproperties) if Method else pd.read_csv(Spindleproperties)
                    Spdllist['toKeep'] = Spdllist['toKeep'].astype(str)
                    dict_Spindleprop[subsession]  = Spdllist[Spdllist['toKeep'].isin(['VRAI', 'True'])]
                    dict_DSprop[subsession]  = pd.read_csv(DownStatesproperties)
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
                SWRlist= pd.read_excel(SWRproperties) if Method else pd.read_csv(SWRproperties)
                SWRlist['toKeep'] = SWRlist['toKeep'].astype(str) if DrugExperiment else 'True'
                dict_SWRprop[session]  =SWRlist[SWRlist['toKeep'].isin(['VRAI', 'True'])]
                Spdllist = pd.read_excel(Spindleproperties) if Method else pd.read_csv(Spindleproperties)
                Spdllist['toKeep'] = Spdllist['toKeep'].astype(str) 
                dict_Spindleprop[session]  = Spdllist[Spdllist['toKeep'].isin(['VRAI', 'True'])]
                dict_DSprop[session]  = pd.read_csv(DownStatesproperties)
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
                                # Cross registration results #
        #######################################################################################

        B = mapping['session']
        if mice == 'Purple' and DrugExperiment==0:
            index = B.columns
            B.columns = index.str.replace('part', 'session2')

        #######################################################################################
        # Distribute Ca2+ intensity & spikes to oscillations for each sessions/subsessions #
        #######################################################################################

        data = {}        
        counter=0
        counter2=0

        norm_freq=20 # final miniscope frequency used for all recordings

        Spindles_GlobalResults= pd.DataFrame(data, columns=['Mice', 'Session','Session_Time','Unique_Unit','UnitNumber','UnitValue','Drug', 'SpdlStatut','SpdlStartLocation',
                                                             'GlobalSpindle', 'SpdlNumber','SpdlDuration','SWR_inside_Spdl','CalciumActivityPreference', 'CalciumActivityBefore',
                                                             'CalciumActivityDuring','CalciumActivityAfter','AUC_calciumBefore','AUC_calciumDuring','AUC_calciumAfter',
                                                             'SpikeActivityPreference','SpikeActivityBefore','SpikeActivityDuring','SpikeActivityAfter'])
        SWR_GlobalResults= pd.DataFrame(data, columns=['Mice', 'Session','Session_Time','Unique_Unit','UnitNumber','UnitValue','Drug','SWRStatut','SpdlLoc', 'SWRNumber','SWRDuration',
                                                       'SWR_inside_Spdl','CalciumActivityPreference', 'CalciumActivityBefore','CalciumActivityDuring','CalciumActivityAfter',
                                                       'AUC_calciumBefore','AUC_calciumDuring','AUC_calciumAfter','SpikeActivityPreference','SpikeActivityBefore','SpikeActivityDuring','SpikeActivityAfter'])
        DS_GlobalResults= pd.DataFrame(data, columns=['Mice', 'Session','Session_Time','Unique_Unit','UnitNumber','UnitValue','Drug',
                                                      'DSStatut','DSLocation','DSNumber','DSDuration','Spdl_in_DS','CalciumActivityPreference','CalciumActivityBefore',
                                                      'CalciumActivityDuring','CalciumActivityAfter','AUC_calciumBefore','AUC_calciumDuring','AUC_calciumAfter',
                                                      'SpikeActivityPreference','SpikeActivityBefore','SpikeActivityDuring','SpikeActivityAfter'])
        for drug in Drugs: 
            for coup in Coupling:
                for ctx in CTX:            
                    locals()[f'dict_All_ActivityCa_{coup}SPDL{ctx}_{drug}']={}
                    locals()[f'dict_All_ActivitySp_{coup}SPDL{ctx}_{drug}']={}
                    locals()[f'dict_All_ActivityCa_{coup}DS{ctx}_{drug}']={}
                    locals()[f'dict_All_ActivitySp_{coup}DS{ctx}_{drug}']={}
                    if coup=='Coupled':
                        locals()[f'dict_All_ActivityCa_{coup}SWR{ctx}_{drug}']={}
                        locals()[f'dict_All_ActivitySp_{coup}SWR{ctx}_{drug}']={}
                    else: 
                        locals()[f'dict_All_ActivityCa_{coup}SWR_{drug}']={}
                        locals()[f'dict_All_ActivitySp_{coup}SWR_{drug}']={}

        previousEndTime=0
        InitialStartTime=0

        for session in list(dict_Stamps.keys()):    
            cCoupled=0
            cUnCoupled=0
            cGlobal=0
            cLocalS1=0
            cLocalPFC=0

            cCoupledSWR=0
            cUnCoupledSWR=0    
            
            cCoupledDS=0
            cUnCoupledDS=0 
            cDSS1=0
            cDSPFC=0
            
            drug=os.path.basename(os.path.dirname(dict_Path[session])) if DrugExperiment else 'Baseline'

            # Start time & freq miniscope

            StartTime = list(dict_Stamps[session][0])[0] # in seconds
            StartTimeO = StartTime
            minian_freq=list(dict_Stamps[session][0])[2] # in Hz
            TimeStamps_miniscope=dict_StampsMiniscope[session]
            TimeStamps_miniscope["Time Stamp (ms)"]=TimeStamps_miniscope["Time Stamp (ms)"] + (StartTimeO*1000)

            minian_freq=round(1/np.mean(np.diff(np.array(TimeStamps_miniscope["Time Stamp (ms)"])/1000)))

            freqLFP=1000

            if minian_freq>=20: # should only remove 1 session                

                # Adjust the StartTime if subsessions

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
                S=dict_Spike[session] 

                Calcium = pd.DataFrame(C, index=C['unit_id'])
                Spike = pd.DataFrame(S, index=S['unit_id'])

                unit_to_drop=dict_TodropFile[session]    
                for u in unit_to_drop: 
                    Calcium=Calcium.drop(index=u) if u in Calcium.index else Calcium #need to know why
                    Spike=Spike.drop(index=u) if u in Spike.index else Spike


                indexMappList=B[session]
                kept_uniq_unit_List=[]
                for unit in Calcium.index:
                    indexMapp = np.where(indexMappList == unit)[0]
                    kept_uniq_unit_List.append(str(indexMapp))
                    
                nb_unit=len(Calcium)
                
              
                Carray=Calcium.values.T.astype(float)
                Sarray=Spike.values.T.astype(float)

                StartFrame_msec=TimeStamps_miniscope['Time Stamp (ms)'][TimeStamps_miniscope['Frame Number'][firstframe]]
                LastFrame_msec=TimeStamps_miniscope['Time Stamp (ms)'][TimeStamps_miniscope['Frame Number'][firstframe+len(Calcium.T)-1]]
                TS_miniscope_sub=TimeStamps_miniscope['Time Stamp (ms)'].iloc[firstframe:firstframe+len(Calcium.T)]
                rec_dur=len(Calcium.T)

                rec_dur_sec= (LastFrame_msec - StartFrame_msec)/1000
                
                nb_of_previousframe=firstframe

                firstframe+=rec_dur

                # Deal with dropped frames (failure to acquire miniscope images)

                list_droppedframes = literal_eval(dict_Stamps[session][0][3])    

                numbdropfr= 0   
                droppedframes_inrec=[]
                for item in list_droppedframes: 
                    if item < (int(StartTimeMiniscope*minian_freq) + rec_dur) and item > int(StartTimeMiniscope*minian_freq):
                        numbdropfr+=1                        

                EndTime = StartTime + rec_dur_sec # (upd_rec_dur/minian_freq) # in seconds
                previousEndTime=EndTime 

                print(session, ': starts at', round(StartTime,1), 's & ends at', round(EndTime,1), 's (', round(rec_dur_sec,1), 's duration, ', numbdropfr, 'dropped frames, minian frequency =', minian_freq, 'Hz, drug = ', drug, ')...') 

                sentence1= f"... kept values = {kept_uniq_unit_List}"
                print(sentence1) 

                # Zscore traces

                #Carray=zscore(Carray, axis=0)
                #Sarray=zscore(Sarray, axis=0)

                if nb_unit==0:
                    continue  # next iteration

                # Align Oscillations to miniscope start 

                SpipropO=dict_Spindleprop[session]
                SpipropM=SpipropO.copy()
                SWRpropO=dict_SWRprop[session]
                SWRpropM=SWRpropO.copy()
                DSpropO=dict_DSprop[session]
                DSpropM=DSpropO.copy()

                SpipropM=SpipropM[SpipropM['start time']> StartFrame_msec]
                SpipropTrunc=SpipropM[SpipropM['end time']< LastFrame_msec]
                SWRpropM=SWRpropM[SWRpropM['start time']> StartFrame_msec]
                SWRpropTrunc=SWRpropM[SWRpropM['end time']< LastFrame_msec]
                DSpropM=DSpropM[DSpropM['start time']> StartFrame_msec]
                DSpropTrunc=DSpropM[DSpropM['end time']< LastFrame_msec]
                
                HalfSpdl = round(durationSpdl*minian_freq)
                HalfSWR = round(durationSWR*minian_freq)
                HalfDS = round(durationDS*minian_freq)

                nb_spindle = SpipropTrunc.shape[0]
                nb_swr = SWRpropTrunc.shape[0]
                nb_ds = DSpropTrunc.shape[0]

                for unit in range(nb_unit): # for each kept units (cause Cseries/Sseries only have kept units)

                    Carray_unit =Carray[:,unit]
                    Sarray_unit =Sarray[:,unit]
                    #peaks, _ = find_peaks(Sarray_unit)#, height=np.std(SpTrace))
                    #Sarray_unit=np.zeros(len(Sarray_unit))
                    #Sarray_unit[peaks]=1

                    #######################################################################################
                                                        # for SPDLs #
                    #######################################################################################
                    for coup in Coupling:
                        for ctx in CTX:            
                            locals()[f'ActivityCa_{coup}Spin{ctx}']=[] #For each unit 
                            locals()[f'ActivitySp_{coup}Spin{ctx}']=[] #For each unit  

                    for Pspin in SpipropTrunc.index: 
                        
                        # Get the calcium and spike trace associated with the spdl
        
                        startSpi=SpipropTrunc.loc[Pspin, "start time"]                
                        endSpi=SpipropTrunc.loc[Pspin, "end time"]    
                        ctxSpi=SpipropTrunc.loc[Pspin, "CTX"]                
                        diffSpi=SpipropTrunc.loc[Pspin, "LocalGlobal"]                
                        StartLocSpi=SpipropTrunc.loc[Pspin, "StartingLoc"]   
                                    

                        TooEarlySpdl=startSpi-durationSpdl*1000<StartFrame_msec # too close to the begining of the recording
                        TooLateSpdl=startSpi+durationSpdl*1000>LastFrame_msec # too close to the end of the recording
                        
                        if TooEarlySpdl or TooLateSpdl:
                            print("/!\ Spindle too close to the begining/end of the recording,", session, ", Spdl n°", Pspin, ", Start Spdl =", round(startSpi/1000,1), "s") if unit==0 else None            
                        else:

                            if ctxSpi=='S1':
                                cLocalS1+=1 if unit==0 else 0
                            elif ctxSpi=='PFC': 
                                cLocalPFC+=1 if unit==0 else 0
                            elif ctxSpi=='S1PFC': 
                                cGlobal+=1 if unit==0 else 0    

                            # Find the index of the closest value in the column
                            Frame_Spindle_start_all = (TS_miniscope_sub - startSpi).abs().idxmin()
                            Frame_Spindle_start=Frame_Spindle_start_all-nb_of_previousframe

                            CaTrace = list(Carray_unit[Frame_Spindle_start-HalfSpdl:Frame_Spindle_start+HalfSpdl])
                            SpTrace = list(Sarray_unit[Frame_Spindle_start-HalfSpdl:Frame_Spindle_start+HalfSpdl]) 

                            ActivityCa_Spin=locals()[f'ActivityCa_Spin{ctxSpi}']
                            ActivitySp_Spin=locals()[f'ActivitySp_Spin{ctxSpi}']
                            ActivityCa_Spin.append(CaTrace)
                            ActivitySp_Spin.append(SpTrace)               

                            # Define if that spindle is coupled with a SWR or not

                            Spdl_statut=[]
                            startSWRList = list(pd.Series(SWRpropTrunc["start time"]))
                            if len(startSWRList)>0:
                                startClosest_SWR_idx = (np.abs(startSWRList - startSpi)).argmin()
                                startClosest_SWR = startSWRList[startClosest_SWR_idx]
                                distance = abs(startClosest_SWR - startSpi)
                                IsTrue=is_between(startSWRList,startSpi, endSpi)
                                if (distance < before) or IsTrue:
                                    Spdl_statut = 'Coupled'
                                    cCoupled+=1 if unit==0 else 0                              
                                else:
                                    Spdl_statut= 'UnCoupled'
                                    cUnCoupled+=1 if unit==0 else 0
                            else:
                                Spdl_statut= 'UnCoupled'
                                cUnCoupled+=1 if unit==0 else 0

                            ActivityCa_SpinCp=locals()[f'ActivityCa_{Spdl_statut}Spin{ctxSpi}']
                            ActivitySp_SpinCp=locals()[f'ActivitySp_{Spdl_statut}Spin{ctxSpi}']
                            ActivityCa_SpinCp.append(CaTrace)
                            ActivitySp_SpinCp.append(SpTrace)
                            
                            # Fill the big summary table Spindles_GlobalResults

                            Spindles_GlobalResults.loc[counter, 'Mice'] = mice
                            Spindles_GlobalResults.loc[counter, 'Session'] = session
                            Spindles_GlobalResults.loc[counter, 'Session_Time'] = None 

                            indexMapp = np.where(B[session] == Calcium.index[unit])[0]
                            Spindles_GlobalResults.loc[counter, 'Unique_Unit'] = indexMapp 
                            Spindles_GlobalResults.loc[counter, 'UnitNumber'] = unit 
                            Spindles_GlobalResults.loc[counter, 'UnitValue'] = Calcium.index[unit] 
                            
                            Spindles_GlobalResults.loc[counter, 'Drug'] =  os.path.basename(os.path.dirname(dict_Path[session])) if DrugExperiment else 'Baseline'

                            Spindles_GlobalResults.loc[counter, 'SpdlStatut'] = Spdl_statut
                            Spindles_GlobalResults.loc[counter, 'SpdlStartLocation'] = StartLocSpi
                            Spindles_GlobalResults.loc[counter, 'GlobalSpindle'] = diffSpi
                            Spindles_GlobalResults.loc[counter, 'SpdlNumber'] = Pspin
                            Spindles_GlobalResults.loc[counter, 'SpdlDuration'] = endSpi- startSpi                        
                            Spindles_GlobalResults.loc[counter, 'SWR_inside_Spdl'] = IsTrue
                            
                            # Activity before/ during/after oscillation

                            durOsc=round((endSpi- startSpi)/1000*minian_freq)
                            durOsc=round(1*minian_freq) # 1 seconds
                            TooEarlySpdl=startSpi/1000<durOsc/minian_freq # too close to the begining of the recording
                            TooLateSpdl=startSpi/1000+(durOsc/minian_freq*2)>LastFrame_msec/1000 # too close to the end of the recording
                            if TooEarlySpdl or TooLateSpdl:
                                print("/!\ Spindle too close to the begining/end of the recording,", session, ", Spdl n°", Pspin, ", Start Spdl =", round(startSpi/1000,1), "s, Spdl duration=", round(durOsc/minian_freq, 1), 's') if unit==0 else None            
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
                    for ctx in CTX: 
                        for coup in Coupling: 
                            # All Ca traces for each spindles per Unique unit (according to cross-registration)
                            ActivityCa = locals()[f'ActivityCa_{coup}Spin{ctx}']
                            dict_All_ActivityCa = locals()[f'dict_All_ActivityCa_{coup}SPDL{ctx}_{drug}']
                            ActivitySp = locals()[f'ActivitySp_{coup}Spin{ctx}']
                            dict_All_ActivitySp = locals()[f'dict_All_ActivitySp_{coup}SPDL{ctx}_{drug}']
                            if len(indexMapp) > 0: #not empty --> cause some units are not in the cross registration..! Need to know why    
                                if len(ActivityCa)>0 :                                
                                    if np.shape(np.array(ActivityCa))[1] == int(norm_freq*durationSpdl*2):  #normalize traces to the same frequency rate         
                                        ActivityCa= np.reshape(np.array(ActivityCa), (-1, len(np.array(ActivityCa)))) if np.ndim(ActivityCa) == 1 else np.array(ActivityCa)    
                                        ActivitySp= np.reshape(np.array(ActivitySp), (-1, len(np.array(ActivitySp)))) if np.ndim(ActivitySp) == 1 else np.array(ActivitySp)    
                                        key=mice + str(indexMapp).replace('[','').replace(']','')
                                        dict_All_ActivityCa[key] = np.append(dict_All_ActivityCa[key], np.array(ActivityCa), axis=0) if key in dict_All_ActivityCa else np.array(ActivityCa)
                                        dict_All_ActivitySp[key] = np.append(dict_All_ActivitySp[key], np.array(ActivitySp), axis=0) if key in dict_All_ActivitySp else np.array(ActivitySp)
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

                                        dataO = np.array(ActivitySp)
                                        data= np.repeat(dataO, 2, axis=0) if dataO.shape[0] == 1 else dataO
                                        x_mesh, y_mesh = np.meshgrid(np.arange(data.shape[1]), np.arange(data.shape[0]))
                                        x_new_mesh, y_new_mesh = np.meshgrid(np.linspace(0, data.shape[1] - 1, int(norm_freq*durationSpdl*2)), np.linspace(0, data.shape[0] - 1, np.shape(data)[0]))
                                        resampled_dataO = griddata((x_mesh.flatten(), y_mesh.flatten()), data.flatten(), (x_new_mesh, y_new_mesh), method='nearest')
                                        resampled_data= resampled_dataO[0,:] if dataO.shape[0] == 1 else resampled_dataO
                                        resampled_data= np.reshape(resampled_data, (-1, len(resampled_data))) if np.ndim(resampled_data) == 1 else resampled_data
                                        key=mice + str(indexMapp).replace('[','').replace(']','')
                                        dict_All_ActivitySp[key] = np.append(dict_All_ActivitySp[key], np.array(resampled_data), axis=0) if key in dict_All_ActivitySp else np.array(resampled_data)
                            #else: 
                                #print(f"/!\ Cell idx {unit} not in the cross registration")
                                    
                    #######################################################################################
                                                        # for SWRs #
                    #######################################################################################
                    for coup in Coupling:
                        for ctx in CTX:   
                            if coup =='Coupled':     
                                locals()[f'ActivityCa_{coup}swr{ctx}']=[] #For each unit 
                                locals()[f'ActivitySp_{coup}swr{ctx}']=[] #For each unit  
                            else:
                                locals()[f'ActivityCa_{coup}swr']=[] #For each unit 
                                locals()[f'ActivitySp_{coup}swr']=[] #For each unit  

                    for Pswr in SWRpropTrunc.index: 

                        # Get the calcium and spike trace associated with the SWR
                        startSwr=SWRpropTrunc.loc[Pswr, "start time"]
                        endSwr=SWRpropTrunc.loc[Pswr, "end time"]
                        
                        TooEarlySWR=startSwr-durationSWR*1000<StartFrame_msec # too close to the begining of the recording
                        TooLateSWR=startSwr+durationSWR*1000>LastFrame_msec # too close to the end of the recording
                        if TooEarlySWR or TooLateSWR:
                            print("/!\ SWR too close to the begining/end of the recording,", session, ", SWR n°", Pswr, ", Start SWR =",  round(startSwr/1000,1), "s") if unit==0 else None 
                        else:

                            Frame_SWR_start_all = (TS_miniscope_sub - startSwr).abs().idxmin()
                            Frame_SWR_start=Frame_SWR_start_all-nb_of_previousframe

                            CaTrace = list(Carray_unit[Frame_SWR_start-HalfSWR:Frame_SWR_start+HalfSWR])
                            SpTrace = list(Sarray_unit[Frame_SWR_start-HalfSWR:Frame_SWR_start+HalfSWR]) 

                            ActivityCa_swr=locals()[f'ActivityCa_swr']
                            ActivitySp_swr=locals()[f'ActivitySp_swr']

                            ActivityCa_swr.append(CaTrace)
                            ActivitySp_swr.append(SpTrace)

                            # Define if that SWR is coupled with a SPDL or not

                            SWR_statut=[]
                            startSpiList = list(pd.Series(SpipropTrunc["start time"]))
                            endSpiList = list(pd.Series(SpipropTrunc["end time"]))
                            ctxSpiList = list(pd.Series(SpipropTrunc["CTX"]))
                            if len(startSpiList)>0:
                                startClosest_Spdl_idx = (np.abs(startSpiList - startSwr)).argmin()
                                startClosest_Spi = startSpiList[startClosest_Spdl_idx]
                                endClosest_Spi=endSpiList[startClosest_Spdl_idx]
                                ctxSpi=ctxSpiList[startClosest_Spdl_idx]
                                distance = abs(startClosest_Spi - startSwr) #  + StartTimeIndexSpi]  
                                IsTrue = startSwr>startClosest_Spi and startSwr<endClosest_Spi #SWR inside the Spindle
                                if distance<before or IsTrue:
                                    SWR_statut = 'Coupled'
                                    cCoupledSWR+=1 if unit==0 else 0
                                else:
                                    SWR_statut= 'UnCoupled'
                                    cUnCoupledSWR+=1 if unit==0 else 0
                                    ctxSpi=''
                            else: 
                                SWR_statut= 'UnCoupled'
                                cUnCoupledSWR+=1 if unit==0 else 0
                                ctxSpi=''

                            ActivityCa_swrCp=locals()[f'ActivityCa_{SWR_statut}swr{ctxSpi}']
                            ActivitySp_swrCp=locals()[f'ActivitySp_{SWR_statut}swr{ctxSpi}']
                            ActivityCa_swrCp.append(CaTrace)
                            ActivitySp_swrCp.append(SpTrace)
                            
                            # Fill the big summary table SWR_GlobalResults

                            SWR_GlobalResults.loc[counter2, 'Mice'] = mice
                            SWR_GlobalResults.loc[counter2, 'Session'] = session
                            SWR_GlobalResults.loc[counter2, 'Session_Time'] = None 
                            indexMapp = np.where(B[session] == Calcium.index[unit])[0]
                            SWR_GlobalResults.loc[counter2, 'Unique_Unit'] = indexMapp 
                            SWR_GlobalResults.loc[counter2, 'UnitNumber'] = unit 
                            SWR_GlobalResults.loc[counter2, 'UnitValue'] = Calcium.index[unit] 
                            
                            SWR_GlobalResults.loc[counter2, 'Drug'] = os.path.basename(os.path.dirname(dict_Path[session])) if DrugExperiment else 'Baseline'

                            SWR_GlobalResults.loc[counter2, 'SWRStatut'] = SWR_statut
                            SWR_GlobalResults.loc[counter2, 'SpdlLoc'] = ctxSpi
                            SWR_GlobalResults.loc[counter2, 'SWRNumber'] = Pswr
                            SWR_GlobalResults.loc[counter2, 'SWRDuration'] = endSwr- startSwr
                            SWR_GlobalResults.loc[counter2, 'SWR_inside_Spdl'] = IsTrue

                            # Activity before/ during/after oscillation

                            durOsc=round((endSwr- startSwr)/1000*minian_freq)
                            durOsc=round(.5*minian_freq) # 1 seconds
                            TooEarlySWR=startSwr/1000<durOsc/minian_freq # too close to the begining of the recording
                            TooLateSWR=startSwr/1000+(durOsc/minian_freq*2)>LastFrame_msec/1000 # too close to the end of the recording
                            #TooEarlySWR=startSwr-durOsc*1000<StartFrame_msec # too close to the begining of the recording
                            #TooLateSWR=startSwr+durOsc*1000>LastFrame_msec # too close to the end of the recording
                            if TooEarlySWR or TooLateSWR:
                                print("/!\ SWR too close to the begining/end of the recording,", session, ", SWR n°", Pswr, ", Start SWR =", round(startSwr/1000,1), "s, SWR duration=", round(durOsc/minian_freq, 1), 's') if unit==0 else None            
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
                    for ctx in CTX: 
                        for coup in Coupling: 
                            if coup=='Coupled':
                                # All Ca traces for each spindles per Unique unit (according to cross-registration)
                                ActivityCa = locals()[f'ActivityCa_{coup}swr{ctx}']
                                dict_All_ActivityCa = locals()[f'dict_All_ActivityCa_{coup}SWR{ctx}_{drug}']
                                ActivitySp = locals()[f'ActivitySp_{coup}swr{ctx}']
                                dict_All_ActivitySp = locals()[f'dict_All_ActivitySp_{coup}SWR{ctx}_{drug}']
                            else:    
                                # All Ca traces for each spindles per Unique unit (according to cross-registration)
                                ActivityCa = locals()[f'ActivityCa_{coup}swr']
                                dict_All_ActivityCa = locals()[f'dict_All_ActivityCa_{coup}SWR_{drug}']
                                ActivitySp = locals()[f'ActivitySp_{coup}swr']
                                dict_All_ActivitySp = locals()[f'dict_All_ActivitySp_{coup}SWR_{drug}']
                            if len(indexMapp) > 0: #not empty --> cause some units are not in the cross registration..! Need to know why 
                                if len(ActivityCa)>0 :                                  
                                    if np.shape(np.array(ActivityCa))[1] == int(norm_freq*durationSWR*2):   #normalize traces to the same frequency rate    
                                        ActivityCa= np.reshape(np.array(ActivityCa), (-1, len(np.array(ActivityCa)))) if np.ndim(ActivityCa) == 1 else np.array(ActivityCa)    
                                        ActivitySp= np.reshape(np.array(ActivitySp), (-1, len(np.array(ActivitySp)))) if np.ndim(ActivitySp) == 1 else np.array(ActivitySp)    
                                        key=mice + str(indexMapp).replace('[','').replace(']','')
                                        dict_All_ActivityCa[key] = np.append(dict_All_ActivityCa[key], np.array(ActivityCa), axis=0) if key in dict_All_ActivityCa else np.array(ActivityCa)
                                        dict_All_ActivitySp[key] = np.append(dict_All_ActivitySp[key], np.array(ActivitySp), axis=0) if key in dict_All_ActivitySp else np.array(ActivitySp)
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
                                                        # for DSs #
                    #######################################################################################
                    for coup in Coupling:
                        for ctx in CTX:            
                            locals()[f'ActivityCa_{coup}DS{ctx}']=[] #For each unit 
                            locals()[f'ActivitySp_{coup}DS{ctx}']=[] #For each unit  

                    for Pds in DSpropTrunc.index: 
                        
                        # Get the calcium and spike trace associated with the spdl
        
                        startDS=DSpropTrunc.loc[Pds, "start time"]                
                        endDS=DSpropTrunc.loc[Pds, "end time"]    
                        ctxDS=DSpropTrunc.loc[Pds, "CTX"]                
                                    
                        TooEarlyDS=startDS-durationDS*1000<StartFrame_msec # too close to the begining of the recording
                        TooLateDS=startDS+durationDS*1000>LastFrame_msec # too close to the end of the recording
                        
                        if TooEarlyDS or TooLateDS:
                            print("/!\ DS too close to the begining/end of the recording,", session, ", DS n°", Pds, ", Start DS =", round(startDS/1000,1), "s") if unit==0 else None            
                        else:

                            if ctxDS=='S1':
                                cDSS1+=1 if unit==0 else 0
                            elif ctxDS=='PFC': 
                                cDSPFC+=1 if unit==0 else 0   

                            # Find the index of the closest value in the column
                            Frame_DS_start_all = (TS_miniscope_sub - startDS).abs().idxmin()
                            Frame_DS_start=Frame_DS_start_all-nb_of_previousframe

                            CaTrace = list(Carray_unit[Frame_DS_start-HalfDS:Frame_DS_start+HalfDS])
                            SpTrace = list(Sarray_unit[Frame_DS_start-HalfDS:Frame_DS_start+HalfDS]) 

                            ActivityCa_DS=locals()[f'ActivityCa_DS{ctxDS}']
                            ActivitySp_DS=locals()[f'ActivitySp_DS{ctxDS}']
                            ActivityCa_DS.append(CaTrace)
                            ActivitySp_DS.append(SpTrace)               

                            # Define if that spindle is coupled with a SWR or not

                            DS_statut=[]
                            startSpiList = list(pd.Series(SpipropTrunc["start time"]))
                            if len(startSpiList)>0:
                                startClosest_Spdl_idx = (np.abs(startSpiList - startDS)).argmin()
                                startClosest_Spdl = startSpiList[startClosest_Spdl_idx]
                                distance = abs(startClosest_Spdl - startDS)
                                IsTrue=is_between(startSpiList,startDS, endDS)
                                if (distance < before) or IsTrue:
                                    DS_statut = 'Coupled'
                                    cCoupledDS+=1 if unit==0 else 0                              
                                else:
                                    DS_statut= 'UnCoupled'
                                    cUnCoupledDS+=1 if unit==0 else 0
                            else:
                                DS_statut= 'UnCoupled'
                                cUnCoupledDS+=1 if unit==0 else 0

                            ActivityCa_DSCp=locals()[f'ActivityCa_{DS_statut}DS{ctxDS}']
                            ActivitySp_DSCp=locals()[f'ActivitySp_{DS_statut}DS{ctxDS}']
                            ActivityCa_DSCp.append(CaTrace)
                            ActivitySp_DSCp.append(SpTrace)
                            
                            # Fill the big summary table DSs_GlobalResults

                            DS_GlobalResults.loc[counter, 'Mice'] = mice
                            DS_GlobalResults.loc[counter, 'Session'] = session
                            DS_GlobalResults.loc[counter, 'Session_Time'] = None 

                            indexMapp = np.where(B[session] == Calcium.index[unit])[0]
                            DS_GlobalResults.loc[counter, 'Unique_Unit'] = indexMapp 
                            DS_GlobalResults.loc[counter, 'UnitNumber'] = unit 
                            DS_GlobalResults.loc[counter, 'UnitValue'] = Calcium.index[unit] 
                            
                            DS_GlobalResults.loc[counter, 'Drug'] =  os.path.basename(os.path.dirname(dict_Path[session])) if DrugExperiment else 'Baseline'

                            DS_GlobalResults.loc[counter, 'DSStatut'] = DS_statut
                            DS_GlobalResults.loc[counter, 'DSLocation'] = ctxDS
                            DS_GlobalResults.loc[counter, 'DSNumber'] = Pds
                            DS_GlobalResults.loc[counter, 'DSDuration'] = endDS- startDS                        
                            DS_GlobalResults.loc[counter, 'Spdl_in_DS'] = IsTrue
                            
                            # Activity before/ during/after oscillation

                            durOsc=round((endDS- startDS)/1000*minian_freq)
                            durOsc=round(.5*minian_freq) # 1 seconds
                            TooEarlyDS=startDS/1000<durOsc/minian_freq # too close to the begining of the recording
                            TooLateDS=startDS/1000+(durOsc/minian_freq*2)>LastFrame_msec/1000 # too close to the end of the recording
                            
                            if TooEarlyDS or TooLateDS:
                                print("/!\ DS too close to the begining/end of the recording,", session, ", DS n°", Pspin, ", Start DS =", round(startDS/1000,1), "s, DS duration=", round(durOsc/minian_freq, 1), 's') if unit==0 else None            
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
                                DS_GlobalResults.loc[counter, 'CalciumActivityPreference'] = pref
                                DS_GlobalResults.loc[counter, 'CalciumActivityBefore'] = ActBefore
                                DS_GlobalResults.loc[counter, 'CalciumActivityDuring'] = ActDuring
                                DS_GlobalResults.loc[counter, 'CalciumActivityAfter'] = ActAfter
                                DS_GlobalResults.loc[counter, 'AUC_calciumBefore'] = np.trapz(CaTrace[:durOsc],np.arange(0,len(CaTrace[:durOsc]),1))
                                DS_GlobalResults.loc[counter, 'AUC_calciumDuring'] = np.trapz(CaTrace[durOsc:durOsc*2],np.arange(0,len(CaTrace[durOsc:durOsc*2]),1))          
                                DS_GlobalResults.loc[counter, 'AUC_calciumAfter'] = np.trapz(CaTrace[durOsc*2:durOsc*3],np.arange(0,len(CaTrace[durOsc*2:durOsc*3]),1))          

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
                                DS_GlobalResults.loc[counter, 'SpikeActivityPreference'] = pref
                                DS_GlobalResults.loc[counter, 'SpikeActivityBefore'] = np.mean(SpTrace[:durOsc],0)
                                DS_GlobalResults.loc[counter, 'SpikeActivityDuring'] = np.mean(SpTrace[durOsc:durOsc*2],0)
                                DS_GlobalResults.loc[counter, 'SpikeActivityAfter'] = np.mean(SpTrace[durOsc*2:durOsc*3],0)                         
                            counter+=1     

                    ## Peristimulus Time Histogram 
                    for ctx in CTX: 
                        for coup in Coupling: 
                            # All Ca traces for each spindles per Unique unit (according to cross-registration)
                            ActivityCa = locals()[f'ActivityCa_{coup}DS{ctx}']
                            dict_All_ActivityCa = locals()[f'dict_All_ActivityCa_{coup}DS{ctx}_{drug}']
                            ActivitySp = locals()[f'ActivitySp_{coup}DS{ctx}']
                            dict_All_ActivitySp = locals()[f'dict_All_ActivitySp_{coup}DS{ctx}_{drug}']
                            if len(indexMapp) > 0: #not empty --> cause some units are not in the cross registration..! Need to know why    
                                if len(ActivityCa)>0 :                                
                                    if np.shape(np.array(ActivityCa))[1] == int(norm_freq*durationDS*2):  #normalize traces to the same frequency rate         
                                        ActivityCa= np.reshape(np.array(ActivityCa), (-1, len(np.array(ActivityCa)))) if np.ndim(ActivityCa) == 1 else np.array(ActivityCa)    
                                        ActivitySp= np.reshape(np.array(ActivitySp), (-1, len(np.array(ActivitySp)))) if np.ndim(ActivitySp) == 1 else np.array(ActivitySp)    
                                        key=mice + str(indexMapp).replace('[','').replace(']','')
                                        dict_All_ActivityCa[key] = np.append(dict_All_ActivityCa[key], np.array(ActivityCa), axis=0) if key in dict_All_ActivityCa else np.array(ActivityCa)
                                        dict_All_ActivitySp[key] = np.append(dict_All_ActivitySp[key], np.array(ActivitySp), axis=0) if key in dict_All_ActivitySp else np.array(ActivitySp)
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

                                        dataO = np.array(ActivitySp)
                                        data= np.repeat(dataO, 2, axis=0) if dataO.shape[0] == 1 else dataO
                                        x_mesh, y_mesh = np.meshgrid(np.arange(data.shape[1]), np.arange(data.shape[0]))
                                        x_new_mesh, y_new_mesh = np.meshgrid(np.linspace(0, data.shape[1] - 1, int(norm_freq*durationDS*2)), np.linspace(0, data.shape[0] - 1, np.shape(data)[0]))
                                        resampled_dataO = griddata((x_mesh.flatten(), y_mesh.flatten()), data.flatten(), (x_new_mesh, y_new_mesh), method='nearest')
                                        resampled_data= resampled_dataO[0,:] if dataO.shape[0] == 1 else resampled_dataO
                                        resampled_data= np.reshape(resampled_data, (-1, len(resampled_data))) if np.ndim(resampled_data) == 1 else resampled_data
                                        key=mice + str(indexMapp).replace('[','').replace(']','')
                                        dict_All_ActivitySp[key] = np.append(dict_All_ActivitySp[key], np.array(resampled_data), axis=0) if key in dict_All_ActivitySp else np.array(resampled_data)
                            #else: 
                                #print(f"/!\ Cell idx {unit} not in the cross registration")
                                    
                sentence2=f"... {nb_spindle} spindles ({cCoupled} Coupled & {cUnCoupled} Uncoupled Spdl // {cGlobal} Global, {cLocalS1} LocalS1 & {cLocalPFC} LocalPFC) and {nb_swr} SWR detected ({cCoupledSWR} Coupled & {cUnCoupledSWR} Uncoupled SWR) and {nb_ds} DS ({cCoupledDS} Coupled & {cUnCoupledDS} Uncoupled DS //{cDSS1} LocalS1 & {cDSPFC} LocalPFC)"
                print(sentence2) 
            else:
                print(f'/!\ {session} not taken into account cause minian frequency = {minian_freq}')
                    

        #######################################################################################
                                # Save Spindles & SWR & DS analysis #
        #######################################################################################
        if saveexcel: 
            # Save the big summary table Spindles_GlobalResults
            filenameOut = folder_to_save / f'Spdl_Global_{mice}.xlsx'
            writer = pd.ExcelWriter(filenameOut)
            Spindles_GlobalResults.to_excel(writer)
            writer.close()

            # Save the big summary table SWR_GlobalResults
            filenameOut = folder_to_save / f'SWR_Global_{mice}.xlsx'
            writer = pd.ExcelWriter(filenameOut)
            SWR_GlobalResults.to_excel(writer)
            writer.close()
                        
            # Save the big summary table SWR_GlobalResults
            filenameOut = folder_to_save / f'DS_Global_{mice}.xlsx'
            writer = pd.ExcelWriter(filenameOut)
            DS_GlobalResults.to_excel(writer)
            writer.close()

        filenameOut = folder_to_save / f'Spdl_Global_{mice}.pkl'
        with open(filenameOut, 'wb') as pickle_file:
            pickle.dump(Spindles_GlobalResults, pickle_file)   
        
        filenameOut = folder_to_save / f'SWR_Global_{mice}.pkl'
        with open(filenameOut, 'wb') as pickle_file:
            pickle.dump(SWR_GlobalResults, pickle_file)

        filenameOut = folder_to_save / f'DS_Global_{mice}.pkl'
        with open(filenameOut, 'wb') as pickle_file:
            pickle.dump(DS_GlobalResults, pickle_file)

        # Do average Calcium & Spike results for Spindles & SWR Peristimulus Time Histogram 

        Data=['Ca', 'Sp']
        for data in Data:
            if saveexcel: 
                filenameOut = folder_to_save / f'Spdl_{data}PSTH_{mice}.xlsx'
                excel_writer = pd.ExcelWriter(filenameOut)
                filenameOut = folder_to_save / f'SWR_{data}PSTH_{mice}.xlsx'
                excel_writer = pd.ExcelWriter(filenameOut)                  
            for ctx in CTX: 
                for coup in Coupling:
                    for drug in Drugs:      
                        dict_All_Activity=locals()[f'dict_All_Activity{data}_{coup}SPDL{ctx}_{drug}']
                        filenameOut = folder_to_save / f'Spdl_{data}PSTH_{coup}{ctx}{drug}_{mice}.pkl' #keep each responses of each cells for all rec Spdl
                        with open(filenameOut, 'wb') as pickle_file:
                            pickle.dump(dict_All_Activity, pickle_file)

                        dict_All_Activity=locals()[f'dict_All_Activity{data}_{coup}DS{ctx}_{drug}']
                        filenameOut = folder_to_save / f'DS_{data}PSTH_{coup}{ctx}{drug}_{mice}.pkl' #keep each responses of each cells for all rec Spdl
                        with open(filenameOut, 'wb') as pickle_file:
                            pickle.dump(dict_All_Activity, pickle_file)

                        if coup=='Coupled' : 
                            dict_All_Activity=locals()[f'dict_All_Activity{data}_{coup}SWR{ctx}_{drug}']
                            filenameOut = folder_to_save / f'SWR_{data}PSTH_{coup}{ctx}{drug}_{mice}.pkl' #keep each responses of each cells for all rec SWR
                            with open(filenameOut, 'wb') as pickle_file:
                                pickle.dump(dict_All_Activity, pickle_file)
                        else:
                            dict_All_Activity=locals()[f'dict_All_Activity{data}_{coup}SWR_{drug}']
                            filenameOut = folder_to_save / f'SWR_{data}PSTH_{coup}{drug}_{mice}.pkl' #keep each responses of each cells for all rec SWR
                            with open(filenameOut, 'wb') as pickle_file:
                                pickle.dump(dict_All_Activity, pickle_file)

        sentence3=f"Nb of unique units for {mice} = {len(dict_All_Activity)}"
        print(sentence3)    