# # Associate Ca2+ signal with spindles for each session & subsessions using crossregistration

#######################################################################################
                            # Define Experiment type #
#######################################################################################

DrugExperiment=0 # 0 if Baseline Experiment / 1 if CGP Experiment

saveexcel=1

AHmethod=0 # 0 if using the method of Aurelie Hay (2025) / 1 if using the method of Audrey Hay (2025)

AnalysisID='_likeAH' if AHmethod else '_pynapple' # '_pynapple' if using the method of Aurelie Hay (2025) / '_minian' if using the method of Audrey Hay (2025)
suffix='no_overlapping_Osc'

CTX=['S1', 'PFC', 'S1PFC']
Coupling=['', 'UnCoupled', 'Coupled']

dir = "//10.69.168.1/crnldata/waking/audrey_hay/L1imaging/Analysed2025_AB/"

before = 500 # Max distance in ms between a SWR and a spindle to be considered as Precoupled
after = 1000 # Max distance in ms between a spindle and a SWR to be considered as Postcoupled
durationSpdl = 1 # number of sec before and after the Spdl onset taken into acount
durationSWR = 1 # number of sec before and after the SWR onset taken into acount

drugs=['baseline', 'CGP'] if DrugExperiment else ['baseline']

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
import time

from itertools import groupby
from ephyviewer import mkQApp, MainViewer, TraceViewer
from IPython.display import display
from ipyfilechooser import FileChooser
from datetime import datetime

import warnings
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

all_expe_types=['baseline', 'preCGP', 'postCGP'] if DrugExperiment else ['baseline', 'preCGP']

# Get the current date and time
FolderNameSave=str(datetime.now())[:19]
FolderNameSave = FolderNameSave.replace(" ", "_").replace(".", "_").replace(":", "_")

destination_folder= f"//10.69.168.1/crnldata/waking/audrey_hay/L1imaging/Analysed2025_AB/_CGP_analysis/Osc_{FolderNameSave}{suffix}{AnalysisID}" if DrugExperiment else f"//10.69.168.1/crnldata/waking/audrey_hay/L1imaging/Analysed2025_AB/_baseline_analysis/Osc_{FolderNameSave}{suffix}{AnalysisID}"
os.makedirs(destination_folder)
folder_to_save=Path(destination_folder)

logfile = open(f"{destination_folder}/output_log.txt", 'w')
sys.stdout = Tee(sys.stdout, logfile)  # print goes to both

# Copy the script file to the destination folder
source_script = "C:/Users/Manip2/SCRIPTS/CodePythonAudrey/CodePythonAurelie/HayLabAnalysis/python/_MINI&OE_6_compute_osc_activity.py"
destination_file_path = f"{destination_folder}/_MINI&OE_6_compute_osc_activity.txt"
shutil.copy(source_script, destination_file_path)

data = {}        
counter=0
counter2=0

norm_freq=20 # final miniscope frequency used for all recordings

Spindles_GlobalResults= pd.DataFrame(data, columns=['Mice', 'NeuronType','Session','Session_Date', 'Session_Time','Unique_Unit','UnitNumber','UnitValue','ExpeType','Drug', 'SpdlStatut','SpdlStartLocation',
                                                        'GlobalSpindle', 'SpdlNumber','SpdlDuration','SWR_inside_Spdl','DistanceSWR_Spdl','DistanceClosestSpdl','CalciumActivityPreference', 'CalciumActivityBefore',
                                                        'CalciumActivityDuring','CalciumActivityAfter', 'AUC_calciumBaseline','AUC_calciumBefore','AUC_calciumDuring','AUC_calciumAfter',
                                                        'SpikeActivityPreference','SpikeActivityBefore','SpikeActivityDuring','SpikeActivityAfter'])
SWR_GlobalResults= pd.DataFrame(data, columns=['Mice', 'NeuronType','Session','Session_Date','Session_Time','Unique_Unit','UnitNumber','UnitValue','ExpeType','Drug','SWRStatut','SpdlLoc', 'SWRNumber','SWRDuration',
                                                'SWR_inside_Spdl','DistanceSWR_Spdl','CalciumActivityPreference', 'CalciumActivityBefore','CalciumActivityDuring','CalciumActivityAfter',
                                                'AUC_calciumBaseline','AUC_calciumBefore','AUC_calciumDuring','AUC_calciumAfter','SpikeActivityPreference','SpikeActivityBefore','SpikeActivityDuring','SpikeActivityAfter'])

for dpath in Path(dir).glob('**/mappingsAB.pkl'):
    
    mappfile = open(dpath.parents[0]/ f'mappingsAB.pkl', 'rb')
    mapping = pickle.load(mappfile)
    mapping_sess = mapping['session']   
        
    centfile = open(dpath.parents[0]/ f'centsAB.pkl', 'rb')
    centroids = pickle.load(centfile) 

    mice = dpath.parents[0].parts[-1]
    NeuronType = dpath.parents[1].parts[-1]
    
    print(f"####################################################################################")
    print(f"################################### {mice} ####################################")
    print(f"####################################################################################")

    nb_minian_total=0

    subsessions = []
    dict_Calcium = {}
    dict_Deconv = {}
    dict_SWRprop = {}
    dict_DSprop={}
    dict_Spindleprop = {}
    dict_Stamps = {}
    dict_StampsMiniscope = {}
    dict_TodropFile = {}
    dict_Path={}
    dict_Scoring = {}

    for drug in drugs: 
        for coup in Coupling:
            for ctx in CTX:            
                locals()[f'dict_All_ActivityCa_{coup}SPDL{ctx}_{drug}']={}
                locals()[f'dict_All_ActivitySp_{coup}SPDL{ctx}_{drug}']={}
                if coup=='Coupled':
                    locals()[f'dict_All_ActivityCa_{coup}SWR{ctx}_{drug}']={}
                    locals()[f'dict_All_ActivitySp_{coup}SWR{ctx}_{drug}']={}
                else: 
                    locals()[f'dict_All_ActivityCa_{coup}SWR_{drug}']={}
                    locals()[f'dict_All_ActivitySp_{coup}SWR_{drug}']={}
    

    previousEndTime=0
    InitialStartTime=0


    minian_folders = [f for f in dpath.parents[0].rglob('minian') if f.is_dir()]

    for minianpath in minian_folders: # for each minian folders found in this mouse

        if any(p in all_expe_types for p in minianpath.parts): # have to be to the expe_types

            start = time.time()

            cCoupled=0
            cUnCoupled=0
            cGlobal=0
            cLocalS1=0
            cLocalPFC=0

            cCoupledSWR=0
            cUnCoupledSWR=0   
            
            session=minianpath.parents[0].name if len(minianpath.parts)==12 else minianpath.parents[1].name.split("_")[-1]
            session_path=minianpath.parents[2] if len(minianpath.parts)==12 else minianpath.parents[1]
            expe_type=minianpath.parents[3].name if len(minianpath.parts)==12 else minianpath.parents[2].name

            session_date= minianpath.parents[2].name.split("_")[0] if len(minianpath.parts)==12 else minianpath.parents[1].name.split("_")[0]
            session_time= minianpath.parents[2].name.split("_")[1] if len(minianpath.parts)==12 else minianpath.parents[1].name.split("_")[1]

            drug='CGP' if expe_type == 'postCGP' else 'baseline'
            dict_Path[session] = session_path
        
            minian_ds = open_minian(minianpath)
            dict_Calcium[session] = minian_ds['C'] # calcium traces 
            dict_Deconv[session] = minian_ds['S'] # estimated spikes deconvolved activity

            dict_Scoring[session]  = pd.read_csv(session_path/ f'OpenEphys/Sleep_Scoring_6Stages_5sEpoch.csv')
            dict_Stamps[session]  = pd.read_excel(session_path/ f'SynchroFileCorrect.xlsx')
            dict_StampsMiniscope[session]  = pd.read_csv(list(session_path.glob('*V4_Miniscope/timeStamps.csv'))[0])

            with open(minianpath / f'TodropFileAB.json', 'r') as f:
                unit_to_drop = json.load(f)
                dict_TodropFile[session]  = unit_to_drop


            SWRlist= pd.read_csv(session_path / f'OpenEphys/SWRproperties.csv' ) if AHmethod else pd.read_csv(session_path / f'OpenEphys/SWR_detection.csv' ) 
            SWRlist['toKeep'] = 'True' # SWRlist['toKeep'].astype(str)
            dict_SWRprop[session]  =SWRlist[SWRlist['toKeep'].isin(['VRAI', 'True'])]

            Spdllist = pd.read_csv(session_path / f'OpenEphys/Spindleproperties_S1&PFC.csv') if AHmethod else pd.read_csv(session_path / f'OpenEphys/SpindlesS1&PFC_detection.csv' ) 
            Spdllist['toKeep'] = 'True' # Spdllist['toKeep'].astype(str)
            dict_Spindleprop[session]  = Spdllist[Spdllist['toKeep'].isin(['VRAI', 'True'])]

            nb_minian_total+=1

            #######################################################################################
            # Distribute Ca2+ intensity & spikes to oscillations for each sessions/subsessions #
            #######################################################################################

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
                    StartTimeMiniscope=0 # start time of miniscope rec of that subsesssions relative to the start of the miniscope recording
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

                Calcium = pd.DataFrame(C_sel, index=C_sel['unit_id'])
                Deconv = pd.DataFrame(D_sel, index=D_sel['unit_id'])

                indexMappList=mapping_sess[session]
                kept_uniq_unit_List=[]
                for unit in Calcium.index:
                    indexMapp = np.where(indexMappList == unit)[0]
                    kept_uniq_unit_List.append(str(indexMapp))

                nb_unit=len(Calcium)
                if nb_unit==0:
                    continue  # next iteration
                
                
                Carray=Calcium.values.T.astype(float)
                Darray=Deconv.values.T.astype(float)

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
                    if item < (round(StartTimeMiniscope*minian_freq) + rec_dur_sec) and item > round(StartTimeMiniscope*minian_freq):
                        numbdropfr+=1                        

                EndTime = StartTime + rec_dur_sec # (upd_rec_dur/minian_freq) # in seconds
                previousEndTime=EndTime 

                print(session, ': starts at', round(StartTime,1), 's & ends at', round(EndTime,1), 's (', round(rec_dur_sec,1), 's duration, ', numbdropfr, 'dropped frames, minian frequency =', minian_freq, 'Hz, experiment type = ', expe_type, ')...') 
                sentence1= f"... kept values = {kept_uniq_unit_List}"
                print(sentence1)

                # Align Oscillations to miniscope start 

                SpipropO=dict_Spindleprop[session]
                SpipropM=SpipropO.copy()
                SWRpropO=dict_SWRprop[session]
                SWRpropM=SWRpropO.copy()

                SpipropM=SpipropM[SpipropM['start time']> StartFrame_msec]
                SpipropTrunc=SpipropM[SpipropM['end time']< LastFrame_msec]
                SWRpropM=SWRpropM[SWRpropM['start time']> StartFrame_msec]
                SWRpropTrunc=SWRpropM[SWRpropM['end time']< LastFrame_msec]
                
                HalfSpdl = round(durationSpdl*minian_freq)
                HalfSWR = round(durationSWR*minian_freq)

                nb_spindle = SpipropTrunc.shape[0]
                nb_swr = SWRpropTrunc.shape[0]

                print(f"... Loading time = {time.time() - start:.2f} seconds")

                start2 = time.time()

                unit_count=0
                for unit in range(nb_unit): # for each kept units (cause Cseries/Sseries only have kept units)
                    
                    indexMapp = np.where(mapping_sess[session] == Calcium.index[unit])[0]
                    
                    if len(indexMapp)>0 : # The neuron needs to be in the cross-registration

                        unit_count+=1
                        
                        Carray_unit =Carray[:,unit]
                        Darray_unit =Darray[:,unit]
                        peaks, _ = find_peaks(Darray_unit)#, height=np.std(SpTrace))
                        Sarray_unit=np.zeros(len(Darray_unit))
                        Sarray_unit[peaks]=1

                        #######################################################################################
                                                            # for SPDLs #
                        #######################################################################################
                        for coup in Coupling:
                            for ctx in CTX:            
                                locals()[f'ActivityCa_{coup}Spin{ctx}']=[] #For each unit 
                                locals()[f'ActivitySp_{coup}Spin{ctx}']=[] #For each unit  

                        start3 = time.time()
                        prevspin=[]

                        for Pspin in SpipropTrunc.index: 
                            
                            # Get the calcium and spike trace associated with the spdl
            
                            startSpi=SpipropTrunc.loc[Pspin, "start time"]    
                            endSpi=SpipropTrunc.loc[Pspin, "end time"]    
                            ctxSpi=SpipropTrunc.loc[Pspin, "CTX"]                
                            diffSpi=SpipropTrunc.loc[Pspin, "LocalGlobal"]                
                            StartLocSpi=SpipropTrunc.loc[Pspin, "StartingLoc"]   
                            closestSpdl=SpipropTrunc.loc[Pspin, "DistanceClosestSpdl"]    

                            endPreviousSpi=SpipropTrunc.loc[prevspin, "end time"] if prevspin else startSpi-durationSpdl*1000 #line and not index cause sometimes, index are not ordered    
                            prevspin=Pspin

                            if startSpi - endPreviousSpi >= durationSpdl*1000 : # if the spindle is not too close from the end of previous one 

                                TooEarlySpdl=startSpi-durationSpdl*1000<StartFrame_msec # too close to the begining of the recording
                                TooLateSpdl=startSpi+durationSpdl*1000>LastFrame_msec # too close to the end of the recording
                        
                                if TooEarlySpdl or TooLateSpdl:
                                    print("/!\ Spindle too close to the begining/end of the recording,", session, ", Spdl n°", Pspin, ", Start Spdl =", round(startSpi/1000,1), "s") if unit_count==1 else None            
                                else:

                                    if ctxSpi=='S1':
                                        cLocalS1+=1 if unit_count==1 else 0
                                    elif ctxSpi=='PFC': 
                                        cLocalPFC+=1 if unit_count==1 else 0
                                    elif ctxSpi=='S1PFC': 
                                        cGlobal+=1 if unit_count==1 else 0    

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
                                    delaiSWRSpdl=None
                                    if len(startSWRList)>0:
                                        startClosest_SWR_idx = (np.abs(startSWRList - startSpi)).argmin()
                                        startClosest_SWR = startSWRList[startClosest_SWR_idx]
                                        distance = abs(startClosest_SWR - startSpi)
                                        IsTrue=is_between(startSWRList,startSpi, endSpi)
                                        if (distance < before) or IsTrue:
                                            Spdl_statut = 'Coupled'
                                            cCoupled+=1 if unit_count==1 else 0   
                                            delaiSWRSpdl=startClosest_SWR - startSpi                           
                                        else:
                                            Spdl_statut= 'UnCoupled'
                                            cUnCoupled+=1 if unit_count==1 else 0
                                    else:
                                        Spdl_statut= 'UnCoupled'
                                        cUnCoupled+=1 if unit_count==1 else 0

                                    ActivityCa_SpinCp=locals()[f'ActivityCa_{Spdl_statut}Spin{ctxSpi}']
                                    ActivitySp_SpinCp=locals()[f'ActivitySp_{Spdl_statut}Spin{ctxSpi}']
                                    ActivityCa_SpinCp.append(CaTrace)
                                    ActivitySp_SpinCp.append(SpTrace)
                                    
                                    # Fill the big summary table Spindles_GlobalResults

                                    Spindles_GlobalResults.loc[counter, 'Mice'] = mice
                                    Spindles_GlobalResults.loc[counter, 'NeuronType'] = NeuronType

                                    Spindles_GlobalResults.loc[counter, 'Session'] = session
                                    Spindles_GlobalResults.loc[counter, 'Session_Date'] = session_date 
                                    Spindles_GlobalResults.loc[counter, 'Session_Time'] = session_time                    

                                    Spindles_GlobalResults.loc[counter, 'Unique_Unit'] = indexMapp 
                                    Spindles_GlobalResults.loc[counter, 'UnitNumber'] = unit 
                                    Spindles_GlobalResults.loc[counter, 'UnitValue'] = Calcium.index[unit] 
                                    
                                    Spindles_GlobalResults.loc[counter, 'ExpeType'] =  expe_type
                                    Spindles_GlobalResults.loc[counter, 'Drug'] =  drug

                                    Spindles_GlobalResults.loc[counter, 'SpdlStatut'] = Spdl_statut
                                    Spindles_GlobalResults.loc[counter, 'SpdlStartLocation'] = StartLocSpi
                                    Spindles_GlobalResults.loc[counter, 'GlobalSpindle'] = diffSpi
                                    Spindles_GlobalResults.loc[counter, 'SpdlNumber'] = Pspin
                                    Spindles_GlobalResults.loc[counter, 'SpdlDuration'] = endSpi- startSpi                        
                                    Spindles_GlobalResults.loc[counter, 'SWR_inside_Spdl'] = IsTrue
                                    Spindles_GlobalResults.loc[counter, 'DistanceSWR_Spdl'] = delaiSWRSpdl 
                                    Spindles_GlobalResults.loc[counter, 'DistanceClosestSpdl'] = closestSpdl 
                
                                    counter+=1    
                            #else: 
                            #    print("/!\ Spindle too close to the end of the previous one,", session, ", Spdl n°", Pspin, ", Start Spdl =", round(startSpi/1000,1), "s") if unit_count==1 else None
                        
                        ## Peristimulus Time Histogram 
                        start4 = time.time()
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

                        print(f'... {len(SpipropTrunc)} Spdl processed in {time.time() - start3:.2f} & PSTH done in {time.time() - start4:.2f} seconds for one cell') if unit_count==1 else None

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

                        start5 = time.time()

                        prevSWR=[]
                        for Pswr in SWRpropTrunc.index: 

                            # Get the calcium and spike trace associated with the SWR
                            startSwr=SWRpropTrunc.loc[Pswr, "start time"]
                            endSwr=SWRpropTrunc.loc[Pswr, "end time"]

                            endPreviousSwr=SWRpropTrunc.loc[prevSWR, "end time"] if prevSWR else startSwr-durationSWR*1000                             
                            prevSWR=Pswr

                            if startSwr - endPreviousSwr >= durationSWR*1000 : # if the spindle is not too close from the end of previous one 
                                
                                TooEarlySWR=startSwr-durationSWR*1000<StartFrame_msec # too close to the begining of the recording
                                TooLateSWR=startSwr+durationSWR*1000>LastFrame_msec # too close to the end of the recording
                                
                                if TooEarlySWR or TooLateSWR:
                                    print("/!\ SWR too close to the begining/end of the recording,", session, ", SWR n°", Pswr, ", Start SWR =",  round(startSwr/1000,1), "s") if unit_count==1 else None 
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
                                    delaiSWRSpdl=None
                                    if len(startSpiList)>0:
                                        startClosest_Spdl_idx = (np.abs(startSpiList - startSwr)).argmin()
                                        startClosest_Spi = startSpiList[startClosest_Spdl_idx]
                                        endClosest_Spi=endSpiList[startClosest_Spdl_idx]
                                        ctxSpi=ctxSpiList[startClosest_Spdl_idx]
                                        distance = abs(startClosest_Spi - startSwr) #  + StartTimeIndexSpi]  
                                        IsTrue = startSwr>startClosest_Spi and startSwr<endClosest_Spi #SWR inside the Spindle
                                        if distance<before or IsTrue:
                                            SWR_statut = 'Coupled'
                                            cCoupledSWR+=1 if unit_count==1 else 0
                                            delaiSWRSpdl=startClosest_Spi - startSwr                           
                                        else:
                                            SWR_statut= 'UnCoupled'
                                            cUnCoupledSWR+=1 if unit_count==1 else 0
                                            ctxSpi=''
                                    else: 
                                        SWR_statut= 'UnCoupled'
                                        cUnCoupledSWR+=1 if unit_count==1 else 0
                                        ctxSpi=''

                                    ActivityCa_swrCp=locals()[f'ActivityCa_{SWR_statut}swr{ctxSpi}']
                                    ActivitySp_swrCp=locals()[f'ActivitySp_{SWR_statut}swr{ctxSpi}']
                                    ActivityCa_swrCp.append(CaTrace)
                                    ActivitySp_swrCp.append(SpTrace)
                                    
                                    # Fill the big summary table SWR_GlobalResults

                                    SWR_GlobalResults.loc[counter2, 'Mice'] = mice
                                    SWR_GlobalResults.loc[counter2, 'NeuronType'] = NeuronType

                                    SWR_GlobalResults.loc[counter2, 'Session'] = session
                                    SWR_GlobalResults.loc[counter2, 'Session_Date'] = session_date 
                                    SWR_GlobalResults.loc[counter2, 'Session_Time'] = session_time                    

                                    SWR_GlobalResults.loc[counter2, 'Unique_Unit'] = indexMapp 
                                    SWR_GlobalResults.loc[counter2, 'UnitNumber'] = unit 
                                    SWR_GlobalResults.loc[counter2, 'UnitValue'] = Calcium.index[unit] 
                                    
                                    SWR_GlobalResults.loc[counter2, 'ExpeType'] =  expe_type
                                    SWR_GlobalResults.loc[counter2, 'Drug'] = drug

                                    SWR_GlobalResults.loc[counter2, 'SWRStatut'] = SWR_statut
                                    SWR_GlobalResults.loc[counter2, 'SpdlLoc'] = ctxSpi
                                    SWR_GlobalResults.loc[counter2, 'SWRNumber'] = Pswr
                                    SWR_GlobalResults.loc[counter2, 'SWRDuration'] = endSwr- startSwr
                                    SWR_GlobalResults.loc[counter2, 'SWR_inside_Spdl'] = IsTrue
                                    SWR_GlobalResults.loc[counter2, 'DistanceSWR_Spdl'] = delaiSWRSpdl 
                                    counter2+=1  
                            #else: 
                            #    print("/!\ SWR too close to the end of the previous one,", session, ", SWR n°", Pswr, ", Start SWR =",  round(startSwr/1000,1), "s") if unit==0 else None

                        ## Peristimulus Time Histogram 
                        start6 = time.time()
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
                        
                        print(f'... {len(SWRpropTrunc)} SWR processed in {time.time() - start5:.2f} seconds & PSTH done in {time.time() - start6:.2f} seconds for one cell') if unit_count==1 else None

                print(f"... The {unit_count} units of {session} analyzed in {time.time() - start2:.2f} seconds")

                sentence2=f"... {cCoupled+cUnCoupled}/{nb_spindle} spindles kept ({cCoupled} Coupled & {cUnCoupled} Uncoupled Spdl // {cGlobal} Global, {cLocalS1} LocalS1 & {cLocalPFC} LocalPFC) and {cCoupledSWR+cUnCoupledSWR}/{nb_swr} SWR kept ({cCoupledSWR} Coupled & {cUnCoupledSWR} Uncoupled SWR)"
                print(sentence2) 
            else:
                print(f'/!\ {session} not taken into account cause minian frequency = {minian_freq}')
                        
    #######################################################################################
                            # Save Spindles & SWR analysis #
    #######################################################################################
    # Do average Calcium & Spike results for Spindles & SWR Peristimulus Time Histogram 
    start7 = time.time()

    Data=['Ca', 'Sp']
    for data in Data:   
        for ctx in CTX: 
            for coup in Coupling:
                for drug in drugs:      
                    dict_All_Activity=locals()[f'dict_All_Activity{data}_{coup}SPDL{ctx}_{drug}']
                    filenameOut = folder_to_save / f'Spdl_{data}PSTH_{coup}{ctx}{drug}_{mice}.pkl' #keep each responses of each cells for all rec Spdl
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

    sentence3=f"{mice} data saved in {time.time() - start2:.2f} seconds"
    print(sentence3)   


start8 = time.time()
if saveexcel: 
    # Save the big summary table Spindles_GlobalResults
    writer = pd.ExcelWriter(folder_to_save / f'Spdl_Global.xlsx')
    Spindles_GlobalResults.to_excel(writer)
    writer.close()

    # Save the big summary table SWR_GlobalResults
    writer = pd.ExcelWriter(folder_to_save / f'SWR_Global.xlsx')
    SWR_GlobalResults.to_excel(writer)
    writer.close()

with open(folder_to_save / f'Spdl_Global.pkl', 'wb') as pickle_file:
    pickle.dump(Spindles_GlobalResults, pickle_file)   

with open(folder_to_save / f'SWR_Global.pkl', 'wb') as pickle_file:
    pickle.dump(SWR_GlobalResults, pickle_file)

print(f"Global matrix saved in {time.time() - start8:.2f} seconds")

sys.stdout = sys.__stdout__
logfile.close()