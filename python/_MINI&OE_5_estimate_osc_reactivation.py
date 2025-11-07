# # Associate Ca2+ signal with sleep stages for each session & subsessions using crossregistration

#######################################################################################
                            # Define Experiment type #
#######################################################################################

DrugExperiment = 0 # = 1 if CGP Experiment // DrugExperiment=0 if Baseline Experiment

AnalysisID = '_calcium_nozscore' 

Cell_Assembly_folder = "2_CellAssemblies_2025-10-25_22_07_22_calcium"

local = True
if local:
    dir = "//10.69.168.1/crnldata/forgetting/Aurelie/MiniscopeOE_data/L2_3_mice/"
else: 
    dir = "/crnldata/forgetting/Aurelie/MiniscopeOE_data/L2_3_mice/"

mapp = {1: 'AW',  2: 'QW', 3: 'NREM',  4: 'IS', 5: 'REM',  6: 'undefined'}
drugs = ['baseline']


permuted = 0 # 0 if not permuted / 1 if permuted

CTX=['S1', 'PFC', 'S1PFC']
CTX2=['CA1', 'RSC', 'CA1RSC']

Coupling=['', 'UnCoupled', 'PreCoupled', 'PostCoupled', 'PrePostCoupled']

before = 500 # Max distance in ms between a SWR and a spindle to be considered as Precoupled
after = 1000 # Max distance in ms between a spindle and a SWR to be considered as Postcoupled
durationSpdl = 1 # number of sec before and after the Spdl onset taken into acount
durationSWR = 1 # number of sec before and after the SWR onset taken into acount


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
import random
import time

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

def is_between(myList, starttime, endtime):
    IsTrue=False
    for ind in range(len(myList)):
        if starttime <= myList[ind] <= endtime:
            IsTrue=True
    return IsTrue

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

destination_folder= f"//10.69.168.1/crnldata/forgetting/Aurelie/MiniscopeOE_analysis/PlaceCells_experiment/4_OscReactivation_{FolderNameSave}{AnalysisID}" if local else f"/crnldata/forgetting/Aurelie/MiniscopeOE_analysis/PlaceCells_experiment/4_OscReactivation_{FolderNameSave}{AnalysisID}"
os.makedirs(destination_folder)
folder_to_save=Path(destination_folder)

logfile = open(f"{destination_folder}/output_log.txt", 'w')
sys.stdout = Tee(sys.stdout, logfile)  # print goes to both


# Copy the script file to the destination folder
source_script = "C:/Users/Manip2/SCRIPTS/HayLabAnalysis/python/_MINI&OE_5_estimate_Osc_reactivation.py" if local else "/home/aurelie.brecier/HayLabAnalysis/python/_MINI&OE_5_estimate_Osc_reactivation.py"
destination_file_path = f"{destination_folder}/_MINI&OE_5_estimate_Osc_reactivation.txt"
shutil.copy(source_script, destination_file_path)

data = {}        
counter=0
counter2=0

norm_freq=20 # final miniscope frequency used for all recordings

Spindles_GlobalResults= pd.DataFrame(data, columns=['Mice', 'NeuronType','Session','Session_Date', 'Session_Time','Assembly_ID', 'Cells_in_Assembly', 'ExpeType','SpdlStatut','SpdlStartLocation',
                                                        'GlobalSpindle', 'SWRLoc', 'SpdlNumber','SpdlDuration','SWR_inside_Spdl','DistanceSWR_Spdl','DistanceClosestSpdl','AUC_1stQuarter','AUC_2ndQuarter','AUC_3rdQuarter','AUC_4thQuarter',
                                                        'FR_1stQuarter','FR_2ndQuarter','FR_3rdQuarter','FR_4thQuarter'])
SWR_GlobalResults= pd.DataFrame(data, columns=['Mice', 'NeuronType','Session','Session_Date','Session_Time','Assembly_ID','Cells_in_Assembly','ExpeType','SWRStatut','SWRStartLocation',
                                               'GlobalSWR', 'SpdlLoc', 'SWRNumber','SWRDuration', 'SWR_inside_Spdl','DistanceSWR_Spdl','DistanceClosestSWR', 'AUC_1stQuarter','AUC_2ndQuarter','AUC_3rdQuarter','AUC_4thQuarter',
                                                        'FR_1stQuarter','FR_2ndQuarter','FR_3rdQuarter','FR_4thQuarter',])


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

    print(f'{len(CellAssembly_dict[mice].T)} cell assemblies detected in {mice}')

    nb_minian_total=0
    dict_Path={}
    dict_Calcium = {}
    dict_Deconv = {}
    dict_Scoring = {}
    dict_Stamps = {}
    dict_TodropFile = {}
    dict_StampsMiniscope = {}
    dict_SWRprop = {}
    dict_Spindleprop = {}

    for drug in drugs: 
        for coup in Coupling:            
            for ctx in CTX:            
                locals()[f'dict_All_ActivityCa_{coup}SPDL{ctx}_{drug}']={}
                if coup=='PreCoupled' or coup == 'PostCoupled' or coup == 'PrePostCoupled': 
                    locals()[f'dict_All_ActivityCa_{coup}SWR{ctx}_{drug}']={}
            for ctx2 in CTX2:       
                locals()[f'dict_All_ActivityCa_{coup}SWR{ctx2}_{drug}']={}     
                if coup=='PreCoupled' or coup == 'PostCoupled' or coup == 'PrePostCoupled': 
                    locals()[f'dict_All_ActivityCa_{coup}SWR{ctx2}_{drug}']={}     

    previousEndTime=0
    InitialStartTime=0

    minian_folders = [f for f in dpath.parents[0].rglob('minian') if f.is_dir()]

    for minianpath in minian_folders: # for each minian folders found in this mouse 
        drug= 'baseline' 

        if 1==1: # any(p in all_expe_types for p in minianpath.parts): # have to be to the expe_types

            cPreCoupled=0
            cPostCoupled=0
            cPrePostCoupled=0
            cUnCoupled=0
            cGlobal=0
            cLocalS1=0
            cLocalPFC=0
            cGlobalSWR=0
            cLocalCA1=0
            cLocalRSC=0

            cPreCoupledSWR=0
            cPostCoupledSWR=0
            cPrePostCoupledSWR=0
            cUnCoupledSWR=0   
            
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


            if session_type == 'SleepAfter':

                start = time.time()

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


                SWRlist= pd.read_csv(session_path.parent / f'OpenEphys/SWRCA1&RSC_finedetection.csv' ) 
                SWRlist['toKeep'] = 'True' # SWRlist['toKeep'].astype(str)
                dict_SWRprop[session]  =SWRlist[SWRlist['toKeep'].isin(['VRAI', 'True'])]

                Spdllist = pd.read_csv(session_path.parent / f'OpenEphys/SpindlesS1&PFC_finedetection.csv' ) 
                Spdllist['toKeep'] = 'True' # Spdllist['toKeep'].astype(str)
                dict_Spindleprop[session]  = Spdllist[Spdllist['toKeep'].isin(['VRAI', 'True'])]

                
                # Start time & freq miniscope
                
                dict_StampsMiniscope[session]  = pd.read_csv(list(session_path.glob('*V4_Miniscope/timeStamps.csv'))[0])

                freqLFP=1000    
                numbdropfr= 0         
                
                nb_minian_total+=1

                #######################################################################################
                # Distribute Ca2+ intensity & spikes to vigilance states for each sessions #
                #######################################################################################
                    
                # Adjust the StartTime if subsessions

                dict_Stamps[session]  = pd.read_excel(session_path/ f'SynchroFile.xlsx')     
                StartTime = list(dict_Stamps[session][0])[0] 
                StartTimeO = StartTime
                minian_freq=list(dict_Stamps[session][0])[2] # in Hz
                TimeStamps_miniscope=dict_StampsMiniscope[session]
                TimeStamps_miniscope["Time Stamp (ms)"]=TimeStamps_miniscope["Time Stamp (ms)"] + (StartTimeO*1000)

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

                Calcium = pd.DataFrame(C_sel, index=C_sel['unit_id'])
                Deconv = pd.DataFrame(D_sel, index=D_sel['unit_id'])

                indexMappList=mapping_sess[session]
                kept_uniq_unit_List=[]
                for unit in Calcium.index:
                    indexMapp = np.where(indexMappList == unit)[0]
                    kept_uniq_unit_List.append(f"{mice}{str(indexMapp).replace('[','').replace(']','')}")

                nb_unit=len(Calcium)
                if nb_unit==0:
                    print(f'no cells kept in the session: {session}')
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

                
                # Align Oscillations to miniscope start 

                SpipropO=dict_Spindleprop[session]
                SpipropM=SpipropO.copy()
                SWRpropO=dict_SWRprop[session]
                SWRpropM=SWRpropO.copy()

                SpipropM=SpipropM[SpipropM['start time']> StartFrame_msec]
                SpipropTrunc=SpipropM[SpipropM['end time']< LastFrame_msec]
                SWRpropM=SWRpropM[SWRpropM['start time']> StartFrame_msec]
                SWRpropTrunc=SWRpropM[SWRpropM['end time']< LastFrame_msec]
                
                

                nb_spindle = SpipropTrunc.shape[0]
                nb_swr = SWRpropTrunc.shape[0]
                

                print(session, ': starts at', round(StartTime,1), 's & ends at', round(EndTime,1), 's (', round(rec_dur_sec,1), 's duration, ', numbdropfr, 'dropped frames, minian frequency =', minian_freq, 'Hz, experiment type = ', session_type, ')...') 
                sentence1= f"... kept values = {kept_uniq_unit_List}"
                print(sentence1)

                # Define cell assemblies
                target_rate = 20 #Hz == 50ms bins
                Array_bin = resample_matrix(Carray, orig_rate=minian_freq, target_rate=target_rate)
                #Array_bin = resample_matrix(Darray, orig_rate=minian_freq, target_rate=target_rate)
                #Array_bin = bin_sum_fractional(Sarray, minian_freq, target_rate)       

                TS_bin = pd.Series(np.linspace(TS_miniscope_sub.iloc[0], TS_miniscope_sub.iloc[-1], len(Array_bin)))

                            
                # Add missing cells from indexMapp
                dN_=  zscore(Array_bin)
                df_array=pd.DataFrame(dN_.T, index=kept_uniq_unit_List)
                dN = df_array.reindex(CellAssembly_dict[mice].index, fill_value=0).to_numpy().T

                patterns = CellAssembly_dict[mice].to_numpy().T
                all_cells = CellAssembly_dict[mice].index.tolist()
                assemblies_ID = CellAssembly_dict[mice].columns.tolist()

                print(f"... Loading time = {time.time() - start:.2f} seconds")
                start2 = time.time()

                HalfSpdl = round(durationSpdl*target_rate)
                HalfSWR = round(durationSWR*target_rate)

                if len(patterns)>0:
                    patterns_th = patterns.copy()                    

                    cellAss_count=0
                    for ass in np.arange(np.shape(patterns)[0]):

                        cellAss_count+=1
                        
                        assembly_ID = assemblies_ID[ass]
                        cells_in_assembly = CellAssembly_members[mice][assembly_ID]

                        template = np.add.outer(abs(patterns[ass]),abs(patterns[ass]))      
                        template = np.nan_to_num(template, nan=0)
                        template = template - np.diag(np.diag(template))
                        tmp = dN @ template 
                        Reactivation_Strength= np.nansum(tmp * dN, axis=1) 
                        #Reactivation_Strength = zscore(Reactivation_Strength)

                        #######################################################################################
                                                            # for SPDLs #
                        #######################################################################################
                        for coup in Coupling:
                            for ctx in CTX:            
                                locals()[f'ActivityCa_{coup}Spin{ctx}']=[] #For each cell assembly 
                        
                        start3 = time.time()
                        prevspin=[]

                        for Pspin in SpipropTrunc.index: 
                            
                            # Get the calcium and spike trace associated with the spdl
                            lag = random.randint(-15*1000, -5*1000) if permuted else 0
    
                            startSpi=SpipropTrunc.loc[Pspin, "start time"] - lag
                            endSpi=SpipropTrunc.loc[Pspin, "end time"]  - lag  
                            ctxSpi=SpipropTrunc.loc[Pspin, "CTX"]                
                            diffSpi=SpipropTrunc.loc[Pspin, "LocalGlobal"]                
                            StartLocSpi=SpipropTrunc.loc[Pspin, "StartingLoc"]   
                            closestSpdl=SpipropTrunc.loc[Pspin, "DistanceClosestSpdl"]    

                            endPreviousSpi=SpipropTrunc.loc[prevspin, "end time"] if prevspin else startSpi-durationSpdl*1000 #line and not index cause sometimes, index are not ordered    
                            prevspin=Pspin

                            if 1==1: # startSpi - endPreviousSpi >= durationSpdl*500 : # if the spindle is not too close from the end of previous one 

                                TooEarlySpdl=startSpi-durationSpdl*1000<StartFrame_msec # too close to the begining of the recording
                                TooLateSpdl=startSpi+durationSpdl*1000>LastFrame_msec # too close to the end of the recording
                                
                            
                                if TooEarlySpdl or TooLateSpdl:
                                    print("/!\ Spindle too close to the begining/end of the recording,", session, ", Spdl n°", Pspin, ", Start Spdl =", round(startSpi/1000,1), "s") if cellAss_count==1 else None            
                                else:

                                    if ctxSpi=='S1':
                                        cLocalS1+=1 if cellAss_count==1 else 0
                                    elif ctxSpi=='PFC': 
                                        cLocalPFC+=1 if cellAss_count==1 else 0
                                    elif ctxSpi=='S1PFC': 
                                        cGlobal+=1 if cellAss_count==1 else 0    

                                    # Find the index of the closest value in the column
                                    Frame_Spindle_start_all = (TS_bin - startSpi).abs().idxmin()
                                    Frame_Spindle_start=Frame_Spindle_start_all#-nb_of_previousframe

                                    Trace = list(Reactivation_Strength[Frame_Spindle_start-HalfSpdl:Frame_Spindle_start+HalfSpdl])

                                    ActivityCa_Spin=locals()[f'ActivityCa_Spin{ctxSpi}']
                                    ActivityCa_Spin.append(Trace)

                                    # Define if that spindle is coupled with a SWR or not

                                    Spdl_statut=[]
                                    startSWRList = list(pd.Series(SWRpropTrunc["start time"]))
                                    ctxSWRList = list(pd.Series(SWRpropTrunc["CTX"]))
                                    delaiSWRSpdl=None
                                    Spdl_statut= 'UnCoupled'
                                    cUnCoupled+=1 if cellAss_count==1 else 0
                                    ctxSWR=''
                                    IsTrue= False
                                    if len(startSWRList)>0:                                        
                                        # if there is a SWR during the Spdl
                                        startClosest_SWR_idx= next((i for i, x in enumerate((startSWRList - startSpi)) if x >= 0), -1) #(np.abs(startSpiList - startSwr)).argmin()
                                        if startClosest_SWR_idx != -1:
                                            startClosest_SWR = startSWRList[startClosest_SWR_idx]
                                            IsTrue = is_between(startSWRList, startSpi, endSpi)
                                            if IsTrue: 
                                                ctxSWR=ctxSWRList[startClosest_SWR_idx]
                                                Spdl_statut = 'PostCoupled'
                                                cPostCoupled+=1 if cellAss_count==1 else 0  
                                                cUnCoupled-=1 if cellAss_count==1 else 0 
                                                delaiSWRSpdl = startClosest_SWR - startSpi 

                                        # if there is a SWR before the Spdl
                                        startClosest_SWR_idx= next((i for i in range(len((startSWRList - startSpi)) - 1, -1, -1) if (startSWRList - startSpi)[i] < 0), -1)
                                        if startClosest_SWR_idx != -1:
                                            startClosest_SWR = startSWRList[startClosest_SWR_idx]
                                            distance = abs(startClosest_SWR - startSpi)
                                            if (distance < before):
                                                ctxSWR=ctxSWRList[startClosest_SWR_idx] 
                                                if Spdl_statut == 'PostCoupled':
                                                    Spdl_statut = 'PrePostCoupled'
                                                    cPrePostCoupled+=1 if cellAss_count==1 else 0
                                                    cPostCoupled-=1 if cellAss_count==1 else 0   
                                                    delaiSWRSpdl = startClosest_SWR - startSpi      
                                                else:   
                                                    Spdl_statut = 'PreCoupled'
                                                    cPreCoupled+=1 if cellAss_count==1 else 0   
                                                    cUnCoupled-=1 if cellAss_count==1 else 0 
                                                    delaiSWRSpdl = startClosest_SWR - startSpi    

                                    ActivityCa_SpinCp=locals()[f'ActivityCa_{Spdl_statut}Spin{ctxSpi}']
                                    ActivityCa_SpinCp.append(Trace)


                                    # Fill the big summary table Spindles_GlobalResults

                                    Spindles_GlobalResults.loc[counter, 'Mice'] = mice
                                    Spindles_GlobalResults.loc[counter, 'NeuronType'] = NeuronType

                                    Spindles_GlobalResults.loc[counter, 'Session'] = session
                                    Spindles_GlobalResults.loc[counter, 'Session_Date'] = session_date 
                                    Spindles_GlobalResults.loc[counter, 'Session_Time'] = session_time                    

                                    Spindles_GlobalResults.loc[counter, 'Assembly_ID'] = assembly_ID
                                    Spindles_GlobalResults.loc[counter, 'Cells_in_Assembly'] = str(cells_in_assembly)                                    
                                    Spindles_GlobalResults.loc[counter, 'ExpeType'] =  session_type

                                    Spindles_GlobalResults.loc[counter, 'SpdlStatut'] = Spdl_statut
                                    Spindles_GlobalResults.loc[counter, 'SpdlStartLocation'] = StartLocSpi
                                    Spindles_GlobalResults.loc[counter, 'GlobalSpindle'] = diffSpi
                                    Spindles_GlobalResults.loc[counter, 'SWRLoc'] = ctxSWR
                                    Spindles_GlobalResults.loc[counter, 'SpdlNumber'] = Pspin
                                    Spindles_GlobalResults.loc[counter, 'SpdlDuration'] = endSpi- startSpi                        
                                    Spindles_GlobalResults.loc[counter, 'SWR_inside_Spdl'] = IsTrue
                                    Spindles_GlobalResults.loc[counter, 'DistanceSWR_Spdl'] = delaiSWRSpdl 
                                    Spindles_GlobalResults.loc[counter, 'DistanceClosestSpdl'] = closestSpdl 
                                    
                                    durOsc=round(durationSpdl*target_rate)  
                                    
                                    Spindles_GlobalResults.loc[counter, 'AUC_1stQuarter'] = 2 * np.trapz(Trace[:round(durOsc/2)],np.arange(0,len(Trace[:round(durOsc/2)]),1)) # *2 cause to be per sec
                                    Spindles_GlobalResults.loc[counter, 'AUC_2ndQuarter'] = 2 * np.trapz(Trace[round(durOsc/2):durOsc],np.arange(0,len(Trace[round(durOsc/2):durOsc]),1))
                                    Spindles_GlobalResults.loc[counter, 'AUC_3rdQuarter'] = 2 * np.trapz(Trace[durOsc:durOsc+round(durOsc/2)],np.arange(0,len(Trace[durOsc:durOsc+round(durOsc/2)]),1))    
                                    Spindles_GlobalResults.loc[counter, 'AUC_4thQuarter'] = 2 * np.trapz(Trace[durOsc+round(durOsc/2):],np.arange(0,len(Trace[durOsc+round(durOsc/2):]),1))          
                                    Spindles_GlobalResults.loc[counter, 'AUC_1stHalf'] = np.trapz(Trace[:durOsc],np.arange(0,len(Trace[:durOsc]),1))          
                                    Spindles_GlobalResults.loc[counter, 'AUC_2ndHalf'] = np.trapz(Trace[durOsc:],np.arange(0,len(Trace[durOsc:]),1))        

                                    Spindles_GlobalResults.loc[counter, 'FR_1stQuarter'] = 2 * np.sum(Trace[:round(durOsc/2)]) # *2 cause to be in Hz
                                    Spindles_GlobalResults.loc[counter, 'FR_2ndQuarter'] = 2 * np.sum(Trace[round(durOsc/2):durOsc])
                                    Spindles_GlobalResults.loc[counter, 'FR_3rdQuarter'] = 2 * np.sum(Trace[durOsc:durOsc+round(durOsc/2)])          
                                    Spindles_GlobalResults.loc[counter, 'FR_4thQuarter'] = 2 * np.sum(Trace[durOsc+round(durOsc/2):])          
                                    Spindles_GlobalResults.loc[counter, 'FR_1stHalf'] = np.sum(Trace[:durOsc])          
                                    Spindles_GlobalResults.loc[counter, 'FR_2ndHalf'] = np.sum(Trace[durOsc:])     

                            
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
                                if len(ActivityCa)>0 :                                
                                    if np.shape(np.array(ActivityCa))[1] == int(norm_freq*durationSpdl*2):  #normalize traces to the same frequency rate         
                                        ActivityCa= np.reshape(np.array(ActivityCa), (-1, len(np.array(ActivityCa)))) if np.ndim(ActivityCa) == 1 else np.array(ActivityCa)    
                                        dict_All_ActivityCa[assembly_ID] = np.append(dict_All_ActivityCa[assembly_ID], np.array(ActivityCa), axis=0) if assembly_ID in dict_All_ActivityCa else np.array(ActivityCa)
                                    else:
                                        dataO = np.array(ActivityCa)
                                        data= np.repeat(dataO, 2, axis=0) if dataO.shape[0] == 1 else dataO
                                        x_mesh, y_mesh = np.meshgrid(np.arange(data.shape[1]), np.arange(data.shape[0]))
                                        x_new_mesh, y_new_mesh = np.meshgrid(np.linspace(0, data.shape[1] - 1, int(norm_freq*durationSpdl*2)), np.linspace(0, data.shape[0] - 1, np.shape(data)[0]))
                                        resampled_dataO = griddata((x_mesh.flatten(), y_mesh.flatten()), data.flatten(), (x_new_mesh, y_new_mesh), method='linear')
                                        resampled_data= resampled_dataO[0,:] if dataO.shape[0] == 1 else resampled_dataO
                                        resampled_data= np.reshape(resampled_data, (-1, len(resampled_data))) if np.ndim(resampled_data) == 1 else resampled_data
                                        dict_All_ActivityCa[assembly_ID] = np.append(dict_All_ActivityCa[assembly_ID], np.array(resampled_data), axis=0) if assembly_ID in dict_All_ActivityCa else np.array(resampled_data)

                        print(f'... {len(SpipropTrunc)} Spdl processed in {time.time() - start3:.2f} & PSTH done in {time.time() - start4:.2f} seconds for one cell assembly') if cellAss_count==1 else None

                        #######################################################################################
                                                            # for SWRs #
                        #######################################################################################
                        
                        for coup in Coupling:
                            for ctx2 in CTX2:                                     
                                locals()[f'ActivityCa_{coup}swr{ctx2}']=[] #For each unit 

                        start5 = time.time()

                        prevSWR=[]
                        for Pswr in SWRpropTrunc.index: 

                            # Get the calcium and spike trace associated with the SWR
                            lag = random.randint(-15*1000, -5*1000) if permuted else 0

                            startSwr=SWRpropTrunc.loc[Pswr, "start time"] - lag
                            endSwr=SWRpropTrunc.loc[Pswr, "end time"] - lag
                            ctxSWR=SWRpropTrunc.loc[Pswr, "CTX"]                
                            diffSWR=SWRpropTrunc.loc[Pswr, "LocalGlobal"]                
                            StartLocSWR=SWRpropTrunc.loc[Pswr, "StartingLoc"] 
                            closestSWR=SWRpropTrunc.loc[Pswr, "DistanceClosestSWR"]     

                            endPreviousSwr=SWRpropTrunc.loc[prevSWR, "end time"] if prevSWR else startSwr-durationSWR*1000                             
                            prevSWR=Pswr

                            if 1==1: #startSwr - endPreviousSwr >= durationSWR*500 : # if the spindle is not too close from the end of previous one 
                                
                                TooEarlySWR=startSwr-durationSWR*1000<StartFrame_msec # too close to the begining of the recording
                                TooLateSWR=startSwr+durationSWR*1000>LastFrame_msec # too close to the end of the recording
                                
                                if TooEarlySWR or TooLateSWR:
                                    print("/!\ SWR too close to the begining/end of the recording,", session, ", SWR n°", Pswr, ", Start SWR =",  round(startSwr/1000,1), "s") if cellAss_count==1 else None 
                                else:

                                    if ctxSWR=='CA1':
                                        cLocalCA1+=1 if cellAss_count==1 else 0
                                    elif ctxSWR=='RSC': 
                                        cLocalRSC+=1 if cellAss_count==1 else 0
                                    elif ctxSWR=='CA1RSC': 
                                        cGlobalSWR+=1 if cellAss_count==1 else 0    

                                    Frame_SWR_start_all = (TS_bin - startSwr).abs().idxmin()
                                    Frame_SWR_start=Frame_SWR_start_all#-nb_of_previousframe

                                    Trace = list(Reactivation_Strength[Frame_SWR_start-HalfSWR:Frame_SWR_start+HalfSWR])

                                    ActivityCa_swr=locals()[f'ActivityCa_swr{ctxSWR}']
                                    ActivityCa_swr.append(Trace)

                                    # Define if that SWR is coupled with a SPDL or not

                                    SWR_statut=[]
                                    startSpiList = list(pd.Series(SpipropTrunc["start time"]))
                                    endSpiList = list(pd.Series(SpipropTrunc["end time"]))
                                    ctxSpiList = list(pd.Series(SpipropTrunc["CTX"]))
                                    delaiSWRSpdl=None
                                    SWR_statut= 'UnCoupled'
                                    cUnCoupledSWR+=1 if cellAss_count==1 else 0
                                    ctxSpi=''
                                    IsTrue= False
                                    if len(startSpiList)>0:
                                        # if there is a spindle start before the SWR
                                        startClosest_Spdl_idx= next((i for i in range(len((startSpiList - startSwr)) - 1, -1, -1) if (startSpiList - startSwr)[i] < 0), -1)
                                        if startClosest_Spdl_idx != -1: 
                                            startClosest_Spi = startSpiList[startClosest_Spdl_idx]
                                            endClosest_Spi=endSpiList[startClosest_Spdl_idx]
                                            IsTrue = startSwr>startClosest_Spi and startSwr<endClosest_Spi #SWR inside the Spindle
                                            if IsTrue: 
                                                ctxSpi=ctxSpiList[startClosest_Spdl_idx]
                                                SWR_statut = 'PostCoupled'
                                                cPostCoupledSWR+=1 if cellAss_count==1 else 0
                                                cUnCoupledSWR-=1 if cellAss_count==1 else 0
                                                delaiSWRSpdl= startClosest_Spi - startSwr 
                                        
                                        # if there is a spindle after the SWR
                                        startClosest_Spdl_idx= next((i for i, x in enumerate((startSpiList - startSwr)) if x >= 0), -1) 
                                        if startClosest_Spdl_idx != -1: 
                                            startClosest_Spi = startSpiList[startClosest_Spdl_idx]
                                            endClosest_Spi=endSpiList[startClosest_Spdl_idx]
                                            distance = startClosest_Spi - startSwr 
                                            if distance<before: 
                                                ctxSpi=ctxSpiList[startClosest_Spdl_idx]
                                                if SWR_statut == 'PostCoupled':
                                                    SWR_statut = 'PrePostCoupled'
                                                    cPrePostCoupledSWR+=1 if cellAss_count==1 else 0
                                                    cPostCoupledSWR-=1 if cellAss_count==1 else 0
                                                    delaiSWRSpdl= distance
                                                else:
                                                    SWR_statut = 'PreCoupled'
                                                    cPreCoupledSWR+=1 if cellAss_count==1 else 0
                                                    cUnCoupledSWR-=1 if cellAss_count==1 else 0
                                                    delaiSWRSpdl=distance                                         

                                    ActivityCa_swrCp=locals()[f'ActivityCa_{SWR_statut}swr{ctxSWR}']
                                    ActivityCa_swrCp.append(Trace)
                                    
                                    # Fill the big summary table SWR_GlobalResults

                                    SWR_GlobalResults.loc[counter2, 'Mice'] = mice
                                    SWR_GlobalResults.loc[counter2, 'NeuronType'] = NeuronType

                                    SWR_GlobalResults.loc[counter2, 'Session'] = session
                                    SWR_GlobalResults.loc[counter2, 'Session_Date'] = session_date 
                                    SWR_GlobalResults.loc[counter2, 'Session_Time'] = session_time                    

                                    SWR_GlobalResults.loc[counter2, 'Assembly_ID'] = assembly_ID                                    
                                    SWR_GlobalResults.loc[counter2, 'Cells_in_Assembly'] = str(cells_in_assembly)                                    
                                    SWR_GlobalResults.loc[counter2, 'ExpeType'] =  session_type


                                    SWR_GlobalResults.loc[counter2, 'SWRStatut'] = SWR_statut
                                    SWR_GlobalResults.loc[counter2, 'SWRStartLocation'] = StartLocSWR
                                    SWR_GlobalResults.loc[counter2, 'GlobalSWR'] = diffSWR
                                    SWR_GlobalResults.loc[counter2, 'SpdlLoc'] = ctxSpi
                                    SWR_GlobalResults.loc[counter2, 'SWRNumber'] = Pswr
                                    SWR_GlobalResults.loc[counter2, 'SWRDuration'] = endSwr- startSwr
                                    SWR_GlobalResults.loc[counter2, 'SWR_inside_Spdl'] = IsTrue
                                    SWR_GlobalResults.loc[counter2, 'DistanceSWR_Spdl'] = delaiSWRSpdl 
                                    SWR_GlobalResults.loc[counter2, 'DistanceClosestSWR'] = closestSWR  

                                    durOsc=round(durationSpdl*target_rate)  
                                    
                                    SWR_GlobalResults.loc[counter2, 'AUC_1stQuarter'] = 2 * np.trapz(Trace[:round(durOsc/2)],np.arange(0,len(Trace[:round(durOsc/2)]),1)) # *2 cause to be per sec
                                    SWR_GlobalResults.loc[counter2, 'AUC_2ndQuarter'] = 2 * np.trapz(Trace[round(durOsc/2):durOsc],np.arange(0,len(Trace[round(durOsc/2):durOsc]),1))
                                    SWR_GlobalResults.loc[counter2, 'AUC_3rdQuarter'] = 2 * np.trapz(Trace[durOsc:durOsc+round(durOsc/2)],np.arange(0,len(Trace[durOsc:durOsc+round(durOsc/2)]),1))    
                                    SWR_GlobalResults.loc[counter2, 'AUC_4thQuarter'] = 2 * np.trapz(Trace[durOsc+round(durOsc/2):],np.arange(0,len(Trace[durOsc+round(durOsc/2):]),1))          
                                    SWR_GlobalResults.loc[counter2, 'AUC_1stHalf'] = np.trapz(Trace[:durOsc],np.arange(0,len(Trace[:durOsc]),1))          
                                    SWR_GlobalResults.loc[counter2, 'AUC_2ndHalf'] = np.trapz(Trace[durOsc:],np.arange(0,len(Trace[durOsc:]),1))          

                                    SWR_GlobalResults.loc[counter2, 'FR_1stQuarter'] = 2 * np.sum(Trace[:round(durOsc/2)]) # *2 cause to be in Hz
                                    SWR_GlobalResults.loc[counter2, 'FR_2ndQuarter'] = 2 * np.sum(Trace[round(durOsc/2):durOsc])
                                    SWR_GlobalResults.loc[counter2, 'FR_3rdQuarter'] = 2 * np.sum(Trace[durOsc:durOsc+round(durOsc/2)])          
                                    SWR_GlobalResults.loc[counter2, 'FR_4thQuarter'] = 2 * np.sum(Trace[durOsc+round(durOsc/2):])          
                                    SWR_GlobalResults.loc[counter2, 'FR_1stHalf'] = np.sum(Trace[:durOsc])          
                                    SWR_GlobalResults.loc[counter2, 'FR_2ndHalf'] = np.sum(Trace[durOsc:])      
 
                                    counter2+=1  
                            #else: 
                            #    print("/!\ SWR too close to the end of the previous one,", session, ", SWR n°", Pswr, ", Start SWR =",  round(startSwr/1000,1), "s") if unit==0 else None

                        ## Peristimulus Time Histogram 
                        start6 = time.time()                        
                        for ctx2 in CTX2: 
                            for coup in Coupling:                                                                                              
                                # All Ca traces for each spindles per Unique unit (according to cross-registration)
                                ActivityCa = locals()[f'ActivityCa_{coup}swr{ctx2}']
                                dict_All_ActivityCa = locals()[f'dict_All_ActivityCa_{coup}SWR{ctx2}_{drug}']
                                if len(ActivityCa)>0 :                                  
                                    if np.shape(np.array(ActivityCa))[1] == int(norm_freq*durationSWR*2):   #normalize traces to the same frequency rate    
                                        ActivityCa= np.reshape(np.array(ActivityCa), (-1, len(np.array(ActivityCa)))) if np.ndim(ActivityCa) == 1 else np.array(ActivityCa)    
                                        dict_All_ActivityCa[assembly_ID] = np.append(dict_All_ActivityCa[assembly_ID], np.array(ActivityCa), axis=0) if assembly_ID in dict_All_ActivityCa else np.array(ActivityCa)
                                    else:
                                        dataO = np.array(ActivityCa)
                                        data= np.repeat(dataO, 2, axis=0) if dataO.shape[0] == 1 else dataO
                                        x_mesh, y_mesh = np.meshgrid(np.arange(data.shape[1]), np.arange(data.shape[0]))
                                        x_new_mesh, y_new_mesh = np.meshgrid(np.linspace(0, data.shape[1] - 1, int(norm_freq*durationSWR*2)), np.linspace(0, data.shape[0] - 1, np.shape(data)[0]))
                                        resampled_dataO = griddata((x_mesh.flatten(), y_mesh.flatten()), data.flatten(), (x_new_mesh, y_new_mesh), method='linear')
                                        resampled_data= resampled_dataO[0,:] if dataO.shape[0] == 1 else resampled_dataO
                                        resampled_data= np.reshape(resampled_data, (-1, len(resampled_data))) if np.ndim(resampled_data) == 1 else resampled_data
                                        dict_All_ActivityCa[assembly_ID] = np.append(dict_All_ActivityCa[assembly_ID], np.array(resampled_data), axis=0) if assembly_ID in dict_All_ActivityCa else np.array(resampled_data)
                
                        print(f'... {len(SWRpropTrunc)} SWR processed in {time.time() - start5:.2f} seconds & PSTH done in {time.time() - start6:.2f} seconds for one cell assembly') if cellAss_count==1 else None

                print(f"... The {cellAss_count} cell assemblies of {session} analyzed in {time.time() - start2:.2f} seconds")

                sentence2=f"... {cPreCoupled+cPostCoupled+cPrePostCoupled+cUnCoupled}/{nb_spindle} spindles kept ({cPreCoupled} PreCoupled & {cPostCoupled} PostCoupled & {cPrePostCoupled} PrePostCoupled & {cUnCoupled} Uncoupled Spdl // {cGlobal} Global, {cLocalS1} LocalS1 & {cLocalPFC} LocalPFC) and {cPreCoupledSWR+cPostCoupledSWR+cPrePostCoupledSWR+cUnCoupledSWR}/{nb_swr} SWR kept ({cPreCoupledSWR} PreCoupled & {cPostCoupledSWR} PostCoupled & {cPrePostCoupledSWR} PrePostcoupled & {cUnCoupledSWR} Uncoupled SWR // {cGlobalSWR} Global, {cLocalCA1} LocalCA1 & {cLocalRSC} LocalRSC)"
                print(sentence2) 

    #######################################################################################
                            # Save Spindles & SWR analysis #
    #######################################################################################
    # Do average Calcium & Spike results for Spindles & SWR Peristimulus Time Histogram 
    start7 = time.time()

    Data=['Ca']
    for data in Data:   
        for ctx in CTX: 
            for coup in Coupling:
                for drug in drugs:      
                    dict_All_Activity=locals()[f'dict_All_Activity{data}_{coup}SPDL{ctx}_{drug}']
                    filenameOut = folder_to_save / f'Spdl_{data}PSTH_{coup}{ctx}{drug}_{mice}.pkl' #keep each responses of each cells for all rec Spdl
                    with open(filenameOut, 'wb') as pickle_file:
                        pickle.dump(dict_All_Activity, pickle_file)
        for ctx2 in CTX2: 
            for coup in Coupling:
                for drug in drugs:      
                    dict_All_Activity=locals()[f'dict_All_Activity{data}_{coup}SWR{ctx2}_{drug}']
                    filenameOut = folder_to_save / f'SWR_{data}PSTH_{coup}{ctx2}{drug}_{mice}.pkl' #keep each responses of each cells for all rec SWR
                    with open(filenameOut, 'wb') as pickle_file:
                        pickle.dump(dict_All_Activity, pickle_file)

    start8 = time.time()
    
    with open(folder_to_save / f'Spdl_Global_{mice}.pkl', 'wb') as pickle_file:
        pickle.dump(Spindles_GlobalResults, pickle_file)   

    with open(folder_to_save / f'SWR_Global_{mice}.pkl', 'wb') as pickle_file:
        pickle.dump(SWR_GlobalResults, pickle_file)

    print(f"Global matrix saved in {time.time() - start8:.2f} seconds")

    sentence3=f"{mice} data saved in {time.time() - start2:.2f} seconds"
    print(sentence3) 


Data=['Ca']
for data in Data:   
    for ctx in CTX: 
        for coup in Coupling:
            for drug in drugs:      
                dict_All_Activity=locals()[f'dict_All_Activity{data}_{coup}SPDL{ctx}_{drug}']
                filenameOut = folder_to_save / f'Spdl_{data}PSTH_{coup}{ctx}{drug}.pkl' #keep each responses of each cells for all rec Spdl
                with open(filenameOut, 'wb') as pickle_file:
                    pickle.dump(dict_All_Activity, pickle_file)
    for ctx2 in CTX2: 
        for coup in Coupling:
            for drug in drugs:                
                dict_All_Activity=locals()[f'dict_All_Activity{data}_{coup}SWR{ctx2}_{drug}']
                filenameOut = folder_to_save / f'SWR_{data}PSTH_{coup}{ctx2}{drug}.pkl' #keep each responses of each cells for all rec SWR
                with open(filenameOut, 'wb') as pickle_file:
                    pickle.dump(dict_All_Activity, pickle_file)

with open(folder_to_save / f'Spdl_Global.pkl', 'wb') as pickle_file:
        pickle.dump(Spindles_GlobalResults, pickle_file)   

with open(folder_to_save / f'SWR_Global.pkl', 'wb') as pickle_file:
    pickle.dump(SWR_GlobalResults, pickle_file)

sys.stdout = sys.__stdout__
logfile.close()