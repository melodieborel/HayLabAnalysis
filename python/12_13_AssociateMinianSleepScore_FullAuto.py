# # Associate Ca2+ signal with sleep stages for each session & subsessions using crossregistration

#######################################################################################
                            # Define Experiment type #
#######################################################################################

#DrugExperiment=1 if CGP Experiment // DrugExperiment=0 if Baseline Experiment
DrugExperiment=1

#Sleep scoring from '_AB' '_AH' or initial ''
suffix='_AB' 
AnalysisID='_ALL_zscored' 

saveexcel=0

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
from scipy.stats import zscore
from itertools import groupby
from IPython.display import display

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
FolderNameSave = FolderNameSave.replace(" ", "_").replace(".", "_").replace(":", "_")

destination_folder= f"//10.69.168.1/crnldata/waking/audrey_hay/L1imaging/AnalysedMarch2023/Gaelle/CGP/AB_Analysis/VigSt_{FolderNameSave}{suffix}{AnalysisID}" if DrugExperiment else f"//10.69.168.1/crnldata/waking/audrey_hay/L1imaging/AnalysedMarch2023/Gaelle/Baseline_recording_ABmodified/AB_Analysis/VigSt_{FolderNameSave}{suffix}{AnalysisID}"
os.makedirs(destination_folder)
folder_to_save=Path(destination_folder)

# Copy the script file to the destination folder
source_script = "C:/Users/Manip2/SCRIPTS/CodePythonAudrey/CodePythonAurelie/HayLabAnalysis/python/12_13_AssociateMinianSleepScore_FullAuto.py"
destination_file_path = f"{destination_folder}/12_13_AssociateMinianSleepScore_FullAuto.txt"
shutil.copy(source_script, destination_file_path)

for mice in MiceList:
    # Load sleep score and Ca2+ time series numpy arrays
    dpath0 = "//10.69.168.1/crnldata/waking/audrey_hay/L1imaging/AnalysedMarch2023/Gaelle/CGP/" if DrugExperiment else "//10.69.168.1/crnldata/waking/audrey_hay/L1imaging/AnalysedMarch2023/Gaelle/Baseline_recording_ABmodified/"
    dpath=dpath0 + mice
    print(f"####################################################################################")
    print(f"################################### {mice} ####################################")
    print(f"####################################################################################")
    print(f"Path to the folder : {dpath}")
    folder_base = Path(dpath)

    mfile = open(folder_base / f'mappingsAB_ALL.pkl', 'rb')
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
    dict_Path={}

    sessions, sessions_path = find_session_folders(folder_base)
    nb_sessions=len(sessions)

    for sess,session in enumerate(sessions):  
        session_path=Path(sessions_path[sess])
        folder_mini = session_path / f'V4_Miniscope'
        nb_subsessions = sum(1 for p in folder_mini.iterdir() if p.is_dir() and p.name.startswith("session"))
        ScoringFile = session_path/ f'OpenEphys/ScoredSleep{suffix}.npy'
        StampsFile = session_path/ f'SynchroFile.xlsx'
        StampsMiniscopeFile = folder_mini / f'timeStamps.csv'

        if nb_subsessions!=0:
            for x in range(1, nb_subsessions+1):            
                subsession= session + str(x)
                subsessions.append(subsession)    
                minian_ds = open_minian(folder_mini / subsession / f'minian')      # OR minianAB
                dict_Path[subsession] = session_path
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
            dict_Path[session] = session_path
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
    if mice == 'Purple' and DrugExperiment==0:
        index = B.columns
        B.columns = index.str.replace('part', 'session2')
    #######################################################################################
      # Distribute Ca2+ intensity & spikes to vigilance states for each sessions/subsessions #
    #######################################################################################
    
    data = {}
    counter=0
    VigilanceState_GlobalResults= pd.DataFrame(data, columns=['Mice','Session', 'Session_Time', 'Unique_Unit','UnitNumber','UnitValue', 'Drug', 'Substate','SubstateNumber','DurationSubstate', 'CalciumActivity', 'Avg_CalciumActivity', 'AUC_calcium','Avg_AUC_calcium', 'DeconvSpikeMeanActivity', 'Avg_DeconvSpikeActivity', 'SpikeActivityHz', 'Avg_SpikeActivityHz'])

    previousEndTime=0
    InitialStartTime=0

    CaCorrWakeMatrixBaseline=[]
    CaCorrNREMMatrixBaseline=[]
    CaCorrN2MatrixBaseline=[]
    CaCorrREMMatrixBaseline=[]

    SpCorrWakeMatrixBaseline=[]
    SpCorrNREMMatrixBaseline=[]
    SpCorrN2MatrixBaseline=[]
    SpCorrREMMatrixBaseline=[]

    RawCaTracesWake_Baseline=[]
    RawCaTracesNREM_Baseline=[]
    RawCaTracesN2_Baseline=[]
    RawCaTracesREM_Baseline=[]

    RawSpTracesWake_Baseline=[]
    RawSpTracesNREM_Baseline=[]
    RawSpTracesN2_Baseline=[]
    RawSpTracesREM_Baseline=[]

    CaCorrWakeMatrixCGP=[]
    CaCorrNREMMatrixCGP=[]
    CaCorrN2MatrixCGP=[]
    CaCorrREMMatrixCGP=[]

    SpCorrWakeMatrixCGP=[]
    SpCorrNREMMatrixCGP=[]
    SpCorrN2MatrixCGP=[]
    SpCorrREMMatrixCGP=[]

    RawCaTracesWake_CGP=[]
    RawCaTracesNREM_CGP=[]
    RawCaTracesN2_CGP=[]
    RawCaTracesREM_CGP=[]

    RawSpTracesWake_CGP=[]
    RawSpTracesNREM_CGP=[]
    RawSpTracesN2_CGP=[]
    RawSpTracesREM_CGP=[]

    Drugs=['Baseline', 'CGP'] if DrugExperiment else ['Baseline']
    if saveexcel: 
        filenameOut = folder_to_save / f'VigSt_CaCorr_{mice}.xlsx'
        excel_writerCa = pd.ExcelWriter(filenameOut)        
        filenameOut = folder_to_save / f'VigSt_SpCorr_{mice}.xlsx'
        excel_writerSp = pd.ExcelWriter(filenameOut)
        filenameOut = folder_to_save / f'VigSt_RawCaTraces_{mice}.xlsx'
        excel_writerRawCa = pd.ExcelWriter(filenameOut)        
        filenameOut = folder_to_save / f'VigSt_RawSpTraces_{mice}.xlsx'
        excel_writerRawSp = pd.ExcelWriter(filenameOut)

    for session in list(dict_Stamps.keys()):

        drug=os.path.basename(os.path.dirname(dict_Path[session])) if DrugExperiment else 'Baseline'

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

        # Normalize traces
        
        Carray=zscore(Carray, axis=0)
        #min=np.min(Carray,axis=0) 
        #Carray=Carray-min
        Sarray=zscore(Sarray, axis=0)
        #min=np.min(Sarray,axis=0) 
        #Sarray=Sarray-min
        
        # Replace dropped frame in calcium and spike traces with the previous value

        for droppedframe in droppedframes_inrec: 
            row_to_repeat = Carray[droppedframe]  
            Carray = np.vstack((Carray[:droppedframe], row_to_repeat, Carray[droppedframe:]))
            row_to_repeat = Sarray[droppedframe]  
            Sarray = np.vstack((Sarray[:droppedframe], row_to_repeat, Sarray[droppedframe:]))

        # Upscale scoring to miniscope frequency

        scale_factor=minian_freq/0.2  #cause scoring was done in 5 seconds bin, ie 0.2 Hz
        SleepScoredTS=dict_Scoring[session]
        SleepScoredTS_upscaled = np.repeat(SleepScoredTS, scale_factor, axis=0)
        StartTime_frame=int(StartTime*minian_freq)
        SleepScoredTS_upscaled_ministart=SleepScoredTS_upscaled[StartTime_frame:StartTime_frame+upd_rec_dur]

        # Remove N2 stage
        SleepScoredTS_upscaled_ministart[SleepScoredTS_upscaled_ministart == 0.5] = 0
        mapp = {1.5: 'Wake', 0: 'NREM',  1: 'REM'}
        # mapp = {1.5: 'Wake', 0.5: 'N2', 0: 'NREM',  1: 'REM'}

        # Determine each substate identity and duration
        array=SleepScoredTS_upscaled_ministart
        substates_duration = [len(list(group)) for key, group in groupby(array)]
        substates_identity = [key for key, _ in groupby(array)]
        substates_end = np.array(substates_duration).cumsum()        
        substates_start =np.append([0],substates_end[:-1]) #substates_start =np.append([1],substates_end+1) create 1 second gap
        substates_identity = [mapp[num] for num in substates_identity]
        substates = pd.DataFrame(list(zip(substates_identity, substates_duration, substates_start, substates_end)), columns=['Identity', 'Duration', 'Start','End'])

        for m in mapp:

            Bool = (SleepScoredTS_upscaled_ministart == m)
            Carray_VigSpe = Carray.copy()
            Carray_VigSpe = Carray_VigSpe[0:np.shape(SleepScoredTS_upscaled_ministart)[0],:] # if Calcium imaging longer than LFP rec
            Carray_VigSpe = Carray_VigSpe[Bool, :]

            RawCaTracesVigStateMatrixName=f'RawCaTraces{mapp[m]}_{drug}'
            RawCaTracesVigStateMatrix = locals()[RawCaTracesVigStateMatrixName]
            RawCaTraces=[]
            RawCaTraces=pd.DataFrame(Carray_VigSpe, columns=[f"{mice}{str(i).replace('[','').replace(']','')}" for i in kept_uniq_unit_List])
            unique_columns = RawCaTraces.columns[RawCaTraces.columns.to_series().duplicated()] # remove units that has an empty unique index '[]'
            RawCaTraces = RawCaTraces.drop(columns=unique_columns)
            RawCaTracesVigStateMatrix.append(RawCaTraces)

            Sarray_VigSpe = Sarray.copy()
            Sarray_VigSpe = Sarray_VigSpe[0:np.shape(SleepScoredTS_upscaled_ministart)[0],:] # if Calcium imaging longer than LFP rec
            Sarray_VigSpe = Sarray_VigSpe[Bool, :]         

            RawSpTracesVigStateMatrixName=f'RawSpTraces{mapp[m]}_{drug}'
            RawSpTracesVigStateMatrix = locals()[RawSpTracesVigStateMatrixName]
            RawSpTraces=[]
            RawSpTraces=pd.DataFrame(Sarray_VigSpe, columns=[f"{mice}{str(i).replace('[','').replace(']','')}" for i in kept_uniq_unit_List])
            unique_columns = RawSpTraces.columns[RawSpTraces.columns.to_series().duplicated()] # remove units that has an empty unique index '[]'
            RawSpTraces = RawSpTraces.drop(columns=unique_columns)
            RawSpTracesVigStateMatrix.append(RawSpTraces) 
            
            CaCorrVigStateMatrixName=f'CaCorr{mapp[m]}Matrix{drug}'
            CaCorrVigStateMatrix = locals()[CaCorrVigStateMatrixName]
            CaCorrMatrix=[]
            CaCorrMatrix = pd.DataFrame(columns=[f'{mice}{str(i)}' for i in range(len(B))], index=[i for i in range(len(B))])
            
            SpCorrVigStateMatrixName=f'SpCorr{mapp[m]}Matrix{drug}'
            SpCorrVigStateMatrix = locals()[SpCorrVigStateMatrixName]
            SpCorrMatrix=[]
            SpCorrMatrix = pd.DataFrame(columns=[f'{mice}{str(i)}' for i in range(len(B))], index=[i for i in range(len(B))])

            for unit in range(nb_unit): 

                Carray_unit =Carray_VigSpe[:,unit]
                Darray_unit =Sarray_VigSpe[:,unit] # on deconv spike not on spike rate                
                #peaks, _ = find_peaks(Darray_unit)
                #Sarray_unit=np.zeros(len(Darray_unit))
                #Sarray_unit[peaks]=1     

                otherunit_range = [x for x in range(nb_unit) if x != unit]

                for unit2 in range(nb_unit):

                    Carray_unit2 =Carray_VigSpe[:,unit2]
                    Darray_unit2 =Sarray_VigSpe[:,unit2]                    
                    #peaks2, _ = find_peaks(Darray_unit2)
                    #Sarray_unit2=np.zeros(len(Darray_unit2))
                    #Sarray_unit2[peaks2]=1

                    indexMapp = str(np.where(B[session] == C_upd_unit_id[unit])[0]).replace('[','').replace(']','')
                    indexMapp2 = np.where(B[session] == C_upd_unit_id[unit2])[0]

                    if any(indexMapp) and len(indexMapp2)>0:      
                                        
                        corr_matrix = np.corrcoef(Carray_unit, Carray_unit2)
                        CaCorrMatrix[f'{mice}{indexMapp}'][indexMapp2]={1: 0.99999, -1: -0.99999}.get(np.round(corr_matrix[0, 1],5), np.round(corr_matrix[0, 1],5))
                        
                        corr_matrix = np.corrcoef(Darray_unit, Darray_unit2)
                        SpCorrMatrix[f'{mice}{indexMapp}'][indexMapp2]={1: 0.99999, -1: -0.99999}.get(np.round(corr_matrix[0, 1],5), np.round(corr_matrix[0, 1],5))
                        
            CaCorrVigStateMatrix.append(CaCorrMatrix)
            SpCorrVigStateMatrix.append(SpCorrMatrix) 

        for unit in range(nb_unit): 

            Carray_unit =Carray[:,unit]
            Darray_unit =Sarray[:,unit]
            peaks, _ = find_peaks(Darray_unit)
            Sarray_unit=np.zeros(len(Darray_unit))
            Sarray_unit[peaks]=1     

            for index in range(len(substates)):
                
                ca_input_sub=Carray_unit[substates.Start[index]:substates.End[index]]
                ds_input_sub=Darray_unit[substates.Start[index]:substates.End[index]]
                sp_input_sub=Sarray_unit[substates.Start[index]:substates.End[index]]

                VigilanceState_GlobalResults.loc[counter, 'Mice'] = mice
                VigilanceState_GlobalResults.loc[counter, 'Session'] = session
                VigilanceState_GlobalResults.loc[counter, 'Session_Time'] = None 
                
                indexMapp = np.where(B[session] == C_upd_unit_id[unit])[0]
                VigilanceState_GlobalResults.loc[counter, 'Unique_Unit'] = indexMapp if len(indexMapp)>0 else None
                VigilanceState_GlobalResults.loc[counter, 'UnitNumber'] = unit
                VigilanceState_GlobalResults.loc[counter, 'UnitValue'] = C_upd_unit_id[unit]

                VigilanceState_GlobalResults.loc[counter, 'Drug'] =  os.path.basename(os.path.dirname(dict_Path[session])) if DrugExperiment else 'Baseline'

                VigilanceState_GlobalResults.loc[counter, 'Substate'] = substates.Identity[index]
                VigilanceState_GlobalResults.loc[counter, 'SubstateNumber'] = substates.index[index]
                VigilanceState_GlobalResults.loc[counter, 'DurationSubstate'] = (substates.Duration[index]/minian_freq)

                VigilanceState_GlobalResults.loc[counter, 'CalciumActivity'] = ca_input_sub.mean()
                VigilanceState_GlobalResults.loc[counter, 'Avg_CalciumActivity'] = Carray_unit.mean()

                VigilanceState_GlobalResults.loc[counter, 'AUC_calcium'] = np.trapz(ca_input_sub,np.arange(0,len(ca_input_sub),1))
                VigilanceState_GlobalResults.loc[counter, 'Avg_AUC_calcium'] = np.trapz(Carray_unit,np.arange(0,len(Carray_unit),1))

                VigilanceState_GlobalResults.loc[counter, 'DeconvSpikeMeanActivity'] = ds_input_sub.mean()
                VigilanceState_GlobalResults.loc[counter, 'Avg_DeconvSpikeActivity'] = Darray_unit.mean()

                VigilanceState_GlobalResults.loc[counter, 'SpikeActivityHz'] = sp_input_sub.sum()/(len(sp_input_sub)/minian_freq)
                VigilanceState_GlobalResults.loc[counter, 'Avg_SpikeActivityHz'] = Sarray_unit.sum()/(len(Sarray_unit)/minian_freq)

                otherunit_range = [x for x in range(nb_unit) if x != unit]

                CaCorrCoeff=[]
                SpCorrCoeff=[]
                Z_CaCorrCoeff=[]
                Z_SpCorrCoeff=[]

                for unit2 in otherunit_range:

                    Carray_unit2 =Carray[:,unit2]
                    Darray_unit2 =Sarray[:,unit2]                    
                    #peaks2, _ = find_peaks(Sarray_unit2)
                    #Sarray_unit2=np.zeros(len(Sarray_unit2))
                    #Sarray_unit2[peaks2]=1
                    ca_input_sub2=Carray_unit2[substates.Start[index]:substates.End[index]]
                    ds_input_sub2=Darray_unit2[substates.Start[index]:substates.End[index]]

                    corr_matrix = np.corrcoef(ca_input_sub, ca_input_sub2)
                    CaCorrCoeff_unit = np.round(corr_matrix[0, 1],5)
                    CaCorrCoeff_unit = {1: 0.99999, -1: -0.99999}.get(CaCorrCoeff_unit, CaCorrCoeff_unit)
                    CaCorrCoeff.append(CaCorrCoeff_unit)
                    Z_CaCorrCoeff_unit=np.arctanh(CaCorrCoeff_unit)
                    Z_CaCorrCoeff.append(Z_CaCorrCoeff_unit)

                    corr_matrix = np.corrcoef(ds_input_sub, ds_input_sub2)
                    SpCorrCoeff_unit =  np.round(corr_matrix[0, 1],5)
                    SpCorrCoeff_unit = {1: 0.99999, -1: -0.99999}.get(SpCorrCoeff_unit, SpCorrCoeff_unit)
                    SpCorrCoeff.append(SpCorrCoeff_unit)
                    Z_SpCorrCoeff_unit=np.arctanh(SpCorrCoeff_unit)
                    Z_SpCorrCoeff.append(Z_SpCorrCoeff_unit)
                
                VigilanceState_GlobalResults.loc[counter, 'CaCorrCoeff'] = np.nanmean(CaCorrCoeff)
                VigilanceState_GlobalResults.loc[counter, 'SpCorrCoeff'] = np.nanmean(SpCorrCoeff)
                
                VigilanceState_GlobalResults.loc[counter, 'Rsquared_CaCorrCoeff'] = np.nanmean([x ** 2 for x in CaCorrCoeff])
                VigilanceState_GlobalResults.loc[counter, 'Rsquared_SpCorrCoeff'] = np.nanmean([x ** 2 for x in SpCorrCoeff])

                VigilanceState_GlobalResults.loc[counter, 'Z_CaCorrCoeff'] = np.nanmean(Z_CaCorrCoeff)
                VigilanceState_GlobalResults.loc[counter, 'Z_SpCorrCoeff'] = np.nanmean(Z_SpCorrCoeff)
                
                VigilanceState_GlobalResults.loc[counter, 'Z_Rsquared_CaCorrCoeff'] = np.nanmean([x ** 2 for x in Z_CaCorrCoeff])
                VigilanceState_GlobalResults.loc[counter, 'Z_Rsquared_SpCorrCoeff'] = np.nanmean([x ** 2 for x in Z_SpCorrCoeff])   

                Carray_Population =np.mean(Carray[:,otherunit_range], axis=1)
                ca_input_sub2=Carray_Population[substates.Start[index]:substates.End[index]]
                corr_matrix = np.corrcoef(ca_input_sub, ca_input_sub2)
                VigilanceState_GlobalResults.loc[counter, 'CaPopCoupling'] = np.round(corr_matrix[0, 1],5)
                corCorrected = {1: 0.99999, -1: -0.99999}.get(np.round(corr_matrix[0, 1],5), np.round(corr_matrix[0, 1],5))
                VigilanceState_GlobalResults.loc[counter, 'Z_CaPopCoupling'] = np.arctanh(corCorrected)

                Sarray_Population =np.mean(Sarray[:,otherunit_range], axis=1)
                ds_input_sub2=Sarray_Population[substates.Start[index]:substates.End[index]]
                corr_matrix = np.corrcoef(ds_input_sub, ds_input_sub2)
                VigilanceState_GlobalResults.loc[counter, 'SpPopCoupling'] = np.round(corr_matrix[0, 1],5)
                corCorrected = {1: 0.99999, -1: -0.99999}.get(np.round(corr_matrix[0, 1],5), np.round(corr_matrix[0, 1],5))
                VigilanceState_GlobalResults.loc[counter, 'Z_SpPopCoupling'] = np.arctanh(corCorrected)

                counter+=1

    dataCaCorr={}
    dataSpCorr={}
    for Drug in Drugs:
        for m in mapp:
            CaCorrVigStateMatrixName=f'CaCorr{mapp[m]}Matrix{Drug}'
            CaCorrVigStateMatrix = locals()[CaCorrVigStateMatrixName]
            if len(CaCorrVigStateMatrix)>0: # cause sometimes no Baseline conditions in CGP experiments
                combined_df = pd.concat(CaCorrVigStateMatrix, ignore_index=False)
                IterationNb=combined_df.groupby(combined_df.index).count()
                combined_df = combined_df.groupby(combined_df.index).sum() #mean
                combined_df.index= [mice + str(idx) for idx in combined_df.index]
                combined_df = combined_df[~(combined_df.fillna(0) == 0).all(axis=1)]
                combined_df = combined_df.loc[:, ~(combined_df.fillna(0) == 0).all(axis=0)]
                IterationNb = IterationNb[~(IterationNb.fillna(0) == 0).all(axis=1)]
                IterationNb = IterationNb.loc[:, ~(IterationNb.fillna(0) == 0).all(axis=0)]
                IterationNb.index=combined_df.index
                if saveexcel: combined_df.to_excel(excel_writerCa, sheet_name=f'{Drug}_{mapp[m]}', index=True, header=True) 
                dataCaCorr[f'{Drug}_{mapp[m]}']=combined_df
                
                combined_df = pd.concat(CaCorrVigStateMatrix, ignore_index=False)
                combined_df = combined_df.applymap(lambda x: np.arctanh(x) if not pd.isna(x) else np.nan)
                combined_df = combined_df.groupby(combined_df.index).sum() #mean
                combined_df.index = [mice + str(idx) for idx in combined_df.index]
                combined_df = combined_df[~(combined_df.fillna(0) == 0).all(axis=1)]
                combined_df = combined_df.loc[:, ~(combined_df.fillna(0) == 0).all(axis=0)]
                if saveexcel: combined_df.to_excel(excel_writerCa, sheet_name=f'Z_{Drug}_{mapp[m]}', index=True, header=True)   
                dataCaCorr[f'Z_{Drug}_{mapp[m]}']=combined_df
                if saveexcel: IterationNb.to_excel(excel_writerCa, sheet_name=f'{Drug}_{mapp[m]}_IterationNb', index=True, header=True) 
                dataCaCorr[f'{Drug}_{mapp[m]}_IterationNb']=IterationNb

                SpCorrVigStateMatrixName=f'SpCorr{mapp[m]}Matrix{Drug}'
                SpCorrVigStateMatrix = locals()[SpCorrVigStateMatrixName]
                combined_df = pd.concat(SpCorrVigStateMatrix, ignore_index=False)
                IterationNb=combined_df.groupby(combined_df.index).count()
                combined_df = combined_df.groupby(combined_df.index).sum()
                combined_df.index = [mice + str(idx) for idx in combined_df.index]
                combined_df = combined_df[~(combined_df.fillna(0) == 0).all(axis=1)]
                combined_df = combined_df.loc[:, ~(combined_df.fillna(0) == 0).all(axis=0)]
                IterationNb = IterationNb[~(IterationNb.fillna(0) == 0).all(axis=1)]
                IterationNb = IterationNb.loc[:, ~(IterationNb.fillna(0) == 0).all(axis=0)]
                IterationNb.index=combined_df.index
                if saveexcel: combined_df.to_excel(excel_writerSp, sheet_name=f'{Drug}_{mapp[m]}', index=True, header=True) 
                dataSpCorr[f'{Drug}_{mapp[m]}']=combined_df
                
                combined_df = pd.concat(SpCorrVigStateMatrix, ignore_index=False)
                combined_df = combined_df.applymap(lambda x: np.arctanh(x) if not pd.isna(x) else np.nan)
                combined_df = combined_df.groupby(combined_df.index).sum()
                combined_df.index = [mice + str(idx) for idx in combined_df.index]
                combined_df = combined_df[~(combined_df.fillna(0) == 0).all(axis=1)]
                combined_df = combined_df.loc[:, ~(combined_df.fillna(0) == 0).all(axis=0)]
                if saveexcel: combined_df.to_excel(excel_writerSp, sheet_name=f'Z_{Drug}_{mapp[m]}', index=True, header=True) 
                dataSpCorr[f'Z_{Drug}_{mapp[m]}']=combined_df                
                if saveexcel: IterationNb.to_excel(excel_writerSp, sheet_name=f'{Drug}_{mapp[m]}_IterationNb', index=True, header=True) 
                dataSpCorr[f'{Drug}_{mapp[m]}_IterationNb']=IterationNb
            
            if saveexcel:
                RawCaTracesVigStateMatrixName=f'RawCaTraces{mapp[m]}_{Drug}'
                RawCaTracesVigStateMatrix= locals()[RawCaTracesVigStateMatrixName]
                if len(RawCaTracesVigStateMatrix)>0: # cause sometimes no Baseline conditions in CGP experiments
                    combined_df = pd.concat(RawCaTracesVigStateMatrix, ignore_index=False)
                    combined_df.index = [mice + str(idx) for idx in combined_df.index]
                    combined_df = combined_df.dropna(axis=0, how='all')
                    combined_df = combined_df.dropna(axis=1, how='all')
                    combined_df.to_excel(excel_writerRawCa, sheet_name=f'{Drug}_{mapp[m]}', index=False, header=True)   

                    RawSpTracesVigStateMatrixName=f'RawSpTraces{mapp[m]}_{Drug}'
                    RawSpTracesVigStateMatrix= locals()[RawSpTracesVigStateMatrixName]
                    combined_df = pd.concat(RawSpTracesVigStateMatrix, ignore_index=False)
                    combined_df.index = [mice + str(idx) for idx in combined_df.index]
                    combined_df = combined_df.dropna(axis=0, how='all')
                    combined_df = combined_df.dropna(axis=1, how='all')
                    combined_df.to_excel(excel_writerRawSp, sheet_name=f'{Drug}_{mapp[m]}', index=False, header=True)   

    if saveexcel: 
        excel_writerCa.close() 
        excel_writerSp.close() 
        excel_writerRawCa.close()
        excel_writerRawSp.close()
        filenameOut = folder_to_save / f'VigSt_Global_{mice}.xlsx'
        writer = pd.ExcelWriter(filenameOut)
        VigilanceState_GlobalResults.to_excel(writer)
        writer.close()

    filenameOut = folder_to_save / f'VigSt_CaCorr_{mice}.pkl'
    with open(filenameOut, 'wb') as pickle_file:
        pickle.dump(dataCaCorr, pickle_file)
    filenameOut = folder_to_save / f'VigSt_SpCorr_{mice}.pkl'
    with open(filenameOut, 'wb') as pickle_file:
        pickle.dump(dataSpCorr, pickle_file)
    filenameOut = folder_to_save / f'VigSt_Global_{mice}.pkl'
    with open(filenameOut, 'wb') as pickle_file:
        pickle.dump(VigilanceState_GlobalResults, pickle_file)

    

