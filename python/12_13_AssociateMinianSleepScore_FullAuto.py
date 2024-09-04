# # Associate Ca2+ signal with sleep stages for each session & subsessions using crossregistration

#######################################################################################
                            # Define Experiment type #
#######################################################################################

#DrugExperiment=1 if CGP Experiment // DrugExperiment=0 if Baseline Experiment

DrugExperimentList=[0,1]

#Sleep scoring from '_AB' '_AH' or initial ''
suffix='_AB' 
AnalysisID='_wRealTS_N2' 

saveexcel=0

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

def filterLFP(lfp, f_lowcut, f_hicut):
    range=int(f_hicut-f_lowcut)
    fs = 1000
    nyq = 0.5 * fs
    N = 3 # Filtre order
    Wn = [f_lowcut/nyq,f_hicut/nyq]  # Nyquist frequency fraction
    b, a = signal.butter(N, Wn, 'band')
    filt= signal.filtfilt(b, a, lfp)
    # Parameter and computation of CWT
    w = 100. #window size
    freq = np.linspace(f_lowcut, f_hicut, range)
    widths = w*fs / (2*freq*np.pi)
    CWT = signal.cwt(filt, signal.morlet2, widths, w=w)
    # Projection calculation
    absCWT = np.absolute(CWT)
    zabsCWT = stats.zscore(absCWT, axis=None)
    proj_CWT = np.max(zabsCWT, axis = 0)/range
    return proj_CWT

def butter_bandpass(lowcut, highcut, fs, order=2):
    nyq = 0.5 * fs  # Nyquist Frequency
    low = lowcut / nyq
    high = highcut / nyq
    b, a = butter(order, [low, high], btype='band')
    return b, a

def bandpass_filter(data, lowcut, highcut, fs, order=2):
    b, a = butter_bandpass(lowcut, highcut, fs, order=order)
    y = filtfilt(b, a, data)
    return y

def moving_average(data, window_size):
    return np.convolve(data, np.ones(window_size)/window_size, mode='valid')

def get_filt_env(lfp_signal, lowcut, highcut, fs):
    filtered_signal = bandpass_filter(lfp_signal, lowcut, highcut, fs)
    analytic_signal = hilbert(filtered_signal)
    envelope = np.abs(analytic_signal)
    window_size = fs*5 #5 sec
    smoothed_envelope = moving_average(envelope, window_size)
    return smoothed_envelope

#######################################################################################
                # Load sleep score and Ca2+ time series numpy arrays #
#######################################################################################

for DrugExperiment in DrugExperimentList: 
    
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

    Channels = '//10.69.168.1/crnldata/waking/audrey_hay/L1imaging/AnalysedMarch2023/LFPChannels_perMice.xlsx' 
    allchannels = pd.read_excel(Channels)

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
        dict_LFP={}

        sessions, sessions_path = find_session_folders(folder_base)
        nb_sessions=len(sessions)

        for sess,session in enumerate(sessions):  
            session_path=Path(sessions_path[sess])
            folder_mini = session_path / f'V4_Miniscope'
            nb_subsessions = sum(1 for p in folder_mini.iterdir() if p.is_dir() and p.name.startswith("session"))
            ScoringFile = session_path/ f'OpenEphys/ScoredSleep{suffix}.npy'
            LFPFile = session_path/ f'OpenEphys/RawDataChannelExtractedDS.npy'
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
                    dict_LFP[subsession] = np.load(LFPFile, mmap_mode= 'r')
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
                dict_LFP[session]  = np.load(LFPFile, mmap_mode= 'r')
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
        VigilanceState_GlobalResults= pd.DataFrame(data, columns=['Mice','Session', 'Session_Time', 'Unique_Unit','UnitNumber','UnitValue', 
                                                                'Drug', 'Substate','SubstateNumber','DurationSubstate', 'CalciumActivity', 
                                                                'Avg_CalciumActivity', 'AUC_calcium','Avg_AUC_calcium', 'DeconvSpikeMeanActivity', 
                                                                'Avg_DeconvSpikeActivity', 'SpikeActivityHz', 'Avg_SpikeActivityHz', 'TotCaPopCoupling', 
                                                                'TotZ_CaPopCoupling', 'TotSpPopCoupling', 'TotZ_SpPopCoupling', 'SigmaPFC_corr',
                                                                'SigmaS1_corr', 'ThetaCA1_corr', 'DeltaS1_corr', 'DeltaPFC_corr','SOS1_corr', 
                                                                'SOPFC_corr','BetaPFC_corr', 'BetaS1_corr', 'Z_SigmaPFC_corr', 'Z_SigmaS1_corr', 
                                                                'Z_ThetaCA1_corr', 'Z_DeltaS1_corr','Z_DeltaPFC_corr','Z_SOS1_corr','Z_SOPFC_corr',
                                                                'Z_BetaPFC_corr', 'Z_BetaS1_corr'])
        
        previousEndTime=0
        InitialStartTime=0

        TotCaCorr=[]
        TotSpCorr=[]

        StatesCaCorrWakeMatrixBaseline= pd.DataFrame(columns=[f'{mice}{str(i)}' for i in range(len(B))], index=[i for i in range(len(B))])
        StatesCaCorrNREMMatrixBaseline=pd.DataFrame(columns=[f'{mice}{str(i)}' for i in range(len(B))], index=[i for i in range(len(B))])
        StatesCaCorrN2MatrixBaseline=pd.DataFrame(columns=[f'{mice}{str(i)}' for i in range(len(B))], index=[i for i in range(len(B))])
        StatesCaCorrREMMatrixBaseline=pd.DataFrame(columns=[f'{mice}{str(i)}' for i in range(len(B))], index=[i for i in range(len(B))])
        StatesCaCorrWakeMatrixCGP=pd.DataFrame(columns=[f'{mice}{str(i)}' for i in range(len(B))], index=[i for i in range(len(B))])
        StatesCaCorrNREMMatrixCGP=pd.DataFrame(columns=[f'{mice}{str(i)}' for i in range(len(B))], index=[i for i in range(len(B))])
        StatesCaCorrN2MatrixCGP=pd.DataFrame(columns=[f'{mice}{str(i)}' for i in range(len(B))], index=[i for i in range(len(B))])
        StatesCaCorrREMMatrixCGP=pd.DataFrame(columns=[f'{mice}{str(i)}' for i in range(len(B))], index=[i for i in range(len(B))])

        ITStatesCaCorrWakeMatrixBaseline=pd.DataFrame(columns=[f'{mice}{str(i)}' for i in range(len(B))], index=[i for i in range(len(B))])
        ITStatesCaCorrNREMMatrixBaseline=pd.DataFrame(columns=[f'{mice}{str(i)}' for i in range(len(B))], index=[i for i in range(len(B))])
        ITStatesCaCorrREMMatrixBaseline=pd.DataFrame(columns=[f'{mice}{str(i)}' for i in range(len(B))], index=[i for i in range(len(B))])
        ITStatesCaCorrN2MatrixBaseline=pd.DataFrame(columns=[f'{mice}{str(i)}' for i in range(len(B))], index=[i for i in range(len(B))])
        ITStatesCaCorrWakeMatrixCGP=pd.DataFrame(columns=[f'{mice}{str(i)}' for i in range(len(B))], index=[i for i in range(len(B))])
        ITStatesCaCorrNREMMatrixCGP=pd.DataFrame(columns=[f'{mice}{str(i)}' for i in range(len(B))], index=[i for i in range(len(B))])
        ITStatesCaCorrN2MatrixCGP=pd.DataFrame(columns=[f'{mice}{str(i)}' for i in range(len(B))], index=[i for i in range(len(B))])
        ITStatesCaCorrREMMatrixCGP=pd.DataFrame(columns=[f'{mice}{str(i)}' for i in range(len(B))], index=[i for i in range(len(B))])

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

        TotCaCorrBaseline=[]
        TotCaCorrCGP=[]
        TotSpCorrBaseline=[]
        TotSpCorrCGP=[]

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
            tsmini=dict_StampsMiniscope[session]['Time Stamp (ms)']
            minian_freq=round(1/np.mean(np.diff(np.array(tsmini)/1000)))

            freqLFP=1000
            

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

            CalciumSub = pd.DataFrame(C, index=C['unit_id'])
            SpikeSub = pd.DataFrame(S, index=S['unit_id'])

            unit_to_drop=dict_TodropFile[session]    
            for u in unit_to_drop: 
                CalciumSub=CalciumSub.drop(index=u) if u in CalciumSub.index else CalciumSub #need to know why
                SpikeSub=SpikeSub.drop(index=u) if u in SpikeSub.index else SpikeSub

            indexMappList=B[session]
            kept_uniq_unit_List=[]
            for unit in CalciumSub.index:
                indexMapp = np.where(indexMappList == unit)[0]
                kept_uniq_unit_List.append(str(indexMapp))

                
            # Realigned to the traces to the recorded timestamps 

            timestamps =  np.array(tsmini[firstframe:firstframe+len(CalciumSub.T)])/freqLFP
            x_values = CalciumSub  # Each row is a feature, each column corresponds to a timestamp
            sample_rate = minian_freq  # Hz
            new_timestamps= np.arange(timestamps[0], timestamps[-1], 1/sample_rate)
            Calcium = pd.DataFrame(index=x_values.index, columns=new_timestamps)
            for feature in x_values.index:
                interpolator = interpolate.interp1d(timestamps, x_values.loc[feature], kind='linear')
                Calcium.loc[feature] = interpolator(new_timestamps)

            x_values = SpikeSub  # Each row is a feature, each column corresponds to a timestamp
            new_timestamps= np.arange(timestamps[0], timestamps[-1], 1/sample_rate)
            Spike = pd.DataFrame(index=x_values.index, columns=new_timestamps)
            for feature in x_values.index:
                interpolator = interpolate.interp1d(timestamps, x_values.loc[feature], kind='linear')
                Spike.loc[feature] = interpolator(new_timestamps)

            Carray=Calcium.values.T.astype(float)
            Sarray=Spike.values.T.astype(float)

            firstframe+=len(CalciumSub.T)
            rec_dur=len(CalciumSub.T)
            rec_dur_sec= timestamps[-1]- timestamps[0]

            print(len(CalciumSub.T), len(Calcium.T), rec_dur_sec)

            # Normalize traces
            
            #Carray=zscore(Carray, axis=0)
            #min=np.min(Carray,axis=0) 
            #Carray=Carray-min
            #Sarray=zscore(Sarray, axis=0)
            #min=np.min(Sarray,axis=0) 
            #Sarray=Sarray-min
            
            
            # Deal with dropped frames (failure to acquire miniscope images) IGNORE CAUSE TIMESTAMPS TAKEN INTO ACCOUNT

            list_droppedframes = literal_eval(dict_Stamps[session][0][3])    

            numbdropfr= 0   
            droppedframes_inrec=[]
            for item in list_droppedframes: 
                if item < (round(StartTimeMiniscope*minian_freq) + rec_dur_sec) and item > round(StartTimeMiniscope*minian_freq):
                    numbdropfr+=1                        

            EndTime = StartTime + rec_dur_sec # (upd_rec_dur/minian_freq) # in seconds
            previousEndTime=EndTime 

            #print(session, ': starts at', round(StartTime,1), 's & ends at', round(EndTime,1), 's (', round(upd_rec_dur/minian_freq,1), 's duration, ', numbdropfr, 'dropped frames, minian frequency =', minian_freq, 'Hz, drug = ', drug, ')...') 
            print(session, ': starts at', round(StartTime,1), 's & ends at', round(EndTime,1), 's (', round(rec_dur_sec,1), 's duration, ', numbdropfr, 'dropped frames, minian frequency =', minian_freq, 'Hz, drug = ', drug, ')...') 
            sentence1= f"... kept values = {kept_uniq_unit_List}"
            print(sentence1)


            nb_unit=len(CalciumSub)
            if nb_unit==0:
                continue  # next iteration


            # Upscale scoring to miniscope frequency

            scale_factor=minian_freq/0.2  #cause scoring was done in 5 seconds bin, ie 0.2 Hz
            SleepScoredTS=dict_Scoring[session]
            SleepScoredTS_upscaled = np.repeat(SleepScoredTS, scale_factor, axis=0)
            StartTime_frame=round(StartTime*minian_freq)
            SleepScoredTS_upscaled_ministart=SleepScoredTS_upscaled[StartTime_frame:StartTime_frame+rec_dur]


            # Remove N2 stage

            #SleepScoredTS_upscaled_ministart[SleepScoredTS_upscaled_ministart == 0.5] = 0
            #mapp = {1.5: 'Wake', 0: 'NREM',  1: 'REM'}
            mapp = {1.5: 'Wake', 0.5: 'N2', 0: 'NREM',  1: 'REM'}

            # Determine each substate identity and duration
            array=SleepScoredTS_upscaled_ministart
            substates_duration = [len(list(group)) for key, group in groupby(array)]
            substates_identity = [key for key, _ in groupby(array)]
            substates_end = np.array(substates_duration).cumsum()        
            substates_start =np.append([0],substates_end[:-1]) #substates_start =np.append([1],substates_end+1) create 1 second gap
            substates_identity = [mapp[num] for num in substates_identity]
            substates = pd.DataFrame(list(zip(substates_identity, substates_duration, substates_start, substates_end)), columns=['Identity', 'Duration', 'Start','End'])

            # Aligned LFP to Calcium traces
            PFCch1=int(allchannels[mice][0].split(',')[0])
            PFCch2=int(allchannels[mice][0].split(',')[1])
            CA1ch1=int(allchannels[mice][2].split(',')[0])
            CA1ch2=int(allchannels[mice][2].split(',')[1])
            S1ch1=int(allchannels[mice][1].split(',')[0])
            S1ch2=int(allchannels[mice][1].split(',')[1])
            EMGch1=int(allchannels[mice][3])
            
            All=dict_LFP[session]
            Allcut=All[round(StartTime*1000): round(EndTime*1000),:]

            PFC  =  Allcut[:, PFCch1]-Allcut[:, PFCch2] 
            CA1  =  Allcut[:, CA1ch1]-Allcut[:, CA1ch2] 
            S1  =  Allcut[:, S1ch1]-Allcut[:, S1ch2] 
            EMG  =  Allcut[:, EMGch1]

            SigmaPFC=get_filt_env(PFC, 10, 18, 1000)
            SigmaS1=get_filt_env(S1, 10, 18, 1000)
            ThetaCA1=get_filt_env(CA1, 5, 9, 1000)
            DeltaS1=get_filt_env(S1, 1, 4, 1000)       
            DeltaPFC=get_filt_env(PFC, 1, 4, 1000)          
            SOS1=get_filt_env(S1, .1, 1, 1000)       
            SOPFC=get_filt_env(PFC, .1, 1, 1000)       
            BetaPFC=get_filt_env(PFC, 16, 35, 1000)
            BetaS1=get_filt_env(S1, 16, 35, 1000)

            # Downscale to Calcium traces
            SigmaPFCds = griddata(np.linspace(0, 1, len(SigmaPFC)), SigmaPFC, np.linspace(0, 1, len(Carray)), method='linear')
            SigmaS1ds = griddata(np.linspace(0, 1, len(SigmaS1)), SigmaS1, np.linspace(0, 1, len(Carray)), method='linear')
            ThetaCA1ds = griddata(np.linspace(0, 1, len(ThetaCA1)), ThetaCA1, np.linspace(0, 1, len(Carray)), method='linear')
            DeltaS1ds = griddata(np.linspace(0, 1, len(DeltaS1)), DeltaS1, np.linspace(0, 1, len(Carray)), method='linear')
            DeltaPFCds = griddata(np.linspace(0, 1, len(DeltaPFC)), DeltaPFC, np.linspace(0, 1, len(Carray)), method='linear')
            BetaPFCds = griddata(np.linspace(0, 1, len(BetaPFC)), BetaPFC, np.linspace(0, 1, len(Carray)), method='linear')
            BetaS1ds = griddata(np.linspace(0, 1, len(BetaS1)), BetaS1, np.linspace(0, 1, len(Carray)), method='linear')
            SOS1ds = griddata(np.linspace(0, 1, len(SOS1)), SOS1, np.linspace(0, 1, len(Carray)), method='linear')
            SOPFCds = griddata(np.linspace(0, 1, len(SOPFC)), SOPFC, np.linspace(0, 1, len(Carray)), method='linear')

            for m in mapp:

                # Correlation between each neurons according to vigilance states 

                Bool = (SleepScoredTS_upscaled_ministart == m)
                Carray_VigSpe = Carray.copy()
                Carray_VigSpe = Carray_VigSpe[0:np.shape(SleepScoredTS_upscaled_ministart)[0],:] # if Calcium imaging longer than LFP rec
                Carray_VigSpe = Carray_VigSpe[Bool, :]
                
                CaCorrVigStateMatrixName=f'CaCorr{mapp[m]}Matrix{drug}'
                CaCorrVigStateMatrix = locals()[CaCorrVigStateMatrixName]
                CaCorrMatrix=[]
                CaCorrMatrix = pd.DataFrame(columns=[f'{mice}{str(i)}' for i in range(len(B))], index=[i for i in range(len(B))])
                
                Sarray_VigSpe = Sarray.copy()
                Sarray_VigSpe = Sarray_VigSpe[0:np.shape(SleepScoredTS_upscaled_ministart)[0],:] # if Calcium imaging longer than LFP rec
                Sarray_VigSpe = Sarray_VigSpe[Bool, :]    

                SpCorrVigStateMatrixName=f'SpCorr{mapp[m]}Matrix{drug}'
                SpCorrVigStateMatrix = locals()[SpCorrVigStateMatrixName]
                SpCorrMatrix=[]
                SpCorrMatrix = pd.DataFrame(columns=[f'{mice}{str(i)}' for i in range(len(B))], index=[i for i in range(len(B))])

                if saveexcel:
                    RawCaTracesVigStateMatrixName=f'RawCaTraces{mapp[m]}_{drug}'
                    RawCaTracesVigStateMatrix = locals()[RawCaTracesVigStateMatrixName]
                    RawCaTraces=[]
                    RawCaTraces=pd.DataFrame(Carray_VigSpe, columns=[f"{mice}{str(i).replace('[','').replace(']','')}" for i in kept_uniq_unit_List])
                    unique_columns = RawCaTraces.columns[RawCaTraces.columns.to_series().duplicated()] # remove units that has an empty unique index '[]'
                    RawCaTraces = RawCaTraces.drop(columns=unique_columns)
                    RawCaTracesVigStateMatrix.append(RawCaTraces)                

                    RawSpTracesVigStateMatrixName=f'RawSpTraces{mapp[m]}_{drug}'
                    RawSpTracesVigStateMatrix = locals()[RawSpTracesVigStateMatrixName]
                    RawSpTraces=[]
                    RawSpTraces=pd.DataFrame(Sarray_VigSpe, columns=[f"{mice}{str(i).replace('[','').replace(']','')}" for i in kept_uniq_unit_List])
                    unique_columns = RawSpTraces.columns[RawSpTraces.columns.to_series().duplicated()] # remove units that has an empty unique index '[]'
                    RawSpTraces = RawSpTraces.drop(columns=unique_columns)
                    RawSpTracesVigStateMatrix.append(RawSpTraces) 
                
                for unit in range(nb_unit): 

                    Carray_unit =Carray_VigSpe[:,unit]
                    Darray_unit =Sarray_VigSpe[:,unit] # on deconv spike not on spike rate 
                    otherunit_range = [x for x in range(nb_unit) if x != unit]

                    for unit2 in range(nb_unit):

                        Carray_unit2 =Carray_VigSpe[:,unit2]
                        Darray_unit2 =Sarray_VigSpe[:,unit2]         

                        indexMapp = str(np.where(B[session] == Calcium.index[unit])[0]).replace('[','').replace(']','')
                        indexMapp2 = np.where(B[session] == Calcium.index[unit2])[0]

                        if any(indexMapp) and len(indexMapp2)>0:      
                                            
                            corr_matrix = np.corrcoef(Carray_unit, Carray_unit2)
                            CaCorrMatrix[f'{mice}{indexMapp}'][indexMapp2]={1: 0.99999, -1: -0.99999}.get(np.round(corr_matrix[0, 1],5), np.round(corr_matrix[0, 1],5))
                            
                            corr_matrix = np.corrcoef(Darray_unit, Darray_unit2)
                            SpCorrMatrix[f'{mice}{indexMapp}'][indexMapp2]={1: 0.99999, -1: -0.99999}.get(np.round(corr_matrix[0, 1],5), np.round(corr_matrix[0, 1],5))
                            
                CaCorrVigStateMatrix.append(CaCorrMatrix)
                SpCorrVigStateMatrix.append(SpCorrMatrix) 

            TotCaCorrName=f'TotCaCorr{drug}'
            TotCaCorr = locals()[TotCaCorrName]
            TotSpCorrName=f'TotSpCorr{drug}'
            TotSpCorr = locals()[TotSpCorrName]
            
            TotCaCorrMatrix=[]
            TotCaCorrMatrix = pd.DataFrame(columns=[f'{mice}{str(i)}' for i in range(len(B))], index=[i for i in range(len(B))])
            TotSpCorrMatrix=[]
            TotSpCorrMatrix = pd.DataFrame(columns=[f'{mice}{str(i)}' for i in range(len(B))], index=[i for i in range(len(B))])

            for unit in range(nb_unit): 

                Carray_unit =Carray[:,unit]
                Darray_unit =Sarray[:,unit]
                peaks, _ = find_peaks(Darray_unit)
                Sarray_unit=np.zeros(len(Darray_unit))
                Sarray_unit[peaks]=1    

                # Correlation with LFP power density
        
                corr_matrix = np.corrcoef(Carray_unit, SigmaPFCds)
                CorrSigmaPFC= np.round(corr_matrix[0, 1],5)
                corCorrected = {1: 0.99999, -1: -0.99999}.get(np.round(corr_matrix[0, 1],5), np.round(corr_matrix[0, 1],5))
                Z_CorrSigmaPFC = np.arctanh(corCorrected)
                corr_matrix = np.corrcoef(Carray_unit, SigmaS1ds)
                CorrSigmaS1= np.round(corr_matrix[0, 1],5)
                corCorrected = {1: 0.99999, -1: -0.99999}.get(np.round(corr_matrix[0, 1],5), np.round(corr_matrix[0, 1],5))
                Z_CorrSigmaS1 = np.arctanh(corCorrected)
                corr_matrix = np.corrcoef(Carray_unit, ThetaCA1ds)
                CorrThetaCA1= np.round(corr_matrix[0, 1],5)
                corCorrected = {1: 0.99999, -1: -0.99999}.get(np.round(corr_matrix[0, 1],5), np.round(corr_matrix[0, 1],5))
                Z_CorrThetaCA1 = np.arctanh(corCorrected)
                corr_matrix = np.corrcoef(Carray_unit, DeltaS1ds)
                CorrDeltaS1= np.round(corr_matrix[0, 1],5)
                corCorrected = {1: 0.99999, -1: -0.99999}.get(np.round(corr_matrix[0, 1],5), np.round(corr_matrix[0, 1],5))
                Z_CorrDeltaS1 = np.arctanh(corCorrected)            
                corr_matrix = np.corrcoef(Carray_unit, DeltaPFCds)
                CorrDeltaPFC= np.round(corr_matrix[0, 1],5)
                corCorrected = {1: 0.99999, -1: -0.99999}.get(np.round(corr_matrix[0, 1],5), np.round(corr_matrix[0, 1],5))
                Z_CorrDeltaPFC = np.arctanh(corCorrected)
                corr_matrix = np.corrcoef(Carray_unit, BetaPFCds)
                CorrBetaPFC= np.round(corr_matrix[0, 1],5)
                corCorrected = {1: 0.99999, -1: -0.99999}.get(np.round(corr_matrix[0, 1],5), np.round(corr_matrix[0, 1],5))
                Z_CorrBetaPFC = np.arctanh(corCorrected)
                corr_matrix = np.corrcoef(Carray_unit, BetaS1ds)
                CorrBetaS1= np.round(corr_matrix[0, 1],5)
                corCorrected = {1: 0.99999, -1: -0.99999}.get(np.round(corr_matrix[0, 1],5), np.round(corr_matrix[0, 1],5))
                Z_CorrBetaS1 = np.arctanh(corCorrected)
                corr_matrix = np.corrcoef(Carray_unit, SOS1ds)
                CorrSOS1= np.round(corr_matrix[0, 1],5)
                corCorrected = {1: 0.99999, -1: -0.99999}.get(np.round(corr_matrix[0, 1],5), np.round(corr_matrix[0, 1],5))
                Z_CorrSOS1 = np.arctanh(corCorrected)            
                corr_matrix = np.corrcoef(Carray_unit, SOPFCds)
                CorrSOPFC= np.round(corr_matrix[0, 1],5)
                corCorrected = {1: 0.99999, -1: -0.99999}.get(np.round(corr_matrix[0, 1],5), np.round(corr_matrix[0, 1],5))
                Z_CorrSOPFC = np.arctanh(corCorrected)
                
                # Population Coupling independent of vigilance states 

                Carray_Population =np.mean(Carray[:,otherunit_range], axis=1)
                corr_matrix = np.corrcoef(Carray_unit, Carray_Population)
                TotCaPopCoupling= np.round(corr_matrix[0, 1],5)
                corCorrected = {1: 0.99999, -1: -0.99999}.get(np.round(corr_matrix[0, 1],5), np.round(corr_matrix[0, 1],5))
                TotZ_CaPopCoupling = np.arctanh(corCorrected)

                Sarray_Population =np.mean(Sarray[:,otherunit_range], axis=1)
                corr_matrix = np.corrcoef(Darray_unit, Sarray_Population)
                TotSpPopCoupling = np.round(corr_matrix[0, 1],5)
                corCorrected = {1: 0.99999, -1: -0.99999}.get(np.round(corr_matrix[0, 1],5), np.round(corr_matrix[0, 1],5))
                TotZ_SpPopCoupling= np.arctanh(corCorrected)

                # Correlation between each neurons independent of vigilance states 

                otherunit_range = [x for x in range(nb_unit) if x != unit]
                for unit2 in otherunit_range:

                    Carray_unit2 =Carray[:,unit2]
                    Darray_unit2 =Sarray[:,unit2]  
                
                    indexMapp = str(np.where(B[session] == Calcium.index[unit])[0]).replace('[','').replace(']','')
                    indexMapp2 = np.where(B[session] == Calcium.index[unit2])[0]

                    if any(indexMapp) and len(indexMapp2)>0:    
                        
                        corr_matrix = np.corrcoef(Carray_unit, Carray_unit2)
                        TotCaCorrMatrix[f'{mice}{indexMapp}'][indexMapp2]={1: 0.99999, -1: -0.99999}.get(np.round(corr_matrix[0, 1],5), np.round(corr_matrix[0, 1],5))
                        
                        corr_matrix = np.corrcoef(Darray_unit, Darray_unit2)
                        TotSpCorrMatrix[f'{mice}{indexMapp}'][indexMapp2]={1: 0.99999, -1: -0.99999}.get(np.round(corr_matrix[0, 1],5), np.round(corr_matrix[0, 1],5))

                # For each substates
                for index in range(len(substates)):
                    
                    ca_input_sub=Carray_unit[substates.Start[index]:substates.End[index]]
                    ds_input_sub=Darray_unit[substates.Start[index]:substates.End[index]]
                    sp_input_sub=Sarray_unit[substates.Start[index]:substates.End[index]]

                    VigilanceState_GlobalResults.loc[counter, 'Mice'] = mice
                    VigilanceState_GlobalResults.loc[counter, 'Session'] = session
                    VigilanceState_GlobalResults.loc[counter, 'Session_Time'] = None 
                    
                    indexMapp = np.where(B[session] == Calcium.index[unit])[0]
                    VigilanceState_GlobalResults.loc[counter, 'Unique_Unit'] = indexMapp if len(indexMapp)>0 else None
                    VigilanceState_GlobalResults.loc[counter, 'UnitNumber'] = unit
                    VigilanceState_GlobalResults.loc[counter, 'UnitValue'] = Calcium.index[unit]

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
                    
                    VigilanceState_GlobalResults.loc[counter, 'TotCaPopCoupling'] = TotCaPopCoupling
                    VigilanceState_GlobalResults.loc[counter, 'TotZ_CaPopCoupling'] = TotZ_CaPopCoupling
                    VigilanceState_GlobalResults.loc[counter, 'TotSpPopCoupling'] = TotSpPopCoupling
                    VigilanceState_GlobalResults.loc[counter, 'TotZ_SpPopCoupling'] = TotZ_SpPopCoupling
                    
                    VigilanceState_GlobalResults.loc[counter, 'SOS1_corr'] = CorrSOS1
                    VigilanceState_GlobalResults.loc[counter, 'SOPFC_corr'] = CorrSOPFC
                    VigilanceState_GlobalResults.loc[counter, 'DeltaS1_corr'] = CorrDeltaS1
                    VigilanceState_GlobalResults.loc[counter, 'DeltaPFC_corr'] = CorrDeltaPFC
                    VigilanceState_GlobalResults.loc[counter, 'ThetaCA1_corr'] = CorrThetaCA1
                    VigilanceState_GlobalResults.loc[counter, 'SigmaS1_corr'] = CorrSigmaS1
                    VigilanceState_GlobalResults.loc[counter, 'SigmaPFC_corr'] = CorrSigmaPFC
                    VigilanceState_GlobalResults.loc[counter, 'BetaS1_corr'] = CorrBetaS1
                    VigilanceState_GlobalResults.loc[counter, 'BetaPFC_corr'] = CorrBetaPFC

                    VigilanceState_GlobalResults.loc[counter, 'Z_SOS1_corr'] = Z_CorrSOS1
                    VigilanceState_GlobalResults.loc[counter, 'Z_SOPFC_corr'] = Z_CorrSOPFC
                    VigilanceState_GlobalResults.loc[counter, 'Z_DeltaS1_corr'] = Z_CorrDeltaS1
                    VigilanceState_GlobalResults.loc[counter, 'Z_DeltaPFC_corr'] = Z_CorrDeltaPFC
                    VigilanceState_GlobalResults.loc[counter, 'Z_ThetaCA1_corr'] = Z_CorrThetaCA1
                    VigilanceState_GlobalResults.loc[counter, 'Z_SigmaS1_corr'] = Z_CorrSigmaS1
                    VigilanceState_GlobalResults.loc[counter, 'Z_SigmaPFC_corr'] = Z_CorrSigmaPFC
                    VigilanceState_GlobalResults.loc[counter, 'Z_BetaS1_corr'] = Z_CorrBetaS1
                    VigilanceState_GlobalResults.loc[counter, 'Z_BetaPFC_corr'] = Z_CorrBetaPFC

                    StatesCaCorrName=f'StatesCaCorr{substates.Identity[index]}Matrix{drug}'
                    StatesCaCorrMatrix = locals()[StatesCaCorrName]
                    IterationMatrixName=f'ITStatesCaCorr{substates.Identity[index]}Matrix{drug}'
                    IterationMatrix = locals()[IterationMatrixName]
                    
                    # Correlation between each neurons dependent of vigilance states 
                    if (substates.Duration[index]/minian_freq)>=20:
                        otherunit_range = [x for x in range(nb_unit) if x != unit]
                        for unit2 in otherunit_range:
                            Carray_unit2 =Carray[:,unit2]
                            ca2_input_sub=Carray_unit2[substates.Start[index]:substates.End[index]]
                        
                            indexMapp = str(np.where(B[session] == Calcium.index[unit])[0]).replace('[','').replace(']','')
                            indexMapp2 = np.where(B[session] == Calcium.index[unit2])[0]

                            if any(indexMapp) and len(indexMapp2)>0:    

                                corr_matrix = np.corrcoef(ca_input_sub, ca2_input_sub)
                                StatesCaCorrMatrix[f'{mice}{indexMapp}'][indexMapp2]=StatesCaCorrMatrix[f'{mice}{indexMapp}'][indexMapp2]+  np.arctanh({1: 0.99999, -1: -0.99999}.get(np.round(corr_matrix[0, 1],5), np.round(corr_matrix[0, 1],5))) if not math.isnan(StatesCaCorrMatrix[f'{mice}{indexMapp}'][indexMapp2]) else  np.arctanh({1: 0.99999, -1: -0.99999}.get(np.round(corr_matrix[0, 1],5), np.round(corr_matrix[0, 1],5)))
                                if not math.isnan(StatesCaCorrMatrix[f'{mice}{indexMapp}'][indexMapp2]):
                                    if not math.isnan(IterationMatrix[f'{mice}{indexMapp}'][indexMapp2]):
                                        IterationMatrix[f'{mice}{indexMapp}'][indexMapp2]= IterationMatrix[f'{mice}{indexMapp}'][indexMapp2]+1 
                                    else:
                                        IterationMatrix[f'{mice}{indexMapp}'][indexMapp2]= 1 
                    counter+=1

            TotCaCorr.append(TotCaCorrMatrix)
            TotSpCorr.append(TotSpCorrMatrix)

        dataCaCorr={}
        dataSpCorr={}
        dataCaCorr2={}
        dataSpCorr2={}
        dataStatesCaCorr={}

        for Drug in Drugs:
            TotCaCorrName=f'TotCaCorr{Drug}'
            TotCaCorr = locals()[TotCaCorrName]
            TotSpCorrName=f'TotSpCorr{Drug}'
            TotSpCorr = locals()[TotSpCorrName]
            if len(TotCaCorr)>0: 
                combined_df = pd.concat(TotCaCorr, ignore_index=False)
                IterationNb=combined_df.groupby(combined_df.index).count()
                combined_df = combined_df.groupby(combined_df.index).sum() #mean
                combined_df.index= [mice + str(idx) for idx in combined_df.index]
                combined_df = combined_df[~(combined_df.fillna(0) == 0).all(axis=1)]
                combined_df = combined_df.loc[:, ~(combined_df.fillna(0) == 0).all(axis=0)]
                IterationNb = IterationNb[~(IterationNb.fillna(0) == 0).all(axis=1)]
                IterationNb = IterationNb.loc[:, ~(IterationNb.fillna(0) == 0).all(axis=0)]
                IterationNb.index=combined_df.index
                dataCaCorr2[f'{Drug}']=combined_df

                combined_df = pd.concat(TotCaCorr, ignore_index=False)
                combined_df = combined_df.applymap(lambda x: np.arctanh(x) if not pd.isna(x) else np.nan)
                combined_df = combined_df.groupby(combined_df.index).sum() #mean
                combined_df.index = [mice + str(idx) for idx in combined_df.index]
                combined_df = combined_df[~(combined_df.fillna(0) == 0).all(axis=1)]
                combined_df = combined_df.loc[:, ~(combined_df.fillna(0) == 0).all(axis=0)]
                dataCaCorr2[f'Z_{Drug}']=combined_df
                dataCaCorr2[f'{Drug}_IterationNb']=IterationNb

                combined_df = pd.concat(TotSpCorr, ignore_index=False)
                IterationNb=combined_df.groupby(combined_df.index).count()
                combined_df = combined_df.groupby(combined_df.index).sum() #mean
                combined_df.index= [mice + str(idx) for idx in combined_df.index]
                combined_df = combined_df[~(combined_df.fillna(0) == 0).all(axis=1)]
                combined_df = combined_df.loc[:, ~(combined_df.fillna(0) == 0).all(axis=0)]
                IterationNb = IterationNb[~(IterationNb.fillna(0) == 0).all(axis=1)]
                IterationNb = IterationNb.loc[:, ~(IterationNb.fillna(0) == 0).all(axis=0)]
                IterationNb.index=combined_df.index
                dataSpCorr2[f'{Drug}']=combined_df

                combined_df = pd.concat(TotSpCorr, ignore_index=False)
                combined_df = combined_df.applymap(lambda x: np.arctanh(x) if not pd.isna(x) else np.nan)
                combined_df = combined_df.groupby(combined_df.index).sum() #mean
                combined_df.index = [mice + str(idx) for idx in combined_df.index]
                combined_df = combined_df[~(combined_df.fillna(0) == 0).all(axis=1)]
                combined_df = combined_df.loc[:, ~(combined_df.fillna(0) == 0).all(axis=0)]
                dataSpCorr2[f'Z_{Drug}']=combined_df
                dataSpCorr2[f'{Drug}_IterationNb']=IterationNb

            for m in mapp:
                IterationMatrixName=f'ITStatesCaCorr{mapp[m]}Matrix{Drug}'
                IterationMatrix = locals()[IterationMatrixName]
                StatesCaCorrMatrixName=f'StatesCaCorr{mapp[m]}Matrix{Drug}'
                StatesCaCorrMatrix = locals()[StatesCaCorrMatrixName]
                if len(StatesCaCorrMatrix)>0: # cause sometimes no Baseline conditions in CGP experiments
                    StatesCaCorrMatrix.index = [mice + str(idx) for idx in StatesCaCorrMatrix.index]
                    IterationMatrix.index = StatesCaCorrMatrix.index
                    dataStatesCaCorr[f'{Drug}_{mapp[m]}']=StatesCaCorrMatrix
                    dataStatesCaCorr[f'Iteration_{Drug}_{mapp[m]}']=IterationMatrix

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
            filenameOut = folder_to_save / f'VigSt_Global_{mice}.xlsx'
            writer = pd.ExcelWriter(filenameOut)
            VigilanceState_GlobalResults.to_excel(writer)
            writer.close()
            excel_writerRawCa.close()
            excel_writerRawSp.close()

        filenameOut = folder_to_save / f'VigSt_CaCorr_{mice}.pkl'
        with open(filenameOut, 'wb') as pickle_file:
            pickle.dump(dataCaCorr, pickle_file)
        filenameOut = folder_to_save / f'VigSt_SpCorr_{mice}.pkl'
        with open(filenameOut, 'wb') as pickle_file:
            pickle.dump(dataSpCorr, pickle_file)
            
        filenameOut = folder_to_save / f'TotCaCorr_{mice}.pkl'
        with open(filenameOut, 'wb') as pickle_file:
            pickle.dump(dataCaCorr2, pickle_file)
        filenameOut = folder_to_save / f'TotSpCorr_{mice}.pkl'
        with open(filenameOut, 'wb') as pickle_file:
            pickle.dump(dataSpCorr2, pickle_file)

        filenameOut = folder_to_save / f'StatesCaCorr_{mice}.pkl'
        with open(filenameOut, 'wb') as pickle_file:
            pickle.dump(dataStatesCaCorr, pickle_file)

        filenameOut = folder_to_save / f'VigSt_Global_{mice}.pkl'
        with open(filenameOut, 'wb') as pickle_file:
            pickle.dump(VigilanceState_GlobalResults, pickle_file)

        