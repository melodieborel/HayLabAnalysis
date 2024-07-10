# Sleep scoring

#######################################################################################
                                # Choose Analysis #
#######################################################################################

#Method=1 # for AB analysis
Method=0 # for AH analysis

DrugExperiment=1 #if CGP Experiment
#DrugExperiment=0 #if Baseline Experiment

#######################################################################################
                                # Load packages #
#######################################################################################

import quantities as pq
import numpy as np
import neo
from pathlib import Path
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
import scipy
from scipy import interpolate
from scipy import fftpack
import os
from scipy import signal
from scipy import fftpack
from datetime import datetime
import shutil
from scipy.signal import find_peaks
from scipy.signal import chirp, find_peaks, peak_widths
from itertools import groupby

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

# Perform analysis for each mouse

MiceList=['BlackLinesOK', 'BlueLinesOK', 'GreenDotsOK','Purple' ,'ThreeColDotsOK'] if DrugExperiment else ['BlackLinesOK', 'BlueLinesOK', 'GreenDotsOK', 'GreenLinesOK', 'Purple', 'RedLinesOK','ThreeColDotsOK', 'ThreeBlueCrossesOK']
dpath0 = "//10.69.168.1/crnldata/waking/audrey_hay/L1imaging/AnalysedMarch2023/Gaelle/CGP/" if DrugExperiment else "//10.69.168.1/crnldata/waking/audrey_hay/L1imaging/AnalysedMarch2023/Gaelle/Baseline_recording_ABmodified/"

# Process

for mice in MiceList:
    
    dpath=Path(dpath0 + mice)
    # Load sleep score and Ca2+ time series numpy arrays
    
    #nb_sessions = sum(1 for p in dpath.iterdir() if p.is_dir() and p.name.startswith("session"))    
    #sessions = [folder.name for folder in dpath.iterdir() if folder.is_dir() and "session" in folder.name]

    sessions, sessions_path = find_session_folders(dpath)
    nb_sessions=len(sessions)

    for sess,session in enumerate(sessions):  
        
        session_path=Path(sessions_path[sess])

        folder_base = Path(session_path) / f'OpenEphys/'
        print(folder_base)

        filename2 = folder_base / f'RawDataChannelExtractedDS.npy'
        if Method: 
            EMGbooleaninput = folder_base / f'EMGframeBoolean_AB.pkl'  
        else:
            try:
                EMGbooleaninput =folder_base / f'EMGframeBoolean.pkl'
                EMGboolean = pd.read_pickle(EMGbooleaninput)
                print('EMGframeBoolean')
            except: 
                EMGbooleaninput =folder_base / f'EMGframeBoolean_AH.pkl'
                EMGboolean = pd.read_pickle(EMGbooleaninput)
                print('EMGframeBoolean_AH')

        Channels = '//10.69.168.1/crnldata/waking/audrey_hay/L1imaging/AnalysedMarch2023/LFPChannels_perMice.xlsx' 

        All = np.load(filename2, mmap_mode= 'r')

        allchannels = pd.read_excel(Channels)
        
        CA1ch1=int(allchannels[mice][2].split(',')[0])
        CA1ch2=int(allchannels[mice][2].split(',')[1])
        CA1  =  All[:, CA1ch1]-All[:, CA1ch2] 

        PFCch1=int(allchannels[mice][0].split(',')[0])
        PFCch2=int(allchannels[mice][0].split(',')[1])
        PFC  =  All[:, PFCch1]-All[:, PFCch2] 

        S1ch1=int(allchannels[mice][1].split(',')[0])
        S1ch2=int(allchannels[mice][1].split(',')[1])
        S1  =  All[:, S1ch1]-All[:, S1ch2] 

        ThetaCh = CA1
        Beta1Ch = PFC
        Beta2Ch = S1

        EMGch=int(allchannels[mice][3])
        EMG  =  All[:, EMGch]

        ScoringVectorLength = len(EMG)
        ScoringVector = np.zeros((ScoringVectorLength))

        WakeStatus = np.zeros((ScoringVectorLength))
        WakeStatusCons = np.zeros((ScoringVectorLength))
        EMGStatusBoolLib = EMGboolean.BooleanLiberal
        EMGStatusBoolCons = EMGboolean.BooleanConservative
        WakeStatus[EMGStatusBoolLib] = 1
        WakeStatusCons[EMGStatusBoolCons] = 1

        #####################################
        ############### Theta ###############
        #####################################    
        
        # Filtre parameter:
        f_lowcut = 5.
        f_hicut = 9.
        fs = 1000
        nyq = 0.5 * fs
        N = 4                 # Filtre order
        Wn = [f_lowcut/nyq,f_hicut/nyq]  # Nyquist frequency fraction

        # Filtering:
        b, a = signal.butter(N, Wn, 'band')

        # Remove Wake activity:
        #ThetaCh[EMGStatusBoolCons]=0

        filt_Theta = signal.filtfilt(b, a, ThetaCh)

        # Parameter and computation of CWT
        w = 30.
        freq = np.linspace(5, 9, 8)
        widths = w*fs / (2*freq*np.pi)
        ThetaCWT = signal.cwt(filt_Theta, signal.morlet2, widths, w=w)

        # Projection calculation
        absThetaCWT = np.absolute(ThetaCWT)
        from scipy import stats

        zabsThetaCWT = stats.zscore(absThetaCWT, axis=None)

        proj_ThetaCWT = np.max(zabsThetaCWT, axis = 0)/8
        sdproj_ThetaCWT = np.std(proj_ThetaCWT)
        meanproj_ThetaCWT = np.mean(proj_ThetaCWT)

        numpnts = EMG.size
        ThetaStatus = np.zeros(numpnts)
        for ind in range(numpnts):
            if proj_ThetaCWT[ind]>(meanproj_ThetaCWT+1.4*sdproj_ThetaCWT):
                ThetaStatus[ind] = 1

        #####################################
        ############### BETA ################
        #####################################     
        # Filtre parameter:
        f_lowcut = 10.
        f_hicut = 18.
        fs = 1000
        nyq = 0.5 * fs
        N = 4                 # Filtre order
        Wn = [f_lowcut/nyq,f_hicut/nyq]  # Nyquist frequency fraction

        # Filtering:
        b, a = signal.butter(N, Wn, 'band')
        filt_Beta1 = signal.filtfilt(b, a, Beta1Ch)
        filt_Beta2 = signal.filtfilt(b, a, Beta2Ch)

        # Parameter and computation of CWT
        w = 10.
        freq = np.linspace(10, 18, 16)
        widths = w*fs / (2*freq*np.pi)
        Beta1CWT = signal.cwt(filt_Beta1, signal.morlet2, widths, w=w)
        Beta2CWT = signal.cwt(filt_Beta2, signal.morlet2, widths, w=w)

        # Projection calculation
        absBeta1CWT = np.absolute(Beta1CWT)
        absBeta2CWT = np.absolute(Beta2CWT)
        from scipy import stats

        zabsBeta1CWT = stats.zscore(absBeta1CWT, axis=None)
        zabsBeta2CWT = stats.zscore(absBeta2CWT, axis=None)

        proj_Beta1CWT = np.max(zabsBeta1CWT, axis = 0)/16
        proj_Beta2CWT = np.max(zabsBeta2CWT, axis = 0)/16
        meanproj_Beta1CWT = np.mean(zabsBeta1CWT)
        meanproj_Beta2CWT = np.mean(zabsBeta2CWT)
        sdproj_Beta1CWT = np.std(proj_Beta1CWT)
        sdproj_Beta2CWT = np.std(proj_Beta2CWT)

        numpnts = EMG.size
        Beta1Status = np.zeros(numpnts)
        Beta2Status = np.zeros(numpnts)
        for ind in range(numpnts):
            if proj_Beta1CWT[ind]>(meanproj_Beta1CWT+3*sdproj_Beta1CWT):
                Beta1Status[ind] = 1
            if proj_Beta2CWT[ind]>(meanproj_Beta2CWT+3*sdproj_Beta2CWT):
                Beta2Status[ind] = 1


        ScoringVector = np.ones((ScoringVectorLength))    
        
        input_arr = ThetaStatus
        R = round(ScoringVectorLength/5000) # 5 sec bins
        split_arr = np.linspace(0, len(input_arr), num=R+1, dtype=int)
        dwnsmpl_subarr = np.split(input_arr, split_arr[1:])
        dwnsmpl_arrT = np.array( list( np.mean(item) for item in dwnsmpl_subarr[:-1] ) )

        th= 0.4 if Method else 0.25

        for i in range(len(dwnsmpl_arrT)):
            if dwnsmpl_arrT[i]<th: 
                dwnsmpl_arrT[i] = 0
            else:
                dwnsmpl_arrT[i] = 1  

        for i in range(len(dwnsmpl_arrT)-3):
            if i > (len(dwnsmpl_arrT)-3):
                break          
            elif (dwnsmpl_arrT[i]>0 and dwnsmpl_arrT[i+2]>0):
                dwnsmpl_arrT[i+1]=1
            elif (dwnsmpl_arrT[i]>0 and dwnsmpl_arrT[i+3]>0):
                dwnsmpl_arrT[i+1]=1
                dwnsmpl_arrT[i+2]=1

        for i in range(len(dwnsmpl_arrT)):
            if i > (len(dwnsmpl_arrT)-3):
                break
            elif (dwnsmpl_arrT[i]<1 and dwnsmpl_arrT[i+1]>0 and dwnsmpl_arrT[i+2]<1):
                dwnsmpl_arrT[i+1]=0


        input_arr1 = proj_Beta1CWT
        input_arr2 = proj_Beta2CWT

        R = round(ScoringVectorLength/5000) # 5 sec bins

        split_arr1 = np.linspace(0, len(input_arr1), num=R+1, dtype=int)
        split_arr2 = np.linspace(0, len(input_arr2), num=R+1, dtype=int)

        dwnsmpl_subarr1 = np.split(input_arr1, split_arr1[1:])
        dwnsmpl_subarr2 = np.split(input_arr2, split_arr2[1:])

        dwnsmpl_arr1B = np.array( list( np.mean(item) for item in dwnsmpl_subarr1[:-1] ) )
        dwnsmpl_arr2B = np.array( list( np.mean(item) for item in dwnsmpl_subarr2[:-1] ) )

        for i in range(len(dwnsmpl_arr1B)):
            if dwnsmpl_arr1B[i]<0.12: # arbitrary set
                dwnsmpl_arr1B[i] = 0
            else:
                dwnsmpl_arr1B[i] = 1  

        for i in range(len(dwnsmpl_arr1B)-3):
            if i > (len(dwnsmpl_arr1B)-3):
                break                 
            elif (dwnsmpl_arr1B[i]>0 and dwnsmpl_arr1B[i+2]>0):
                dwnsmpl_arr1B[i+1]=1
            elif (dwnsmpl_arr1B[i]>0 and dwnsmpl_arr1B[i+3]>0):
                dwnsmpl_arr1B[i+1]=1
                dwnsmpl_arr1B[i+2]=1

        for i in range(len(dwnsmpl_arr2B)):
            if dwnsmpl_arr2B[i]<0.12: # arbitrary set
                dwnsmpl_arr2B[i] = 0
            else:
                dwnsmpl_arr2B[i] = 1  

        for i in range(len(dwnsmpl_arr2B)-3):          
            if (dwnsmpl_arr2B[i]>0 and dwnsmpl_arr2B[i+2]>0):
                dwnsmpl_arr2B[i+1]=1
            elif (dwnsmpl_arr2B[i]>0 and dwnsmpl_arr2B[i+3]>0):
                dwnsmpl_arr2B[i+1]=1
                dwnsmpl_arr2B[i+2]=1


        input_arrW = WakeStatus
        R = round(ScoringVectorLength/5000) # 5 sec bins
        split_arr = np.linspace(0, len(input_arr), num=R+1, dtype=int)
        dwnsmpl_subarr = np.split(input_arrW, split_arr[1:])
        dwnsmpl_arrW = np.array( list( np.mean(item) for item in dwnsmpl_subarr[:-1] ) )

        for i in range(len(dwnsmpl_arrW)):
            if dwnsmpl_arrW[i]<0.4: # arbitrary set #0.2
                dwnsmpl_arrW[i] = 0
            else:
                dwnsmpl_arrW[i] = 1 

        dwnsmpl_arrT = dwnsmpl_arrT * 1
        dwnsmpl_arrW = dwnsmpl_arrW * 1.5
        dwnsmpl_arr1B = dwnsmpl_arr1B * 0.5
        dwnsmpl_arr2B = dwnsmpl_arr2B * 0.5

        ScoringVectorS = np.zeros((len(dwnsmpl_arrW)))
        if Method: 
            for ind in range(len(dwnsmpl_arrW)):
                if dwnsmpl_arr1B[ind]>0:
                    ScoringVectorS[ind] = 0.5
                if dwnsmpl_arr2B[ind]>0:
                    ScoringVectorS[ind] = 0.5
                if dwnsmpl_arrT[ind]>0:
                    ScoringVectorS[ind] = 1
                if dwnsmpl_arrW[ind]>0:
                    ScoringVectorS[ind] = 1.5

            array=ScoringVectorS
            substates_duration = [len(list(group)) for key, group in groupby(array)]
            substates_identity = [key for key, _ in groupby(array)]
            substates_end = np.array(substates_duration).cumsum()        
            substates_start =np.append([0],substates_end[:-1]) #substates_start =np.append([1],substates_end+1) create 1 second gap
            mapp = {0: 'NREM', 0.5: 'N2', 1: 'REM', 1.5: 'Wake'}
            substates_identity = [mapp[num] for num in substates_identity]
            substates = pd.DataFrame(list(zip(substates_identity, substates_duration, substates_start, substates_end)), columns=['Identity', 'Duration', 'Start','End'])

            ScoringVectorS2=ScoringVectorS.copy()
            for index in substates.index[substates.Identity == 'REM'].tolist():
                if index+2<=len(substates):
                    if substates.Identity[index+2]=='REM' and substates.Duration[index+1]<10: #if duration of a state between 2 REM inferior at 10bins(50s)
                        start=substates.Start[index+1]
                        end=substates.End[index+1]
                        ScoringVectorS2[start:end]=1
            ScoringVectorS=ScoringVectorS2
            #for i in range(len(ScoringVectorS)):
            #    if i > (len(ScoringVectorS)-4):
            #        break       
            #    #elif not ScoringVectorS[i]==1.5 and ScoringVectorS[i+1]==1.5 and not ScoringVectorS[i+2]==1.5:
            #    #    ScoringVectorS[i+1]=ScoringVectorS[i] # need at least 2consecutive bins / 10sec of Wake
            #    elif (ScoringVectorS[i]==0.5 and ScoringVectorS[i+1]==0 and ScoringVectorS[i+2]>0):
            #        ScoringVectorS[i+1]=0.5
            #    elif (ScoringVectorS[i]==1 and ScoringVectorS[i+1]==0):
            #        ScoringVectorS[i+1]=1 #can't be NREM after REM
            #    elif (ScoringVectorS[i]==1 and ScoringVectorS[i+1]==0.5):
            #        ScoringVectorS[i+1]=1 #can't be N2 after REM

        else: 
                
            ScoringVectorS = np.zeros((len(dwnsmpl_arrW)))
            for ind in range(len(dwnsmpl_arrW)):
                if dwnsmpl_arr2B[ind]>0:
                    ScoringVectorS[ind] = 0.5
                if dwnsmpl_arrT[ind]>0:
                    ScoringVectorS[ind] = 1
                if dwnsmpl_arrW[ind]>0:
                    ScoringVectorS[ind] = 1.5

            for i in range(len(ScoringVectorS)):
                if i > (len(ScoringVectorS)-4):
                    break          
                elif (ScoringVectorS[i]==0.5 and ScoringVectorS[i+1]==0 and ScoringVectorS[i+2]>0):
                    ScoringVectorS[i+1]=0.5
                elif (ScoringVectorS[i]==1 and ScoringVectorS[i+1]==0 and ScoringVectorS[i+2]>0):
                    ScoringVectorS[i+1]=1
                elif (ScoringVectorS[i]==1.5 and ScoringVectorS[i+1]<1.5 and (ScoringVectorS[i+2]==1.5 or ScoringVectorS[i+3]==1.5)):
                    ScoringVectorS[i+1]=1.5

        suffix='_AB' if Method else '_AH'

        filenameOut = folder_base / f'ScoredSleep{suffix}.npy'
        np.save(filenameOut, ScoringVectorS)