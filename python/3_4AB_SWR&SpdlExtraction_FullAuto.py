
# ## Load LFP and packages

from scipy import signal
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, Cursor
from scipy import fftpack
import pandas as pd
from pathlib import Path
import os
from IPython.display import display
from ipyfilechooser import FileChooser
from datetime import datetime
import shutil
from scipy.signal import find_peaks
from scipy.signal import chirp, find_peaks, peak_widths

# Perform analysis for each mouse

MiceList=['BlackLinesOK', 'BlueLinesOK', 'GreenDotsOK', 'GreenLinesOK', 'Purple', 'RedLinesOK','ThreeColDotsOK', 'ThreeBlueCrossesOK']
MiceList=['GreenLinesOK', 'Purple', 'RedLinesOK','ThreeColDotsOK', 'ThreeBlueCrossesOK']
dpath0 = "//10.69.168.1/crnldata/waking/audrey_hay/L1imaging/AnalysedMarch2023/Gaelle/Baseline_recording_ABmodified/"

for micename in MiceList:
    
    dpath=Path(dpath0 + micename)
    # Load sleep score and Ca2+ time series numpy arrays
    nb_sessions = sum(1 for p in dpath.iterdir() if p.is_dir() and p.name.startswith("session"))    
    sessions = [folder.name for folder in dpath.iterdir() if folder.is_dir() and "session" in folder.name]

    for session in sessions:  
        folder_base = Path(dpath) / session / f'OpenEphys/'
        print(folder_base)

        filename = folder_base / f'LFPwake0.npy'
        filename3 = folder_base / f'LFPwakeremoved.npy'
        filename2 = folder_base / f'RawDataChannelExtractedDS.npy'
        EMGbooleaninput = folder_base / f'EMGframeBoolean.pkl'
        Channels = '//10.69.168.1/crnldata/waking/audrey_hay/L1imaging/AnalysedMarch2023/LFPChannels_perMice.xlsx' 

        EMGboolean = pd.read_pickle(EMGbooleaninput)
        LFPwakeremoved = np.load(filename3, mmap_mode= 'r')
        All = np.load(filename2, mmap_mode= 'r')

        mice = os.path.basename(os.path.dirname(os.path.dirname(folder_base)))
        allchannels = pd.read_excel(Channels)
        
        CA1ch1=int(allchannels[mice][2].split(',')[0])
        CA1ch2=int(allchannels[mice][2].split(',')[1])
        CA1  =  All[:, CA1ch1]-All[:, CA1ch2] 
        CA1wakeremoved = LFPwakeremoved[:,CA1ch1]-LFPwakeremoved[:,CA1ch2] 

        PFCch1=int(allchannels[mice][0].split(',')[0])
        PFCch2=int(allchannels[mice][0].split(',')[1])
        PFC  =  All[:, PFCch1]-All[:, PFCch2] 
        PFCwakeremoved = LFPwakeremoved[:,PFCch1]-LFPwakeremoved[:,PFCch2] 

        S1ch1=int(allchannels[mice][1].split(',')[0])
        S1ch2=int(allchannels[mice][1].split(',')[1])
        S1  =  All[:, S1ch1]-All[:, S1ch2] 
        S1wakeremoved = LFPwakeremoved[:,S1ch1]-LFPwakeremoved[:,S1ch2]

        # SWR: 120-200 Hz

        # Filtre parameter:
        f_lowcut = 120.
        f_hicut = 200.
        fs = 1000
        nyq = 0.5 * fs
        N = 6                 # Filtre order
        Wn = [f_lowcut/nyq,f_hicut/nyq]  # Nyquist frequency fraction

        # Filtering:
        b, a = signal.butter(N, Wn, 'band')
        filt_CA1 = signal.filtfilt(b, a, CA1)
        filt_CA1wakeremoved = signal.filtfilt(b, a, CA1wakeremoved)

        # Parameter and computation of CWT
        w = 10.
        freq = np.linspace(120, 200, 80)
        widths = w*fs / (2*freq*np.pi)
        CA1NWcwt = signal.cwt(filt_CA1wakeremoved, signal.morlet2, widths, w=w)

        # Projection calculation
        absCA1NWcwt = np.absolute(CA1NWcwt)
        proj_CA1NWcwt = np.sum(absCA1NWcwt, axis = 0)/80
        sdproj_CA1cwt = np.std(proj_CA1NWcwt)
        #sd3proj_CA1cwt = sdproj_CA1cwt*3
        #sd5proj_CA1cwt = sdproj_CA1cwt*5
        sd6proj_CA1cwt = sdproj_CA1cwt*6
        #sd7proj_CA1cwt = sdproj_CA1cwt*7
        #sd8proj_CA1cwt = sdproj_CA1cwt*8

        # Conservative boolean filtering of CA1 filtered signal
        BooleanCons = EMGboolean['BooleanConservative']
        fCA1wake0C = filt_CA1.copy()
        fCA1wake0C[BooleanCons] = 0
        CA1wake0C = CA1.copy()
        CA1wake0C[BooleanCons] = 0
        # Liberal boolean filtering of CA1 filtered signal
        BooleanLib = EMGboolean['BooleanLiberal']
        fCA1wake0L = filt_CA1.copy()
        fCA1wake0L[BooleanLib] = 0
        CA1wake0L = CA1.copy()
        CA1wake0L[BooleanLib] = 0

        # Computation of CWT
        CA1cwtWake0cons = signal.cwt(fCA1wake0C, signal.morlet2, widths, w=w)
        CA1cwtWake0lib = signal.cwt(fCA1wake0L, signal.morlet2, widths, w=w)

        # Projection calculation
        absCA1W0Ccwt = np.absolute(CA1cwtWake0cons)
        proj_CA1W0Ccwt = np.sum(absCA1W0Ccwt, axis = 0)/80
        absCA1W0Lcwt = np.absolute(CA1cwtWake0lib)
        proj_CA1W0Lcwt = np.sum(absCA1W0Lcwt, axis = 0)/80

        combined = np.stack([CA1, filt_CA1, proj_CA1W0Ccwt, proj_CA1W0Lcwt], axis = 1)

        sample_rate = 1000.
        t_start = 0.

        # First extraction of SWR peaks, initiation, end and width
        peaks, properties = find_peaks(proj_CA1W0Lcwt, prominence=1, width=20, height=sd6proj_CA1cwt) #2AB detection with 6*SD #AB detection with 8*SD // Audrey's detection=3*SD
        properties["prominences"], properties["widths"]

        # SWR boundaries taken at 70% from peak of intensity. This means that the SWRs with small amplitude will be longer than the big ones.
        results_width = peak_widths(proj_CA1W0Lcwt, peaks, rel_height=0.7)

        # Organise results in numpy array
        peaks2 = peaks.reshape(len(peaks),1)
        npresults_width = np.array(results_width).reshape(4,-1)
        SWR_prop = np.append(peaks2, results_width).reshape(5,len(peaks2)).round()

        # Second extraction of main frequency and power 

        projMaxP_cwtmg = np.max(CA1cwtWake0lib, axis = 0)
        projMaxF_cwtmg = np.argmax(CA1cwtWake0lib, axis = 0) + 120
        nb_SWR = np.arange(0,len(peaks),1)
        data = np.zeros((len(peaks),4))

        for tt in nb_SWR:
            SWR_start = int(SWR_prop[3,tt])
            SWR_stop = int(SWR_prop[4,tt])
            SWR_MaxP = projMaxP_cwtmg[SWR_start:SWR_stop]
            SWR_MaxP = [val.real for val in SWR_MaxP]  # Convert to real numbers
            SWR_MaxF = projMaxF_cwtmg[SWR_start:SWR_stop]
            SWR_MaxF = [val.real for val in SWR_MaxF]  # Convert to real numbers
            data[tt, 0] = max(SWR_MaxF).round()
            data[tt, 1] = max(SWR_MaxP).round()
            data[tt, 2] = round(sum(SWR_MaxF)/len(SWR_MaxF))
            data[tt, 3] = round(sum(SWR_MaxP)/len(SWR_MaxP))

        param_SWR = pd.DataFrame(data, columns = ['Max freq', 'Max int', 'Avg freq', 'Avg int'])
        tSWR_prop = SWR_prop.transpose()
        pd_prop_SWR = pd.DataFrame(tSWR_prop, columns = ['peak time', 'Duration', 'peak amp', 'start time', 'end time'])
        All_SWR = pd.concat([pd_prop_SWR, param_SWR], axis=1)

        SWR_peak = peaks
        SWR_start = SWR_prop[3,:].astype(int)
        SWR_end = SWR_prop[4,:].astype(int)

        # Store the results in All_SWR_prop pd dataframe and save as pkl/csv for post processing.
        filename3 = folder_base / f'SWRpropertiesInitial_sd6_AB.csv'
        All_SWR.to_pickle(filename3)
        All_SWR.to_csv(filename3, sep = ',')

        All_SWR2=All_SWR.copy()
        #All_SWR2.loc[:, :] = np.nan
        filename3 = folder_base / f'SWRproperties_sd6_AB.csv'
        All_SWR2.to_pickle(filename3)
        All_SWR2.to_csv(filename3, sep = ',')


        # Filter parameter :
        f_lowcut = 10.
        f_hicut = 16.
        N = 4
        fs = 1000
        nyq = 0.5 * fs
        Wn = [f_lowcut/nyq,f_hicut/nyq]  # Nyquist frequency fraction

        # Filter creation :
        b, a = signal.butter(N, Wn, 'band')

        filt_S1 = signal.filtfilt(b, a, S1)
        filt_S1wakeremoved = signal.filtfilt(b, a, S1wakeremoved)

        filt_PFC = signal.filtfilt(b, a, PFC)
        filt_PFCwakeremoved = signal.filtfilt(b, a, PFCwakeremoved)

        # Parameter and computation of CWT
        w = 10.
        freq = np.linspace(10, 16, 6)#18)
        widths = w*fs / (2*freq*np.pi)
        PFCNWcwt = signal.cwt(filt_PFCwakeremoved, signal.morlet2, widths, w=w)
        S1NWcwt = signal.cwt(filt_S1wakeremoved, signal.morlet2, widths, w=w)

        # Projection calculation PFC
        absPFCNWcwt = np.absolute(PFCNWcwt)
        proj_PFCNWcwt = np.sum(absPFCNWcwt, axis = 0)/24
        sdproj_PFCcwt = np.std(proj_PFCNWcwt)
        sd5proj_PFCcwt = sdproj_PFCcwt*5
        #sd7proj_PFCcwt = sdproj_PFCcwt*7

        # Projection calculation S1
        absS1NWcwt = np.absolute(S1NWcwt)
        proj_S1NWcwt = np.sum(absS1NWcwt, axis = 0)/24
        sdproj_S1cwt = np.std(proj_S1NWcwt)
        sd5proj_S1cwt = sdproj_S1cwt*5
        #sd7proj_S1cwt = sdproj_S1cwt*7

        #####################################
        ########         PFC         #########
        #####################################
        # Conservative boolean filtering of PFC filtered signal
        BooleanCons = EMGboolean['BooleanConservative']
        fPFCwake0C = filt_PFC.copy()
        fPFCwake0C[BooleanCons] = 0
        PFCwake0C = PFC.copy()
        PFCwake0C[BooleanCons] = 0
        # Liberal boolean filtering of PFC filtered signal
        BooleanLib = EMGboolean['BooleanLiberal']
        fPFCwake0L = filt_PFC.copy()
        fPFCwake0L[BooleanLib] = 0
        PFCwake0L = PFC.copy()
        PFCwake0L[BooleanLib] = 0

        # Computation of CWT
        PFCcwtWake0cons = signal.cwt(fPFCwake0C, signal.morlet2, widths, w=w)
        PFCcwtWake0lib = signal.cwt(fPFCwake0L, signal.morlet2, widths, w=w)

        # Projection calculation
        absPFCW0Ccwt = np.absolute(PFCcwtWake0cons)
        proj_PFCW0Ccwt = np.sum(absPFCW0Ccwt, axis = 0)/24
        absPFCW0Lcwt = np.absolute(PFCcwtWake0lib)
        proj_PFCW0Lcwt = np.sum(absPFCW0Lcwt, axis = 0)/24


        #####################################
        ########         S1         #########
        #####################################
        # Conservative boolean filtering of S1 filtered signal
        BooleanCons = EMGboolean['BooleanConservative']
        fS1wake0C = filt_S1.copy()
        fS1wake0C[BooleanCons] = 0
        S1wake0C = S1.copy()
        S1wake0C[BooleanCons] = 0
        # Liberal boolean filtering of S1 filtered signal
        BooleanLib = EMGboolean['BooleanLiberal']
        fS1wake0L = filt_S1.copy()
        fS1wake0L[BooleanLib] = 0
        S1wake0L = S1.copy()
        S1wake0L[BooleanLib] = 0

        # Computation of CWT
        S1cwtWake0cons = signal.cwt(fS1wake0C, signal.morlet2, widths, w=w)
        S1cwtWake0lib = signal.cwt(fS1wake0L, signal.morlet2, widths, w=w)

        # Projection calculation
        absS1W0Ccwt = np.absolute(S1cwtWake0cons)
        proj_S1W0Ccwt = np.sum(absS1W0Ccwt, axis = 0)/24
        absS1W0Lcwt = np.absolute(S1cwtWake0lib)
        proj_S1W0Lcwt = np.sum(absS1W0Lcwt, axis = 0)/24

        # 7 sd threshold
        peaks, properties = find_peaks(proj_PFCW0Lcwt, width=200, height=sd5proj_PFCcwt) #AB detection first sd=7
        properties["prominences"], properties["widths"]

        # Spindles boundaries taken at 70% from peak of intensity. This means that the spindles with small amplitude will be longer than the big ones.
        results_width = peak_widths(proj_PFCW0Lcwt, peaks, rel_height=0.7)

        # Organise results in numpy array
        peaks2 = peaks.reshape(len(peaks),1)
        npresults_width = np.array(results_width).reshape(4,-1)
        Spindle_prop = np.append(peaks2, results_width).reshape(5,len(peaks2)).round()
        
        projMaxP_cwtmg = np.max(PFCcwtWake0lib, axis = 0)
        projMaxF_cwtmg = np.argmax(PFCcwtWake0lib, axis = 0)+ 10 #/2 + 8
        projMaxP_cwtmg.shape

        nb_Spindles = np.arange(0,len(peaks),1)
        data = np.zeros((len(peaks),4))

        for tt in nb_Spindles:
            Spindle_start = int(Spindle_prop[3,tt])
            Spindle_stop = int(Spindle_prop[4,tt])
            Spindle_MaxP = projMaxP_cwtmg[Spindle_start:Spindle_stop]
            Spindle_MaxF = projMaxF_cwtmg[Spindle_start:Spindle_stop]
            Spindle_MaxP = [val.real for val in Spindle_MaxP]  # Convert to real numbers
            Spindle_MaxF = [val.real for val in Spindle_MaxF]  # Convert to real numbers
            data[tt, 0] = max(Spindle_MaxF).round()
            data[tt, 1] = max(Spindle_MaxP).round()
            data[tt, 2] = round(sum(Spindle_MaxF)/len(Spindle_MaxF))
            data[tt, 3] = round(sum(Spindle_MaxP)/len(Spindle_MaxP))

        param_Spindle = pd.DataFrame(data, columns = ['Max freq', 'Max int', 'Avg freq', 'Avg int'])
        tSpindle_prop = Spindle_prop.transpose()
        pd_prop_Spindle = pd.DataFrame(tSpindle_prop, columns = ['peak time', 'Duration', 'peak amp', 'start time', 'end time'])
        All_Spindle = pd.concat([pd_prop_Spindle, param_Spindle], axis=1)

        nb_spindle = All_Spindle.shape[0]
        listtodrop = []
        for tt in range(nb_spindle-1):
            if(All_Spindle['end time'][tt]>All_Spindle['start time'][tt + 1]):
                if(All_Spindle['Duration'][tt]<All_Spindle['Duration'][tt + 1]):
                    if(All_Spindle['start time'][tt]<All_Spindle['start time'][tt + 1]):
                        All_Spindle['start time'][tt+1] = All_Spindle['start time'][tt]
                        listtodrop.append(tt)
                    else:
                        listtodrop.append(tt)
                if(All_Spindle['Duration'][tt]>All_Spindle['Duration'][tt + 1]):
                    if(All_Spindle['end time'][tt]<All_Spindle['end time'][tt + 1]):
                        All_Spindle['end time'][tt] = All_Spindle['end time'][tt + 1]
                        listtodrop.append(tt+1)
                    else:
                        listtodrop.append(tt+1)

        for tt in range(nb_spindle-1):
            if((All_Spindle['start time'][tt + 1] - All_Spindle['end time'][tt])<200):
                if((All_Spindle['Duration'][tt]+300)<All_Spindle['Duration'][tt + 1]):
                    All_Spindle['start time'][tt + 1] = All_Spindle['start time'][tt]
                    listtodrop.append(tt)
                if((All_Spindle['Duration'][tt+1]+300)<All_Spindle['Duration'][tt]):
                    All_Spindle['end time'][tt] = All_Spindle['start time'][tt + 1]
                    listtodrop.append(tt+1)

        for tt in range(nb_spindle):
            All_Spindle['Duration'][tt]=All_Spindle['end time'][tt]-All_Spindle['start time'][tt]
        All_Spindle = All_Spindle.drop(listtodrop) 

        filename3 = folder_base / f'Spindlesproperties_PFCInitial_sd5_AB.csv'
        All_Spindle.to_pickle(filename3)
        All_Spindle.to_csv(filename3, sep = ',')

        filename3 = folder_base / f'Spindlesproperties_PFC_sd5_AB.csv'
        All_Spindle.to_pickle(filename3)
        All_Spindle.to_csv(filename3, sep = ',')
        
            
        # 5 sd threshold
        peaks, properties = find_peaks(proj_S1W0Lcwt, width=200, height=sd5proj_S1cwt)
        properties["prominences"], properties["widths"]

        # Spindles boundaries taken at 70% from peak of intensity. This means that the spindles with small amplitude will be longer than the big ones.
        results_width = peak_widths(proj_S1W0Lcwt, peaks, rel_height=0.7)

        # Organise results in numpy array
        peaks2 = peaks.reshape(len(peaks),1)
        npresults_width = np.array(results_width).reshape(4,-1)
        Spindle_prop = np.append(peaks2, results_width).reshape(5,len(peaks2)).round()

        projMaxP_cwtmg = np.max(S1cwtWake0lib, axis = 0)
        projMaxF_cwtmg = np.argmax(S1cwtWake0lib, axis = 0)+ 10 #/2 + 8

        nb_Spindles = np.arange(0,len(peaks),1)
        data = np.zeros((len(peaks),4))

        for tt in nb_Spindles:
            Spindle_start = int(Spindle_prop[3,tt])
            Spindle_stop = int(Spindle_prop[4,tt])
            Spindle_MaxP = projMaxP_cwtmg[Spindle_start:Spindle_stop]
            Spindle_MaxF = projMaxF_cwtmg[Spindle_start:Spindle_stop]
            Spindle_MaxP = [val.real for val in Spindle_MaxP]  # Convert to real numbers
            Spindle_MaxF = [val.real for val in Spindle_MaxF]  # Convert to real numbers
            data[tt, 0] = max(Spindle_MaxF).round()
            data[tt, 1] = max(Spindle_MaxP).round()
            data[tt, 2] = round(sum(Spindle_MaxF)/len(Spindle_MaxF))
            data[tt, 3] = round(sum(Spindle_MaxP)/len(Spindle_MaxP))

        param_Spindle = pd.DataFrame(data, columns = ['Max freq', 'Max int', 'Avg freq', 'Avg int'])
        tSpindle_prop = Spindle_prop.transpose()
        pd_prop_Spindle = pd.DataFrame(tSpindle_prop, columns = ['peak time', 'Duration', 'peak amp', 'start time', 'end time'])
        All_Spindle = pd.concat([pd_prop_Spindle, param_Spindle], axis=1)

        nb_spindle = All_Spindle.shape[0]
        listtodrop = []
        for tt in range(nb_spindle-1):
            if(All_Spindle['end time'][tt]>All_Spindle['start time'][tt + 1]):
                if(All_Spindle['Duration'][tt]<All_Spindle['Duration'][tt + 1]):
                    if(All_Spindle['start time'][tt]<All_Spindle['start time'][tt + 1]):
                        All_Spindle['start time'][tt+1] = All_Spindle['start time'][tt]
                        listtodrop.append(tt)
                    else:
                        listtodrop.append(tt)
                if(All_Spindle['Duration'][tt]>All_Spindle['Duration'][tt + 1]):
                    if(All_Spindle['end time'][tt]<All_Spindle['end time'][tt + 1]):
                        All_Spindle['end time'][tt] = All_Spindle['end time'][tt + 1]
                        listtodrop.append(tt+1)
                    else:
                        listtodrop.append(tt+1)

        for tt in range(nb_spindle-1):
            if((All_Spindle['start time'][tt + 1] - All_Spindle['end time'][tt])<200):
                if((All_Spindle['Duration'][tt]+300)<All_Spindle['Duration'][tt + 1]):
                    All_Spindle['start time'][tt + 1] = All_Spindle['start time'][tt]
                    listtodrop.append(tt)
                if((All_Spindle['Duration'][tt+1]+300)<All_Spindle['Duration'][tt]):
                    All_Spindle['end time'][tt] = All_Spindle['start time'][tt + 1]
                    listtodrop.append(tt+1)

        for tt in range(nb_spindle):
            All_Spindle['Duration'][tt]=All_Spindle['end time'][tt]-All_Spindle['start time'][tt]
        All_Spindle = All_Spindle.drop(listtodrop) 
        All_Spindle.shape[0]

        filename3 = folder_base / f'Spindlesproperties_S1Initial_sd5_AB.csv'
        All_Spindle.to_pickle(filename3)
        All_Spindle.to_csv(filename3, sep = ',')

        filename3 = folder_base / f'Spindlesproperties_S1_sd5_AB.csv'
        All_Spindle.to_pickle(filename3)
        All_Spindle.to_csv(filename3, sep = ',')