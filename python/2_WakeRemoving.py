# Extracting wake and filtering it out of data

#######################################################################################
                                # Load packages #
#######################################################################################

import scipy
from scipy import signal
from scipy import interpolate
from scipy import fftpack
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, Cursor
import pandas as pd
import quantities as pq
import neo
from pathlib import Path
import xarray as xr
from IPython.display import display
from ipyfilechooser import FileChooser
import os
from ephyviewer import mkQApp, MainViewer, TraceViewer

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

MiceList=['BlackLinesOK', 'BlueLinesOK', 'GreenDotsOK', 'GreenLinesOK', 'Purple', 'RedLinesOK','ThreeColDotsOK', 'ThreeBlueCrossesOK']
MiceList=['ThreeColDotsOK']

dpath0 = "//10.69.168.1/crnldata/waking/audrey_hay/L1imaging/AnalysedMarch2023/Gaelle/CGP/"

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
        folder_base = session_path / f'OpenEphys/'
        print(folder_base)
                
        filename2 = folder_base / f'RawDataChannelExtractedDS.npy'
        All = np.load(filename2, mmap_mode= 'r')

        Channels = '//10.69.168.1/crnldata/waking/audrey_hay/L1imaging/AnalysedMarch2023/LFPChannels_perMice.xlsx' 
        allchannels = pd.read_excel(Channels)

        EMGch=int(allchannels[mice][3])
        EMG  =  All[:, EMGch]

        # Filter parameter :
        f_lowcut = 200.
        f_hicut = 400.
        N = 4
        fs = 1000
        nyq = 0.5 * fs
        Wn = [f_lowcut/nyq,f_hicut/nyq]  # Nyquist frequency fraction

        # Filter creation :
        b, a = signal.butter(N, Wn, 'band')
        filt_EMG = signal.filtfilt(b, a, EMG)

        # Parameter and computation of CWT
        w = 4.
        freq = np.linspace(200, 400, 50)
        widths = w*fs / (2*freq*np.pi)
        EMGcwt = signal.cwt(EMG, signal.morlet2, widths, w=w)

        # Projection calculation
        absEMGcwt = np.absolute(EMGcwt)
        proj_EMGcwt = np.sum(absEMGcwt, axis = 0)/50
        mproj_EMGcwt = np.mean(proj_EMGcwt)
        sdproj_EMGcwt = np.std(proj_EMGcwt)
        sd3proj_EMGcwt = mproj_EMGcwt + sdproj_EMGcwt*3
        sd05proj_EMGcwt = mproj_EMGcwt + sdproj_EMGcwt*0.5
        sd1proj_EMGcwt = mproj_EMGcwt + sdproj_EMGcwt

        # Assigning values wake (1, 2) and sleep (0)
        numpnts = EMG.size
        EMGstatusRaw = np.zeros(numpnts)
        for ind in range(numpnts):
            if proj_EMGcwt[ind]<sd1proj_EMGcwt:
                EMGstatusRaw[ind] = 0
            elif proj_EMGcwt[ind]>sd3proj_EMGcwt:
                EMGstatusRaw[ind] = 2
            else:
                EMGstatusRaw[ind] = 1

        # Expanding borders for wake (1, 2) and sleep (0) to ±1 s around detected muscular activity
        EMGstatusRaw2 = np.zeros(numpnts)
        for ind in range(numpnts):
            if EMGstatusRaw[ind]>1:
                EMGstatusRaw2[ind-1000:ind+1000] = 2
            elif EMGstatusRaw[ind]==1:
                for ind2 in range(ind-1000, ind+1000):
                    if ind2==numpnts:
                        break
                    elif EMGstatusRaw2[ind2]<2:
                        EMGstatusRaw2[ind2] = 1


        EMGStatusBoolLib = (EMGstatusRaw2>1)
        EMGStatusBoolCons = (EMGstatusRaw2>0)

        # Save files
        LFP = All[:,:]
        LFPwake0 = LFP.copy()
        LFPwake0[EMGStatusBoolLib] = 0
        filename = folder_base/ f'LFPwake0_AB.npy'
        np.save(filename, LFPwake0)

        LFPwakeremoved = LFP.copy()
        LFPwakeremoved = LFPwakeremoved[~EMGStatusBoolLib, :]
        filename = folder_base/ f'LFPwakeremoved_AB.npy'
        np.save(filename, LFPwakeremoved)
        data = {
            'EMGstatus': EMGstatusRaw2,
            'BooleanLiberal' : EMGStatusBoolLib,
            'BooleanConservative' : EMGStatusBoolCons
        }
        WakeFrame = pd.DataFrame(data, columns=['EMGstatus', 'BooleanLiberal', 'BooleanConservative'])
        filename = folder_base/ f'EMGframeBoolean_AB.pkl'
        WakeFrame.to_pickle(filename)