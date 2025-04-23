#######################################################################################
                            # Define Directory #
#######################################################################################

dir = "//10.69.168.1/crnldata/waking/audrey_hay/L1imaging/Analysed2025_AB/L1NDNF_mice/BlackLines/"

#######################################################################################
                                # Load Signals #
#######################################################################################

# Load packages
import os
from scipy import signal
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, Cursor
from scipy import fftpack
import pandas as pd
from pathlib import Path
from IPython.display import display
from datetime import datetime
import shutil
from scipy.signal import find_peaks
from scipy.signal import chirp, find_peaks, peak_widths

#Load LFP coordinates 
Channels = f'{os.getcwd()}\python\_LFP_coordinates_of_all_mice.csv'
all_LFPcoordinates = pd.read_csv(Channels, index_col=0)

for dpath in Path(dir).glob('**/DataFrame_rawdataDS.pkl'):
    #Load signals
    folder_base = Path(dpath).parent
    print(folder_base)
    LFPfile = Path(f'{folder_base}\DataFrame_rawdataDS.pkl')
    LFPs_df = pd.read_pickle(LFPfile)
    samplerate = 1000 
    numchannel = LFPs_df.shape[1]
    rec_ch_list = LFPs_df.columns.values
    # Load LFPs timestamps 
    for file_pathTS in folder_base.parent.parent.glob('**/continuous/*/timeStampsDS.npy'):
        print('LFPs timestamps file found')
        LFPtimestamps = np.load(file_pathTS)  
    print(round(LFPs_df.shape[0]/samplerate/60), 'min of recording')

    # Identify mouse & choose threshold for detection
    mouse = []
    pos_mice = []
    for mouse_name in all_LFPcoordinates.index:
        if mouse_name in LFPfile.__str__():
            mouse.append(mouse_name)
            pos_mice.append(LFPfile.__str__().find(mouse_name)) 
    mouse = [x for _, x in sorted(zip(pos_mice, mouse))] # sort mouse in the same order as they appear in the path
    mouse=mouse[0]   

    """
    SWRfactor={'BlackLines': 4, 'BlueLines':4, 'GreenDots': 5, 'GreenLines':4, 'PurpleSquare':4,  'RedLines':4, 'ThreeColDots':4,'ThreeBlueCrosses':4}
    SpdlfactorPFC={'BlackLines': 6, 'BlueLines':4, 'GreenDots': 3, 'GreenLines':4, 'PurpleSquare':3,  'RedLines':4, 'ThreeColDots':4,'ThreeBlueCrosses':4}
    SpdlfactorS1={'BlackLines': 6, 'BlueLines':4, 'GreenDots': 3, 'GreenLines':4, 'PurpleSquare':3,  'RedLines':4, 'ThreeColDots':4,'ThreeBlueCrosses':4}
    """
    
    SWRfactor={'BlackLines': 4, 'BlueLines':4, 'GreenDots': 5, 'GreenLines':4, 'PurpleSquare':4,  'RedLines':4, 'ThreeColDots':4,'ThreeBlueCrosses':4}
    SpdlfactorPFC={'BlackLines': 5.5, 'BlueLines':4, 'GreenDots': 3, 'GreenLines':4, 'PurpleSquare':3,  'RedLines':4, 'ThreeColDots':4,'ThreeBlueCrosses':4}
    SpdlfactorS1={'BlackLines': 3, 'BlueLines':4, 'GreenDots': 3, 'GreenLines':4, 'PurpleSquare':3,  'RedLines':4, 'ThreeColDots':4,'ThreeBlueCrosses':4}

    SWRfactor=SWRfactor[mouse]
    SpdlfactorPFC=SpdlfactorPFC[mouse]
    SpdlfactorS1=SpdlfactorS1[mouse]

    # Identify electrodes & create differential LFPs
    all_LFPcoordinates = all_LFPcoordinates.astype(str)
    for region in all_LFPcoordinates.loc[mouse].index:
        locals()[region] = []
        locals()[f'{region}_ch'] = []
    ID=0
    rec_ch_list_mouse = [value for value in rec_ch_list if 0+(ID*32) <= value <= 31+(ID*32)]
    for rec_ch in rec_ch_list_mouse:
        for idx, LFPcoord_str in enumerate(all_LFPcoordinates.loc[mouse]):
            region = all_LFPcoordinates.loc[mouse].index[idx]
            if LFPcoord_str != 'nan':
                LFPcoord = LFPcoord_str.split('_')[:2] # only take into account the 2 first of electrode of that region 
                num_ch = np.where(str(rec_ch-(ID*32)) == np.array(LFPcoord))[0]
                if len(num_ch) > 0:
                    region = all_LFPcoordinates.loc[mouse].index[idx]
                    LFP = locals()[region]
                    LFP = LFP-np.array(LFPs_df[(rec_ch)]) if len(LFP) > 0 else np.array(LFPs_df[(rec_ch)])
                    locals()[region] = LFP
                    locals()[f'{region}_ch'].append(rec_ch)
                    break
                continue    
    for region in all_LFPcoordinates.loc[mouse].index:
        LFP = locals()[region]
        LFP_ch = locals()[f'{region}_ch']

    # Load Sleep Scoring    
    ScoringFile= folder_base / f'Sleep_Scoring_6Stages_5sEpoch.csv'
    SleepScored = pd.read_csv(ScoringFile)
    SleepScored['label']= SleepScored['label'].str.extract(r'(\d+)', expand=False)
    SleepScoredTS=np.array(SleepScored['label'])
    scale_factor=samplerate/0.2  #cause scoring was done in 5 seconds bin, ie 0.2 Hz              
    SleepScoredTS_upscaled0 = np.repeat(SleepScoredTS, scale_factor, axis=0)
    SleepScoredTS_upscaled0=SleepScoredTS_upscaled0.astype(int)
    SleepScoredTS_upscaled=SleepScoredTS_upscaled0.copy()

    SleepScoredTS_upscaled[SleepScoredTS_upscaled==2]=1 # Transform QW in AW

    if len(CA1) > len(SleepScoredTS_upscaled):
        SleepScoredTS_upscaled = np.pad(SleepScoredTS_upscaled, (0, (len(CA1) - len(SleepScoredTS_upscaled))), constant_values=SleepScoredTS_upscaled[-1])
    elif len(CA1) < len(SleepScoredTS_upscaled):
        SleepScoredTS_upscaled=SleepScoredTS_upscaled[:len(CA1)]

    WakeBool = SleepScoredTS_upscaled==1 

    CA1wakeremoved = CA1.copy()
    PFCwakeremoved = PFC.copy()
    S1wakeremoved = S1.copy()

    CA1wakeremoved=CA1wakeremoved[~WakeBool]  #CA1wakeremoved[WakeBool] = np.nan
    PFCwakeremoved=PFCwakeremoved[~WakeBool]
    S1wakeremoved=S1wakeremoved[~WakeBool]

    ###########################################
            # SWR in CA1: 120-200 Hz #
    ###########################################

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
    Fsdproj_CA1cwt = sdproj_CA1cwt*SWRfactor

    fCA1wake0L = filt_CA1.copy()
    fCA1wake0L[WakeBool] = 0
    CA1wake0L = CA1.copy()
    CA1wake0L[WakeBool] = 0

    # Computation of CWT
    CA1cwtWake0lib = signal.cwt(fCA1wake0L, signal.morlet2, widths, w=w)

    # Projection calculation
    absCA1W0Lcwt = np.absolute(CA1cwtWake0lib)
    proj_CA1W0Lcwt = np.sum(absCA1W0Lcwt, axis = 0)/80

    # First extraction of SWR peaks, initiation, end and width
    peaks, properties = find_peaks(proj_CA1W0Lcwt, prominence=1, width=20, height=Fsdproj_CA1cwt) #2AB detection with 6*SD #AB detection with 8*SD // Audrey's detection=3*SD
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
    filename3 = folder_base / f'SWR_detection_initial.xlsx'
    All_SWR.to_excel(filename3)
    filename = folder_base / f'SWR_detection.xlsx'
    All_SWR.to_excel(filename)

    print(len(All_SWR), 'SWR detected in CA1')

    ########################################
                # Spdl: 10-16 Hz #
    ########################################

    # Filter parameter :
    f_lowcut = 10.
    f_hicut = 16.
    N = 4
    fs = 1000
    nyq = 0.5 * fs
    Wn = [f_lowcut/nyq,f_hicut/nyq]  # Nyquist frequency fraction

    # Filter creation :
    b, a = signal.butter(N, Wn, 'band')

    # Parameter and computation of CWT
    w = 10.
    freq = np.linspace(10, 16, 6)#18)
    widths = w*fs / (2*freq*np.pi)

    ##################################
        ##         PFC         ##
    
    filt_PFC = signal.filtfilt(b, a, PFC)
    filt_PFCwakeremoved = signal.filtfilt(b, a, PFCwakeremoved)
    PFCNWcwt = signal.cwt(filt_PFCwakeremoved, signal.morlet2, widths, w=w)

    # Projection calculation PFC
    absPFCNWcwt = np.absolute(PFCNWcwt)
    proj_PFCNWcwt = np.sum(absPFCNWcwt, axis = 0)/24
    sdproj_PFCcwt = np.std(proj_PFCNWcwt)
    #sd5proj_PFCcwt = sdproj_PFCcwt*5
    Fsdproj_PFCcwt = sdproj_PFCcwt*SpdlfactorPFC

    # Conservative boolean filtering of PFC filtered signal
    fPFCwake0L = filt_PFC.copy()
    fPFCwake0L[WakeBool] = 0
    PFCwake0L = PFC.copy()
    PFCwake0L[WakeBool] = 0

    # Computation of CWT
    PFCcwtWake0lib = signal.cwt(fPFCwake0L, signal.morlet2, widths, w=w)

    # Projection calculation
    absPFCW0Lcwt = np.absolute(PFCcwtWake0lib)
    proj_PFCW0Lcwt = np.sum(absPFCW0Lcwt, axis = 0)/24

    # SD threshold
    peaks, properties = find_peaks(proj_PFCW0Lcwt, width=200, height=Fsdproj_PFCcwt) #AB detection second sd=5 #AB detection first sd=7
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

    filename = folder_base / f'SpindlesPFC_detection_intial.xlsx'
    All_Spindle.to_excel(filename)
    filename = folder_base / f'SpindlesPFC_detection.xlsx'
    All_Spindle.to_excel(filename)

    print(len(All_Spindle), 'Spdl detected in PFC')

    #####################################
            ##         S1         ##

    filt_S1 = signal.filtfilt(b, a, S1)
    filt_S1wakeremoved = signal.filtfilt(b, a, S1wakeremoved)
    S1NWcwt = signal.cwt(filt_S1wakeremoved, signal.morlet2, widths, w=w)

    # Projection calculation S1
    absS1NWcwt = np.absolute(S1NWcwt)
    proj_S1NWcwt = np.sum(absS1NWcwt, axis = 0)/24
    sdproj_S1cwt = np.std(proj_S1NWcwt)
    #sd5proj_S1cwt = sdproj_S1cwt*5
    Fsdproj_S1cwt = sdproj_S1cwt*SpdlfactorS1

    # Liberal boolean filtering of S1 filtered signal
    fS1wake0L = filt_S1.copy()
    fS1wake0L[WakeBool] = 0
    S1wake0L = S1.copy()
    S1wake0L[WakeBool] = 0

    # Computation of CWT
    S1cwtWake0lib = signal.cwt(fS1wake0L, signal.morlet2, widths, w=w)

    # Projection calculation
    absS1W0Lcwt = np.absolute(S1cwtWake0lib)
    proj_S1W0Lcwt = np.sum(absS1W0Lcwt, axis = 0)/24
        
    # Sd threshold
    peaks, properties = find_peaks(proj_S1W0Lcwt, width=200, height=Fsdproj_S1cwt) #AB detection second sd=5 #AB detection first sd=7
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

    filename = folder_base / f'SpindlesS1_detection_intial.xlsx'
    All_Spindle.to_excel(filename)
    filename = folder_base / f'SpindlesS1_detection.xlsx'
    All_Spindle.to_excel(filename)
    
    print(len(All_Spindle), 'Spdl detected in S1')
