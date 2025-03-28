##### Called by convert_dat_to_pkl.sh  #####

##################################
        # SETTING-UP #
##################################

import numpy as np
from pathlib import Path
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, Cursor
from scipy import signal
import quantities as pq
import os
import json 
import time

st = time.time()

##################################
        # LOAD RECORDINGS #
##################################
dpath='/mnt/data/AurelieB_other/'
folder_base = Path(dpath) 

# Load LFPs data 
if Path(f'{folder_base}/RawDataChannelExtractedDS.npy').exists():
    print('RawDataChannelExtractedDS.npy file')
    LFPfile = Path(f'{folder_base}/RawDataChannelExtractedDS.npy')
    DataRec = np.load(LFPfile, mmap_mode= 'r')
    samplerate = 1000 
    numchannel = np.shape(DataRec)[1]
    rec_ch_list = np.arange(0,numchannel,1)
    LFPs_df=pd.DataFrame(DataRec, columns=rec_ch_list)

elif Path(f'{folder_base}/continuous.dat').exists():
    LFPfile = Path(f'{folder_base}/continuous.dat')
    print('continuous.dat file')
    LFPfile = Path(f'{folder_base}/continuous.dat')
    DataRec = np.fromfile(LFPfile, dtype="int16")
    filepath = Path(os.path.join(folder_base, f'structure.oebin'))
    with open(filepath) as f:
        metadata = json.load(f)
    samplerate = metadata['continuous'][0]['sample_rate']  
    numchannel = metadata['continuous'][0]['num_channels'] 
    rec_ch_list = np.array([int(''.join(c for c in metadata['continuous'][0]['channels'][x]['channel_name'] if c.isdigit()))-1 for x in range(0, len(metadata['continuous'][0]['channels']))])
    DataRec = DataRec.reshape(-1,numchannel)
    LFPs_df=pd.DataFrame(DataRec, columns=rec_ch_list)
    """
    # Load LFPs timestamps 
    file_pathTS = Path(os.path.join(folder_base, f'timeStamps.npy'))
    print('LFPs timestamps file found')
    LFPtimestamps = np.load(file_pathTS) 
    """

print('sample rate =', samplerate, 'Hz')
print(numchannel, 'channels recorded')
print(round(LFPs_df.shape[0]/samplerate/60), 'min of recording')

##################################
    # DOWNSAMPLE RECORDINGS #
##################################

if samplerate > 1000:
    new_sampling_rate = 1000 # Hz
    Nmber_points = int(np.shape(LFPs_df)[0] * new_sampling_rate / samplerate)
    LFPs_df_DS = pd.DataFrame(signal.resample(LFPs_df, Nmber_points, axis = 0), columns=LFPs_df.columns.values)
    #LFPtimestampsDS = LFPtimestamps[::int(samplerate/new_sampling_rate)][:-1]
    LFPs_df = LFPs_df_DS
    #LFPtimestamps = LFPtimestampsDS
    samplerate = new_sampling_rate

##################################
        # SAVE RECORDINGS #
##################################
    
LFPs_df.to_pickle(f'{LFPfile.parent}/DataFrame_rawdataDS.pkl')
#np.save(f'{file_pathTS.parent}/timeStampsDS.npy', LFPtimestamps)

elapsed_time = time.time() - st
print('Execution time:', time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))