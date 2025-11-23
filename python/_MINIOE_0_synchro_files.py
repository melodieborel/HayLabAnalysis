
##################################
        # SETTING-UP #
##################################

import numpy as np
import csv
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import xarray as xr

##################################
        # DEFINE PATH #
##################################

local = True
if local:
    DIR= Path("//10.69.168.1/crnldata/forgetting/Aurelie/MiniscopeOE_data/L2_3_mice/")
else: 
    DIR= Path("/crnldata/forgetting/Aurelie/MiniscopeOE_data/L2_3_mice/")

##################################
    # CREATE SYNCHRO FILE #
##################################

for dpath in DIR.glob('**/*V4_Miniscope*/'):

    folder = Path(dpath.parent)
    SynchroFile = next(folder.rglob('SynchroFile.xlsx'), None)
    
# REMOVE CONDTION AND REDO CAUSE SOME OPENEPHYS FOLDERS HAVE BEEN CREATED AFTER
    if not SynchroFile: # process only if there are no SynchroFile already 
# REMOVE CONDTION AND REDO CAUSE SOME OPENEPHYS FOLDERS HAVE BEEN CREATED AFTER
        for dir_path in folder.glob('**/*V4_Miniscope*/'):

            print("####################################################################################################################################################")
            print(folder)
            print("####################################################################################################################################################")

            file_path= f'{dir_path}/timeStamps.csv'
            stamps_miniscope = pd.read_csv(file_path)  

            stamps_miniscope_time = stamps_miniscope['Time Stamp (ms)']
            delay_stamps = []
            number_frames = stamps_miniscope['Time Stamp (ms)'].count()
            for i in range(number_frames -1):
                delay_stamps.append(stamps_miniscope_time[i+1] - stamps_miniscope_time[i])


            freq_acq = round(1000/(sum(delay_stamps)/len(delay_stamps)))
            print("The calculated frame rate is : {} Hz".format(freq_acq) )

            delaystampMini = (1000/freq_acq)*1.5
            print('Threshold dropped frame=',delaystampMini)

            delay_stamps = []
            dropped_frames = []
            number_frames = stamps_miniscope['Time Stamp (ms)'].count()
            for i in range(number_frames -1):
                delay_stamps.append(stamps_miniscope_time[i+1] - stamps_miniscope_time[i])
                if delay_stamps[i] > delaystampMini:
                    dropped_frames.append(i)

            print("{} frame(s) were dropped : {}".format(len(dropped_frames),dropped_frames))
        
            if any(folder.parent.rglob('OpenEphys/*')):
                
                try: # in case no sync file            

                    recordings = ['TTL_full_words', 'TTL_sample_numbers', 'TTL_states', 'TTL_timestamps']
                    Allstamps = {}
                    for file_path in folder.parent.glob('**/TTL*.npy'):
                        filename = file_path.stem
                        print(f"Loading {filename}")
                        Allstamps[filename] = np.load(file_path)

                    datalen = max(len(arr) for arr in Allstamps.values())
                    fullwords      = Allstamps.get('TTL_full_words', np.full(datalen, np.nan))
                    timestamps     = Allstamps.get('TTL_timestamps', np.full(datalen, np.nan))
                    channelstates  = Allstamps.get('TTL_states', np.full(datalen, np.nan))
                    channels       = Allstamps.get('TTL_sample_numbers', np.full(datalen, np.nan))

                    OE_stamps_miniscope = []
                    for i in range(datalen):
                        if channelstates[i] == 8: # TTL 8 is Miniscope
                            OE_stamps_miniscope.append(timestamps[i])

                    for file_path in folder.parent.glob('**/sync_messages.txt'):
                        with open(file_path, "r") as f:
                            lines = f.readlines()
                        # Loop through all lines, ignoring header.
                            for l in lines[1:]:
                                initial_OE_start = int(l.split()[-1] )
                                acqFreqOE = int(l.split()[-3] )

                    # transform in pd series for easier manipulation
                    B = pd.Series(OE_stamps_miniscope)

                    # normalise to ms
                    OE_stamps_miniscope_n = B - (initial_OE_start/acqFreqOE)

                    # _n is from 0
                    acquisition_mini_start_n = OE_stamps_miniscope_n[0]
                    # _a is from acquisition time software
                    acquisition_mini_start_a = B[0]

                    outSumm = pd.Series([acquisition_mini_start_n, acquisition_mini_start_a, freq_acq, dropped_frames],
                                index=['Miniscope start from 0', 'Miniscope start from Acq time soft', 'mini acq freq', 'dropped frames'])
                    print(outSumm)

                    filenameOut = folder / f'SynchroFile.xlsx'
                    writer = pd.ExcelWriter(filenameOut)
                    outSumm.to_excel(writer)
                    writer.close()

                except:
                    print('NO sync_messages.txt FILE SAVED. Skipping...')
                    acquisition_mini_start_n = np.nan
                    acquisition_mini_start_a = np.nan

                    outSumm = pd.Series([acquisition_mini_start_n, acquisition_mini_start_a, freq_acq, dropped_frames],
                                index=['Miniscope start from 0', 'Miniscope start from Acq time soft', 'mini acq freq', 'dropped frames'])
                    print(outSumm)

                    filenameOut = folder / f'SynchroFile.xlsx'
                    writer = pd.ExcelWriter(filenameOut)
                    outSumm.to_excel(writer)
                    writer.close()
                    continue
            else:
                print('NO OpenEphys FOLDER. ...')
                acquisition_mini_start_n = np.nan
                acquisition_mini_start_a = np.nan

                outSumm = pd.Series([acquisition_mini_start_n, acquisition_mini_start_a, freq_acq, dropped_frames],
                            index=['Miniscope start from 0', 'Miniscope start from Acq time soft', 'mini acq freq', 'dropped frames'])
                print(outSumm)

                filenameOut = folder / f'SynchroFile.xlsx'
                writer = pd.ExcelWriter(filenameOut)
                outSumm.to_excel(writer)
                writer.close()
                
                continue
        