
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
DIR= Path("//10.69.168.1/crnldata/waking/audrey_hay/L1imaging/Analysed2025_AB/L1NDNF_mice/RedLines/")

##################################
    # CREATE SYNCHRO FILE #
##################################

for dpath in DIR.glob('**/*V4_Miniscope*/'):

    folder = Path(dpath.parent)
    file = next(folder.rglob('SynchroFileCorrect.xlsx'), None)

    if not file: # process only if there are no SynchroFile already 
        
        try: # in case no sync file

            print("####################################################################################################################################################")
            print(folder)
            print("####################################################################################################################################################")

            for dir_path in folder.glob('**/*V4_Miniscope*/'):
                file_path= f'{dir_path}/timeStamps.csv'
                stamps_miniscope = pd.read_csv(file_path)

            for file_path in folder.glob('**/*.npy'):
                subfolder = file_path.parents[0].stem
                if subfolder == 'TTL':
                    file = file_path.stem
                    print(file_path)
                    stamps_OEmini = np.load(file_path)
                    datalen = len(stamps_OEmini)
                    coords = {
                        'recordings' : np.array(['full_words', 'timestamps', 'channel_states', 'channels']),
                        'duration_rec' : np.arange(datalen)
                    }
                    Allstamps = xr.DataArray(coords=coords, dims=['recordings', 'duration_rec'])

            for file_path in folder.glob('**/*.npy'):
                subfolder = file_path.parents[0].stem
                if subfolder == 'TTL':
                    file = file_path.stem
                    stamps_OEmini = np.load(file_path)
                    Allstamps.loc[file,:] = stamps_OEmini


            time = range(datalen)
            fullwords = Allstamps.loc['full_words',:].values
            timestamps = Allstamps.loc['timestamps',:].values
            channelstates = Allstamps.loc['channel_states',:]
            channels = Allstamps.loc['channels',:]

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


            OE_stamps_miniscope = []
            OE_stamps_laser = []
            for i in range(datalen):
                if channels[i] == 2:
                    OE_stamps_miniscope.append(timestamps[i])
                elif channels[i] == 1:
                    OE_stamps_laser.append(timestamps[i])


            A = []
            for file_path in folder.glob('**/sync_messages.txt'):
                with open(file_path, "r") as f:
                    lines = f.readlines()
                # Loop through all lines, ignoring header.
                    for l in lines[1:]:
                        A.append(l.split()[-1])# take last element to list (i.e. the process name)

            # remove acquisition frequency that is normally always 25 kHz
            initial_OE_start = int(' '.join([x.split('@')[0] for x in A]))
            acqFreqOE = int(' '.join([x[:-2].split('@')[1] for x in A]))

            # transform in pd series for easier manipulation
            B = pd.Series(OE_stamps_miniscope)
            C = pd.Series(OE_stamps_laser)

            # normalise to ms
            OE_stamps_miniscope_n = (B - initial_OE_start)/acqFreqOE
            OE_stamps_laser_inter = (C - initial_OE_start)/acqFreqOE

            # take only the middle of the laser pulse
            OE_stamps_laser_n = [] 
            for i in range(len(OE_stamps_laser_n) - 1):
                if (OE_stamps_laser_n[i+1] - OE_stamps_laser_n[i]) < freq_acq:
                    interm = OE_stamps_laser_n[i] + 10
                    OE_stamps_laser_n.append(interm)


            # _n is from 0
            acquisition_mini_start_n = OE_stamps_miniscope_n[0]
            # _a is from acquisition time software
            acquisition_mini_start_a = B[0]

            outSumm = pd.Series([acquisition_mini_start_n, acquisition_mini_start_a, freq_acq, dropped_frames],
                        index=['Miniscope start from 0', 'Miniscope start from Acq time soft', 'mini acq freq', 'dropped frames'])
            print(outSumm)

            filenameOut = folder / f'SynchroFileCorrect.xlsx'
            writer = pd.ExcelWriter(filenameOut)
            outSumm.to_excel(writer)
            writer.close()

        except:
            print('NO SYNC FILE SAVED. Skipping...')
            continue
