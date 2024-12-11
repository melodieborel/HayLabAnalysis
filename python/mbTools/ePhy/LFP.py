from .. import experiment
from . import ePhy
from datetime import datetime, timedelta
import os
import re
from scipy import signal, ndimage

import numpy as np
import matplotlib.pyplot as plt


def find_nearest(array, value):
   idx = (np.abs(array - value)).argmin()
   return idx

class IntanLFP(ePhy):
   def __init__(self, parent: experiment, files_list, numChannels = 32, recSyst = 'Bonsai') -> None:
      super().__init__(parent, numChannels = numChannels)
      self.files_list = files_list
      self.recSyst = recSyst
      self.signal = None
      self.bufferCount=1024
      self.times=None
      if recSyst=='Bonsai':
         print('data recorded with Bonsai')
         self.fileType = 'IntanLFP'
         self.dtype=np.uint16
         self.offset = int(np.iinfo(self.dtype).max/2)
      else:
         print('data not recorded with Bonsai')
         self.fileType = 'OE32channels.bin'
         self.dtype=np.int16
         self.offset = 0
      self.sampling_rate=20000
      self.signal = self.loadData()
      self.loadMetaData()

   def loadData(self):
      return super().loadData()
   

   def reAlignTimes(self):
      self.times=np.linspace(0,self.signal.shape[0]/self.sampling_rate,self.signal.shape[0])+self.start


   def resetAlign(self):
      freq=20000
      self.times=np.linspace(0,self.signal.shape[0]/freq,self.signal.shape[0])

   def loadMetaData(self):
      super().loadMetaData()
      self.reAlignTimes()
      
   def updateParser(self,key,value):
      expeConfigFN=self.files_list[0].with_name('expeConfig.ini')
      self.parser['OE_LFP'][key]=str(value)
      with open(expeConfigFN, 'w') as configfile:
            self.parser.write(configfile)



   def loadTimeStamps(self):
      if len(self.files_list)>1:
         raise Exception(f"Multiple files not implemented yet, please contact MB if you are interested by this option")
      for f in self.files_list:
         fTS=f.replace('OE_32ch_data','OE_32ch_timestamps').replace('.bin','.csv')
         seps=[m.start() for m in re.finditer('_',f.name)]
         datestr = fn[seps[2]+1:-4]
         launch_start = datetime.strptime(datestr, '%Y-%m-%dT%H_%M_%S').astimezone()
         
         import pandas as pd
         df = pd.read_csv(fTS, sep=',', header=None, names=['ts'])
         df.ts = pd.to_datetime(df.ts, format='ISO8601')
         df['delays']=(df.ts-df.iloc[0]['ts']).dt.total_seconds()
         self.times=df.delays.values
         grid=np.arange(1024)
         instFreq=1024/np.diff(self.times)
         self.times=np.concatenate([grid/instFreq[i]+self.times[i] for i in range(self.times.shape[0]-1)])

         self.sampling_rate = self.bufferCount/(df.ts.diff().mean().total_seconds())

         self.start = (df.loc[0]['ts']-launch_start).total_seconds()
         #self.times+=self.start
         #launch_start = df.loc[0]['ts'] - timedelta(seconds=self.start)
         print(f"the calculated sampling rate is {self.sampling_rate} Hz")
         print(f"the recording tool {self.start} s to start")
         print(f"the first timestamp is {df.iloc[0]['ts']}")
         print(f"the before last timestamp is {df.iloc[-2]['ts']}")
         print(f"the last timestamp is {df.iloc[-1]['ts']}")
         print(f"the calculated launch start is {launch_start}")

         return df
   
   def convertTime2Index(self,t):
      self.times
      idx = 0
      return idx
   
   def combineStructures(self, structures=None, start=0, end=None):
      return super().combineStructures(structures, start, end)
   
class NPX(ePhy):
   def __init__(self, parent: experiment, files_list, numChannels = 384) -> None:
      super().__init__(parent, numChannels = numChannels)
      self.files_list = files_list
      self.fileType = 'NPX'
      self.signal = {}
      self.dtype=np.uint16
      self.sampling_rate=30000
      self.acquisitionClockHz=250000000
      self.signal = self.loadData()

   def loadData(self):
      spikesPrefix = 'NP_spikes_'
      if len(self.files_list)>1:
         raise Exception(f"Multiple files not implemented yet, please contact MB if you are interested by this option")
      for filename in self.files_list:
         # Neuropixels V1 Probe
         npix = {}
         #npix['frame-counter'] = np.fromfile(filename.replace(spikesPrefix,'NP_FrameCounter_'), dtype=np.int32)

         self.loadTimeStamps()
         import spikeinterface.full as si
         import pickle
         npix['spike'] = si.read_binary(filename, dtype='uint16', num_channels=384, sampling_frequency=self.sampling_rate)
         try:
            with open('//10.69.168.1/crnldata/waking/audrey_hay/NPX/NPXprobe.pkl', 'rb') as outp: 
               probe = pickle.load(outp)
         except Exception:
            with open('/Volumes/waking/audrey_hay/NPX/NPXprobe.pkl', 'rb') as outp: 
               probe = pickle.load(outp)
         probe.set_device_channel_indices(np.arange(384))

         npix['spike'] = npix['spike'].set_probe(probe)

         #npix['lfp'] = np.fromfile(filename.replace(spikesPrefix,'NP_LFPdata_'), dtype=np.uint16).reshape(-1, self.numChannels)
         # npix['lfp-clock'] = np.fromfile(filename.replace(spikesPrefix,'NP_FrameCounter_'), dtype=np.uint64)
      
      if False: #self.times.shape[0]<npix['spike'].get_num_frames():
         rab=npix['spike'].get_num_frames()-self.times.shape[0]
         print(f"there are {rab} less timestamps than samples so extrapoling rest")
         extraTS = (np.arange(rab) + 1) * self.sampling_rate + self.times[-1]
         self.times = np.append(self.times, extraTS)
      print(self.times.shape)
      print(npix['spike'].get_num_frames())
      npix['spike'].set_times(self.times[:npix['spike'].get_num_frames()]) #rq not working if set before probe
      self.signal = npix
      self.channelLabels = [i for i in range(384)]
      return self.signal
   
   def loadTimeStamps(self):
      spikesPrefix = 'NP_spikes_'
      if len(self.files_list)>1:
         raise Exception(f"Multiple files not implemented yet, please contact MB if you are interested by this option")
      for filename in self.files_list:
         seps=[m.start() for m in re.finditer('_',filename.name)]
         datestr = filename.name[seps[1]+1:-4]
         launch_start = datetime.strptime(datestr, '%Y-%m-%dT%H_%M_%S').astimezone()
         offset = 0.896598400
         launch_start+= timedelta(seconds=offset)
         self.signal['spike-clock'] = np.fromfile(filename.with_name(filename.name.replace(spikesPrefix,'NP_timestamps_')), dtype=np.uint64)
      self.times = self.signal['spike-clock']/self.acquisitionClockHz
      self.sampling_rate = self.acquisitionClockHz/np.diff(self.signal['spike-clock']).mean()
      print(f"the calculated sampling rate is {self.sampling_rate} Hz")
      print(f"launch start would be {launch_start}")
      self.start = (timedelta(seconds=self.signal['spike-clock'][0]/self.acquisitionClockHz)).total_seconds()
      print(f"the interval to first clock is {self.start}")
      print(f"the first timestamp for {self.signal['spike-clock'][0]} samples, corresponding to {(timedelta(seconds=self.signal['spike-clock'][0]/self.acquisitionClockHz)).total_seconds()} s, would be {launch_start + timedelta(seconds=self.signal['spike-clock'][0]/self.acquisitionClockHz)}")
      #print(f"the before last timestamp {self.signal['spike-clock'][222587890]} would be {launch_start + timedelta(seconds=self.signal['spike-clock'][-2]/self.acquisitionClockHz)}")
      #print(f"the last timestamp {self.signal['spike-clock'][222587891]} would be {launch_start + timedelta(seconds=self.signal['spike-clock'][-1]/self.acquisitionClockHz)}")
      print(f"there are {len(self.signal['spike-clock'])} timestamps")


class LFP_DS(ePhy):
   def __init__(self, parent: experiment, files_list, numChannels = 32) -> None:
      super().__init__(parent, numChannels=numChannels)
      self.files_list = files_list
      self.fileType = 'NPY'
      self.signal = {}
      self.dtype=np.int16
      self.offset = 0
      self.sampling_rate=1000
      self.loadMetaData()
      self.signal = self.loadData()
      self.reAlignData(realSR=self.sampling_rate)

   def loadData(self):
      if len(self.files_list) == 1:
         file = self.files_list[0]
         signal = np.load(file, mmap_mode= 'r', allow_pickle=True)
         if signal.shape[0]<=32*4:
            signal = np.transpose(signal)
         self.numChannels = signal.shape[1]
         self.channelsMap = self.expe.channelsMap
      else:
         raise Exception(f"Several npy files ({len(self.files_list)}) to merge but this is not in place yet")
      return signal
   
   def combineStructures(self, structures=None, start=0, end=None):
      return super().combineStructures(structures, start, end)
   
   def reAlignData(self,realSR=20000):
      freqInitTheoric=20000
      freqDS=1000
      realignFactor=freqInitTheoric/realSR
      print(self.sampling_rate)
      self.sampling_rate=freqDS*realignFactor
      print(realignFactor)
      print(self.sampling_rate)

   def loadMetaData(self):
      return super().loadMetaData()
   
   def filterEMG(self, f_lowcut = 200., f_hicut = 400., N = 4):
      if 'EMG' in self.channelsMap.keys():
         EMG = self.combineStructures(['EMG'])[:,0]
      
         nyq = 0.5 * self.sampling_rate
         Wn = [f_lowcut/nyq,f_hicut/nyq]  # Nyquist frequency fraction

         # Filter creation :
         b, a = signal.butter(N, Wn, 'band')
         filt_EMG = signal.filtfilt(b, a, EMG)

         # Plot
         timesmin = np.linspace(0, EMG.size/self.sampling_rate/60, num=EMG.size)

         plt.close()
         fig, ax = plt.subplots()
         ax.plot(timesmin, EMG)
         ax.plot(timesmin, filt_EMG)
         return filt_EMG
      else:
         print("Couldn't find an EMG signal")
         return None

   def sleep_scoring(self, display = False, low_cwt_factor=0.5, high_cwt_factor=3):
      if 'EMG' in self.channelsMap.keys():
         EMG = self.combineStructures(['EMG'])[:,0]

         # Parameter and computation of CWT
         w = 4.
         freq = np.linspace(200, 400, 50)
         widths = w*self.sampling_rate / (2*freq*np.pi)
         EMGcwt = signal.cwt(EMG, signal.morlet2, widths, w=w)

         # Projection calculation
         absEMGcwt = np.absolute(EMGcwt)
         proj_EMGcwt = np.sum(absEMGcwt, axis = 0)/50
         sdproj_EMGcwt = np.std(proj_EMGcwt)
         mproj_EMGcwt = np.mean(proj_EMGcwt)

         sd_low_proj_EMGcwt = mproj_EMGcwt + sdproj_EMGcwt * low_cwt_factor
         sd_high_proj_EMGcwt = mproj_EMGcwt + sdproj_EMGcwt * high_cwt_factor

         if display:

            # Defining subset
            start = 00000
            end = 400000

            times = np.arange(0, EMG.size/self.sampling_rate, 1./self.sampling_rate)
            tt = times[start:end]
            EMGt = EMG[start:end]
            EMGcwtt = EMGcwt[:, start:end]
            proj_EMGcwtt = proj_EMGcwt[start:end]

            plt.close()
            plt.axhline(sdproj_EMGcwt, color='r') # horizontal
            plt.axhline(sd_high_proj_EMGcwt, color='g') # horizontal
            plt.axhline(sd_low_proj_EMGcwt, color='b') # horizontal
            plt.plot(tt, EMGt)
            plt.plot(tt, proj_EMGcwtt)
            plt.show()

         EMGstatusRaw2 = self.assignSleep(proj_EMGcwt, sd_low_proj_EMGcwt, sd_high_proj_EMGcwt)

         return EMGstatusRaw2
      else:
         print("Couldn't find an EMG signal")
         return None
      

   def assignSleep(self, proj_EMGcwt, sd_low_proj_EMGcwt, sd_high_proj_EMGcwt):
      numpnts = self.signal.shape[0]
      EMGstatusRaw = np.zeros(numpnts)

      # Assigning values wake (1, 2) and sleep (0)
      EMGstatusRaw[proj_EMGcwt>sd_low_proj_EMGcwt]=1
      EMGstatusRaw[proj_EMGcwt>sd_high_proj_EMGcwt]=2

      # Expanding borders for wake (1, 2) and sleep (0) to ±1 s around detected muscular activity
      lim=1000
      kernel = np.ones(2*lim)
      EMGstatusRaw2=EMGstatusRaw.copy()
      mask2 = ndimage.binary_dilation(EMGstatusRaw2==2, kernel)
      EMGstatusRaw2[mask2]=2
      mask1 = ndimage.binary_dilation(EMGstatusRaw2==1, kernel)
      EMGstatusRaw2[mask1]=np.maximum(EMGstatusRaw2[mask1],np.ones(EMGstatusRaw2[mask1].shape[0]))

      plt.close()
      plt.plot(EMGstatusRaw)
      plt.plot(EMGstatusRaw2)
      plt.show()
      
      return EMGstatusRaw2