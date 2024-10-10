import mbTools.mbTools
from . import ePhy
from datetime import datetime, timezone, timedelta
import os
import re
import configparser
import ast

import numpy as np

def find_nearest(array, value):
    idx = (np.abs(array - value)).argmin()
    return idx

class IntanLFP(ePhy):
   def __init__(self, parent: mbTools.mbTools.experiment, files_list, numChannels = 32, recSyst = 'Bonsai') -> None:
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
      #self.loadTimeStamps()

   def loadData(self):
      return super().loadData()
   

   def reAlignTimes(self):
      self.times=np.linspace(0,self.signal.shape[0]/self.sampling_rate,self.signal.shape[0])+self.start


   def resetAlign(self):
      freq=20000
      self.times=np.linspace(0,self.signal.shape[0]/freq,self.signal.shape[0])

   def loadMetaData(self):
      bn=os.path.split(self.files_list[0])[0]
      expeConfigFN=os.path.sep.join([bn,'expeConfig.ini'])
      self.parser = configparser.ConfigParser()
      self.parser.read(expeConfigFN)

      if os.path.isfile(expeConfigFN):
         print('mapping exists so loading it')
         self.channelsMap = ast.literal_eval(self.parser['OE_LFP']['channelsMap'])
         self.start=ast.literal_eval(self.parser['OE_LFP']['start'])
         self.sampling_rate=ast.literal_eval(self.parser['OE_LFP']['freq'])
         NPX = ast.literal_eval(self.parser['OE_LFP']['NPX'])
         timesreset = ast.literal_eval(self.parser['OE_LFP']['timesreset'])
      else:
         print("mapping doesn't exist so generating it")
         self.channelsMap = dict( \
                  M1 = [dict(canal = 17, status=1),
                     dict(canal = 16, status=2)],
            )

         self.parser['OE_LFP'] = {'channelsMap': self.channelsMap}


         artefacts=[]
         self.parser['OE_LFP']['NPX']=str(artefacts)
         self.parser['OE_LFP']['timesreset']=str(artefacts)


         self.parser['OE_LFP']['start']=str(52)
         self.parser['OE_LFP']['freq']=str(2000)

         with open(expeConfigFN, 'w') as configfile:
            self.parser.write(configfile)

      print("the mapping:", self.channelsMap)
      print("the offset: ", self.start)
      print("the sampling rate: ", self.sampling_rate)

      self.reAlignTimes()
        
   def updateParser(self,key,value):
      bn=os.path.split(self.files_list[0])[0]
      expeConfigFN=os.path.sep.join([bn,'expeConfig.ini'])
      self.parser['OE_LFP'][key]=str(value)
      with open(expeConfigFN, 'w') as configfile:
            self.parser.write(configfile)



   def loadTimeStamps(self):
      if len(self.files_list)>1:
         raise Exception(f"Multiple files not implemented yet, please contact MB if you are interested by this option")
      for f in self.files_list:
         fTS=f.replace('OE_32ch_data','OE_32ch_timestamps').replace('.bin','.csv')
         
         fn=os.path.split(f)[1]
         seps=[m.start() for m in re.finditer('_',fn)]
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
      return idx
   
   def combineStructures(self,structures, start = 0, end = None):
      if end is None:
         end = self.signal.shape[0]
      combined = np.empty((end-start,0),np.int16)
      self.channelLabels = []
      if structures is None:
         #self.generateChannelsMap()
         combined = self.signal
         self.channelLabels = [i for i in range(self.signal.shape[0])]
      else:
         for region in structures:
            print(region, "->", self.channelsMap[region])
            if len([canal["canal"] for canal in self.channelsMap[region] if canal["status"]==2])>0:
               c2 = [canal["canal"] for canal in self.channelsMap[region] if canal["status"]==2][0]
               c1 = [canal["canal"] for canal in self.channelsMap[region] if canal["status"]==1][0]
               print("Getting differential signal of channel {} - channel {} for {}".format(c2,c1,region))
               self.channelLabels.append(region)
               combined = np.append(combined, self.signal[start:end, c2, np.newaxis] - self.signal[:, c1, np.newaxis], axis=1)
            elif len([canal["canal"] for canal in self.channelsMap[region] if canal["status"]==1])>0:
               c = [canal["canal"] for canal in self.channelsMap[region] if canal["status"]==1][0]
               print("Getting floating signal of channel {} for {}".format(c,region))
               combined = np.append(combined, self.signal[start:end,c, np.newaxis], axis=1)
               self.channelLabels.append(region)
      return combined
   
class NPX(ePhy):
   def __init__(self, parent: mbTools.mbTools.experiment, files_list, numChannels = 384) -> None:
      super().__init__(parent, numChannels = numChannels)
      self.files_list = files_list
      self.fileType = 'NPX'
      self.signal = {}
      self.dtype=np.uint16
      self.sampling_rate=30000
      self.acquisitionClockHz=250000000
      self.signal = self.loadData()
      self.loadTimeStamps()

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
      
      self.signal = npix
      self.channelLabels = [i for i in range(384)]
      return self.signal
   
   def loadTimeStamps(self):
      spikesPrefix = 'NP_spikes_'
      if len(self.files_list)>1:
         raise Exception(f"Multiple files not implemented yet, please contact MB if you are interested by this option")
      for filename in self.files_list:
         fn=os.path.split(filename)[1]
         seps=[m.start() for m in re.finditer('_',fn)]
         datestr = fn[seps[1]+1:-4]
         launch_start = datetime.strptime(datestr, '%Y-%m-%dT%H_%M_%S').astimezone()
         offset = 0.896598400
         launch_start+= timedelta(seconds=offset)
         self.signal['spike-clock'] = np.fromfile(filename.replace(spikesPrefix,'NP_timestamps_'), dtype=np.uint64)
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