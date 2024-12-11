
import numpy as np
import os
import configparser
import ast
from pathlib import Path

class ePhy():
   def __init__(self, parent, numChannels = None) -> None:
      #super().__init__()
      self.numChannels = int(numChannels)
      self.recSyst = None
      self.fileType = None
      self.dtype = None
      self.channelsMap = dict()
      self.channelLabels = []
      self.files_list = []
      self.signal = np.zeros((0,0))
      self.offset = 0
      self.expe = parent
      self.start = 0
      self.sampling_rate = 1000


   def loadData(self):
      signal = np.zeros([0], dtype=self.dtype)
      for file in self.files_list:
         print("importing {}".format(file))
         file_data = np.fromfile(file, dtype=self.dtype)
         signal=np.append(signal, file_data, axis=0)
      if self.offset != 0:
         print('applying offset')
         signal = signal - self.offset
      if signal.dtype is not np.dtype(np.int16):
         print('converting to int16')
         signal = signal.astype(np.int16)
      signal = signal.reshape(-1,self.numChannels)

      print(f'{self.fileType} file loaded, with {self.numChannels} channels and {signal.shape[0]} datapoint')
      self.signal = signal
      return self.signal

   def loadMetaData(self):
      self.channelsMap = self.expe.channelsMap
      self.start=ast.literal_eval(self.expe.parser['OE_LFP']['start'])
      self.sampling_rate=ast.literal_eval(self.expe.parser['OE_LFP']['freq'])
      NPX = ast.literal_eval(self.expe.parser['OE_LFP']['NPX'])
      timesreset = ast.literal_eval(self.expe.parser['OE_LFP']['timesreset'])

      print("the mapping:", self.channelsMap)
      print("the offset: ", self.start)
      print("the sampling rate: ", self.sampling_rate)
   
   def combineStructures(self, structures=None, start = 0, end = None):
      """Retrieve a combined array with either all cannals (if structures is None (default)), or the differential signals correspondin to the mapped structures

      Args:
         structures (None, "All", or array of structures, optional): indicates what data to combine. Defaults to None.
         start (int, optional): if only part of the data to display. Defaults to 0.
         end (optional): if only part of the data to display. Defaults to None.

      Returns:
         array: combined numpy array ready to visualize
      """
      if end is None:
         end = self.signal.shape[0]
      combined = np.empty((end-start,0),np.int16)
      self.channelLabels = []
      if structures is None:
         #self.generateChannelsMap()
         combined = self.signal
         self.channelLabels = [i for i in range(self.signal.shape[0])]
      else:
         if structures=='All':
            structures = self.channelsMap.keys()
            print(structures)
         for region in structures:
            print(region, "->", self.channelsMap[region])
            if len([canal["canal"] for canal in self.channelsMap[region] if canal["status"]==2])>0:
               c2 = int([canal["canal"] for canal in self.channelsMap[region] if canal["status"]==2][0])
               c1 = int([canal["canal"] for canal in self.channelsMap[region] if canal["status"]==1][0])
               print("Getting differential signal of channel {} - channel {} for {}".format(c2,c1,region))
               self.channelLabels.append(region)
               combined = np.append(combined, self.signal[start:end, c2, np.newaxis] - self.signal[:, c1, np.newaxis], axis=1)
            elif len([canal["canal"] for canal in self.channelsMap[region] if canal["status"]==1])>0:
               c = int([canal["canal"] for canal in self.channelsMap[region] if canal["status"]==1][0])
               print("Getting floating signal of channel {} for {}".format(c,region))
               combined = np.append(combined, self.signal[start:end,c, np.newaxis], axis=1)
               self.channelLabels.append(region)
      return combined