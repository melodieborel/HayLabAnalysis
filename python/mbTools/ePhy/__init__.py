
import numpy as np

class ePhy():
   def __init__(self, parent, numChannels = None) -> None:
      #super().__init__()
      self.numChannels = numChannels
      self.recSyst = None
      self.fileType = None
      self.dtype = None
      self.channelsMap = dict()
      self.channelLabels = []
      self.files_list = []
      self.signal = None
      self.offset = 0
      self.expe = parent


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

   def combineStructures(self,structures, start = 0, end = None):
      if end is None:
         end = self.All.shape[0]
      combined = np.empty((end-start,0),np.int16)
      self.channelLabels = []
      if structures is None:
         #self.generateChannelsMap()
         combined = self.All
         self.channelLabels = [i for i in range(self.All.shape[0])]
      else:
         for region in structures:
            print(region, "->", self.channelsMap[region])
            if len([canal["canal"] for canal in self.channelsMap[region] if canal["status"]==2])>0:
               c2 = [canal["canal"] for canal in self.channelsMap[region] if canal["status"]==2][0]
               c1 = [canal["canal"] for canal in self.channelsMap[region] if canal["status"]==1][0]
               print("Getting differential signal of channel {} - channel {} for {}".format(c2,c1,region))
               self.channelLabels.append(region)
               combined = np.append(combined, self.All[start:end, c2, np.newaxis] - self.All[:, c1, np.newaxis], axis=1)
            elif len([canal["canal"] for canal in self.channelsMap[region] if canal["status"]==1])>0:
               c = [canal["canal"] for canal in self.channelsMap[region] if canal["status"]==1][0]
               print("Getting floating signal of channel {} for {}".format(c,region))
               combined = np.append(combined, self.All[start:end,c, np.newaxis], axis=1)
               self.channelLabels.append(region)
      return combined