import mbTools.mbTools
from . import ePhy

import numpy as np

class IntanLFP(ePhy):
   def __init__(self, parent: mbTools.mbTools.experiment, files_list, numChannels = 32, recSyst = 'Bonsai') -> None:
      super().__init__(parent, numChannels = numChannels)
      self.files_list = files_list
      self.recSyst = recSyst
      self.fileType = 'IntanLFP'
      self.signal = None
      self.dtype=np.uint16
      self.offset = int(np.iinfo(self.dtype).max/2)
      self.signal = self.loadData()


   def loadData(self):
      return super().loadData()
   
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
      self.signal = None
      self.dtype=np.uint16
      self.signal = self.loadData()

   def loadData(self):
      spikesPrefix = 'NP_spikes_'
      if len(self.files_list)>1:
         raise Exception(f"Multiple files not yet in place, please contact MB if you are interested by this option")
      for filename in self.files_list:
         # Neuropixels V1 Probe
         npix = {}
         npix['frame-counter'] = np.fromfile(filename.replace(spikesPrefix,'NP_FrameCounter_'), dtype=np.int32)


         npix['spike'] = np.fromfile(filename, dtype=np.uint16).reshape(-1, self.numChannels)
         npix['spike-clock'] = np.fromfile(filename.replace(spikesPrefix,'NP_timestamps_'), dtype=np.uint64)

         npix['lfp'] = np.fromfile(filename.replace(spikesPrefix,'NP_LFPdata_'), dtype=np.uint16).reshape(-1, self.numChannels)
         # npix['lfp-clock'] = np.fromfile(filename.replace(spikesPrefix,'NP_FrameCounter_'), dtype=np.uint64)
      
      self.signal = npix
      self.channelLabels = [i for i in range(384)]
      return self.signal