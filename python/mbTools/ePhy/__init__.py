
import numpy as np
import os
import pandas as pd

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

   
   def loadSpindles(self, relativePath='', structure = "M1", suffix=''):
      base = os.path.normpath(self.expe.expePath)
      alltracesFN = os.path.sep.join([base,relativePath,'RawDataChannelExtractedDS.npy'])
      filename2 = os.path.sep.join([base,relativePath,f'Spindlesproperties_{structure}{suffix}.pkl'])
      filename3 = os.path.sep.join([base,relativePath,f'Spindlesproperties_{structure}{suffix}.csv'])
      filenameData =os.path.sep.join([base,relativePath,f'Signal{structure}{suffix}.npy'])
      if os.path.isfile(filename3):
         All_Spindle = pd.read_csv(filename3, sep=',', header=0, index_col=0)
         if 'toKeep' not in All_Spindle:
            All_Spindle['toKeep'] = True
         print(f"file {filename3} was found so loading it")
      else:
         print(f"file {filename3} doesn't exist yet which is probably because you haven't done the analysis or the relative path is wrong")
         All_Spindle = None
      if os.path.isfile(filenameData):
         signal=np.load(filenameData, mmap_mode= 'r')[:,0]
         print(f"file {filenameData} was found so loading it")
      else:
         print(f"no file {filenameData} was found, please make sure it exists")
         Signal = None
      return All_Spindle, signal