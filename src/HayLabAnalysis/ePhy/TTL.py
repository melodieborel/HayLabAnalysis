from .. import experiment
from . import ePhy
from datetime import datetime, timedelta
import os
import re
from scipy import signal, ndimage

import numpy as np
import matplotlib.pyplot as plt



class TTLin():
   def __init__(self, parent: experiment, files_list) -> None:
      self.files_list = files_list
      self.expe = parent
      self.ttlIndexes = None
      self.signal = self.loadData()
      self.getOnsets()
      self.times = self.getTTLTimes(self.expe.data['OE_LFP'].sampling_rate)

   def loadData(self):
      if len(self.files_list)>1:
         raise Exception(f"Multiple files not implemented yet, please contact MB if you are interested by this option")
      d = np.fromfile(self.files_list[0], dtype=np.uint16)
      d=(d-np.iinfo(np.uint16).max/2)/np.power(2,7)
      return d.astype(int)
   
   def getOnsets(self):
      self.signal = self.loadData()#[:150000]
      threshold_crossings = np.diff(self.signal > 0.7, prepend=self.signal[0])
      #threshold_crossings[threshold_crossings<0]=0 #keep only rising edges
      indexes = np.argwhere(threshold_crossings)[:,0]
      d = np.diff(indexes,prepend=[0])
      findburst = np.argwhere(d>200) #plus de 200 samples entre 2 ttlIns
      bigOnsets = indexes[findburst].flatten()
      split_lim=np.where(np.isin(indexes,bigOnsets))[0]
      self.ttlIndexes = np.split(indexes, split_lim)[1:]
      return bigOnsets

   def getTTLTimes(self,freq):
      times = [x / freq for x in self.ttlIndexes]
      return times
