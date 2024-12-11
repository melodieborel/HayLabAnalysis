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
      self.signal = None
      self.expe = parent
      self.signal = self.loadData()
      self.getOnsets()

   def loadData(self):
      if len(self.files_list)>1:
         raise Exception(f"Multiple files not implemented yet, please contact MB if you are interested by this option")
      #return np.fromfile("C:/Users/Manip7/Documents/data/Bonsai_data/Lou/Habituation/2024_12_11/SleepBefore_10_03_20/Lou_OE_ttlin_2024-12-11T10_04_04.bin", dtype=np.uint8)
      d = np.fromfile(self.files_list[0], dtype=np.uint16)
      d=(d-np.iinfo(np.uint16).max/2)/np.power(2,7)
      return d.astype(int)
   
   def getOnsets(self):
      threshold_crossings = np.diff(self.signal > 0.7, prepend=False)
      self.ttlIndexes = np.argwhere(threshold_crossings)[:,0]
      d = np.diff(self.ttlIndexes,prepend=[0])
      bigOnsets =  d[d > 200] #plus de 200 samples entre 2 ttlIns
      return bigOnsets
