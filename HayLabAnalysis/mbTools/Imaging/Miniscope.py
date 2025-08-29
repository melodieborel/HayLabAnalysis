from .. import experiment, tools
from . import Imaging
from datetime import datetime, timedelta
import os
import re
from scipy import signal, ndimage

import numpy as np
import matplotlib.pyplot as plt


class Miniscope(Imaging):
   def __init__(self, parent: experiment, files_list) -> None:
      super().__init__(parent)
      self.files_list = files_list
      self.times = []
      self.signal = self.loadData()
      self.loadMetaData()
      self.loadTimeStamps()
      self.alignMiniscopeWithLFP(self.expe.data['OE_LFP'].TS_times[0])

   def loadData(self):
      return super().loadData()
   
   def loadMetaData(self):
      super().loadMetaData()

   def loadTimeStamps(self):
      super().loadTimeStamps("miniscope","miniscopeTS")

   def alignMiniscopeWithLFP(self, deltaT):
      self.times = []
      for t in self.TS_times:
         self.times.append(t - deltaT)
      print("miniscope realigned")


class Webcam(Imaging):
   def __init__(self, parent: experiment, files_list) -> None:
      super().__init__(parent)
      self.files_list = files_list
      self.times = []
      self.signal = self.loadData()
      self.loadMetaData()
      self.loadTimeStamps()
      self.alignWebcamWithMiniscope()

   def loadData(self):
      return super().loadData()
   
   def loadMetaData(self):
      super().loadMetaData()

   def loadTimeStamps(self):
      super().loadTimeStamps("webcam","webcamTS")

   def alignWebcamWithMiniscope(self):
      from scipy import interpolate
      self.times = []
      for trial in range(len(self.expe.data['miniscope'].TS_times)):
         y = self.expe.data['miniscope'].TS_times[trial]
         x = self.expe.data['ttl'].times[trial][:y.shape[0]]
         arrayToAlign = self.TS_times[trial]
         
         #newArray = np.interp(arrayToAlign, y, x)
         f = interpolate.interp1d(y, x, fill_value = "extrapolate", kind='linear')
         newArray = f(arrayToAlign)
         print("new array is", newArray)
         if False:
            plt.close()
            plt.plot(newArray,arrayToAlign)
            plt.plot(x,y)
            #plt.plot(newArray)
            plt.show()
         if True:
            print("first miniscope image ts :", y[0])
            print("corresponding real time based on ttls", x[0])
            b = tools.find_nearest(arrayToAlign,y[0])
            print("miniscope onset at webcam image number: ",b)
            print("realigned times aroud that frame are",newArray[b-2:b+2])
         self.times.append(newArray)
      print("webcam realigned")