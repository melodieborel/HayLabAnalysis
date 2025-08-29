
import numpy as np
import os
import configparser
import ast
from pathlib import Path

class Imaging():
   def __init__(self, parent) -> None:
      #super().__init__()
      self.files_list = []
      self.expe = parent
      self.signal = None

   def loadData(self):
      pass

   def loadMetaData(self):
      pass

   
   def loadTimeStamps(self, data_pattern, ts_pattern ):
      self.TS_times = []
      for f in self.files_list:
         fTS=f.with_name(f.name.replace(data_pattern,ts_pattern).replace('.avi','.csv'))

         import pandas as pd
         df = pd.read_csv(fTS, sep=',', header=None, names=['ts'],usecols=[0])
         df.ts = pd.to_datetime(df.ts, format='ISO8601')
         if self.expe.refTime is None:
            refT = df.iloc[0]['ts']
         else:
            refT = self.expe.refTime
         df['delays']=(df.ts-refT).dt.total_seconds()
         times=df.delays.values
         #times-=self.start
         #self.TS_sampling_rate = 1/(df.ts.diff().mean().total_seconds())

         self.TS_times.append(times)
