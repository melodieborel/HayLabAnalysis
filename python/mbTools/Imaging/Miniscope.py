from .. import experiment
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
      self.signal = self.loadData()
      self.loadMetaData()

   def loadData(self):
      return super().loadData()
   
   def loadMetaData(self):
      super().loadMetaData()

   def loadTimeStamps(self):
      pass
