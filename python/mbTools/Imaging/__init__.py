
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