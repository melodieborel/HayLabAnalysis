import os
import configparser
import ast
import fnmatch
from pathlib import Path

import numpy as np
import pandas as pd

import re

from ipyfilechooser import FileChooser
import ipywidgets as widgets
from IPython.display import display

from .tools import expeConfigDict, magicstore

class experiment():
   """experiment is a class for a whole experiment including all its component (NPX, miniscope, intan...)
   """
   def __init__(self, expe = None, numChannels = 32) -> None:
      if isinstance(expe,expeConfigDict):
         self.expe = expe
         self.expePath = self.expe.rawDataPath
         self.interimAnalysisPath = self.expe.interimAnalysisPath
      else:
         self.expe = None
         self.expePath = expe
         self.interimAnalysisPath = ''
      self.numChannels = numChannels
      self.recSyst = None
      self.fileType = None
      self.channelsMap = dict()
      self.All = None
      self.channelLabels = []
      self.data = dict()

      try:
         fc1 = FileChooser(self.expePath,select_default=True, show_only_dirs = True, title = "<b>OpenEphys Folder</b>")
         display(fc1)
         fc1.register_callback(self.update_my_expe_choice)
      except Exception as error:
         print(f"something went wrong, make sure the experiment path ({self.expePath}) is a folder")

   # Function to find files containing a specific string
   def find_files_with_string(self, folder_path, search_string):
      matching_files = []
      # Traverse the folder to find files
      for root, _, files in os.walk(folder_path):
         for file in files:
               if fnmatch.fnmatch(file, f"*{search_string}*"):
                  matching_files.append(os.path.join(root, file))
      return matching_files
   
   def update_my_expe_choice(self,chooser):
      dpath = chooser.selected
      magicstore('dpath', dpath)
   

   def analyseExpe_findData(self, fullSampling=False, spindleBN='Spindlesproperties',suffix=''):
      """findData: function that analyse the content of the raw data folder and detects component of the experiment to load all of them
      """
      from .ePhy.LFP import IntanLFP, NPX, LFP_DS
      self.loadAnimalInfo()
      folderpath = Path(self.expePath)
      print(folderpath)
      DSdata=False

      if self.find_files_with_string(self.interimAnalysisPath,  "RawDataChannelExtractedDS.npy"): # pre-analysed data
         print('found some RawDataChannelExtractedDS.npy files')
         matching_files = self.find_files_with_string(self.interimAnalysisPath, "RawDataChannelExtractedDS.npy")
         self.data['LFP_DS'] = LFP_DS(self, matching_files)
         DSdata=True

      if self.find_files_with_string(self.interimAnalysisPath,  f"{spindleBN}_*.csv"): #NPX's Data
         print('found some Spindles files')
         matching_files = self.find_files_with_string(self.interimAnalysisPath, f"{spindleBN}_*.csv")
         All_Spindles = self.loadSpindles(matching_files,spindleBN,suffix)
         self.data['Spindles'] = All_Spindles

      if not DSdata and fullSampling:
         if self.find_files_with_string(folderpath,  ".bin"): #Bonsai or IgorPro
            print('found some .bin files')
            matching_files = self.find_files_with_string(folderpath, ".bin")
            self.data['OE_LFP'] = IntanLFP(self, matching_files)
      
         if self.find_files_with_string(folderpath,  "continuous.dat"): #OpenEphys
            print('found some continuous.dat files')
            matching_files = self.find_files_with_string(folderpath, "continuous.dat")
            print('carrefull, to match my case, numChannels is set to 64')
            self.data['OE_LFP'] = IntanLFP(self, matching_files, recSyst = 'OpenEphys', numChannels=64)


      if self.find_files_with_string(folderpath,  "NP_spikes_*.raw"): #NPX's Data
         print('found some NPX files')
         matching_files = self.find_files_with_string(folderpath, "NP_spikes_*.raw")
         self.data['NPX'] = NPX(self, matching_files)

   def loadAnimalInfo(self):
      print(self.expe.interimAnalysisPath)
      return
      bn=os.path.split(self.files_list[0])[0]
      expeConfigFN=os.path.sep.join([bn,'expeConfig.ini'])
      self.parser = configparser.ConfigParser()
      self.parser.read(expeConfigFN)

      if os.path.isfile(expeConfigFN):
         print('mapping exists so loading it')
         self.channelsMap = ast.literal_eval(self.parser['OE_LFP']['channelsMap'])
         self.start=ast.literal_eval(self.parser['OE_LFP']['start'])
         self.sampling_rate=ast.literal_eval(self.parser['OE_LFP']['freq'])
         NPX = ast.literal_eval(self.parser['OE_LFP']['NPX'])
         timesreset = ast.literal_eval(self.parser['OE_LFP']['timesreset'])
      else:
         print("mapping doesn't exist so generating it")
         self.channelsMap = dict( \
                  M1 = [dict(canal = 17, status=1),
                     dict(canal = 16, status=2)],
            )
         self.start=52
         self.sampling_rate=20046

         self.parser['OE_LFP'] = {'channelsMap': self.channelsMap}


         artefacts=[]
         self.parser['OE_LFP']['NPX']=str(artefacts)
         self.parser['OE_LFP']['timesreset']=str(artefacts)


         self.parser['OE_LFP']['start']=str(self.start)
         self.parser['OE_LFP']['freq']=str(self.sampling_rate)

         with open(expeConfigFN, 'w') as configfile:
            self.parser.write(configfile)

      print("the mapping:", self.channelsMap)
      print("the offset: ", self.start)
      print("the sampling rate: ", self.sampling_rate)


   def loadSpindles(self,matching_files,spindleBN,suffix):
      All_Spindle = dict()
      for f in matching_files:
         try:
            print(f)
            structure = re.findall(f'{spindleBN}_(.*){suffix}.csv', f)[0]
            print(structure)
            spindles = pd.read_csv(f, sep=',', header=0, index_col=0)
            if 'toKeep' not in spindles:
               spindles['toKeep'] = True
            print(f"file {f} was found so loading it")
            All_Spindle[structure]=spindles
         except Exception as error:
            print(error)
      return All_Spindle
   