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

from .tools import getPathComponent
from .localConfigurations import localConf

class experiment():
   """experiment is a class for a whole experiment including all its component (NPX, miniscope, intan...)
   """
   def __init__(self) -> None:
      self.config = localConf()
      self.expe_path = Path("")
      self.raw_data_path = ""
      self.interim_analysis_path = ""
      self.parser_fn = 'expeConfig1.ini'
      self.parser = configparser.ConfigParser()
      self.parser.read('defaultExpeConfig.ini') # check if there are modifs to load
      self.data = dict()
      self.expe_info = dict()
      self.num_lfp_channels = 32

      self.remote_prefix = Path(self.config.get('GENERAL','remote_prefix'))
      currentFolder = self.remote_prefix / self.config.get('GENERAL','currentfolder')
      self.loadCurrentFolder(currentFolder)
      
      try:
         fc1 = FileChooser(self.expe_path, select_default=True, show_only_dirs = True, title = "<b>OpenEphys Folder</b>")
         display(fc1)
         fc1.register_callback(self.update_my_expe_choice)
      except Exception as error:
         print(f"something went wrong, make sure the experiment path ({self.expe_path}) is a folder")

   def loadCurrentFolder(self,currentFolder):
      parserName = currentFolder / self.parserFN
      if currentFolder == Path(""):
         print("no folder currently selected")
         if self.config.getboolean('DATA','is_remote'):
            path=self.config.get('DATA','remote_path')
         else:
            path=self.config.get('DATA','local_path')
         self.expe_path = Path(path)
      elif parserName.is_file():
         print(f'current folder {currentFolder} contains a config file')
         self.expe_path = currentFolder
         self.parser.read(parserName)
         self.__dict__.update(dict(self.parser.items('ALL')))
         #TODO: soon remove next line, ,it is for compatibility only
         try:
            self.__dict__.update(ast.literal_eval(self.parser.get('ALL',"expe_info")))
         except:
            pass
         self.updateExpeConfigFile()
         
         self.expe_path = self.remote_prefix / self.interim_analysis_path
         self.config.set('GENERAL','current_folder', self.interim_analysis_path)
         self.config.updateConf()
      else:
         print(f'current folder {currentFolder} does not contain a config file, it must be the raw data folder')
         self.raw_data_path = currentFolder

         self.project_type = self.config.get('ANALYSIS','project_type')
         self.__dict__.update(getPathComponent(self.raw_data_path,self.project_type))

         if self.project_type == 1:
            self.interim_analysis_path = self.remote_prefix / self.analysis_path_root / self.project_id / self.sub_project_id / self.config['ANALYSIS']['interim_analysis_path'] / self.condition_id / self.animal_id / self.recording_id
         else:
            self.interim_analysis_path = self.remote_prefix / self.analysis_path_root / self.project_id / self.sub_project_id / self.config['ANALYSIS']['interim_analysis_path'] / self.animal_id / self.condition_id / self.recording_id
         os.makedirs(self.interim_analysis_path, exist_ok=True)

         self.updateExpeConfigFile()

         self.expe_path = self.interim_analysis_path
         self.config.set('GENERAL','current_folder', self.expe_path)
         self.config.updateConf()


   def updateExpeConfigFile(self):
      if self.interim_analysis_path == '': # for compatibility
         self.interim_analysis_path = self.interimAnalysisPath
      configFN = self.remote_prefix / self.interim_analysis_path / self.parser_fn
      for item in self.__dict__.keys():
         self.parser.set('ALL',item,str(self.__dict__[item]))
      item_to_ignore = ["config","parser","data"]
      for item in item_to_ignore:
         self.parser.remove_option('ALL',item)
      with open(configFN, 'w') as configfile:
         self.parser.write(configfile)
         #print(f"{configFN} saved")

   def setnum_lfp_channels(self,num_lfp_channels):
      self.num_lfp_channels = num_lfp_channels
      self.updateExpeConfigFile()

   # Function to find files containing a specific string
   def find_files_with_string(self, folder_path, search_string):
      matching_files = []
      # Traverse the folder to find files
      for path_object in folder_path.rglob(f"*{search_string}*"):
         matching_files.append(path_object)
      return matching_files
   
   def update_my_expe_choice(self,chooser):
      selection = Path(chooser.selected)
      self.loadCurrentFolder(selection)


   def analyseExpe_findData(self, fullSampling=False):
      """findData: function that analyse the content of the raw data folder and detects component of the experiment to load all of them
      """
      from .ePhy.LFP import IntanLFP, NPX, LFP_DS
      self.loadAnimalInfo()
      DSdata=False
      
      interim_analysis_folder = self.remote_prefix / self.interim_analysis_path
      raw_data_folder = self.remote_prefix / self.raw_data_path

      if self.find_files_with_string(interim_analysis_folder,  "RawDataChannelExtractedDS.npy"): # pre-analysed data
         print('********found some RawDataChannelExtractedDS.npy files********')
         matching_files = self.find_files_with_string(interim_analysis_folder, "RawDataChannelExtractedDS.npy")
         self.data['LFP_DS'] = LFP_DS(self, matching_files)
         DSdata=True

      if self.find_files_with_string(interim_analysis_folder,  f"{self.spindle_bn}_*.csv"): #NPX's Data
         print('********found some Spindles files********')
         matching_files = self.find_files_with_string(interim_analysis_folder, f"{self.spindle_bn}_*.csv")
         All_Spindles = self.loadSpindles(matching_files,self.spindle_bn,self.suffix)
         self.data['Spindles'] = All_Spindles

      if fullSampling or not DSdata:
         if self.find_files_with_string(raw_data_folder,  ".bin"): #Bonsai or IgorPro
            print('********found some .bin files********')
            matching_files = self.find_files_with_string(raw_data_folder, ".bin")
            self.data['OE_LFP'] = IntanLFP(self, matching_files, numChannels=self.num_lfp_channels)
      
         if self.find_files_with_string(raw_data_folder,  "continuous.dat"): #OpenEphys
            print('********found some continuous.dat files********')
            matching_files = self.find_files_with_string(raw_data_folder, "continuous.dat")
            self.data['OE_LFP'] = IntanLFP(self, matching_files, recSyst = 'OpenEphys', numChannels=self.num_lfp_channels)


      if self.find_files_with_string(raw_data_folder,  "NP_spikes_*.raw"): #NPX's Data
         print('********found some NPX files********')
         matching_files = self.find_files_with_string(raw_data_folder, "NP_spikes_*.raw")
         self.data['NPX'] = NPX(self, matching_files)

   def loadAnimalInfo(self):
      """ loads channelMaps.ini located at the animal level of the interimAnalysis tree that contains animal metadata.
      Notably, channelsMap is a dict of brain structures with the corresponding electrode numbers.
      They are labeled with dinstinct status: 0 for unused, 1 for floating electrode, 1 and 2 for differential
      """
      animalConfBN='channelMaps.ini'
      if int(self.project_type) == 0:
         animalConf = self.remote_prefix / self.analysis_path_root / self.project_id / self.sub_project_id / self.config['ANALYSIS']['interim_analysis_path'] / self.condition_id / self.animal_id / animalConfBN
      else:
         animalConf = self.remote_prefix / self.analysis_path_root / self.project_id / self.sub_project_id / self.config['ANALYSIS']['interim_analysis_path'] / self.animal_id / animalConfBN
      animalParser = configparser.ConfigParser()
      animalParser.read(animalConf)
      self.channelsMap = ast.literal_eval(animalParser[self.animal_id]['channelsMap'])
      print("Mapping found and loaded")
      print(self.channelsMap)


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
   