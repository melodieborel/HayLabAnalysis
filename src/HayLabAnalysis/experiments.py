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
      self.setInstanceVars()

      currentFolder = self.remote_prefix / self.config.get('GENERAL','current_folder')
      self.loadCurrentFolder(currentFolder)
      
      try:
         fc1 = FileChooser(self.expe_path, select_default=True, show_only_dirs = True, title = "<b>OpenEphys Folder</b>")
         display(fc1)
         fc1.register_callback(self.update_my_expe_choice)
      except Exception as error:
         print(f"something went wrong, make sure the experiment path ({self.expe_path}) is a folder")
         
   def setInstanceVars(self):
      print("reseting vars")
      self.expe_path = Path("")
      self.raw_data_path = ""
      self.interim_analysis_path = ""
      self.parser_fn = 'expeConfig1.ini'
      self.parser = configparser.ConfigParser()
      
      self.parser.read('defaultExpeConfig.ini') # check if there are modifs to load
      self.data = dict()
      self.expe_info = dict()
      self.num_lfp_channels = 32
      self.refTime = None
      
      self.remote_prefix = Path(self.config.get('GENERAL','remote_prefix'))

      
      
   def loadCurrentFolder(self,currentFolder):
      parserName = currentFolder / self.parser_fn
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
         self.raw_data_path = str(currentFolder.relative_to(self.remote_prefix))
         
         self.project_type = self.config.get('ANALYSIS','project_type')
         self.__dict__.update(getPathComponent(currentFolder,self.project_type))
                  
         if "iterim_location" in self.config['ANALYSIS']:
            print('a location for interim analysis was provided so using it')
            if self.project_type == 1:
               self.expe_path = self.remote_prefix / self.config['ANALYSIS']['iterim_location'] / self.config['ANALYSIS']['interim_path'] / self.condition_id / self.animal_id / self.recording_id
            else:
               self.expe_path = self.remote_prefix / self.config['ANALYSIS']['iterim_location'] / self.config['ANALYSIS']['interim_path'] / self.animal_id / self.condition_id / self.recording_id
         else:
            print('no location for interim analysis was provided so using default')
            if self.project_type == 1:
               self.expe_path = self.remote_prefix / self.analysis_path_root / self.project_id / self.sub_project_id / self.config['ANALYSIS']['interim_path'] / self.condition_id / self.animal_id / self.recording_id
            else:
               self.expe_path = self.remote_prefix / self.analysis_path_root / self.project_id / self.sub_project_id / self.config['ANALYSIS']['interim_path'] / self.animal_id / self.condition_id / self.recording_id


         self.interim_analysis_path = str(self.expe_path.relative_to(self.remote_prefix))
         print(self.expe_path)
         print(self.interim_analysis_path)
         
         self.expe_path.mkdir(exist_ok=True, parents=True)
         
         self.updateExpeConfigFile()
         
         self.config.set('GENERAL','current_folder', self.interim_analysis_path)
         self.config.updateConf()


   def updateExpeConfigFile(self):
      if type(self.remote_prefix) == str:
         self.remote_prefix = Path(self.remote_prefix)
      configFN = self.remote_prefix / self.interim_analysis_path / self.parser_fn
      for item in self.__dict__.keys():
         self.parser.set('ALL',item,str(self.__dict__[item]))
      item_to_ignore = ["config","parser","data","expe_path","expe_info","channelsMap"]
      for item in item_to_ignore:
         self.parser.remove_option('ALL',item)
      configFN.parent.mkdir(exist_ok=True, parents=True)
      with open(configFN, 'w') as configfile:
         self.parser.write(configfile)
         print(f"{configFN} saved")

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
      self.setInstanceVars()
      self.loadCurrentFolder(selection)


   def analyseExpe_findData(self, fullSampling=False):
      """analyses the content of the raw data folder, detects component of the experiment and add them as entries of the data dictionary

      This function will look for specific file patterns in the raw data folder, and initialize the corresponding data objects based on the data type
      and its acquisition context. Each data object is then added to the data dictionary and contains specific methods for data manipulation and analysis.
      All data objects are based on either ePhy (LFP, NPX) or Imaging (Miniscope, webcam) data types.

      Args:
          fullSampling (bool, optional): imports raw data at their original sampling rate regardless of the presence of downsampled data. Defaults to False.
      """
      from .ePhy.LFP import IntanLFP, NPX, LFP_DS
      from .ePhy.TTL import TTLin
      from .Imaging.Miniscope import Miniscope, Webcam
      self.loadAnimalInfo()
      DSdata=False
      
      interim_analysis_folder = self.remote_prefix / self.interim_analysis_path
      raw_data_folder = self.remote_prefix / self.raw_data_path

      if self.find_files_with_string(raw_data_folder,  "recTS*.csv"): # pre-analysed data
         print('********found recording start timestamp********')
         self.refTime = self.loadTimeStamp(self.find_files_with_string(raw_data_folder,  "recTS*.csv"))
         print(f"recording started on {self.refTime}")

      if self.find_files_with_string(interim_analysis_folder,  "RawDataChannelExtractedDS.npy"): # pre-analysed data
         print('********found some RawDataChannelExtractedDS.npy files********')
         matching_files = self.find_files_with_string(interim_analysis_folder, "RawDataChannelExtractedDS.npy")
         self.data['LFP_DS'] = LFP_DS(self, matching_files, numChannels=self.num_lfp_channels)
         DSdata=True

      if self.find_files_with_string(interim_analysis_folder,  f"{self.spindle_bn}_*.csv"): #NPX's Data
         print('********found some Spindles files********')
         matching_files = self.find_files_with_string(interim_analysis_folder, f"{self.spindle_bn}_*.csv")
         All_Spindles = self.loadSpindles(matching_files,self.spindle_bn,self.suffix)
         self.data['Spindles'] = All_Spindles

      if fullSampling or not DSdata:
         if self.find_files_with_string(raw_data_folder,  "OE_*data*.bin"): #Bonsai or IgorPro
            print('********found some .bin files********')
            matching_files = self.find_files_with_string(raw_data_folder, "OE_*data*.bin")
            print(matching_files)
            self.data['OE_LFP'] = IntanLFP(self, matching_files, numChannels=self.num_lfp_channels)
      
         if self.find_files_with_string(raw_data_folder,  "continuous.dat"): #OpenEphys
            print('********found some continuous.dat files********')
            matching_files = self.find_files_with_string(raw_data_folder, "continuous.dat")
            self.data['OE_LFP'] = IntanLFP(self, matching_files, recSyst = 'OpenEphys', numChannels=self.num_lfp_channels)


      if self.find_files_with_string(raw_data_folder,  "NP_spikes_*.raw"): #NPX's Data
         print('********found some NPX files********')
         matching_files = self.find_files_with_string(raw_data_folder, "NP_spikes_*.raw")
         self.data['NPX'] = NPX(self, matching_files)

      if self.find_files_with_string(raw_data_folder,  "_ttlin*.bin"): #NPX's Data
         print('********found some OE TTL in files********')
         matching_files = self.find_files_with_string(raw_data_folder, "_ttlin*.bin")
         self.data['ttl'] = TTLin(self, matching_files)

      if self.find_files_with_string(raw_data_folder,  "_miniscope_*.avi"): #NPX's Data
         print('********found some Miniscope files********')
         matching_files = self.find_files_with_string(raw_data_folder, "_miniscope_*.avi")
         self.data['miniscope'] = Miniscope(self, matching_files)

      if self.find_files_with_string(raw_data_folder,  "_webcam_*.avi"): #NPX's Data
         print('********found some webcam files********')
         matching_files = self.find_files_with_string(raw_data_folder, "_webcam_*.avi")
         self.data['webcam'] = Webcam(self, matching_files)

   def loadAnimalInfo(self):
      """ loads channelMaps.ini located at the animal level of the interimAnalysis tree that contains animal metadata.
      Notably, channelsMap is a dict of brain structures with the corresponding electrode numbers.
      They are labeled with dinstinct status: 0 for unused, 1 for floating electrode, 1 and 2 for differential
      """
      animalConfBN='channelMaps.ini'
      if int(self.project_type) == 0:
         animalConf = self.remote_prefix / self.analysis_path_root / self.project_id / self.sub_project_id / self.config['ANALYSIS']['interim_path'] / self.condition_id / self.animal_id / animalConfBN
      else:
         animalConf = self.remote_prefix / self.analysis_path_root / self.project_id / self.sub_project_id / self.config['ANALYSIS']['interim_path'] / self.animal_id / animalConfBN
      try:
         animalParser = configparser.ConfigParser()
         animalParser.read(animalConf)
         self.channelsMap = ast.literal_eval(animalParser[self.animal_id]['channelsMap'])
         print("Mapping found and loaded")
         print(self.channelsMap)
      except:
         print(f"could not load channel map. Please make sure the animalID {self.animal_id} is mapped in the file {animalConf}")
         self.channelsMap = None

   def loadSpindles(self,matching_files,spindleBN,suffix):
      All_Spindle = dict()
      for f in matching_files:
         try:
            print(f)
            structure = re.findall(f'{spindleBN}_(.*){suffix}.csv', str(f))[0]
            print(structure)
            spindles = pd.read_csv(f, sep=',', header=0, index_col=0)
            if 'toKeep' not in spindles:
               spindles['toKeep'] = True
            print(f"file {f} was found so loading it")
            All_Spindle[structure]=spindles
         except Exception as error:
            print(error)
      return All_Spindle
   
   def loadTimeStamp(self, file):
      import pandas as pd
      df = pd.read_csv(file[0], sep=',', header=None, names=['ts'])
      df.ts = pd.to_datetime(df.ts, format='ISO8601')
      return df.ts[0]