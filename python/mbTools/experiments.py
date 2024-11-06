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
      self.expePath = ""
      self.rawDataPath = ""
      self.interimAnalysisPath = ""
      self.parser = configparser.ConfigParser()
      self.parserFN = 'expeConfig1.ini'
      self.parser.read('defaultExpeConfig.ini') # check if there are modifs to load
      self.data = dict()
      self.expeInfo = dict()
      self.numLFPchannels = 32

      currentFolder = self.config.get('GENERAL','currentfolder')
      self.loadCurrentFolder(currentFolder)
      
      try:
         fc1 = FileChooser(self.expePath, select_default=True, show_only_dirs = True, title = "<b>OpenEphys Folder</b>")
         display(fc1)
         fc1.register_callback(self.update_my_expe_choice)
      except Exception as error:
         print(f"something went wrong, make sure the experiment path ({self.expePath}) is a folder")

   def loadCurrentFolder(self,currentFolder):
      if currentFolder == "":
         print("no folder currently selected")
         if self.config.getboolean('DATA','isremote'):
             path=self.config.get('DATA','remotepath')
         else:
             path=self.config.get('DATA','localpath')
         self.expePath = path
      elif os.path.isfile(os.path.join(currentFolder,self.parserFN)):
         print(f'current folder {currentFolder} contains a config file')
         self.expePath = currentFolder
         parserName = os.path.join(currentFolder,self.parserFN)
         print(parserName)
         self.parser.read(parserName)
         self.rawDataPath = self.parser.get('ALL','rawdatapath')
         self.interimAnalysisPath = self.parser.get('ALL','interimanalysispath')
         self.expeInfo = ast.literal_eval(self.parser.get('ALL','expeinfo'))
         self.numLFPchannels = int(self.parser.get('ALL','numLFPchannels'))
         self.updateExpeConfigFile()
      else:
         print(f'current folder {currentFolder} does not contain a config file, it must be the raw data folder')
         self.rawDataPath = currentFolder

         self.projectType = int(self.config['ANALYSIS']['projectType'])
         self.expeInfo = getPathComponent(self.rawDataPath,self.projectType)

         if self.projectType == 1:
            self.interimAnalysisPath = os.path.join(self.expeInfo['analysisPath'], self.expeInfo['ProjectID'], self.expeInfo['subProjectID'], self.config['ANALYSIS']['interimpath'], self.expeInfo['conditionID'], self.expeInfo['AnimalID'], self.expeInfo['recordingID'])
         else:
            self.interimAnalysisPath = os.path.join(self.expeInfo['analysisPath'], self.expeInfo['ProjectID'], self.expeInfo['subProjectID'], self.config['ANALYSIS']['interimpath'], self.expeInfo['AnimalID'], self.expeInfo['conditionID'], self.expeInfo['recordingID'])
         os.makedirs(self.interimAnalysisPath, exist_ok=True)

         self.parser.set('ALL','rawdatapath', self.rawDataPath)
         self.parser.set('ALL','interimanalysispath', self.interimAnalysisPath)
         self.parser.set('ALL','expeinfo', str(self.expeInfo))

         self.updateExpeConfigFile()

         self.expePath = self.interimAnalysisPath
         self.config.set('GENERAL','currentFolder', self.expePath)
         self.config.updateConf()


   def updateExpeConfigFile(self):
      configFN = os.path.join(self.interimAnalysisPath,self.parserFN)
      with open(configFN, 'w') as configfile:
         self.parser.write(configfile)
         #print(f"{configFN} saved")

   def setNumLFPchannels(self,numLFPchannels):
      self.numLFPchannels = numLFPchannels
      self.parser.set('ALL','numLFPchannels',self.numLFPchannels)
      self.updateExpeConfigFile()

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
      selection = chooser.selected
      self.loadCurrentFolder(selection)


   def analyseExpe_findData(self, fullSampling=False, spindleBN='Spindlesproperties',suffix=''):
      """findData: function that analyse the content of the raw data folder and detects component of the experiment to load all of them
      """
      from .ePhy.LFP import IntanLFP, NPX, LFP_DS
      self.loadAnimalInfo()
      DSdata=False

      if self.find_files_with_string(self.interimAnalysisPath,  "RawDataChannelExtractedDS.npy"): # pre-analysed data
         print('********found some RawDataChannelExtractedDS.npy files********')
         matching_files = self.find_files_with_string(self.interimAnalysisPath, "RawDataChannelExtractedDS.npy")
         self.data['LFP_DS'] = LFP_DS(self, matching_files)
         DSdata=True

      if self.find_files_with_string(self.interimAnalysisPath,  f"{spindleBN}_*.csv"): #NPX's Data
         print('********found some Spindles files********')
         matching_files = self.find_files_with_string(self.interimAnalysisPath, f"{spindleBN}_*.csv")
         All_Spindles = self.loadSpindles(matching_files,spindleBN,suffix)
         self.data['Spindles'] = All_Spindles

      if fullSampling or not DSdata:
         if self.find_files_with_string(self.rawDataPath,  ".bin"): #Bonsai or IgorPro
            print('********found some .bin files********')
            matching_files = self.find_files_with_string(self.rawDataPath, ".bin")
            self.data['OE_LFP'] = IntanLFP(self, matching_files, numChannels=self.numLFPchannels)
      
         if self.find_files_with_string(self.rawDataPath,  "continuous.dat"): #OpenEphys
            print('********found some continuous.dat files********')
            matching_files = self.find_files_with_string(self.rawDataPath, "continuous.dat")
            self.data['OE_LFP'] = IntanLFP(self, matching_files, recSyst = 'OpenEphys', numChannels=self.numLFPchannels)


      if self.find_files_with_string(self.rawDataPath,  "NP_spikes_*.raw"): #NPX's Data
         print('********found some NPX files********')
         matching_files = self.find_files_with_string(self.rawDataPath, "NP_spikes_*.raw")
         self.data['NPX'] = NPX(self, matching_files)

   def loadAnimalInfo(self):
      """ loads channelMaps.ini located at the animal level of the interimAnalysis tree that contains animal metadata.
      Notably, channelsMap is a dict of brain structures with the corresponding electrode numbers.
      They are labeled with dinstinct status: 0 for unused, 1 for floating electrode, 1 and 2 for differential
      """
      animalConfBN='channelMaps.ini'
      if int(self.expeInfo['projectType']) == 0:
         animalConf = os.path.join(self.expeInfo['analysisPath'], self.expeInfo['ProjectID'], self.expeInfo['subProjectID'], self.config['ANALYSIS']['interimpath'], self.expeInfo['conditionID'], self.expeInfo['AnimalID'], animalConfBN)
      else:
         animalConf = os.path.join(self.expeInfo['analysisPath'], self.expeInfo['ProjectID'], self.expeInfo['subProjectID'], self.config['ANALYSIS']['interimpath'], self.expeInfo['AnimalID'], animalConfBN)
      animalParser = configparser.ConfigParser()
      animalParser.read(animalConf)
      self.channelsMap = ast.literal_eval(animalParser[self.expeInfo['AnimalID']]['channelsMap'])
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
   