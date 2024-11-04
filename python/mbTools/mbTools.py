import os
import configparser
import pickle
import ast
import fnmatch
from pathlib import Path

import numpy as np
import pandas as pd

import re

from ipyfilechooser import FileChooser
import ipywidgets as widgets
from IPython.display import display
from IPython import get_ipython
import IPython

import pprint

class color:
   PURPLE = '\033[95m'
   CYAN = '\033[96m'
   DARKCYAN = '\033[36m'
   BLUE = '\033[94m'
   GREEN = '\033[92m'
   YELLOW = '\033[93m'
   RED = '\033[91m'
   BOLD = '\033[1m'
   UNDERLINE = '\033[4m'
   END = '\033[0m'

def superCleanPlot(lfp, npx, canauxLFP=None, structureLFP=None, canauxNPX=[0,1], time=0, pre=1, post=4, offset=0, scaleNPX=1, scaleLFP=1):
   """
   superCleanPlot plots very precisely aligned NPX and LFP data

   :lfp: the lfp object containing signal, times, sampling_rate...
   :npx: the npx object containing all infos as well
   :canauxLFP: (facultative, None by default) array of int that indicates the channels to display. StructureLFP will be ignored if canauxLFP is not None
   :structureLFP: (facultative, None by default) array of strings that indicates the brain structures defined by mapping to display. Will be ignored if canauxLFP is not None
   :canauxNPX: ((facultative, default is [0,1]) array of ints that indicates the channels to display. For ann interval, please enter np.arange(start, stop).
   :time: the time to center on display in seconds
   :pre: duration (float in s) before time of interest to display 
   :post: duration (float in s) before time of interest to display 
   """
   import matplotlib.pyplot as plt
   plt.close()

   #offset=51.51#51.4576900#t_start['LFP']#52.6734#52.68
   # perfect align manual offset=51.5146156977
   #offset=51.52262754

   idx=find_nearest(lfp.times, time)
   print(lfp.times[idx])
   print(idx)

   x=lfp.times[idx-int(pre*lfp.sampling_rate):idx+int(post*lfp.sampling_rate)]
   if canauxLFP is not None:
      y=lfp.signal[idx-int(pre*lfp.sampling_rate):idx+int(post*lfp.sampling_rate),canauxLFP]
   elif structureLFP is not None:
      y=lfp.combineStructures(structureLFP)[idx-int(pre*lfp.sampling_rate):idx+int(post*lfp.sampling_rate),:]
   else:
      y=lfp.signal[idx-int(pre*lfp.sampling_rate):idx+int(post*lfp.sampling_rate),:]


   print(npx.times)
   idx2=find_nearest(npx.times-npx.times[0], time)
   print(idx2)
   x2=npx.times[idx2-int(pre*npx.sampling_rate):idx2+int(post*npx.sampling_rate)]-npx.times[0]
   y2=npx.signal['spike'].select_channels(canauxNPX).get_traces(start_frame=idx2-int(pre*npx.sampling_rate), end_frame=idx2+int(post*npx.sampling_rate), return_scaled=False)

   plt.plot(x, y*scaleLFP,'-')
   plt.plot(x2, y2*scaleNPX+offset,'-')
   plt.show()

def convertTheoricIndex2realTime(thIdx,realFreq=1, offset=0):
    realTime=thIdx/realFreq + offset
    return realTime

def find_nearest(array, value):
    idx = (np.abs(array - value)).argmin()
    return idx

class localConf(configparser.ConfigParser):
   def __init__(self, configFN = 'localConfig.ini') -> None:
      super().__init__()
      self.configFN = configFN
      self.read('defaultLocalConfig.ini') # check if there are modifs to load

      if os.path.isfile(self.configFN):
         self.read(self.configFN)
         print(f'Local config file loaded from {configFN}')
      else:
         self.set('DATA', 'localPath', os.path.expanduser("~"))
         self.set('ANALYSIS', 'interimPath', 'interimAnalysis')
         print(f'Local config file did not exist, it was successfully created at {configFN}')
      
      with open(self.configFN, 'w') as configfile:
         self.write(configfile)
      print(f'Local config file updated')

   def completeConf(self):
      """maybe should add the possibility to ensure all parts of the config is there"""
      pass

   def updateConf(self):
      with open(self.configFN, 'w') as configfile:
         self.write(configfile)
   
   def getProjects(self):
         return [p.split('.')[0] for p in self.sections() if p not in ['DATA','ANALYSIS']]
   
   def getSubProjects(self, projectID):
      return [p.split('.')[1] for p in self.sections() if p.split('.')[0]==projectID]
   
   def getSubProjectsWidget(self):
      return {p.split('.')[0]: widgets.Dropdown(
         options=[p.split('.')[1]],
         description='Sub-project (you can update the list in your localConfig.ini file):',
         ) for p in self.sections() if p not in ['DATA','ANALYSIS','GENERAL']} 
   
   def printAll(self):
      for section in self.sections():
         print(section)
         pprint.pp(self.items(section))

class expeConfigDict(dict):
   def __init__(self, config: localConf = None) -> None:
      super().__init__()
      if config is None:
         self.config = localConf()
      else:
         self.config = config
      self.expePath = self.config.get('GENERAL','currentFile', fallback="")
      self.rawDataPath = ""
      self.interimAnalysisPath = ""
      self.projectType = None
      self.expeInfo = dict()
      self.iWidget = None
      self.parser = configparser.ConfigParser()

      if self.expePath is not None and os.path.isfile(self.expePath): # a file is currently being used
         print(f"the file is {self.expePath}")
         self.pathName, self.fileName = os.path.split(self.expePath)
         self.loadExpeConfigDict()
      else:
         print(f"the file {self.expePath} was not found")
         self.pathName = self.rawDataPath
         self.fileName = ""

      fc = FileChooser(path=self.pathName, filename='', select_default=True, show_only_dirs = False, title = "<b>Select file</b>")
      display(fc)
         
      # Register callback function
      fc.register_callback(self.update_my_expe_choice)

      self.analysisPath = self.config['ANALYSIS']['interimPath']
      self.projectType = int(self.config['ANALYSIS']['projectType'])
      self.ProjectID = self.config['ANALYSIS']['ProjectID']
      self.subProjectID = self.config['ANALYSIS']['subProjectID']
      self.conditionID = self.config['ANALYSIS']['conditionID']
      self.AnimalID = int(self.config['ANALYSIS']['AnimalID'])
      self.recordingID = int(self.config['ANALYSIS']['recordingID'])

      self.constructWidgets()


   def constructWidgets(self):
      self.wProject = widgets.Dropdown(
         options=self.config.getProjects(),
         value=self.ProjectID,
         description='Project (you can update the list in your localConfig.ini file):',
         disabled=False,
      )
      self.wProject.observe(self.updateProject, 'value')

      #wSubProject.observe(updateSubProject)

      designs = ['independant groups', 'within subject']
      self.wDesign = widgets.RadioButtons(
         options=designs,
         value=designs[self.projectType], # Defaults to 'independant groups'
         description='Experiment design:'
      )
      self.wDesign.observe(self.update_design, names=['index'])

      self.wAnimal = widgets.BoundedIntText(
         value=self.AnimalID,
         min=0,
         #max=10,
         step=1,
         description='Animal ID:'
      )

      conditions = ast.literal_eval(self.config["{}.{}".format(self.ProjectID, self.subProjectID)]['conditions'])
      self.wCondition = widgets.Dropdown(
         options=conditions,
         value=self.conditionID,
         description='Condition:',
      )

      self.wRec = widgets.BoundedIntText(
         value=self.recordingID,
         min=0,
         #max=10,
         step=1,
         description='Recording ID:'
      )

      self.wValidateBtn = widgets.Button(description="Validate")
      self.wValidateBtn.on_click(self.validate_dict)

      self.wOutput = widgets.Output()

      if self.projectType == 0:
         self.iWidget = self.iWidgetConstructor()
      else:
         self.iWidget = self.iWidgetConstructor()

   def generateExpeConfigParser(self, expePath, rawDataPath = None):
      print("generating")
      print(self.expePath)
      self.expePath = expePath
      self.rawDataPath = rawDataPath

      self.projectType = int(self.config['ANALYSIS']['projectType'])
      self.expeInfo = getPathComponent(self.rawDataPath,self.projectType)
      if 'ALL' not in self.parser.sections():
         self.parser.add_section('ALL')
      self.parser.set('ALL','rawDataPath',self.rawDataPath)
      self.parser.set('ALL','expeInfo',str(self.expeInfo))
      self.parser.set('ALL','interimAnalysisPath',os.split(self.expePath[0]))
      
      with open(self.expePath, 'w') as configfile:
         self.parser.write(configfile)

   def loadExpeConfigDict(self, expePath = None):
      if expePath is None:
         expePath = self.expePath
      self.parser.read(expePath)

      self.rawDataPath = self.parser.get('ALL','rawDataPath')
      if 'expeInfo' in self.parser:
         self.expeInfo = self.parser.get('ALL','expeInfo')
      else:
         self.expeInfo = getPathComponent(expePath,self.projectType)
         self.parser.set('ALL','expeInfo',str(self.expeInfo))
         self.updateExpeConfigDict(expePath)
      if 'interimAnalysisPath' in self.parser:
         self.interimAnalysisPath = self.parser.get('ALL','interimAnalysisPath')
      else:
         self.interimAnalysisPath = os.path.split(expePath)[0]
         self.parser.set('ALL','interimAnalysisPath',self.interimAnalysisPath)
         self.updateExpeConfigDict(expePath)


   def updateExpeConfigDict(self, configFN):
      with open(configFN, 'w') as configfile:
         self.parser.write(configfile)

   def printExpeInfo(self, **func_kwargs):
      pass#print(expeInfo)

   def updateProject(self, widget):
      self.ProjectID = widget.new
      self.expeInfo['ProjectID'] = self.ProjectID
      if self.subProjectID not in self.config.getSubProjects(self.ProjectID):
         self.subProjectID = self.config.getSubProjects(self.ProjectID)[0]
      new_i = self.iWidgetConstructor()
      self.iWidget.children = new_i.children

   def iWidgetConstructor(self):
      if self.projectType == 0:
         iWidget = widgets.interactive(self.printExpeInfo, project=self.wProject, subProject=self.config.getSubProjectsWidget()[self.ProjectID], design=self.wDesign, condition=self.wCondition, animal=self.wAnimal, rec=self.wRec)
      else:
         iWidget = widgets.interactive(self.printExpeInfo, project=self.wProject, subProject=self.config.getSubProjectsWidget()[self.ProjectID], design=self.wDesign, animal=self.wAnimal, condition=self.wCondition, rec=self.wRec)
      return iWidget
   

   def validate_dict(self,b):
      self.AnimalID = self.wAnimal.value
      self.recordingID = self.wRec.value
      self.conditionID = self.wCondition.value
      with self.wOutput:
         print("Button clicked.")
         print(self.AnimalID)
      if self.projectType == 0:
         path = os.path.join(self.config['ANALYSIS']['path'], self.ProjectID, self.subProjectID, self.conditionID, str(self.AnimalID), str(self.recordingID))
      else:
         path = os.path.join(self.config['ANALYSIS']['path'], self.ProjectID, self.subProjectID, str(self.AnimalID), self.conditionID, str(self.recordingID))
      os.makedirs(path, exist_ok=True)
      currentFile = os.path.join(path,'saved_dictionary.ini')
      self.generateExpeConfigParser(currentFile, rawDataPath = self.rawDataPath)
      self.loadExpeConfigDict(expePath = currentFile)
      self.config.set('GENERAL','currentFile', currentFile)
      self.config.updateConf()

   def saveInterimAnalysisFolder(self):
      self.projectType = self.expeInfo['projectType']
      self.ProjectID = self.expeInfo['ProjectID']
      self.subProjectID = self.expeInfo['subProjectID']
      self.AnimalID = self.expeInfo['AnimalID']
      self.recordingID = self.expeInfo['recordingID']
      self.conditionID = self.expeInfo['conditionID']
      print(self.expePath)
      if self.projectType == 0:
         self.interimAnalysisPath = os.path.join(self.expeInfo['analysisPath'], self.ProjectID, self.subProjectID, self.config['ANALYSIS']['interimpath'], self.conditionID, str(self.AnimalID), str(self.recordingID))
      else:
         self.interimAnalysisPath = os.path.join(self.expeInfo['analysisPath'], self.ProjectID, self.subProjectID, self.config['ANALYSIS']['interimpath'], str(self.AnimalID), self.conditionID, str(self.recordingID))
      os.makedirs(self.interimAnalysisPath, exist_ok=True)
      currentFile = os.path.join(self.interimAnalysisPath,'saved_dictionary.ini')
      self.generateExpeConfigParser(currentFile, rawDataPath = self.rawDataPath)
      #self.loadExpeConfigDict(expePath = currentFile)
      self.config.set('GENERAL','currentFile', currentFile)
      self.config.updateConf()

   def updateSubProject(self, widget):
      if widget['type'] == 'change' and widget['name'] == 'value':
         self.expeInfo['subProjectID'] = widget.new

   def update_design(self, widget):
      projectType = widget.new
      new_i = self.iWidgetConstructor()
      self.iWidget.children = new_i.children
      magicstore('projectType', projectType)

   def update_my_expe_choice(self,chooser):
      print("this is a config file so we are loading it")
      selection = chooser.selected
      if selection.endswith("ini") and os.path.isfile(selection):
         currentFile = str(selection)
         self.interimAnalysisPath,_ = os.path.split(selection)
         self.loadExpeConfigDict(expePath = selection)
      else:
         print("this is not a config file and we should deal with that")
         self.rawDataPath = selection
         if True:
            self.expeInfo = getPathComponent(self.rawDataPath,self.projectType)
            self.saveInterimAnalysisFolder()
            print("get path automatically")
         else:
            display(self.iWidget)
            display(self.wValidateBtn, self.wOutput)

   def rawDataSelector(self):
      if self.rawDataPath is not None and os.path.isdir(self.rawDataPath):
         rfc = FileChooser(path=self.rawDataPath,select_default=True, show_only_dirs = True, title = "<b>Expe data folder</b>")
      else:
         rfc = FileChooser(show_only_dirs = True, title = "<b>Expe data folder</b>")
      display(rfc)
      # Register callback function
      rfc.register_callback(self.update_rawDataPath)

   def update_rawDataPath(self,chooser):
      self.rawDataPath = chooser.selected
      self.updateExpeConfigDict('rawDataPath', self.rawDataPath)

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

   def defineMap(self,structure,channels: dict):
      self.channelsMap[structure]= channels
   
   def generateChannelsMap(self):
      self.channelsMap = dict( \
         EMG = [dict(canal = 6, status=1)],
         PFC = [dict(canal = 5, status=1),
            dict(canal = 4, status=2)
            ],
         CA1 = [dict(canal = 8, status=1),
            dict(canal = 0, status=0),
            dict(canal = 1, status=0),
            ],
         TTL = [dict(canal = 10, status=1)],
      )
   
   def combineStructures(self,structures, start = 0, end = None):
      if end is None:
         end = self.All.shape[0]
      combined = np.empty((end-start,0),np.int16)
      self.channelLabels = []
      if structures is None:
         #self.generateChannelsMap()
         combined = self.All
         self.channelLabels = [i for i in range(self.All.shape[0])]
      else:
         for region in structures:
            print(region, "->", self.channelsMap[region])
            if len([canal["canal"] for canal in self.channelsMap[region] if canal["status"]==2])>0:
               c2 = [canal["canal"] for canal in self.channelsMap[region] if canal["status"]==2][0]
               c1 = [canal["canal"] for canal in self.channelsMap[region] if canal["status"]==1][0]
               print("Getting differential signal of channel {} - channel {} for {}".format(c2,c1,region))
               self.channelLabels.append(region)
               combined = np.append(combined, self.All[start:end, c2, np.newaxis] - self.All[:, c1, np.newaxis], axis=1)
            elif len([canal["canal"] for canal in self.channelsMap[region] if canal["status"]==1])>0:
               c = [canal["canal"] for canal in self.channelsMap[region] if canal["status"]==1][0]
               print("Getting floating signal of channel {} for {}".format(c,region))
               combined = np.append(combined, self.All[start:end,c, np.newaxis], axis=1)
               self.channelLabels.append(region)
      return combined

   def analyseExpe_findData(self, fullSampling=False, spindleBN='Spindlesproperties',suffix=''):
      """findData: function that analyse the content of the raw data folder and detects component of the experiment to load all of them
      """
      from .ePhy.LFP import IntanLFP, NPX, LFP_DS
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


   def loadSpindles(self,matching_files,spindleBN,suffix):
      All_Spindle = dict()
      for f in matching_files:
         try:
            print(f)
            structure = re.findall(f'{spindleBN}_(.*).csv', f)[0]
            print(structure)
            spindles = pd.read_csv(f, sep=',', header=0, index_col=0)
            if 'toKeep' not in spindles:
               spindles['toKeep'] = True
            print(f"file {f} was found so loading it")
            All_Spindle[structure]=spindles
         except Exception as error:
            print(error)
      return All_Spindle
   
   def loadLFP(self):
      folderpath = Path(self.expePath)
      isNPY = False

      if self.find_files_with_string(folderpath,  ".bin"): #Bonsai or IgorPro
         self.recSyst = "Bonsai" #"IgorPro"
         matching_files = self.find_files_with_string(folderpath, ".bin")
         self.fileType="OE32channels.bin"
         print('found some .bin files')
      elif self.find_files_with_string(folderpath,  "continuous.dat"): #OpenEphys
         self.recSyst = "OpenEphys" #"IgorPro"
         self.fileType="OE32channels.bin"
         matching_files = self.find_files_with_string(folderpath, "continuous.dat")
         print('found some .dat files')
      elif self.find_files_with_string(folderpath,  "RawDataChannelExtractedDS.npy"): #OpenEphys Gaelle's Data
         isNPY = True
         self.fileType="NPY"
         matching_files = self.find_files_with_string(folderpath, "RawDataChannelExtractedDS.npy")
      else:
         raise Exception(f"Couldn't find any .bin or .dat file. Please check your path : {folderpath}")

      if 'sommeil' in self.expePath:
         self.numChannels=64

      match self.recSyst:
         case "Bonsai":
            dtype=np.uint16
            offset=int(np.iinfo(dtype).max/2)
         case _:
            dtype=np.int16
            offset=0

      if isNPY==True:
         if len(matching_files) == 1:
            file = matching_files[0]
            All = np.load(file, mmap_mode= 'r')
            self.numChannels = All.shape[1]
         else:
            raise Exception(f"Several npy files ({len(matching_files)}) to merge but this is not in place yet")
      else:
         All = np.zeros([0], dtype=dtype)
         for file in matching_files:
            print("importing {}".format(file))
            file_data = np.fromfile(file, dtype=dtype)
            All=np.append(All, file_data, axis=0)
         if offset != 0:
            All = All - offset
         if All.dtype is not np.dtype(np.int16):
            All = All.astype(np.int16)
         All = All.reshape(-1,self.numChannels)

      print(f'{self.fileType} file loaded, with {self.numChannels} channels and {All.shape[0]} datapoint')
      self.All = All
      return self.All

   def loadRecording_TimeStamps(self):
      print('warning: this function is not confirmed yet. Please see with MB, that it is porperly working when you want to use it')
      if True:
         folder = Path('.').absolute()
         print(folder)
         path_list_ERS = []
         Ephys_rec_stamps = {}

         for file in folder.glob('**/*continuous/Rhythm_FPGA-112.0/timestamps.npy'):
            path_list_ERS.append(file)

         for file_path in folder.glob('**/*continuous/Rhythm_FPGA-112.0/timestamps.npy'):
            recording = file_path.parents[2].stem
            arr = np.load(file_path)
            Ephys_rec_stamps[recording] = arr
         ## Here the timestamps are stored in a dict which is not necessarily what I want, maybe will need to amend that
      else:
         # Not working yet as synchronised timestamps and timestamps for recording 2 are of different size.
         TTL_stamp2 = []
         for file_path in folder.glob('**/*.npy'):
            subfolder = file_path.parents[1].stem
            if subfolder == 'continuous':
               recording = file_path.parents[2].stem.replace('recording','')
               print(recording)
               file = file_path.stem
               print(recording, file)
               np_arr = np.load(file_path)
               datalen = len(np_arr)
               print(file, datalen)
               if recording == 1: #not in TTL_stamp2:
                     TTL_stamp2.append(recording)
                     coords = {
                        'channels' : np.array(['synchronized_timestamps', 'timestamps']),
                        'duration_rec' : np.arange(datalen)
                     }
                     globals()[f"StampsCont_{recording}"] = xr.DataArray(coords=coords, dims=['channels', 'duration_rec'])
               globals()[f"StampsCont_{recording}"].loc[file,:] = np_arr   

   def loadTTL_TimeStamps(self):
      print('warning: this function is not confirmed yet. Please see with MB, that it is porperly working when you want to use it')
      folder = Path('.').absolute()
      print(folder)
      TTL_stamps = []
      list_recordings = []
      for file_path in folder.glob('**/*.npy'):
         subfolder = file_path.parents[0].stem
         if subfolder == 'TTL_1':
            recording = file_path.parents[3].stem.replace('recording','')
            file = file_path.stem
            np_arr = np.load(file_path)
            datalen = len(np_arr)
            if recording not in TTL_stamps:
                  TTL_stamps.append(recording)
                  list_recordings.append(file_path.parents[3].stem)
      return TTL_stamps

def magicstore(stored_var, value):
   # myvar will contain the variable previously stored with "%store test"
   myvar_filename = get_ipython().ipython_dir + '/profile_default/db/autorestore/' + stored_var
   with open(myvar_filename, 'wb') as f:
      pickle.dump(value,f)

def getPathComponent(filename,projectType):
   if not os.path.isdir(filename):
      filename = os.path.split(os.path.normpath(filename))[0]
   dirPathComponents = os.path.normpath(filename).split(os.sep)
   expeInfo = dict()

   expeInfo['analysisPath'] = os.path.sep.join([*dirPathComponents[0:-5]])
   expeInfo['ProjectID'] = dirPathComponents[-5]
   expeInfo['subProjectID'] = dirPathComponents[-4]

   projectConfig = os.path.sep.join([*dirPathComponents[0:-3],'projectConfig.ini'])
   projParser = configparser.ConfigParser()
   if os.path.isfile(projectConfig):
      projParser.read(projectConfig)
      expeInfo['projectType'] = projParser.get('ALL','projectType')
   else:
      projParser.add_section('ALL')
      projParser.set('ALL','projectType',str(projectType))
      with open(projectConfig, 'w') as configfile:
         projParser.write(configfile)

   if expeInfo['projectType'] == 0:
      expeInfo['conditionID'] = dirPathComponents[-3]
      expeInfo['AnimalID'] = dirPathComponents[-2]
   else:
      expeInfo['AnimalID'] = dirPathComponents[-3]
      expeInfo['conditionID'] = dirPathComponents[-2]
      
   expeInfo['recordingID'] = dirPathComponents[-1]

   return expeInfo