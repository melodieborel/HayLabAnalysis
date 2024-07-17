import os
import configparser
import pickle
import ast
import fnmatch
from pathlib import Path

import numpy as np

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


class localConf(configparser.ConfigParser):
   def __init__(self, configFN = 'localConfig.ini') -> None:
      super().__init__()
      self.configFN = configFN

      if os.path.isfile(self.configFN):
         self.read(self.configFN)
         print(f'Local config file loaded from {configFN}')
      else:
         self.generateLocalConfigFile(self.configFN)
         print(f'Local config file did not exist, it was successfully created at {configFN}')

   def generateLocalConfigFile(self,configFN):
      self['DATA'] = {'path': os.path.expanduser("~")}
      self['ANALYSIS'] = {
         'path': os.path.join(os.path.expanduser("~"),'Analysis'),
         'projecttype': 0,
         'animalid': 0,
         'projectid': 'AProject',
         'subprojectid': 'OneOfItsSubProject',
         'conditionid': 'control',
         'recordingID': 0,
         'suffix': ''
         }
      self['AProject.OneOfItsSubProject'] = {
         'design': 0,
         'nAnimal': 6,
         'conditions': ["control"],
         'nrecordings': 1
         }
      with open(configFN, 'w') as configfile:
         self.write(configfile)

   def updateConf(self):
      with open(self.configFN, 'w') as configfile:
         self.write(configfile)
   
   def getProjects(self):
         return [p.split('.')[0] for p in self.sections() if p not in ['DATA','ANALYSIS']]
   
   def getSubProjects(self):
      return {p.split('.')[0]: widgets.Dropdown(
         options=[p.split('.')[1]],
         description='Sub-project (you can update the list in your localConfig.ini file):',
         ) for p in self.sections() if p not in ['DATA','ANALYSIS']} 
   
   def printAll(self):
      for section in self.sections():
         print(section)
         pprint.pp(self.items(section))

class expeConfigDict(dict):
   def __init__(self, expePath = None) -> None:
      super().__init__()
      self.config = localConf()
      self.expePath = expePath
      self.rawDataPath = ""
      self.numChanels = None
      self.channelsMap = dict()
      self.projectType = None
      self.expeInfo = dict()
      self.iWidget = None

      if self.expePath is not None and os.path.isfile(self.expePath): # a file is currently being used
         self.pathName, self.fileName = os.path.split(self.expePath)
         self.loadExpeConfigDict()
      else:
         self.pathName = self.rawDataPath
         self.fileName = ""

      fc = FileChooser(path=self.pathName, filename=self.fileName, select_default=True, show_only_dirs = False, title = "<b>Select file</b>")
      display(fc)
         
      # Register callback function
      fc.register_callback(self.update_my_expe_choice)

      self.analysisPath = self.config['ANALYSIS']['path']
      self.projectType = int(self.config['ANALYSIS']['projectType'])
      self.ProjectID = self.config['ANALYSIS']['ProjectID']
      self.subProjectID = self.config['ANALYSIS']['subProjectID']
      self.conditionID = self.config['ANALYSIS']['conditionID']
      self.AnimalID = int(self.config['ANALYSIS']['AnimalID'])
      self.recordingID = int(self.config['ANALYSIS']['recordingID'])

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

      if self.projectType == 0:
         self.iWidget = widgets.interactive(self.printExpeInfo, project=self.wProject, subProject=self.config.getSubProjects()[self.ProjectID], design=self.wDesign, condition=self.wCondition, animal=self.wAnimal, rec=self.wRec)
      else:
         self.iWidget = widgets.interactive(self.printExpeInfo, project=self.wProject, subProject=self.config.getSubProjects()[self.ProjectID], design=self.wDesign, animal=self.wAnimal, condition=self.wCondition, rec=self.wRec)

   def generateExpeConfigDict(self, expePath, rawDataPath = None):
    self.expePath = expePath
    self.rawDataPath = rawDataPath
    # TODO: idéallement, les infos devraient être demandées au fur et à mesur ou bien récupérées depuis le fichier d'aurélie
    self.numChanels=64

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
    self.projectType = int(self.config['ANALYSIS']['projectType'])
    self.expeInfo = getPathComponent(self.expePath,self.projectType)

    allParamsDict = dict(channelsMap = self.channelsMap, numChanels = self.numChanels, rawDataPath = self.rawDataPath, expeInfo = self.expeInfo)

    with open(self.expePath, 'wb') as f:
        pickle.dump(allParamsDict, f)

   def loadExpeConfigDict(self, expePath = None):
    if expePath is None:
       expePath = self.expePath
    with open(expePath, 'rb') as f:
        loaded_dict = pickle.load(f, encoding='UTF8')
        self.numChanels = loaded_dict['numChanels']
        self.rawDataPath = loaded_dict['rawDataPath']
        if 'expeInfo' in loaded_dict:
            self.expeInfo = loaded_dict['expeInfo']
        else:
            self.expeInfo = getPathComponent(expePath,self.projectType)
            self.updateExpeConfigDict(expePath,'expeInfo',self.expeInfo)
        self.channelsMap = loaded_dict['channelsMap']

   def updateExpeConfigDict(self, key, value):
      with open(self.expePath, 'rb') as f:
         loaded_dict = pickle.load(f)
      with open(self.expePath, 'wb') as f:
         loaded_dict[key] = value
         pickle.dump(loaded_dict,f)


   def printExpeInfo(self, **func_kwargs):
      pass#print(expeInfo)

   def updateProject(self, widget):
      ProjectID = widget.new
      self.expeInfo['ProjectID'] = ProjectID
      new_i = widgets.interactive(self.printExpeInfo, project=self.wProject, subProject=self.config.getSubProjects()[ProjectID])
      self.iWidget.children = new_i.children


   def updateSubProject(self, widget):
      if widget['type'] == 'change' and widget['name'] == 'value':
         self.expeInfo['subProjectID'] = widget.new

   def update_design(self, widget):
      projectType = widget.new
      if projectType == 0:
         new_i = widgets.interactive(self.printExpeInfo, project=self.wProject, subProject=self.config.getSubProjects()[self.ProjectID], design=self.wDesign, condition=self.wCondition, animal=self.wAnimal, rec=self.wRec)
      else:
         new_i = widgets.interactive(self.printExpeInfo, project=self.wProject, subProject=self.config.getSubProjects()[self.ProjectID], design=self.wDesign, animal=self.wAnimal, condition=self.wCondition, rec=self.wRec)
      self.iWidget.children = new_i.children
      #%store projectType

   def update_my_expe_choice(self,chooser):
    selection = chooser.selected
    if selection.endswith("pkl"):
        currentFile = str(selection)
        self.loadExpeConfigDict(expePath = selection)
    else:
        print("this is not a config file and we should deal with that")
        display(self.iWidget)
        if self.projectType == 0:
            path = os.path.join(self.config['ANALYSIS']['path'], self.ProjectID, self.subProjectID, self.conditionID, str(self.AnimalID), str(self.recordingID))
            
        else:
            path = os.path.join(self.config['ANALYSIS']['path'], self.ProjectID, self.subProjectID, str(self.AnimalID), self.conditionID, str(self.recordingID))
        os.makedirs(path, exist_ok=True)
        currentFile = os.path.join(os.path.split(path)[0],'saved_dictionary.pkl')
        self.generateExpeConfigDict(currentFile, rawDataPath = selection)
        self.loadExpeConfigDict(expePath = currentFile)
    magicstore('currentFile', currentFile)

   def rawDataSelector(self):
      #print(rawDataPath)
      if self.rawDataPath is not None:
         rawDirname, rawFN = os.path.split(self.rawDataPath)
         rfc = FileChooser(path=rawDirname, filename=rawFN,select_default=True, show_only_dirs = False, title = "<b>ePhys data</b>")
      else:
         rfc = FileChooser(show_only_dirs = False, title = "<b>ePhys data</b>")
      display(rfc)
      # Register callback function
      rfc.register_callback(self.update_rawDataPath)

   def update_rawDataPath(self,chooser):
      self.rawDataPath = chooser.selected
      self.updateExpeConfigDict('rawDataPath', self.rawDataPath)

class experiment():
   def __init__(self, expePath = None) -> None:
      self.expePath = expePath
      self.numChanels = None
      self.channelsMap = dict()

      fc1 = FileChooser(expePath,select_default=True, show_only_dirs = True, title = "<b>OpenEphys Folder</b>")
      display(fc1)
      fc1.register_callback(self.update_my_expe_choice)

   # Function to find files containing a specific string
   def find_files_with_string(self, folder_path, search_string):
      matching_file = []
      # Traverse the folder to find files
      for root, _, files in os.walk(folder_path):
         for file in files:
               if fnmatch.fnmatch(file, f"*{search_string}*"):
                  matching_file=os.path.join(root, file)
      return matching_file
   
   def update_my_expe_choice(self,chooser):
    dpath = chooser.selected
    magicstore('dpath', dpath)
   
   def loadLFP(self):
      folderpath = Path(self.expePath)
      if self.find_files_with_string(folderpath,  ".bin"): #Bonsai
         matching_file = self.find_files_with_string(folderpath, ".bin")
         All = np.fromfile(matching_file, dtype=np.uint16)
         All = All - int(65535/2)
         All = All.astype(np.int16)
         if 'sommeil' in self.expePath:
            self.numchannels=64
         All = All.reshape(-1,self.numchannels)
         print(f'File loaded: OE32channels.bin, {self.numchannels} channels')
         
      elif self.find_files_with_string(folderpath,  "RawDataChannelExtractedDS.npy"): #OpenEphys Gaelle's Data
         matching_file = self.find_files_with_string(folderpath, "RawDataChannelExtractedDS.npy")
         All = np.load(matching_file, mmap_mode= 'r')
         self.numchannels=32
         #All = All.reshape(-1,numchannels)
         print(f'File loaded: RawDataChannelExtractedDS.npy, {self.numchannels} channels')


      elif self.find_files_with_string(folderpath,  "continuous.dat"): #OpenEphys
         matching_file = self.find_files_with_string(folderpath, "continuous.dat")
         All = np.fromfile(matching_file, dtype=np.int16)

         #All =np.memmap(matching_file, dtype='int16', mode='c')
         
         #All = All - int(65535/2)
         #All = All.astype(np.int16)    
         if 'sommeil' in self.expePath:
            self.numchannels=64
         All = All.reshape(-1,self.numchannels)
         print(f'File loaded: continuous.dat, {self.numchannels} channels')

      return All


def magicstore(stored_var, value):
   # myvar will contain the variable previously stored with "%store test"
   myvar_filename = get_ipython().ipython_dir + '/profile_default/db/autorestore/' + stored_var
   with open(myvar_filename, 'wb') as f:
      pickle.dump(value,f)

def getPathComponent(filename,projectType):
    
    dirPathComponents = os.path.normpath(filename).split(os.sep)
    expeInfo = dict()

    expeInfo['analysisPath'] = os.path.join('/',*dirPathComponents[0:-5])
    expeInfo['ProjectID'] = dirPathComponents[-5]
    expeInfo['subProjectID'] = dirPathComponents[-4]

    projectConfig = os.path.join('/',*dirPathComponents[0:-3],'projectConfig.pkl')
    if os.path.isfile(projectConfig):
        with open(projectConfig, 'rb') as f:
            loaded_dict = pickle.load(f)
            expeInfo['projectType'] = loaded_dict['projectType']
    else:
        with open(projectConfig, 'wb') as f:
            projDict = dict(projectType = projectType)
            pickle.dump(projDict, f)
            print('Project config dict created')

    if projectType == 0:
        expeInfo['conditionID'] = dirPathComponents[-3]
        expeInfo['AnimalID'] = dirPathComponents[-2]
    else:
        expeInfo['AnimalID'] = dirPathComponents[-3]
        expeInfo['conditionID'] = dirPathComponents[-2]
        
    expeInfo['recordingID'] = dirPathComponents[-1]

    return expeInfo