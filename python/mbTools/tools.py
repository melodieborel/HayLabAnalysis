import os
import configparser
import pickle
import ast

import numpy as np

from ipyfilechooser import FileChooser
import ipywidgets as widgets
from IPython.display import display
from IPython import get_ipython

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