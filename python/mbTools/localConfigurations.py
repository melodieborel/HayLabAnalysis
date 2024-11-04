import os
import configparser

import ipywidgets as widgets

import pprint

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