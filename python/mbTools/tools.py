import os
import configparser
import pickle
import ast

import numpy as np

from ipyfilechooser import FileChooser
import ipywidgets as widgets
from IPython.display import display
from IPython import get_ipython

from .localConfigurations import localConf

from pathlib import Path

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
   
   offsetsLFP = np.arange(y.shape[1])
   y = np.add(y,offsetsLFP*offset)
   
   offsetsNPX = np.arange(y2.shape[1])
   print(y2.shape)
   y2 = np.add(y2,offsetsNPX*offset)
   
   plt.plot(x, y*scaleLFP,'-')
   plt.plot(x2, y2*scaleNPX+offset,'-')
   plt.show()

def convertTheoricIndex2realTime(thIdx,realFreq=1, offset=0):
   realTime=thIdx/realFreq + offset
   return realTime

def find_nearest(array, value):
   idx = (np.abs(array - value)).argmin()
   return idx

def getPathComponent(filename,project_type):
   if not filename.is_dir():
      filename = filename.parent
   dirPathComponents = filename.parts
   expeInfo = dict()
   remote_prefix = Path(filename.anchor)
   print(remote_prefix)
   expeInfo['analysis_path_root'] = remote_prefix.joinpath(*dirPathComponents[:-5]).relative_to(remote_prefix)
   expeInfo['project_id'] = dirPathComponents[-5]
   expeInfo['sub_project_id'] = dirPathComponents[-4]

   projectConfig = remote_prefix.joinpath(*dirPathComponents[0:-3],'projectConfig.ini')
   projParser = configparser.ConfigParser()
   try:
      if projectConfig.is_file():
         projParser.read(projectConfig)
         try:
            expeInfo['project_type'] = projParser.get('ALL','project_type')
         except:
            #TODO: soon remove, it was only for copatibility
            expeInfo['project_type'] = projParser.get('ALL','projecttype')
            projParser.set('ALL','project_type',expeInfo['project_type'])
            projParser.remove_option('ALL','projecttype')
            with open(projectConfig, 'w') as configfile:
               projParser.write(configfile)
      else:
         projParser.add_section('ALL')
         projParser.set('ALL','project_type',str(project_type))
         with open(projectConfig, 'w') as configfile:
            projParser.write(configfile)
   except Exception as e:
      print(f"there was an error: {e}. Possibly, make sure this is not related to the configFile name that should be {projectConfig}. If another with a different name exists, try to duplicate it with the expected name.")

   if expeInfo['project_type'] == 0:
      expeInfo['condition_id'] = dirPathComponents[-3]
      expeInfo['animal_id'] = dirPathComponents[-2]
   else:
      expeInfo['animal_id'] = dirPathComponents[-3]
      expeInfo['condition_id'] = dirPathComponents[-2]
      
   expeInfo['recording_id'] = dirPathComponents[-1]

   return expeInfo


def replacePathComponent(original_path: Path,old_component,new_component):
   try:
      path_parts =  list(Path(original_path).parts)
      index = path_parts.index(old_component)
      new_component = Path(new_component)
      path_parts[index:index+len(new_component.parts)-1] = new_component.parts
      new_path = Path().joinpath(*tuple(path_parts))
      print(f"path {original_path} manipulated to {new_path}")
      return new_path
   except ValueError as e:
      print(f"component {old_component} not found in path {original_path}. Error: {e}")
      return original_path