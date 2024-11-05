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