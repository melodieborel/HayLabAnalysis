# # Test

#######################################################################################
                            # Define Experiment type #
#######################################################################################

DrugExperiment=0 # =1 if CGP Experiment // DrugExperiment=0 if Baseline Experiment

AnalysisID='' 

saveexcel=1

local = False
if local:
    dir= "//10.69.168.1/crnldata/forgetting/Aurelie/MiniscopeOE_data/L1NDNF_mice/"
else: 
    dir= "/crnldata/forgetting/Aurelie/MiniscopeOE_data/L1NDNF_mice/OW"

mapp = {1: 'AW',  2: 'QW', 3: 'NREM',  4: 'IS', 5: 'REM',  6: 'undefined'}
drugs=['baseline']

#######################################################################################
                                # Load packages #
#######################################################################################

import os
import numpy as np
from scipy import signal
import quantities as pq
import math 
import neo
import json
from pathlib import Path
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, Cursor
import pickle
import sys 
from datetime import datetime
import shutil
from ast import literal_eval
from scipy.signal import find_peaks
from scipy.stats import zscore
from scipy import stats
from itertools import groupby
from IPython.display import display
from scipy.interpolate import griddata
import numpy as np
from scipy.signal import butter, filtfilt, hilbert
import matplotlib.pyplot as plt
from scipy import interpolate
import numpy.matlib
from sklearn.decomposition import PCA
from scipy import stats
import numpy as np
import pickle
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import resample
from scipy.signal import resample_poly
from math import gcd
import warnings
import numpy as np
from scipy.stats import zscore
import sys

warnings.filterwarnings("ignore")


minian_path = os.path.join(os.path.abspath('..'),'minian')
print("The folder used for minian procedures is : {}".format(minian_path))
sys.path.append(minian_path)

from minian.utilities import (
    TaskAnnotation,
    get_optimal_chk,
    load_videos,
    open_minian,
    save_minian,
)

#######################################################################################
                                # Define functions #
#######################################################################################
class Tee:
    def __init__(self, *files):
        self.files = files
    def write(self, obj):
        for f in self.files:
            f.write(obj)
            f.flush()
    def flush(self):
        for f in self.files:
            f.flush()

#######################################################################################
                # Load sleep score and Ca2+ time series numpy arrays #
#######################################################################################

all_expe_types=['baseline','preCGP', 'postCGP'] if DrugExperiment else ['baseline', 'preCGP']

# Get the current date and time
FolderNameSave=str(datetime.now())[:19]
FolderNameSave = FolderNameSave.replace(" ", "_").replace(".", "_").replace(":", "_")

destination_folder= f"//10.69.168.1/crnldata/forgetting/Théa/MiniscopeOE_analysis/Exploration_task/1_VigSt_{FolderNameSave}{AnalysisID}" if local else f"/crnldata/forgetting/Théa/MiniscopeOE_analysis/Exploration_task/1_VigSt_{FolderNameSave}{AnalysisID}"
os.makedirs(destination_folder)
folder_to_save=Path(destination_folder)


logfile = open(f"{destination_folder}/output_log.txt", 'w')
sys.stdout = Tee(sys.stdout, logfile)  # print goes to both


# Copy the script file to the destination folder
source_script = "C:/Users/Acquisition/HayLabAnalysis/python/_MINIOE_1_compute_vigstates_activity.py" if local else "/home/thea.michel/HayLabAnalysis/python/_MINIOE_1_compute_vig_activity.py"
destination_file_path = f"{destination_folder}/_MINIOE_1_compute_vig_activity.txt"
shutil.copy(source_script, destination_file_path)


data = {}
counter=0
VigilanceState_GlobalResults= pd.DataFrame(data, columns=['Mice','NeuronType','Session', 'Session_Date', 'Session_Time', 'Unique_Unit','UnitNumber','UnitValue','Unit_ID', 'UnitLocation',
                                                        'ExpeType', 'Drug', 'Substate','Substate_ID','Session_ID','SubstateNumber','DurationSubstate', 'CalciumActivity', 
                                                        'Avg_CalciumActivity', 'AUC_calcium','Avg_AUC_calcium', 'NormalizedAUC_calcium', 'DeconvSpikeMeanActivity', 
                                                        'Avg_DeconvSpikeActivity', 'SpikeActivityHz', 'Avg_SpikeActivityHz', 'TotCaPopCoupling', 
                                                        'TotZ_CaPopCoupling', 'TotSpPopCoupling', 'TotZ_SpPopCoupling'])

print("Testing path:", dir)
print("Exists:", os.path.exists(dir))

files = list(Path(dir).glob('**/mappingsAB.pkl'))
print("Nb files found:", len(files))

for dpath in Path(dir).glob('**/mappingsAB.pkl'):

    mappfile = open(dpath.parents[0]/ f'mappingsAB.pkl', 'rb')
    mapping = pickle.load(mappfile)
    mapping_sess = mapping['session']   

    centroidfile = open(dpath.parents[0]/ f'centsAB.pkl', 'rb')
    centroids = pickle.load(centroidfile)
    centroids_sess = centroids['session']   

    mice = dpath.parents[1].parts[-1]
    NeuronType = dpath.parents[2].parts[-1]
    
    print(f"####################################################################################")
    print(f"################################### {mice} ####################################")
    print(f"####################################################################################")

    if saveexcel: 
        filenameOut = folder_to_save / f'VigSt_CaCorr_{mice}.xlsx'
        excel_writerCa = pd.ExcelWriter(filenameOut)        
        filenameOut = folder_to_save / f'VigSt_SpCorr_{mice}.xlsx'
        excel_writerSp = pd.ExcelWriter(filenameOut)
        filenameOut = folder_to_save / f'VigSt_RawCaTraces_{mice}.xlsx'
        excel_writerRawCa = pd.ExcelWriter(filenameOut)        
        filenameOut = folder_to_save / f'VigSt_RawSpTraces_{mice}.xlsx'
        excel_writerRawSp = pd.ExcelWriter(filenameOut)

    nb_minian_total=0
    dict_Path={}
    dict_Calcium = {}
    dict_Deconv = {}
    dict_Scoring = {}
    dict_Stamps = {}
    dict_TodropFile = {}
    dict_StampsMiniscope = {}

    TotCaCorr=[]
    TotSpCorr=[]

    for drug in drugs: 
        for m in mapp: 
            locals()[f'CaCorr{mapp[m]}Matrix{drug}']=[]
            locals()[f'SpCorr{mapp[m]}Matrix{drug}']=[]
            locals()[f'RawCaTraces{mapp[m]}_{drug}']=[]
            locals()[f'RawSpTraces{mapp[m]}_{drug}']=[]    
            locals()[f'StatesCaCorr{mapp[m]}Matrix{drug}']=pd.DataFrame(columns=[f'{mice}{str(i)}' for i in range(len(mapping_sess))], index=[i for i in range(len(mapping_sess))])  
            locals()[f'ITStatesCaCorr{mapp[m]}Matrix{drug}']=pd.DataFrame(columns=[f'{mice}{str(i)}' for i in range(len(mapping_sess))], index=[i for i in range(len(mapping_sess))])
            locals()[f'StatesSpCorr{mapp[m]}Matrix{drug}']=pd.DataFrame(columns=[f'{mice}{str(i)}' for i in range(len(mapping_sess))], index=[i for i in range(len(mapping_sess))])  
            locals()[f'ITStatesSpCorr{mapp[m]}Matrix{drug}']=pd.DataFrame(columns=[f'{mice}{str(i)}' for i in range(len(mapping_sess))], index=[i for i in range(len(mapping_sess))])
        locals()[f'TotCaCorr{drug}']=[]    
        locals()[f'TotSpCorr{drug}']=[]    

    previousEndTime=0
    InitialStartTime=0

    assembly_nb=0

    counter=0
    counterProbe=0
    day=1
    trial=1
    previousmice=0
    previous_session_date=0
    previous_session_time=0
    
    minian_folders = [f for f in dpath.parents[0].rglob('minian') if f.is_dir()]

    for minianpath in minian_folders: # for each minian folders found in this mouse
    
        # garder uniquement Training
        if "training" not in [p.lower() for p in minianpath.parts]:
            continue

        drug= 'baseline' 

        if 1==1: # any(p in all_expe_types for p in minianpath.parts): # have to be to the expe_types

            if any(minianpath.parents[1].glob('*V4_Miniscope/timeStamps.csv')): 
                session_path=minianpath.parents[1]        
                tsmini= pd.read_csv(list(session_path.glob('*V4_Miniscope/timeStamps.csv'))[0])['Time Stamp (ms)']
                V4subfolder=False               
                if any(session_path.parent.rglob('OpenEphys/EphyViewer_ManualScor5s_6Stages.csv')):
                    session_type=minianpath.parents[3].name
                    session_date=minianpath.parents[4].name
                    session_time=minianpath.parents[1].name
                else:
                    session_type=minianpath.parents[2].name
                    session_date=minianpath.parents[3].name
                    session_time=minianpath.parents[1].name
                session=session_time        
            elif any(minianpath.parents[2].glob('*V4_Miniscope/timeStamps.csv')):
                session_path=minianpath.parents[2]      
                tsmini= pd.read_csv(list(session_path.glob('*V4_Miniscope/timeStamps.csv'))[0])['Time Stamp (ms)']                      
                V4subfolder=True
                session_type=minianpath.parents[4].name
                session_date=minianpath.parents[5].name
                session_time=minianpath.parents[0].name
                session=session_time




            if mice == previousmice:   
                if session_date == previous_session_date : 
                    if session_time==previous_session_time:
                        trial+=1
                    else:
                        day+=1
                        trial=1
                else: 
                    day=1
                    trial=1
            else:
                day=1
                trial=1



            print(f"Processing {session_type} session: {session} on the {session_date}, subfolders = {V4subfolder}")