# # Associate Ca2+ signal with sleep stages for each session & subsessions using crossregistration

#######################################################################################
                            # Define Experiment type #
#######################################################################################

DrugExperiment=0 # =1 if CGP Experiment // DrugExperiment=0 if Baseline Experiment

AnalysisID='' 

saveexcel=0

local = True
if local:
    dir= "//10.69.168.1/crnldata/forgetting/Aurelie/MiniscopeOE_data/L2_3_mice/"
else: 
    dir= "/crnldata/forgetting/Aurelie/MiniscopeOE_data/L2_3_mice/"

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

destination_folder= f"//10.69.168.1/crnldata/forgetting/Aurelie/MiniscopeOE_analysis/PlaceCells_experiment/1_VigSt_{FolderNameSave}{AnalysisID}" if local else f"/crnldata/forgetting/Aurelie/MiniscopeOE_analysis/PlaceCells_experiment/1_VigSt_{FolderNameSave}{AnalysisID}"
os.makedirs(destination_folder)
folder_to_save=Path(destination_folder)

logfile = open(f"{destination_folder}/output_log.txt", 'w')
sys.stdout = Tee(sys.stdout, logfile)  # print goes to both


# Copy the script file to the destination folder
source_script = "C:/Users/Manip2/SCRIPTS/HayLabAnalysis/python/_MINI&OE_1_compute_vigstates_activity.py" if local else "/home/aurelie.brecier/HayLabAnalysis/python/_MINI&OE_1_compute_vigstates_activity.py"
destination_file_path = f"{destination_folder}/_MINI&OE_1_compute_vigstates_activity.txt"
shutil.copy(source_script, destination_file_path)


data = {}
counter=0
VigilanceState_GlobalResults= pd.DataFrame(data, columns=['Mice','NeuronType','Session', 'Session_Date', 'Session_Time', 'Unique_Unit','UnitNumber','UnitValue','Unit_ID', 'UnitLocation',
                                                        'ExpeType', 'Drug', 'Substate','Substate_ID','Session_ID','SubstateNumber','DurationSubstate', 'CalciumActivity', 
                                                        'Avg_CalciumActivity', 'AUC_calcium','Avg_AUC_calcium', 'NormalizedAUC_calcium', 'DeconvSpikeMeanActivity', 
                                                        'Avg_DeconvSpikeActivity', 'SpikeActivityHz', 'Avg_SpikeActivityHz', 'TotCaPopCoupling', 
                                                        'TotZ_CaPopCoupling', 'TotSpPopCoupling', 'TotZ_SpPopCoupling'])
for dpath in Path(dir).glob('**/PlaceCells_experiment/mappingsAB.pkl'):

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

    minian_folders = [f for f in dpath.parents[0].rglob('minian') if f.is_dir()]

    for minianpath in minian_folders: # for each minian folders found in this mouse
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
            print(f"Processing {session_type} session: {session} on the {session_date}, subfolders = {V4subfolder}")

            minian_ds = open_minian(minianpath)
            dict_Calcium[session] = minian_ds['C'] # calcium traces 
            dict_Deconv[session] = minian_ds['S'] # estimated spikes deconvolved activity
            
            try: 
                try: 
                    with open(minianpath.parent / f'TodropFileAB.json', 'r') as f:
                        unit_to_drop = json.load(f)
                        dict_TodropFile[session]  = unit_to_drop
                except: 
                    with open(minianpath/ f'TodropFileAB.json', 'r') as f:
                        unit_to_drop = json.load(f)
                        dict_TodropFile[session]  = unit_to_drop
                nb_minian_total+=1
            except:
                print(' !!!! Minian session not validated !!!!')
                continue
            print(session_path)

            
            # Start time & freq miniscope
            
            dict_StampsMiniscope[session]  = pd.read_csv(list(session_path.glob('*V4_Miniscope/timeStamps.csv'))[0])
            tsmini=dict_StampsMiniscope[session]['Time Stamp (ms)']
            minian_freq=round(1/np.mean(np.diff(np.array(tsmini)/1000)))  

            freqLFP=1000    
            numbdropfr= 0         
            
            nb_minian_total+=1

            #######################################################################################
            # Distribute Ca2+ intensity & spikes to vigilance states for each sessions #
            #######################################################################################
                
            # Adjust the StartTime if subsessions

            dict_Stamps[session]  = pd.read_excel(session_path/ f'SynchroFile.xlsx')     
            StartTime = list(dict_Stamps[session][0])[0] 
            if math.isnan(StartTime):   
                StartTime = 0 # No OpenEphys file found, so Miniscope start relative to OpenEphys start equal 0

            if InitialStartTime==0:
                InitialStartTime=StartTime    
                firstframe=0
                StartTimeMiniscope=0 # start time of miniscope rec of that subsesssions relative to the start of the mniscope recording
            else:
                if StartTime == InitialStartTime: # just a subsession
                    StartTime = previousEndTime + 1/minian_freq #  +1 frame in seconds
                    StartTimeMiniscope= StartTime-InitialStartTime
                else:  
                    InitialStartTime=StartTime # this is a new session
                    firstframe=0
                    StartTimeMiniscope=0 
                    

            # Remove bad units from recordings

            C=dict_Calcium[session]
            D=dict_Deconv[session] 
            unit_to_drop=dict_TodropFile[session]    

            C_sel=C.drop_sel(unit_id=unit_to_drop)
            D_sel=D.drop_sel(unit_id=unit_to_drop)

            CalciumSub = pd.DataFrame(C_sel, index=C_sel['unit_id'])
            DeconvSub = pd.DataFrame(D_sel, index=D_sel['unit_id'])

            indexMappList=mapping_sess[session]
            kept_uniq_unit_List=[]
            for unit in CalciumSub.index:
                indexMapp = np.where(indexMappList == unit)[0]
                kept_uniq_unit_List.append(f"{mice}{str(indexMapp).replace('[','').replace(']','')}")

            nb_unit=len(CalciumSub)
            if nb_unit==0:
                print(f'no cells kept in the session: {session}')
                continue  # next iteration
            

            # Realigned the traces to the recorded timestamps 

            timestamps =  np.array(tsmini[firstframe:firstframe+len(CalciumSub.T)])/freqLFP
            new_timestamps= np.arange(timestamps[0], timestamps[-1], 1/minian_freq)
            
            Calcium = pd.DataFrame(index=CalciumSub.index, columns=new_timestamps)
            for feature in CalciumSub.index:
                interpolator = interpolate.interp1d(timestamps, CalciumSub.loc[feature], kind='linear')
                Calcium.loc[feature] = interpolator(new_timestamps)
            Carray=Calcium.values.T.astype(float) # Calcium activity

            Deconv = pd.DataFrame(index=DeconvSub.index, columns=new_timestamps)
            for feature in DeconvSub.index:
                interpolator = interpolate.interp1d(timestamps, DeconvSub.loc[feature], kind='linear')
                Deconv.loc[feature] = interpolator(new_timestamps)
            Darray=Deconv.values.T.astype(float) # Deconvolved activity

            Sarray= np.zeros((np.shape(Darray)[0], np.shape(Darray)[1])) # Spike activity
            for i in np.arange(np.shape(Darray)[1]):
                Darray_unit =Darray[:,i]
                peaks, _ = find_peaks(Darray_unit)
                Sarray_unit=np.zeros(len(Darray_unit))
                Sarray_unit[peaks]=1   
                Sarray[:,i]=Sarray_unit
            Spike= pd.DataFrame(Sarray.T, index=CalciumSub.index, columns=new_timestamps)

            firstframe+=len(CalciumSub.T)
            rec_dur=len(CalciumSub.T)
            rec_dur_sec= timestamps[-1]- timestamps[0] 
            EndTime = StartTime + rec_dur_sec # (upd_rec_dur/minian_freq) # in seconds
            previousEndTime=EndTime 

            if any(session_path.parent.rglob('OpenEphys/EphyViewer_ManualScor5s_6Stages.csv')):
                dict_Scoring[session]  = pd.read_csv(session_path.parent/ f'OpenEphys/EphyViewer_ManualScor5s_6Stages.csv')
                               
                # Upscale scoring to miniscope frequency
                SleepScored=dict_Scoring[session]
                Bool= np.sum(SleepScored['duration'])/len(SleepScored) == 5

                if not Bool:
                    print('/!\ Scoring freq=', np.sum(SleepScored['duration'])/len(SleepScored), 'Hz !!!!')
                    continue 

                SleepScored['label']= SleepScored['label'].str.extract(r'(\d+)', expand=False)
                SleepScoredTS=np.array(SleepScored['label'])

                scale_factor=minian_freq/0.2  #cause scoring was done in 5 seconds bin, ie 0.2 Hz              
                SleepScoredTS_upscaled = np.repeat(SleepScoredTS, scale_factor, axis=0)
                StartTime_frame=round(StartTime*minian_freq)
                SleepScoredTS_upscaled_ministart=SleepScoredTS_upscaled[StartTime_frame:StartTime_frame+rec_dur]
                #SleepScoredTS_upscaled_ministart[SleepScoredTS_upscaled_ministart == '2'] = '1'

                # Determine each substate identity and duration
                SleepScoredTS_upscaled_ministart=SleepScoredTS_upscaled_ministart.astype(int)
                substates_duration = [len(list(group)) for key, group in groupby(SleepScoredTS_upscaled_ministart)]
                substates_identity = [key for key, _ in groupby(SleepScoredTS_upscaled_ministart)]
                substates_end = np.array(substates_duration).cumsum()        
                substates_start =np.append([0],substates_end[:-1]) #substates_start =np.append([1],substates_end+1) create 1 second gap
                substates_identity = [mapp[num] for num in substates_identity]
                substates = pd.DataFrame(list(zip(substates_identity, substates_duration, substates_start, substates_end)), columns=['Identity', 'Duration', 'Start','End'])
    
            else:
                if session_type == 'Cheeseboard':
                    print(f"No scoring file found, assuming it's all Wake cause it's a Cheeseboard session")
                    substates = pd.DataFrame(data={'Identity': ['AW'], 'Duration': [Carray.shape[0]], 'Start': [StartTime], 'End': [Carray.shape[0]]})
                    SleepScoredTS_upscaled_ministart = np.ones(Carray.shape[0]).astype(int) # all awake
                    SleepScoredTS = np.ones(int(Carray.shape[0]*1.5)).astype(int) # all awake
                else: 
                    print(f'!!!! No scoring file found in a {session_type} session !!!!')
                    continue

            print(session, ': starts at', round(StartTime,1), 's & ends at', round(EndTime,1), 's (', round(rec_dur_sec,1), 's duration, ', numbdropfr, 'dropped frames, minian frequency =', minian_freq, 'Hz, experiment type = ', session_type, ')...') 
            sentence1= f"... kept values = {kept_uniq_unit_List}"
            print(sentence1)

            for m in mapp:

                # Correlation between each neurons according to vigilance states 
                Bool = (SleepScoredTS_upscaled_ministart == m)
                Carray_VigSpe = Carray.copy()
                Carray_VigSpe = Carray_VigSpe[0:np.shape(SleepScoredTS_upscaled_ministart)[0],:] # if Calcium imaging longer than LFP rec
                Carray_VigSpe = Carray_VigSpe[Bool, :]
                
                CaCorrVigStateMatrix = locals()[f'CaCorr{mapp[m]}Matrix{drug}']
                CaCorrMatrix=[]
                CaCorrMatrix = pd.DataFrame(columns=[f'{mice}{str(i)}' for i in range(len(mapping_sess))], index=[i for i in range(len(mapping_sess))])
                
                Sarray_VigSpe = Sarray.copy()
                Sarray_VigSpe = Sarray_VigSpe[0:np.shape(SleepScoredTS_upscaled_ministart)[0],:] # if Calcium imaging longer than LFP rec
                Sarray_VigSpe = Sarray_VigSpe[Bool, :]    

                SpCorrVigStateMatrix = locals()[f'SpCorr{mapp[m]}Matrix{drug}']
                SpCorrMatrix=[]
                SpCorrMatrix = pd.DataFrame(columns=[f'{mice}{str(i)}' for i in range(len(mapping_sess))], index=[i for i in range(len(mapping_sess))])

                if saveexcel:
                    RawCaTracesVigStateMatrix = locals()[f'RawCaTraces{mapp[m]}_{drug}']
                    RawCaTraces=[]
                    RawCaTraces=pd.DataFrame(Carray_VigSpe, columns=kept_uniq_unit_List)
                    unique_columns = RawCaTraces.columns[RawCaTraces.columns.to_series().duplicated()] # remove units that has an empty unique index '[]'
                    RawCaTraces = RawCaTraces.drop(columns=unique_columns)
                    RawCaTracesVigStateMatrix.append(RawCaTraces)                

                    RawSpTracesVigStateMatrix = locals()[f'RawSpTraces{mapp[m]}_{drug}']
                    RawSpTraces=[]
                    RawSpTraces=pd.DataFrame(Sarray_VigSpe, columns= kept_uniq_unit_List)
                    unique_columns = RawSpTraces.columns[RawSpTraces.columns.to_series().duplicated()] # remove units that has an empty unique index '[]'
                    RawSpTraces = RawSpTraces.drop(columns=unique_columns)
                    RawSpTracesVigStateMatrix.append(RawSpTraces) 
                
                for unit in range(nb_unit): 

                    Carray_unit =Carray_VigSpe[:,unit]
                    Sarray_unit =Sarray_VigSpe[:,unit] 
                    otherunit_range = [x for x in range(nb_unit) if x != unit]

                    for unit2 in range(nb_unit):

                        Carray_unit2 =Carray_VigSpe[:,unit2]
                        Sarray_unit2 =Sarray_VigSpe[:,unit2]         

                        indexMapp = str(np.where(mapping_sess[session] == Calcium.index[unit])[0]).replace('[','').replace(']','')
                        indexMapp2 = np.where(mapping_sess[session] == Calcium.index[unit2])[0]

                        if any(indexMapp) and len(indexMapp2)>0:      
                                            
                            corr_matrix = np.corrcoef(Carray_unit, Carray_unit2)
                            CaCorrMatrix[f'{mice}{indexMapp}'][indexMapp2]={1: 0.99999, -1: -0.99999}.get(np.round(corr_matrix[0, 1],5), np.round(corr_matrix[0, 1],5))
                            
                            corr_matrix = np.corrcoef(Sarray_unit, Sarray_unit2)
                            SpCorrMatrix[f'{mice}{indexMapp}'][indexMapp2]={1: 0.99999, -1: -0.99999}.get(np.round(corr_matrix[0, 1],5), np.round(corr_matrix[0, 1],5))
                            
                CaCorrVigStateMatrix.append(CaCorrMatrix)
                SpCorrVigStateMatrix.append(SpCorrMatrix) 
            
            TotCaCorr = locals()[f'TotCaCorr{drug}']
            TotSpCorr = locals()[f'TotSpCorr{drug}']
            
            TotCaCorrMatrix=[]
            TotCaCorrMatrix = pd.DataFrame(columns=[f'{mice}{str(i)}' for i in range(len(mapping_sess))], index=[i for i in range(len(mapping_sess))])
            TotSpCorrMatrix=[]
            TotSpCorrMatrix = pd.DataFrame(columns=[f'{mice}{str(i)}' for i in range(len(mapping_sess))], index=[i for i in range(len(mapping_sess))])

            for unit in range(nb_unit): 
                indexMapp = np.where(mapping_sess[session] == Calcium.index[unit])[0]
                
                if len(indexMapp)>0 : # The neuron needs to be in the cross-registration
                    
                    
                    Carray_unit =Carray[:,unit]
                    Darray_unit =Darray[:,unit]
                    Sarray_unit =Sarray[:,unit]
                    
                    # Population Coupling independent of vigilance states 

                    Carray_Population =np.mean(Carray[:,otherunit_range], axis=1)
                    corr_matrix = np.corrcoef(Carray_unit, Carray_Population)
                    TotCaPopCoupling= np.round(corr_matrix[0, 1],5)
                    corCorrected = {1: 0.99999, -1: -0.99999}.get(np.round(corr_matrix[0, 1],5), np.round(corr_matrix[0, 1],5))
                    TotZ_CaPopCoupling = np.arctanh(corCorrected)

                    Sarray_Population =np.mean(Sarray[:,otherunit_range], axis=1)
                    corr_matrix = np.corrcoef(Sarray_unit, Sarray_Population)
                    TotSpPopCoupling = np.round(corr_matrix[0, 1],5)
                    corCorrected = {1: 0.99999, -1: -0.99999}.get(np.round(corr_matrix[0, 1],5), np.round(corr_matrix[0, 1],5))
                    TotZ_SpPopCoupling= np.arctanh(corCorrected)

                    # Correlation between each neurons independent of vigilance states 

                    otherunit_range = [x for x in range(nb_unit) if x != unit]
                    for unit2 in otherunit_range:

                        Carray_unit2 =Carray[:,unit2]
                        Sarray_unit2 =Sarray[:,unit2]  
                    
                        indexMapp = str(np.where(mapping_sess[session] == Calcium.index[unit])[0]).replace('[','').replace(']','')
                        indexMapp2 = np.where(mapping_sess[session] == Calcium.index[unit2])[0]

                        if any(indexMapp) and len(indexMapp2)>0:    
                            
                            corr_matrix = np.corrcoef(Carray_unit, Carray_unit2)
                            TotCaCorrMatrix[f'{mice}{indexMapp}'][indexMapp2]={1: 0.99999, -1: -0.99999}.get(np.round(corr_matrix[0, 1],5), np.round(corr_matrix[0, 1],5))
                            
                            corr_matrix = np.corrcoef(Sarray_unit, Sarray_unit2)
                            TotSpCorrMatrix[f'{mice}{indexMapp}'][indexMapp2]={1: 0.99999, -1: -0.99999}.get(np.round(corr_matrix[0, 1],5), np.round(corr_matrix[0, 1],5))
                    
                    # For each substates
                    for index in range(len(substates)):
                        
                        ca_input_sub=Carray_unit[substates.Start[index]:substates.End[index]]
                        ds_input_sub=Darray_unit[substates.Start[index]:substates.End[index]]
                        sp_input_sub=Sarray_unit[substates.Start[index]:substates.End[index]]

                        VigilanceState_GlobalResults.loc[counter, 'Mice'] = mice
                        VigilanceState_GlobalResults.loc[counter, 'NeuronType'] = NeuronType
                        
                        VigilanceState_GlobalResults.loc[counter, 'Session'] = session
                        VigilanceState_GlobalResults.loc[counter, 'Session_Date'] = session_date 
                        VigilanceState_GlobalResults.loc[counter, 'Session_Time'] = session_time                    

                        VigilanceState_GlobalResults.loc[counter, 'Unique_Unit'] = indexMapp 
                        VigilanceState_GlobalResults.loc[counter, 'UnitNumber'] = unit
                        VigilanceState_GlobalResults.loc[counter, 'UnitValue'] = Calcium.index[unit]                      
                        VigilanceState_GlobalResults.loc[counter, 'Unit_ID'] = mice  + indexMapp
                        
                        centroids_sess=centroids[centroids['session']==session]
                        VigilanceState_GlobalResults.loc[counter, 'UnitLocation'] = [centroids_sess[centroids_sess['unit_id']==Calcium.index[unit]]['height'].tolist(), centroids_sess[centroids_sess['unit_id']==Calcium.index[unit]]['width'].tolist()] 

                        VigilanceState_GlobalResults.loc[counter, 'ExpeType'] =  session_type
                        VigilanceState_GlobalResults.loc[counter, 'Drug'] =  drug

                        VigilanceState_GlobalResults.loc[counter, 'Substate'] = substates.Identity[index]
                        VigilanceState_GlobalResults.loc[counter, 'Substate_ID'] = mice + session + str(substates.Identity[index]) + str(substates.index[index])
                        VigilanceState_GlobalResults.loc[counter, 'Session_ID'] = session_date + '_' + session_time
                        VigilanceState_GlobalResults.loc[counter, 'SubstateNumber'] = substates.index[index]
                        VigilanceState_GlobalResults.loc[counter, 'DurationSubstate'] = (substates.Duration[index]/minian_freq)

                        VigilanceState_GlobalResults.loc[counter, 'CalciumActivity'] = ca_input_sub.mean()
                        VigilanceState_GlobalResults.loc[counter, 'Avg_CalciumActivity'] = Carray_unit.mean()

                        VigilanceState_GlobalResults.loc[counter, 'AUC_calcium'] = np.trapz(ca_input_sub,np.arange(0,len(ca_input_sub),1))
                        VigilanceState_GlobalResults.loc[counter, 'Avg_AUC_calcium'] = np.trapz(Carray_unit,np.arange(0,len(Carray_unit),1))
                        VigilanceState_GlobalResults.loc[counter, 'NormalizedAUC_calcium'] = np.trapz(ca_input_sub,np.arange(0,len(ca_input_sub),1))/ (substates.Duration[index]/minian_freq)

                        VigilanceState_GlobalResults.loc[counter, 'DeconvSpikeMeanActivity'] = ds_input_sub.mean()
                        VigilanceState_GlobalResults.loc[counter, 'Avg_DeconvSpikeActivity'] = Darray_unit.mean()

                        VigilanceState_GlobalResults.loc[counter, 'SpikeActivityHz'] = sp_input_sub.sum()/(len(sp_input_sub)/minian_freq)
                        VigilanceState_GlobalResults.loc[counter, 'Avg_SpikeActivityHz'] = Sarray_unit.sum()/(len(Sarray_unit)/minian_freq)

                        otherunit_range = [x for x in range(nb_unit) if x != unit]
                        Carray_Population =np.mean(Carray[:,otherunit_range], axis=1)
                        ca_input_sub2=Carray_Population[substates.Start[index]:substates.End[index]]
                        corr_matrix = np.corrcoef(ca_input_sub, ca_input_sub2)

                        VigilanceState_GlobalResults.loc[counter, 'CaPopCoupling'] = np.round(corr_matrix[0, 1],5)
                        corCorrected = {1: 0.99999, -1: -0.99999}.get(np.round(corr_matrix[0, 1],5), np.round(corr_matrix[0, 1],5))
                        VigilanceState_GlobalResults.loc[counter, 'Z_CaPopCoupling'] = np.arctanh(corCorrected)

                        Sarray_Population =np.mean(Sarray[:,otherunit_range], axis=1)
                        ds_input_sub2=Sarray_Population[substates.Start[index]:substates.End[index]]
                        corr_matrix = np.corrcoef(ds_input_sub, ds_input_sub2)

                        VigilanceState_GlobalResults.loc[counter, 'SpPopCoupling'] = np.round(corr_matrix[0, 1],5)
                        corCorrected = {1: 0.99999, -1: -0.99999}.get(np.round(corr_matrix[0, 1],5), np.round(corr_matrix[0, 1],5))
                        VigilanceState_GlobalResults.loc[counter, 'Z_SpPopCoupling'] = np.arctanh(corCorrected)
                        
                        
                        VigilanceState_GlobalResults.loc[counter, 'TotCaPopCoupling'] = TotCaPopCoupling
                        VigilanceState_GlobalResults.loc[counter, 'TotZ_CaPopCoupling'] = TotZ_CaPopCoupling
                        VigilanceState_GlobalResults.loc[counter, 'TotSpPopCoupling'] = TotSpPopCoupling
                        VigilanceState_GlobalResults.loc[counter, 'TotZ_SpPopCoupling'] = TotZ_SpPopCoupling
                        

                        StatesCaCorrName=f'StatesCaCorr{substates.Identity[index]}Matrix{drug}'
                        StatesCaCorrMatrix = locals()[StatesCaCorrName]
                        IterationMatrixName=f'ITStatesCaCorr{substates.Identity[index]}Matrix{drug}'
                        IterationMatrix = locals()[IterationMatrixName]
                        
                        # Correlation between each neurons dependent of vigilance states 
                        if (substates.Duration[index]/minian_freq)>=20:
                            otherunit_range = [x for x in range(nb_unit) if x != unit]
                            for unit2 in otherunit_range:
                                Carray_unit2 =Carray[:,unit2]
                                ca2_input_sub=Carray_unit2[substates.Start[index]:substates.End[index]]
                            
                                indexMapp = str(np.where(mapping_sess[session] == Calcium.index[unit])[0]).replace('[','').replace(']','')
                                indexMapp2 = np.where(mapping_sess[session] == Calcium.index[unit2])[0]

                                if any(indexMapp) and len(indexMapp2)>0:    

                                    corr_matrix = np.corrcoef(ca_input_sub, ca2_input_sub)
                                    StatesCaCorrMatrix[f'{mice}{indexMapp}'][indexMapp2]=StatesCaCorrMatrix[f'{mice}{indexMapp}'][indexMapp2]+  np.arctanh({1: 0.99999, -1: -0.99999}.get(np.round(corr_matrix[0, 1],5), np.round(corr_matrix[0, 1],5))) if not math.isnan(StatesCaCorrMatrix[f'{mice}{indexMapp}'][indexMapp2]) else  np.arctanh({1: 0.99999, -1: -0.99999}.get(np.round(corr_matrix[0, 1],5), np.round(corr_matrix[0, 1],5)))
                                    if not math.isnan(StatesCaCorrMatrix[f'{mice}{indexMapp}'][indexMapp2]):
                                        if not math.isnan(IterationMatrix[f'{mice}{indexMapp}'][indexMapp2]):
                                            IterationMatrix[f'{mice}{indexMapp}'][indexMapp2]= IterationMatrix[f'{mice}{indexMapp}'][indexMapp2]+1 
                                        else:
                                            IterationMatrix[f'{mice}{indexMapp}'][indexMapp2]= 1 
                        counter+=1

            TotCaCorr.append(TotCaCorrMatrix)
            TotSpCorr.append(TotSpCorrMatrix)

    dataCaCorr={}
    dataSpCorr={}
    dataCaCorr2={}
    dataSpCorr2={}
    dataStatesCaCorr={}

    for d in drugs:
        TotCaCorrName=f'TotCaCorr{d}'
        TotCaCorr = locals()[TotCaCorrName]
        TotSpCorrName=f'TotSpCorr{d}'
        TotSpCorr = locals()[TotSpCorrName]
        if len(TotCaCorr)>0: 
            combined_df = pd.concat(TotCaCorr, ignore_index=False)
            IterationNb=combined_df.groupby(combined_df.index).count()
            combined_df = combined_df.groupby(combined_df.index).sum() #mean
            combined_df.index= [mice + str(idx) for idx in combined_df.index]
            combined_df = combined_df[~(combined_df.fillna(0) == 0).all(axis=1)]
            combined_df = combined_df.loc[:, ~(combined_df.fillna(0) == 0).all(axis=0)]
            IterationNb = IterationNb[~(IterationNb.fillna(0) == 0).all(axis=1)]
            IterationNb = IterationNb.loc[:, ~(IterationNb.fillna(0) == 0).all(axis=0)]
            IterationNb.index=combined_df.index
            dataCaCorr2[f'{d}']=combined_df

            combined_df = pd.concat(TotCaCorr, ignore_index=False)
            combined_df = combined_df.applymap(lambda x: np.arctanh(x) if not pd.isna(x) else np.nan)
            combined_df = combined_df.groupby(combined_df.index).sum() #mean
            combined_df.index = [mice + str(idx) for idx in combined_df.index]
            combined_df = combined_df[~(combined_df.fillna(0) == 0).all(axis=1)]
            combined_df = combined_df.loc[:, ~(combined_df.fillna(0) == 0).all(axis=0)]
            dataCaCorr2[f'Z_{d}']=combined_df
            dataCaCorr2[f'{d}_IterationNb']=IterationNb

            combined_df = pd.concat(TotSpCorr, ignore_index=False)
            IterationNb=combined_df.groupby(combined_df.index).count()
            combined_df = combined_df.groupby(combined_df.index).sum() #mean
            combined_df.index= [mice + str(idx) for idx in combined_df.index]
            combined_df = combined_df[~(combined_df.fillna(0) == 0).all(axis=1)]
            combined_df = combined_df.loc[:, ~(combined_df.fillna(0) == 0).all(axis=0)]
            IterationNb = IterationNb[~(IterationNb.fillna(0) == 0).all(axis=1)]
            IterationNb = IterationNb.loc[:, ~(IterationNb.fillna(0) == 0).all(axis=0)]
            IterationNb.index=combined_df.index
            dataSpCorr2[f'{d}']=combined_df

            combined_df = pd.concat(TotSpCorr, ignore_index=False)
            combined_df = combined_df.applymap(lambda x: np.arctanh(x) if not pd.isna(x) else np.nan)
            combined_df = combined_df.groupby(combined_df.index).sum() #mean
            combined_df.index = [mice + str(idx) for idx in combined_df.index]
            combined_df = combined_df[~(combined_df.fillna(0) == 0).all(axis=1)]
            combined_df = combined_df.loc[:, ~(combined_df.fillna(0) == 0).all(axis=0)]
            dataSpCorr2[f'Z_{d}']=combined_df
            dataSpCorr2[f'{d}_IterationNb']=IterationNb

        for m in mapp:
            IterationMatrixName=f'ITStatesCaCorr{mapp[m]}Matrix{d}'
            IterationMatrix = locals()[IterationMatrixName]
            StatesCaCorrMatrixName=f'StatesCaCorr{mapp[m]}Matrix{d}'
            StatesCaCorrMatrix = locals()[StatesCaCorrMatrixName]
            if len(StatesCaCorrMatrix)>0: # cause sometimes no Baseline conditions in CGP experiments
                StatesCaCorrMatrix.index = [mice + str(idx) for idx in StatesCaCorrMatrix.index]
                IterationMatrix.index = StatesCaCorrMatrix.index
                dataStatesCaCorr[f'{d}_{mapp[m]}']=StatesCaCorrMatrix
                dataStatesCaCorr[f'Iteration_{d}_{mapp[m]}']=IterationMatrix

            CaCorrVigStateMatrixName=f'CaCorr{mapp[m]}Matrix{d}'
            CaCorrVigStateMatrix = locals()[CaCorrVigStateMatrixName]
            if len(CaCorrVigStateMatrix)>0: # cause sometimes no Baseline conditions in CGP experiments
                
                combined_df = pd.concat(CaCorrVigStateMatrix, ignore_index=False)
                IterationNb=combined_df.groupby(combined_df.index).count()
                combined_df = combined_df.groupby(combined_df.index).sum() #mean
                combined_df.index= [mice + str(idx) for idx in combined_df.index]
                combined_df = combined_df[~(combined_df.fillna(0) == 0).all(axis=1)]
                combined_df = combined_df.loc[:, ~(combined_df.fillna(0) == 0).all(axis=0)]
                IterationNb = IterationNb[~(IterationNb.fillna(0) == 0).all(axis=1)]
                IterationNb = IterationNb.loc[:, ~(IterationNb.fillna(0) == 0).all(axis=0)]
                IterationNb.index=combined_df.index
                if saveexcel: combined_df.to_excel(excel_writerCa, sheet_name=f'{d}_{mapp[m]}', index=True, header=True) 
                dataCaCorr[f'{d}_{mapp[m]}']=combined_df
                
                combined_df = pd.concat(CaCorrVigStateMatrix, ignore_index=False)
                combined_df = combined_df.applymap(lambda x: np.arctanh(x) if not pd.isna(x) else np.nan)
                combined_df = combined_df.groupby(combined_df.index).sum() #mean
                combined_df.index = [mice + str(idx) for idx in combined_df.index]
                combined_df = combined_df[~(combined_df.fillna(0) == 0).all(axis=1)]
                combined_df = combined_df.loc[:, ~(combined_df.fillna(0) == 0).all(axis=0)]
                if saveexcel: combined_df.to_excel(excel_writerCa, sheet_name=f'Z_{d}_{mapp[m]}', index=True, header=True)   
                dataCaCorr[f'Z_{d}_{mapp[m]}']=combined_df
                if saveexcel: IterationNb.to_excel(excel_writerCa, sheet_name=f'{d}_{mapp[m]}_IterationNb', index=True, header=True) 
                dataCaCorr[f'{d}_{mapp[m]}_IterationNb']=IterationNb
                

                SpCorrVigStateMatrixName=f'SpCorr{mapp[m]}Matrix{d}'
                SpCorrVigStateMatrix = locals()[SpCorrVigStateMatrixName]
                combined_df = pd.concat(SpCorrVigStateMatrix, ignore_index=False)
                IterationNb=combined_df.groupby(combined_df.index).count()
                combined_df = combined_df.groupby(combined_df.index).sum()
                combined_df.index = [mice + str(idx) for idx in combined_df.index]
                combined_df = combined_df[~(combined_df.fillna(0) == 0).all(axis=1)]
                combined_df = combined_df.loc[:, ~(combined_df.fillna(0) == 0).all(axis=0)]
                IterationNb = IterationNb[~(IterationNb.fillna(0) == 0).all(axis=1)]
                IterationNb = IterationNb.loc[:, ~(IterationNb.fillna(0) == 0).all(axis=0)]
                IterationNb.index=combined_df.index
                if saveexcel: combined_df.to_excel(excel_writerSp, sheet_name=f'{d}_{mapp[m]}', index=True, header=True) 
                dataSpCorr[f'{d}_{mapp[m]}']=combined_df
                
                combined_df = pd.concat(SpCorrVigStateMatrix, ignore_index=False)
                combined_df = combined_df.applymap(lambda x: np.arctanh(x) if not pd.isna(x) else np.nan)
                combined_df = combined_df.groupby(combined_df.index).sum()
                combined_df.index = [mice + str(idx) for idx in combined_df.index]
                combined_df = combined_df[~(combined_df.fillna(0) == 0).all(axis=1)]
                combined_df = combined_df.loc[:, ~(combined_df.fillna(0) == 0).all(axis=0)]
                if saveexcel: combined_df.to_excel(excel_writerSp, sheet_name=f'Z_{d}_{mapp[m]}', index=True, header=True) 
                dataSpCorr[f'Z_{d}_{mapp[m]}']=combined_df                
                if saveexcel: IterationNb.to_excel(excel_writerSp, sheet_name=f'{d}_{mapp[m]}_IterationNb', index=True, header=True) 
                dataSpCorr[f'{d}_{mapp[m]}_IterationNb']=IterationNb
            
            if saveexcel:
                RawCaTracesVigStateMatrixName=f'RawCaTraces{mapp[m]}_{d}'
                RawCaTracesVigStateMatrix= locals()[RawCaTracesVigStateMatrixName]
                if len(RawCaTracesVigStateMatrix)>0: # cause sometimes no Baseline conditions in CGP experiments
                    combined_df = pd.concat(RawCaTracesVigStateMatrix, ignore_index=False)
                    combined_df.index = [mice + str(idx) for idx in combined_df.index]
                    combined_df = combined_df.dropna(axis=0, how='all')
                    combined_df = combined_df.dropna(axis=1, how='all')
                    combined_df.to_excel(excel_writerRawCa, sheet_name=f'{d}_{mapp[m]}', index=False, header=True)   

                    RawSpTracesVigStateMatrixName=f'RawSpTraces{mapp[m]}_{d}'
                    RawSpTracesVigStateMatrix= locals()[RawSpTracesVigStateMatrixName]
                    combined_df = pd.concat(RawSpTracesVigStateMatrix, ignore_index=False)
                    combined_df.index = [mice + str(idx) for idx in combined_df.index]
                    combined_df = combined_df.dropna(axis=0, how='all')
                    combined_df = combined_df.dropna(axis=1, how='all')
                    combined_df.to_excel(excel_writerRawSp, sheet_name=f'{d}_{mapp[m]}', index=False, header=True)   
        
    if saveexcel: 
        excel_writerCa.close() 
        excel_writerSp.close() 
        excel_writerRawCa.close()
        excel_writerRawSp.close()

    filenameOut = folder_to_save / f'VigSt_CaCorr_{mice}.pkl'
    with open(filenameOut, 'wb') as pickle_file:
        pickle.dump(dataCaCorr, pickle_file)
    filenameOut = folder_to_save / f'VigSt_SpCorr_{mice}.pkl'
    with open(filenameOut, 'wb') as pickle_file:
        pickle.dump(dataSpCorr, pickle_file)
        
    filenameOut = folder_to_save / f'TotCaCorr_{mice}.pkl'
    with open(filenameOut, 'wb') as pickle_file:
        pickle.dump(dataCaCorr2, pickle_file)
    filenameOut = folder_to_save / f'TotSpCorr_{mice}.pkl'
    with open(filenameOut, 'wb') as pickle_file:
        pickle.dump(dataSpCorr2, pickle_file)

    filenameOut = folder_to_save / f'StatesCaCorr_{mice}.pkl'
    with open(filenameOut, 'wb') as pickle_file:
        pickle.dump(dataStatesCaCorr, pickle_file)    
    
    filenameOut = folder_to_save / f'VigStates_Global_{mice}.pkl'
    with open(filenameOut, 'wb') as pickle_file:
        pickle.dump(VigilanceState_GlobalResults, pickle_file)

filenameOut = folder_to_save / f'VigStates_Global.pkl'
with open(filenameOut, 'wb') as pickle_file:
    pickle.dump(VigilanceState_GlobalResults, pickle_file)

if saveexcel: 
    filenameOut = folder_to_save / f'VigStates_Global.xlsx'
    writer = pd.ExcelWriter(filenameOut)
    VigilanceState_GlobalResults.to_excel(writer)
    writer.close()

sys.stdout = sys.__stdout__
logfile.close()