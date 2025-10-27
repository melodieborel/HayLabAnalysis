# # Associate Ca2+ signal with sleep stages for each session & subsessions using crossregistration

#######################################################################################
                            # Define Experiment type #
#######################################################################################

DrugExperiment = 0 # =1 if CGP Experiment // DrugExperiment=0 if Baseline Experiment

AnalysisID = '_spikes_200ms' 

saveexcel = 0

local = True
if local:
    dir = "//10.69.168.1/crnldata/forgetting/Aurelie/MiniscopeOE_data/L2_3_mice/"
else: 
    dir = "/crnldata/forgetting/Aurelie/MiniscopeOE_data/L2_3_mice/"

drugs = ['baseline']

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

warnings.filterwarnings("ignore")
import sys
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

def resample_matrix(data, orig_rate, target_rate, axis=0):
    if orig_rate == target_rate:
        return data.copy()
    # Compute integer up/down factors using GCD
    up = int(target_rate)
    down = int(orig_rate)
    factor = gcd(up, down)
    up //= factor
    down //= factor
    return resample_poly(data, up=up, down=down, axis=axis)


def bin_sum_fractional(x, Fs1, Fs2):
    N = len(x)
    bin_width = Fs1 / Fs2
    y = []
    start = 0
    while start < N:
        end = start + bin_width
        i_start = int(np.floor(start))
        i_end = int(np.floor(end))
        if i_end >= N:
            i_end = N - 1
        total = 0.0
        total += x[i_start] * (1 - (start - i_start))
        for i in range(i_start + 1, i_end):
            total += x[i]
        if i_end > i_start:
            total += x[i_end] * (end - i_end)
        y.append(total)
        start = end
    return np.array(y)

def marcenkopastur(significance):
    nbins = significance.nbins
    nneurons = significance.nneurons
    tracywidom = significance.tracywidom
    q = float(nbins)/float(nneurons)
    lambdaMax = pow((1+np.sqrt(1/q)),2)
    lambdaMax += tracywidom*pow(nneurons,-2./3)
    return lambdaMax

def getlambdacontrol(zactmat_):
    significance_ = PCA()
    significance_.fit(zactmat_.T)
    lambdamax_ = np.max(significance_.explained_variance_)
    return lambdamax_

def binshuffling(zactmat,significance):
    np.random.seed()
    lambdamax_ = np.zeros(significance.nshu)
    for shui in range(significance.nshu):
        zactmat_ = np.copy(zactmat)
        for (neuroni,activity) in enumerate(zactmat_):
            randomorder = np.argsort(np.random.rand(significance.nbins))
            zactmat_[neuroni,:] = activity[randomorder]
        lambdamax_[shui] = getlambdacontrol(zactmat_)
    lambdaMax = np.percentile(lambdamax_,significance.percentile)
    return lambdaMax

def circshuffling(zactmat,significance):
    np.random.seed()
    lambdamax_ = np.zeros(significance.nshu)
    for shui in range(significance.nshu):
        zactmat_ = np.copy(zactmat)
        for (neuroni,activity) in enumerate(zactmat_):
            cut = int(np.random.randint(significance.nbins*2))
            zactmat_[neuroni,:] = np.roll(activity,cut)
        lambdamax_[shui] = getlambdacontrol(zactmat_)
    lambdaMax = np.percentile(lambdamax_,significance.percentile)
    return lambdaMax

def runSignificance(zactmat,significance):
    if significance.nullhyp == 'mp':
        lambdaMax = marcenkopastur(significance)
    elif significance.nullhyp == 'bin':
        lambdaMax = binshuffling(zactmat,significance)
    elif significance.nullhyp == 'circ':
        lambdaMax = circshuffling(zactmat,significance)
    else:
        print('ERROR !')
        print('    nyll hypothesis method '+str(nullhyp)+' not understood')
        significance.nassemblies = np.nan
    nassemblies = np.sum(significance.explained_variance_>lambdaMax)
    significance.nassemblies = nassemblies
    return significance

def extractPatterns(actmat,significance,method):
    nassemblies = significance.nassemblies
    if method == 'pca':
        idxs = np.argsort(-significance.explained_variance_)[0:nassemblies]
        patterns = significance.components_[idxs,:]
    elif method == 'ica':
        from sklearn.decomposition import FastICA
        ica = FastICA(n_components=nassemblies, max_iter=1000)
        ica.fit(actmat.T)
        patterns = ica.components_
    else:
        print('ERROR !')
        print('    assembly extraction method '+str(method)+' not understood')
        patterns = np.nan
    if patterns is not np.nan:
        patterns = patterns.reshape(nassemblies,-1)
        norms = np.linalg.norm(patterns,axis=1)
        patterns /= np.matlib.repmat(norms,np.size(patterns,1),1).T
    return patterns

def runPatterns(actmat, method='ica', nullhyp = 'mp', nshu = 1000, percentile = 99, tracywidom = False):
    nneurons = np.size(actmat,0)
    nbins = np.size(actmat,1)
    silentneurons = np.var(actmat,axis=1)==0
    actmat_ = actmat[~silentneurons,:]
    zactmat_ = stats.zscore(actmat_,axis=1)
    significance = PCA()
    significance.fit(zactmat_.T)
    significance.nneurons = nneurons
    significance.nbins = nbins
    significance.nshu = nshu
    significance.percentile = percentile
    significance.tracywidom = tracywidom
    significance.nullhyp = nullhyp
    significance = runSignificance(zactmat_,significance)
    if np.isnan(significance.nassemblies):
        patterns = []
        zactmat = []
        significance = []
        #return
    if significance.nassemblies<1:
        print('WARNING 1!')
        print('    no assembly detected!')
        patterns = []
        zactmat = []
        significance = []
    else:
        patterns_ = extractPatterns(zactmat_,significance,method)
        if patterns_ is np.nan:
            patterns = []
            zactmat = []
            significance = []
            #return
        patterns = np.zeros((np.size(patterns_,0),nneurons))
        patterns[:,~silentneurons] = patterns_
        zactmat = np.copy(actmat)
        zactmat[~silentneurons,:] = zactmat_
    return patterns,significance,zactmat

def computeAssemblyActivity(patterns,zactmat,zerodiag = True):
    if len(patterns) == 0:
        print('WARNING 2!')
        print('    no assembly detecded!')
        assemblyAct = []
    else:
        nassemblies = len(patterns)
        nbins = np.size(zactmat,1)
        assemblyAct = np.zeros((nassemblies,nbins))
        for (assemblyi,pattern) in enumerate(patterns):
            projMat = np.outer(pattern,pattern)
            projMat -= zerodiag*np.diag(np.diag(projMat))
            for bini in range(nbins):
                assemblyAct[assemblyi,bini] = np.dot(np.dot(zactmat[:,bini],projMat),zactmat[:,bini])
    return assemblyAct

#######################################################################################
                # Load sleep score and Ca2+ time series numpy arrays #
#######################################################################################

all_expe_types=['baseline','preCGP', 'postCGP'] if DrugExperiment else ['baseline', 'preCGP']

# Get the current date and time
FolderNameSave=str(datetime.now())[:19]
FolderNameSave = FolderNameSave.replace(" ", "_").replace(".", "_").replace(":", "_")

destination_folder= f"//10.69.168.1/crnldata/forgetting/Aurelie/MiniscopeOE_analysis/PlaceCells_experiment/2_CellAssemblies_{FolderNameSave}{AnalysisID}" if local else f"/crnldata/forgetting/Aurelie/MiniscopeOE_analysis/PlaceCells_experiment/2_CellAssemblies_{FolderNameSave}{AnalysisID}"
os.makedirs(destination_folder)
folder_to_save=Path(destination_folder)

logfile = open(f"{destination_folder}/output_log.txt", 'w')
sys.stdout = Tee(sys.stdout, logfile)  # print goes to both


# Copy the script file to the destination folder
source_script = "C:/Users/Manip2/SCRIPTS/HayLabAnalysis/python/_MINI&OE_2_detect_cell_assembly.py" if local else "/home/aurelie.brecier/HayLabAnalysis/python/_MINI&OE_2_detect_cell_assembly.py"
destination_file_path = f"{destination_folder}/_MINI&OE_2_detect_cell_assembly.txt"
shutil.copy(source_script, destination_file_path)

CellAssembly_dict={}
CellAssembly_members={}

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

    nb_minian_total=0
    dict_Path={}
    dict_Calcium = {}
    dict_Deconv = {}
    dict_Scoring = {}
    dict_Stamps = {}
    dict_TodropFile = {}
    dict_StampsMiniscope = {}

    CellAssembly_patterns=pd.DataFrame()  
    CellAssembly_members[mice]={}

    previousEndTime=0
    InitialStartTime=0

    assembly_nb=0

    minian_folders = [f for f in dpath.parents[0].rglob('minian') if f.is_dir()]

    for minianpath in minian_folders: # for each minian folders found in this mouse 

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

            if session_type == 'Cheeseboard':

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

                print(session, ': starts at', round(StartTime,1), 's & ends at', round(EndTime,1), 's (', round(rec_dur_sec,1), 's duration, ', numbdropfr, 'dropped frames, minian frequency =', minian_freq, 'Hz, experiment type = ', session_type, ')...') 
                sentence1= f"... kept values = {kept_uniq_unit_List}"
                print(sentence1)

                # Define cell assemblies
            
                target_rate = 5 #Hz == 50ms bins
                #Array_bin = resample_matrix(Carray, orig_rate=minian_freq, target_rate=target_rate)
                #Array_bin = resample_matrix(Darray, orig_rate=minian_freq, target_rate=target_rate)
                Array_bin = bin_sum_fractional(Sarray, minian_freq, target_rate)            
                patterns,significance,zactmat= runPatterns(Array_bin.T, method='ica', nullhyp = 'mp', nshu = 1000, percentile = 99, tracywidom = False)       
                all_patterns = pd.DataFrame({ass: patterns[ass].tolist() for ass in np.arange(np.shape(patterns)[0])}, index=kept_uniq_unit_List).add_prefix(f"{session_time}_CellAss")
                CellAssembly_patterns = CellAssembly_patterns.join(all_patterns, how="outer") if not CellAssembly_patterns.empty else all_patterns
                
                if len(patterns)>0:
                    patterns_th = patterns.copy()
                    for ass in np.arange(np.shape(patterns)[0]):
                        thresh = np.mean(patterns_th[ass]) + 3 * np.std(patterns_th[ass])
                        patterns_th[ass][abs(patterns_th[ass])<thresh]=np.nan
                        non_nan_indices = np.where(~np.isnan(patterns_th[ass]))[0] 
                        assembly_ID = all_patterns.columns.tolist()[ass]
                        cells_in_assembly=np.array(kept_uniq_unit_List)[non_nan_indices].tolist()  
                        CellAssembly_members[mice][assembly_ID]=cells_in_assembly
                
    CellAssembly_dict[mice]= CellAssembly_patterns

filenameOut = folder_to_save / f'CellAssembly_dict.pkl'
with open(filenameOut, 'wb') as pickle_file:
    pickle.dump(CellAssembly_dict, pickle_file)

filenameOut = folder_to_save / f'CellAssembly_members.pkl'
with open(filenameOut, 'wb') as pickle_file:
    pickle.dump(CellAssembly_members, pickle_file)

sys.stdout = sys.__stdout__
logfile.close()