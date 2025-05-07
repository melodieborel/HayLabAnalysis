# # Associate Ca2+ signal with sleep stages for each session & subsessions using crossregistration

#######################################################################################
                            # Define Experiment type #
#######################################################################################

AnalysisID='' 

saveexcel=1

dir = "//10.69.168.1/crnldata/forgetting/Clementine/CheeseboardExperiment/L2_3/"

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
import re
import warnings
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


minian_path = os.path.join(os.path.abspath('.'),'minian')
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

def Convert(string):
            li = list(string.split(", "))
            li2 = len(li)
            return li2

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
        #print('WARNING 1!')
        #print('    no assembly detected!')
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
        #print('WARNING 2!')
        #print('    no assembly detecded!')
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

#all_expe_types=['Habituation','Training', 'Test'] 
all_expe_types=['Training'] 

# Get the current date and time
FolderNameSave=str(datetime.now())[:19]
FolderNameSave = FolderNameSave.replace(" ", "_").replace(".", "_").replace(":", "_")

destination_folder= f"//10.69.168.1/crnldata/forgetting/Clementine/Neuron_analysis/{FolderNameSave}{AnalysisID}"
os.makedirs(destination_folder)
folder_to_save=Path(destination_folder)

logfile = open(f"{destination_folder}/output_log.txt", 'w')
sys.stdout = Tee(sys.stdout, logfile)  # print goes to both

# Copy the script file to the destination folder
source_script = "C:/Users/Manip2/SCRIPTS/CodePythonAudrey/CodePythonAurelie/HayLabAnalysis/python/_DLC&MINI_1_average_activity.py"
destination_file_path = f"{destination_folder}/_DLC&MINI_1_average_activity.txt"
shutil.copy(source_script, destination_file_path)

data = {}
counter=0
NeuronActivity_GlobalResults= pd.DataFrame(data, columns=['Mice','NeuronType', 'SessionType', 'TrialType', 'Session_Date', 'Day#','Trial_Time','Trial#', 'Trial_Duration_in_sec',
                                                          'Unit#_in_crossreg','Unit#_in_session','UnitLocation',
                                                            'NormalizedAUC_calcium', 'DeconvSpikeMeanActivity','SpikeActivityHz'])
counter2=0

CellAssembly_GlobalResults= pd.DataFrame(data, columns=['Mice','NeuronType', 'SessionType', 'TrialType', 'Session_Date', 'Day#', 'Trial_Time','Trial#', 'Trial_Duration_in_sec',
                                                        'Assembly_ID', 'Assembly_size', 'Cells_in_Assembly',
                                                            'Avg_Activity', 'EventFreq', 'EventTime' ])

previousmice=0
previous_session_time=0

for dpath in Path(dir).glob('**/mappingsAB.pkl'):
    
    mappfile = open(dpath.parents[0]/ f'mappingsAB.pkl', 'rb')
    mapping = pickle.load(mappfile)
    mapping_sess = mapping['session']   
        
    centfile = open(dpath.parents[0]/ f'centsAB.pkl', 'rb')
    centroids = pickle.load(centfile) 

    mice = dpath.parents[0].parts[-1]
    NeuronType = dpath.parents[1].parts[-1]
    
    print(f"####################################################################################")
    print(f"################################### {mice} ####################################")
    print(f"####################################################################################")

    if saveexcel: 
        filenameOut = folder_to_save / f'RawCaTraces_{mice}.xlsx'
        excel_writerRawCa = pd.ExcelWriter(filenameOut)        
        filenameOut = folder_to_save / f'RawSpTraces_{mice}.xlsx'
        excel_writerRawSp = pd.ExcelWriter(filenameOut)

    nb_minian_total=0    
    day=0
    trial=0  
    prevsession_date=''
    prevtrial_type=''


    minian_folders = [f for f in dpath.parents[0].rglob('minian') if f.is_dir()]

    for minianpath in minian_folders: # for each minian folders found in this mouse 

        if any(p in all_expe_types for p in minianpath.parts): # have to be to the expe_types

            trial_time= minianpath.parents[1].name
            trial_type=minianpath.parents[2].name 
            session_date= minianpath.parents[3].name
            session_type=minianpath.parents[4].name
            print(minianpath)

            if prevsession_date == session_date:
                if prevtrial_type == trial_type:
                    trial+=1
                else: 
                    trial=1
            else: 
                day+=1
                trial=1
            
            prevsession_date = session_date
            prevtrial_type = trial_type


            assembly_nb=0

        
            minian_ds = open_minian(minianpath)
            C = minian_ds['C'] # calcium traces 
            D = minian_ds['S'] # estimated spikes deconvolved activity

            try: 
                with open(minianpath.parents[0] / f'TodropFileAB.json', 'r') as f:
                    unit_to_drop = json.load(f)
            except:
                unit_to_drop = []
                print(f"!!!!!! Need to define unit to drop. By default, all units are kept in this session !!!!!!!")

            nb_minian_total+=1

            #######################################################################################
            # Distribute Ca2+ intensity & spikes for each sessions #
            #######################################################################################

            minian_freq=30

            # Remove bad units from recordings
            C_sel=C.drop_sel(unit_id=unit_to_drop)
            D_sel=D.drop_sel(unit_id=unit_to_drop)

            Calcium = pd.DataFrame(C_sel, index=C_sel['unit_id'])
            Deconv = pd.DataFrame(D_sel, index=D_sel['unit_id'])

            indexMappList=mapping_sess[trial_time]
            kept_uniq_unit_List=[]
            for unit in Calcium.index:
                indexMapp = np.where(indexMappList == unit)[0]
                kept_uniq_unit_List.append(str(indexMapp))

            nb_unit=len(Calcium)
            if nb_unit==0:
                print(f'no cells kept in the session: {minian_path}')
                continue  # next iteration            
            
            Carray=Calcium.values.T.astype(float) # Calcium activity
            Darray=Deconv.values.T.astype(float) # Deconvolved activity
            Sarray= np.zeros((np.shape(Darray)[0], np.shape(Darray)[1])) # Spike activity
            for i in np.arange(np.shape(Darray)[1]):
                Darray_unit =Darray[:,i]
                peaks, _ = find_peaks(Darray_unit)
                Sarray_unit=np.zeros(len(Darray_unit))
                Sarray_unit[peaks]=1   
                Sarray[:,i]=Sarray_unit 

            sentence1= f"... number of units = {nb_unit}"
            print(sentence1)        

            # Define cell assemblies
            patterns,significance,zactmat= runPatterns(Carray.T, method='ica', nullhyp = 'mp', nshu = 1000, percentile = 99, tracywidom = False)
            if len(patterns)>0:
                thresh = np.mean(patterns)+2*np.std(patterns)
                patterns_th=patterns.copy()
                patterns_th[patterns_th<thresh]=np.nan
                for ass in np.arange(np.shape(patterns)[0]):
                    non_nan_indices = np.where(~np.isnan(patterns_th[ass]))[0] 
                    if len(non_nan_indices)>1:
                        cells_in_assembly=[]
                        for unit in non_nan_indices:
                            indexMapp = np.where(indexMappList == Calcium.index[unit])[0]
                            if len(indexMapp)>0:
                                cells_in_assembly.append(f"{mice}{str(indexMapp).replace('[','').replace(']','')}")
                        assembly_activity = zscore(Calcium.iloc[non_nan_indices].values.astype(float).T).mean(axis=1)
                        assembly_nb+=1
                        assembly_ID=f'Assembly_{mice}_{assembly_nb}'

                        mean = np.mean(assembly_activity)
                        std = np.std(assembly_activity)
                        prominence_threshold = mean + 2 * std

                        mean_act_ass = np.nanmean(assembly_activity)
                        
                        peaks, properties = find_peaks(assembly_activity, prominence=prominence_threshold)

                        CellAssembly_GlobalResults.loc[counter2, 'Mice'] = mice
                        CellAssembly_GlobalResults.loc[counter2, 'NeuronType'] = NeuronType
                            
                        CellAssembly_GlobalResults.loc[counter2, 'SessionType'] =  session_type                        
                        CellAssembly_GlobalResults.loc[counter2, 'TrialType'] = trial_type
                        CellAssembly_GlobalResults.loc[counter2, 'Session_Date'] = session_date 
                        CellAssembly_GlobalResults.loc[counter2, 'Day#'] = day 
                        
                        CellAssembly_GlobalResults.loc[counter2, 'Trial_Time'] = trial_time 
                        CellAssembly_GlobalResults.loc[counter2, 'Trial#'] = trial                    
                        CellAssembly_GlobalResults.loc[counter2, 'Trial_Duration_in_sec'] = len(Carray)/minian_freq                    

                        CellAssembly_GlobalResults.loc[counter2, 'Assembly_ID'] = assembly_ID
                        CellAssembly_GlobalResults.loc[counter2, 'Assembly_size'] = len(non_nan_indices)
                        CellAssembly_GlobalResults.loc[counter2, 'Cells_in_Assembly'] = cells_in_assembly                    

                        CellAssembly_GlobalResults.loc[counter2, 'Avg_Activity'] = mean_act_ass
                        CellAssembly_GlobalResults.loc[counter2, 'EventFreq'] = len(peaks)/(len(Carray)/minian_freq)
                        CellAssembly_GlobalResults.loc[counter2, 'EventTime'] = str(np.round(peaks/minian_freq, 2))

                        counter2+=1
            print(f"... number of assembly = {assembly_nb}")

            for unit in range(nb_unit): 
                indexMapp = np.where(mapping_sess[trial_time] == Calcium.index[unit])[0]
                
                if len(indexMapp)>0 : # The neuron needs to be in the cross-registration
                    
                    Carray_unit =Carray[:,unit]
                    Darray_unit =Darray[:,unit]
                    Sarray_unit =Sarray[:,unit]                                   

                    NeuronActivity_GlobalResults.loc[counter, 'Mice'] = mice
                    NeuronActivity_GlobalResults.loc[counter, 'NeuronType'] = NeuronType
                    
                    NeuronActivity_GlobalResults.loc[counter, 'SessionType'] =  session_type                        
                    NeuronActivity_GlobalResults.loc[counter, 'TrialType'] = trial_type
                    NeuronActivity_GlobalResults.loc[counter, 'Session_Date'] = session_date 
                    NeuronActivity_GlobalResults.loc[counter, 'Day#'] = day 
                    NeuronActivity_GlobalResults.loc[counter, 'Trial_Time'] = trial_time  
                    NeuronActivity_GlobalResults.loc[counter, 'Trial#'] = trial
                    NeuronActivity_GlobalResults.loc[counter, 'Trial_Duration_in_sec'] = len(Carray)/minian_freq  
                    indexMapp = np.where(indexMappList == Calcium.index[unit])[0]
                    NeuronActivity_GlobalResults.loc[counter, 'Unit#_in_crossreg'] = f"{mice}{str(indexMapp).replace('[','').replace(']','')}"
                    NeuronActivity_GlobalResults.loc[counter, 'Unit#_in_session'] = Calcium.index[unit]

                    centroids_sess=centroids[centroids['session']==trial_time]
                    NeuronActivity_GlobalResults.loc[counter, 'UnitLocation'] = [centroids_sess[centroids_sess['unit_id']==Calcium.index[unit]]['height'].tolist(), centroids_sess[centroids_sess['unit_id']==Calcium.index[unit]]['width'].tolist()] 

                    NeuronActivity_GlobalResults.loc[counter, 'NormalizedAUC_calcium'] = np.trapz(Carray_unit,np.arange(0,len(Carray_unit),1))/len(Carray_unit)
                    NeuronActivity_GlobalResults.loc[counter, 'DeconvSpikeMeanActivity'] = Darray_unit.mean()
                    NeuronActivity_GlobalResults.loc[counter, 'SpikeActivityHz'] = Sarray_unit.sum()/(len(Sarray_unit)/minian_freq)
                    
                    counter+=1
           
    if saveexcel:  
        excel_writerRawCa.close()
        excel_writerRawSp.close()

filenameOut = folder_to_save / f'NeuronActivity_Global.pkl'
with open(filenameOut, 'wb') as pickle_file:
    pickle.dump(NeuronActivity_GlobalResults, pickle_file)

filenameOut = folder_to_save / f'CellAssembly_Global.pkl'
with open(filenameOut, 'wb') as pickle_file:
    pickle.dump(CellAssembly_GlobalResults, pickle_file)

if saveexcel: 
    filenameOut = folder_to_save / f'NeuronActivity_Global.xlsx'
    writer = pd.ExcelWriter(filenameOut)
    NeuronActivity_GlobalResults.to_excel(writer)
    writer.close()
    filenameOut = folder_to_save / f'CellAssembly_Global.xlsx'
    writer = pd.ExcelWriter(filenameOut)
    CellAssembly_GlobalResults.to_excel(writer)
    writer.close()

sys.stdout = sys.__stdout__
logfile.close()
