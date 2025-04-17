# # Associate Ca2+ signal with sleep stages for each session & subsessions using crossregistration

#######################################################################################
                            # Define Experiment type #
#######################################################################################

DrugExperiment=0 # =1 if CGP Experiment // DrugExperiment=0 if Baseline Experiment

suffix='' 
AnalysisID='_CellAssembly' 

saveexcel=1

dir = "//10.69.168.1/crnldata/waking/audrey_hay/L1imaging/Analysed2025_AB/"

mapp = {1: 'AW',  2: 'QW', 3: 'NREM',  4: 'IS', 5: 'REM',  6: 'undefined'}
drugs=['baseline', 'CGP']

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

import warnings
warnings.filterwarnings("ignore")

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
        ica = FastICA(n_components=nassemblies)
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
        return
    if significance.nassemblies<1:
        print('WARNING !')
        print('    no assembly detecded!')
        patterns = []
        zactmat = []
    else:
        patterns_ = extractPatterns(zactmat_,significance,method)
        if patterns_ is np.nan:
            return
        patterns = np.zeros((np.size(patterns_,0),nneurons))
        patterns[:,~silentneurons] = patterns_
        zactmat = np.copy(actmat)
        zactmat[~silentneurons,:] = zactmat_
    return patterns,significance,zactmat

def computeAssemblyActivity(patterns,zactmat,zerodiag = True):
    if len(patterns) == 0:
        print('WARNING !')
        print('    no assembly detecded!')
        assemblyAct = []
        return assemblyAct
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

all_expe_types=['preCGP', 'postCGP'] if DrugExperiment else ['baseline', 'preCGP']

# Get the current date and time
FolderNameSave=str(datetime.now())[:19]
FolderNameSave = FolderNameSave.replace(" ", "_").replace(".", "_").replace(":", "_")

destination_folder= f"//10.69.168.1/crnldata/waking/audrey_hay/L1imaging/Analysed2025_AB/_CGP_analysis/VigSt_{FolderNameSave}{suffix}{AnalysisID}" if DrugExperiment else f"//10.69.168.1/crnldata/waking/audrey_hay/L1imaging/Analysed2025_AB/_baseline_analysis/VigSt_{FolderNameSave}{suffix}{AnalysisID}"
os.makedirs(destination_folder)
folder_to_save=Path(destination_folder)

# Copy the script file to the destination folder
source_script = "C:/Users/Manip2/SCRIPTS/CodePythonAudrey/CodePythonAurelie/HayLabAnalysis/python/_MINI&OE_1_compute_vigstates_activity.py"
destination_file_path = f"{destination_folder}/_MINI&OE_1_compute_vigstates_activity.txt"
shutil.copy(source_script, destination_file_path)


data = {}
counter=0
VigilanceState_GlobalResults= pd.DataFrame(data, columns=['Mice','NeuronType','Session', 'Session_Date', 'Session_Time', 'Unique_Unit','UnitNumber','UnitValue', 'UnitLocation',
                                                        'ExpeType', 'Drug', 'Substate','SubstateNumber','DurationSubstate', 'CalciumActivity', 
                                                        'Avg_CalciumActivity', 'AUC_calcium','Avg_AUC_calcium', 'DeconvSpikeMeanActivity', 
                                                        'Avg_DeconvSpikeActivity', 'SpikeActivityHz', 'Avg_SpikeActivityHz', 'TotCaPopCoupling', 
                                                        'TotZ_CaPopCoupling', 'TotSpPopCoupling', 'TotZ_SpPopCoupling'])
counter2=0
CellAssembly_GlobalResults= pd.DataFrame(data, columns=['Mice','NeuronType','Session', 'Session_Date', 'Session_Time', 'Assembly_ID', 'Assembly_size', 'Cells_in_Assembly',
                                                            'ExpeType', 'Drug', 'Substate', 'Avg_Activity' ])

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

        if any(p in all_expe_types for p in minianpath.parts): # have to be to the expe_types

            session=minianpath.parents[0].name if len(minianpath.parts)==12 else minianpath.parents[1].name.split("_")[-1]
            session_path=minianpath.parents[2] if len(minianpath.parts)==12 else minianpath.parents[1]
            expe_type=minianpath.parents[3].name if len(minianpath.parts)==12 else minianpath.parents[2].name

            session_date= minianpath.parents[2].name.split("_")[0] if len(minianpath.parts)==12 else minianpath.parents[1].name.split("_")[0]
            session_time= minianpath.parents[2].name.split("_")[1] if len(minianpath.parts)==12 else minianpath.parents[1].name.split("_")[1]

            drug='CGP' if expe_type == 'postCGP' else 'baseline'
            dict_Path[session] = session_path
        
            minian_ds = open_minian(minianpath)
            dict_Calcium[session] = minian_ds['C'] # calcium traces 
            dict_Deconv[session] = minian_ds['S'] # estimated spikes deconvolved activity

            dict_Scoring[session]  = pd.read_csv(session_path/ f'OpenEphys/Sleep_Scoring_6Stages_5sEpoch.csv')
            dict_Stamps[session]  = pd.read_excel(session_path/ f'SynchroFileCorrect.xlsx')
            dict_StampsMiniscope[session]  = pd.read_csv(list(session_path.glob('*V4_Miniscope/timeStamps.csv'))[0])

            with open(minianpath / f'TodropFileAB.json', 'r') as f:
                unit_to_drop = json.load(f)
                dict_TodropFile[session]  = unit_to_drop

            nb_minian_total+=1

            #######################################################################################
            # Distribute Ca2+ intensity & spikes to vigilance states for each sessions #
            #######################################################################################

            # Start time & freq miniscope

            StartTime = list(dict_Stamps[session][0])[0]
            tsmini=dict_StampsMiniscope[session]['Time Stamp (ms)']
            minian_freq=round(1/np.mean(np.diff(np.array(tsmini)/1000)))

            freqLFP=1000
            
            # Adjust the StartTime if subsessions

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
                kept_uniq_unit_List.append(str(indexMapp))

            nb_unit=len(CalciumSub)
            if nb_unit==0:
                continue  # next iteration
            
            # Realigned the traces to the recorded timestamps 

            timestamps =  np.array(tsmini[firstframe:firstframe+len(CalciumSub.T)])/freqLFP
            new_timestamps= np.arange(timestamps[0], timestamps[-1], 1/minian_freq)
            Calcium = pd.DataFrame(index=CalciumSub.index, columns=new_timestamps)
            for feature in CalciumSub.index:
                interpolator = interpolate.interp1d(timestamps, CalciumSub.loc[feature], kind='linear')
                Calcium.loc[feature] = interpolator(new_timestamps)

            new_timestamps= np.arange(timestamps[0], timestamps[-1], 1/minian_freq)
            Deconv = pd.DataFrame(index=DeconvSub.index, columns=new_timestamps)
            for feature in DeconvSub.index:
                interpolator = interpolate.interp1d(timestamps, DeconvSub.loc[feature], kind='linear')
                Deconv.loc[feature] = interpolator(new_timestamps)

            Carray=Calcium.values.T.astype(float) # Calcium activity
            Darray=Deconv.values.T.astype(float) # Deconvolved activity
            Sarray= np.zeros((np.shape(Darray)[0], np.shape(Darray)[1])) # Spike activity
            for i in np.arange(np.shape(Darray)[1]):
                Darray_unit =Darray[:,i]
                peaks, _ = find_peaks(Darray_unit)
                Sarray_unit=np.zeros(len(Darray_unit))
                Sarray_unit[peaks]=1   
                Sarray[:,i]=Sarray_unit 

            firstframe+=len(CalciumSub.T)
            rec_dur=len(CalciumSub.T)
            rec_dur_sec= timestamps[-1]- timestamps[0]            


            # Deal with dropped frames (failure to acquire miniscope images) IGNORE CAUSE TIMESTAMPS TAKEN INTO ACCOUNT

            list_droppedframes = literal_eval(dict_Stamps[session][0][3])    
            numbdropfr= 0   
            droppedframes_inrec=[]
            for item in list_droppedframes: 
                if item < (round(StartTimeMiniscope*minian_freq) + rec_dur_sec) and item > round(StartTimeMiniscope*minian_freq):
                    numbdropfr+=1                        

            EndTime = StartTime + rec_dur_sec # (upd_rec_dur/minian_freq) # in seconds
            previousEndTime=EndTime 

            print(session, ': starts at', round(StartTime,1), 's & ends at', round(EndTime,1), 's (', round(rec_dur_sec,1), 's duration, ', numbdropfr, 'dropped frames, minian frequency =', minian_freq, 'Hz, experiment type = ', expe_type, ')...') 
            sentence1= f"... kept values = {kept_uniq_unit_List}"
            print(sentence1)
        

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
            #SleepScoredTS_upscaled_ministart[SleepScoredTS_upscaled_ministart == 0.5] = 0

            # Determine each substate identity and duration
            SleepScoredTS_upscaled_ministart=SleepScoredTS_upscaled_ministart.astype(int)
            substates_duration = [len(list(group)) for key, group in groupby(SleepScoredTS_upscaled_ministart)]
            substates_identity = [key for key, _ in groupby(SleepScoredTS_upscaled_ministart)]
            substates_end = np.array(substates_duration).cumsum()        
            substates_start =np.append([0],substates_end[:-1]) #substates_start =np.append([1],substates_end+1) create 1 second gap
            substates_identity = [mapp[num] for num in substates_identity]
            substates = pd.DataFrame(list(zip(substates_identity, substates_duration, substates_start, substates_end)), columns=['Identity', 'Duration', 'Start','End'])


            # Define cell assemblies

            patterns,significance,zactmat= runPatterns(Carray.T, method='ica', nullhyp = 'mp', nshu = 1000, percentile = 99, tracywidom = False)
            print(len(patterns), 'assemblies found')
            if len(patterns)>0:
                thresh = np.mean(patterns)+2*np.std(patterns)
                patterns_th=patterns.copy()
                patterns_th[patterns_th<thresh]=np.nan
                for ass in np.arange(np.shape(patterns)[0]):
                    non_nan_indices = np.where(~np.isnan(patterns_th[ass]))[0] 
                    if len(non_nan_indices)>1:
                        indexMappList=mapping_sess[session]
                        cells_in_assembly=[]
                        for unit in non_nan_indices:
                            indexMapp = np.where(indexMappList == Calcium.index[unit])[0]
                            if len(indexMapp)>0:
                                cells_in_assembly.append(f"{mice}{str(indexMapp).replace('[','').replace(']','')}")
                        assembly_activity = zscore(Calcium.iloc[non_nan_indices].values.astype(float).T).mean(axis=1)
                        assembly_nb+=1
                        assembly_ID=f'Assembly_{mice}_{assembly_nb}'
                        print(len(cells_in_assembly))

                        for m in mapp:                   
                            Bool = (SleepScoredTS_upscaled_ministart == m)
                            assembly_activity_VigSpe = assembly_activity.copy()
                            assembly_activity_VigSpe = assembly_activity_VigSpe[0:np.shape(SleepScoredTS_upscaled_ministart)[0]] # if Calcium imaging longer than LFP rec
                            mean_act_ass = np.mean(assembly_activity_VigSpe[Bool])

                            CellAssembly_GlobalResults.loc[counter2, 'Mice'] = mice
                            CellAssembly_GlobalResults.loc[counter2, 'NeuronType'] = NeuronType
                            
                            CellAssembly_GlobalResults.loc[counter2, 'Session'] = session
                            CellAssembly_GlobalResults.loc[counter2, 'Session_Date'] = session_date 
                            CellAssembly_GlobalResults.loc[counter2, 'Session_Time'] = session_time                    

                            CellAssembly_GlobalResults.loc[counter2, 'Assembly_ID'] = assembly_ID
                            CellAssembly_GlobalResults.loc[counter2, 'Assembly_size'] = len(non_nan_indices)
                            CellAssembly_GlobalResults.loc[counter2, 'Cells_in_Assembly'] = cells_in_assembly
                            
                            CellAssembly_GlobalResults.loc[counter2, 'ExpeType'] =  expe_type
                            CellAssembly_GlobalResults.loc[counter2, 'Drug'] =  drug

                            CellAssembly_GlobalResults.loc[counter2, 'Substate'] = mapp[m]
                            CellAssembly_GlobalResults.loc[counter2, 'Avg_Activity'] = mean_act_ass

                            counter2+=1

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
                    RawCaTraces=pd.DataFrame(Carray_VigSpe, columns=[f"{mice}{str(i).replace('[','').replace(']','')}" for i in kept_uniq_unit_List])
                    unique_columns = RawCaTraces.columns[RawCaTraces.columns.to_series().duplicated()] # remove units that has an empty unique index '[]'
                    RawCaTraces = RawCaTraces.drop(columns=unique_columns)
                    RawCaTracesVigStateMatrix.append(RawCaTraces)                

                    RawSpTracesVigStateMatrix = locals()[f'RawSpTraces{mapp[m]}_{drug}']
                    RawSpTraces=[]
                    RawSpTraces=pd.DataFrame(Sarray_VigSpe, columns=[f"{mice}{str(i).replace('[','').replace(']','')}" for i in kept_uniq_unit_List])
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

                        centroids_sess=centroids[centroids['session']==session]
                        VigilanceState_GlobalResults.loc[counter, 'UnitLocation'] = [centroids_sess[centroids_sess['unit_id']==Calcium.index[unit]]['height'].tolist(), centroids_sess[centroids_sess['unit_id']==Calcium.index[unit]]['width'].tolist()] 

                        VigilanceState_GlobalResults.loc[counter, 'ExpeType'] =  expe_type
                        VigilanceState_GlobalResults.loc[counter, 'Drug'] =  drug

                        VigilanceState_GlobalResults.loc[counter, 'Substate'] = substates.Identity[index]
                        VigilanceState_GlobalResults.loc[counter, 'SubstateNumber'] = substates.index[index]
                        VigilanceState_GlobalResults.loc[counter, 'DurationSubstate'] = (substates.Duration[index]/minian_freq)

                        VigilanceState_GlobalResults.loc[counter, 'CalciumActivity'] = ca_input_sub.mean()
                        VigilanceState_GlobalResults.loc[counter, 'Avg_CalciumActivity'] = Carray_unit.mean()

                        VigilanceState_GlobalResults.loc[counter, 'AUC_calcium'] = np.trapz(ca_input_sub,np.arange(0,len(ca_input_sub),1))
                        VigilanceState_GlobalResults.loc[counter, 'Avg_AUC_calcium'] = np.trapz(Carray_unit,np.arange(0,len(Carray_unit),1))

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
                dataCaCorr[f'Z_{drug}_{mapp[m]}']=combined_df
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
                    combined_df.to_excel(excel_writerRawSp, sheet_name=f'{drug}_{mapp[m]}', index=False, header=True)   
        
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



filenameOut = folder_to_save / f'VigStates_Global.pkl'
with open(filenameOut, 'wb') as pickle_file:
    pickle.dump(VigilanceState_GlobalResults, pickle_file)

filenameOut = folder_to_save / f'CellAssembly_Global.pkl'
with open(filenameOut, 'wb') as pickle_file:
    pickle.dump(CellAssembly_GlobalResults, pickle_file)

if saveexcel: 
    filenameOut = folder_to_save / f'VigStates_Global.xlsx'
    writer = pd.ExcelWriter(filenameOut)
    VigilanceState_GlobalResults.to_excel(writer)
    writer.close()
    filenameOut = folder_to_save / f'CellAssembly_Global.xlsx'
    writer = pd.ExcelWriter(filenameOut)
    CellAssembly_GlobalResults.to_excel(writer)
    writer.close()