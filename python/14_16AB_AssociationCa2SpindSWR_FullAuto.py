# # Associate Ca2+ signal with spindles for each session & subsessions using crossregistration

# Load packages

import os

import quantities as pq
import numpy as np
import math 
import neo
import json
from pathlib import Path
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, Cursor
from scipy.interpolate import interp2d
from scipy.signal import find_peaks
import pickle
import os
from scipy.interpolate import griddata
import logging
import sys 

from itertools import groupby
from ephyviewer import mkQApp, MainViewer, TraceViewer
from IPython.display import display
from ipyfilechooser import FileChooser


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

# Perform analysis for each mouse

MiceList=['BlackLinesOK', 'BlueLinesOK', 'GreenDotsOK', 'GreenLinesOK', 'Purple', 'RedLinesOK','ThreeColDotsOK', 'ThreeBlueCrossesOK']

for micename in MiceList:

    # Load sleep score and Ca2+ time series numpy arrays

    dpath0 = "//10.69.168.1/crnldata/waking/audrey_hay/L1imaging/AnalysedMarch2023/Gaelle/Baseline_recording/"
    dpath=dpath0 + micename
    print(dpath)
    folder_base = Path(dpath)

    nb_sessions = sum(1 for p in folder_base.iterdir() if p.is_dir() and p.name.startswith("session"))
    try:
        mfile = open(folder_base / f'mappingsAB.pkl', 'rb')
        mapping = pickle.load(mfile)
        print('mappingsAB.pkl opened')
    except:
        mfile = open(folder_base / f'mappings.pkl', 'rb')
        mapping = pickle.load(mfile)
        print('mappings.pkl opened')

    sessions = []
    subsessions = []
    nb_minian_total=0
    dict_Calcium = {}
    dict_Spike = {}
    dict_SWRprop = {}
    dict_Spindleprop_PFC = {}
    dict_Spindleprop_S1 = {}
    dict_Stamps = {}
    dict_StampsMiniscope = {}
    dict_TodropFile = {}

    for y in range(1, nb_sessions+1):
        session= 'session' + str(y)
        print(session)
        sessions.append(session)
        folder_mini = folder_base / f'session{y}/V4_Miniscope'
        nb_subsessions = sum(1 for p in folder_mini.iterdir() if p.is_dir() and p.name.startswith("session"))
        SWRproperties = folder_base / f'session{y}/OpenEphys/SWRproperties_AB.csv'
        Spindleproperties_PFC = folder_base / f'session{y}/OpenEphys/Spindlesproperties_PFC_AB.csv'
        Spindleproperties_S1 = folder_base / f'session{y}/OpenEphys/Spindlesproperties_S1_AB.csv'
        StampsFile = folder_base / f'session{y}/SynchroFile.xlsx'
        StampsMiniscopeFile = folder_mini / f'timeStamps.csv'
        if nb_subsessions!=0:
            for x in range(1, nb_subsessions+1):            
                subsession= "session"  + str(y) + str(x)
                subsessions.append(subsession)    
                minian_ds = open_minian(folder_mini / subsession / f'minian')      # OR minianAB
                dict_Calcium[subsession] = minian_ds['C'] # calcium traces 
                dict_Spike[subsession] = minian_ds['S'] # estimated spikes
                dict_SWRprop[subsession]  = pd.read_csv(SWRproperties)
                dict_Spindleprop_PFC[subsession]  = pd.read_csv(Spindleproperties_PFC)
                dict_Spindleprop_S1[subsession]  = pd.read_csv(Spindleproperties_S1)
                dict_Stamps[subsession]  = pd.read_excel(StampsFile)
                dict_StampsMiniscope[subsession]  = pd.read_csv(StampsMiniscopeFile)
                try:
                    TodropFile = folder_mini / subsession / f'minian/TodropFileAB.json'
                    with open(TodropFile, 'r') as f:
                        unit_to_drop = json.load(f)
                        dict_TodropFile[subsession]  = unit_to_drop
                except:
                    TodropFile = folder_mini / subsession / f'minian/TodropFile.json'
                    with open(TodropFile, 'r') as f:
                        unit_to_drop = json.load(f)
                        dict_TodropFile[subsession]  = unit_to_drop
                nb_minian_total+=1
                print(nb_minian_total)
        else:
            minian_ds = open_minian(folder_mini / f'minian')            # OR minianAB
            dict_Calcium[session] = minian_ds['C'] # calcium traces 
            dict_Spike[session] = minian_ds['S'] # estimated spikes
            dict_SWRprop[session]  = pd.read_csv(SWRproperties)
            dict_Spindleprop_PFC[session]  = pd.read_csv(Spindleproperties_PFC)
            dict_Spindleprop_S1[session]  = pd.read_csv(Spindleproperties_S1)
            dict_Stamps[session]  = pd.read_excel(StampsFile)
            dict_StampsMiniscope[session]  = pd.read_csv(StampsMiniscopeFile)
            try:
                TodropFile = folder_mini / f'minian/TodropFileAB.json'
                with open(TodropFile, 'r') as f:
                    unit_to_drop = json.load(f)
                    dict_TodropFile[session]  = unit_to_drop
            except:
                TodropFile = folder_mini / f'minian/TodropFile.json'
                with open(TodropFile, 'r') as f:
                    unit_to_drop = json.load(f)
                    dict_TodropFile[session]  = unit_to_drop
            nb_minian_total+=1  
            print(nb_minian_total)

    # Cross registration results

    B = mapping['session']
    if os.path.basename(folder_base) == 'Purple':
        index = B.columns
        B.columns = index.str.replace('part', 'session2')

    # Define functions

    from bisect import bisect_left

    def take_closest(myList, myNumber):
        """
        Assumes myList is sorted. Returns closest value to myNumber.
        If two numbers are equally close, return the smallest number.
        """
        pos = bisect_left(myList, myNumber)
        if pos == 0:
            return 0
        if pos == len(myList):
            return len(myList)
        before = myList[pos - 1]
        after = myList[pos]
        if after - myNumber < myNumber - before:
            return after
        else:
            return before

    def take_closest2(myList, myNumber):
        value2 = 10000000
        for ind in range(len(myList)):
            value = abs(myList[ind]-myNumber)
            if value < value2:
                value2 = value
                index = myList[ind]
        return index


    def take_closest3(myList, myNumber):
        pos = bisect_left(myList, myNumber)
        if pos == 0:
            return myList[0]
        if pos == len(myList):
            return myList[-1]
        before = myList[pos - 1]
        after = myList[pos]
        if after - myNumber < myNumber - before:
            dummy = myList.index(after)
            return dummy
        else:
            dummy = myList.index(before)
            return dummy
        

    def is_between(myList, starttime, endtime):
        IsTrue='False'
        for ind in range(len(myList)):
            if starttime <= myList[ind] <= endtime:
                IsTrue='True'
        return IsTrue

    # Distribute Ca2+ intensity & spikes to vigilance state for each sessions/subsessions

    CortexList= ['PFC', 'S1']

    for Cortex in CortexList:

        data = {}
        
        before = 1000 # Max distance in ms between a SWR and a spindle to be considered as Precoupled
        after = 1000 # Max distance in ms between a spindle and a SWR to be considered as Postcoupled
        durationSpdl = 2 # number of sec before and after the Spdl onset taken into acount
        durationSWR = 0.5 # number of sec before and after the Spdl onset taken into acount
        counter=0
        counter2=0

        norm_freq=10 # final miniscope frequency used for all recordings

        Spindles_GlobalResults= pd.DataFrame(data, columns=['Mice', 'Session','Session_Time','Unique_Unit','UnitNumber','UnitValue','SpdlStatut','SpdlNumber','SpdlDuration (ms)','SWR inside Spdl','CalciumActivityPreference', 'CalciumActivityBefore','CalciumActivityAfter','AUC_calciumBefore','AUC_calciumAfter','SpikeActivityPreference','SpikeActivityBefore','SpikeActivityAfter','AUC_spikeBefore', 'AUC_spikeAfter'])
        SWR_GlobalResults= pd.DataFrame(data, columns=['Mice', 'Session','Session_Time','Unique_Unit','UnitNumber','UnitValue','SWRStatut','SWRNumber','SWRDuration (ms)','SWR inside Spdl','CalciumActivityPreference', 'CalciumActivityBefore','CalciumActivityAfter','AUC_calciumBefore','AUC_calciumAfter','SpikeActivityPreference','SpikeActivityBefore','SpikeActivityAfter','AUC_spikeBefore', 'AUC_spikeAfter'])

        dict_All_ActivityCa_spin={}
        dict_All_ActivityCa_spin_Precoupled={}
        dict_All_ActivityCa_spin_Postcoupled={}
        dict_All_ActivityCa_spin_Uncoupled={}

        dict_All_ActivityCa_swr={}
        dict_All_ActivityCa_swr_Precoupled={}
        dict_All_ActivityCa_swr_Postcoupled={}
        dict_All_ActivityCa_swr_Uncoupled={}

        dict_All_ActivitySp_spin={}
        dict_All_ActivitySp_spin_Precoupled={}
        dict_All_ActivitySp_spin_Postcoupled={}
        dict_All_ActivitySp_spin_Uncoupled={}

        dict_All_ActivitySp_swr={}
        dict_All_ActivitySp_swr_Precoupled={}
        dict_All_ActivitySp_swr_Postcoupled={}
        dict_All_ActivitySp_swr_Uncoupled={}

        previousEndTime=0
        InitialStartTime=0

        for i in list(dict_Stamps.keys()):    
            cPreCoupled=0
            cPostCoupled=0
            cUnCoupled=0

            cPreCoupledSWR=0
            cPostCoupledSWR=0
            cUnCoupledSWR=0
            
            # Start time & freq miniscope
            StartTime = list(dict_Stamps[i][0])[0] # in seconds
            minian_freq=list(dict_Stamps[i][0])[2] # in Hz

            # start time session 2
            def Convert(string):
                li = list(string.split(", "))
                li2 = len(li)
                return li2
            stri = dict_Stamps[i][0][3]
            numbdropfr = Convert(stri)

            from ast import literal_eval
            list_droppedframes = literal_eval(dict_Stamps[i][0][3])

            First_frame = 0 #StartTime*minian_freq 
            C=dict_Calcium[i]
            Cupd = C.loc[:, :] #C.loc[:, First_frame:]
            rec_dur = Cupd.shape[1]
            S=dict_Spike[i] 
            Supd = S.loc[:, :] #S.loc[:, :]

            if InitialStartTime==0:
                InitialStartTime=StartTime    
            else:
                if StartTime == InitialStartTime:
                    StartTime = previousEndTime + 1/minian_freq #  +1 frame in seconds 
                else:  
                    InitialStartTime=StartTime   

            if len(list_droppedframes) > 0:
                numbdropfr = sum(1 for item in list_droppedframes if item < (int(StartTime*minian_freq) + rec_dur) and item > int(StartTime*minian_freq))
            else:
                numbdropfr = 0   

            EndTime = StartTime + ((rec_dur + numbdropfr)/minian_freq) # in seconds
            previousEndTime=EndTime     

            print(i, ': starts at', round(StartTime,1), 's & ends at', round(EndTime,1), 's (', round((rec_dur + numbdropfr)/minian_freq,1), 's duration, ', numbdropfr, 'dropped frames, minian frequency =', minian_freq, ')...') 
            
            AA = C['unit_id']
            copyAA = list(AA.copy())
            unit_to_drop=dict_TodropFile[i]    
            for u in unit_to_drop:
                copyAA.remove(u)
            unit_to_keep = copyAA
            Cupd = Cupd.loc[unit_to_keep,:]
            Supd = Supd.loc[unit_to_keep,:]
            nb_unit = Cupd.shape[0]
            units = range(nb_unit)
            
            Cseries = Cupd.to_series()
            C_upd_unit_id = Cupd['unit_id'].values
            Sseries = Supd.to_series()
            #sentence1= f"... kept values = {C_upd_unit_id}"

            kept_uniq_unit_List=[]
            for unit in units:
                indexMapp = np.where(B[i] == C_upd_unit_id[unit])[0]
                kept_uniq_unit_List.append(str(indexMapp))

            sentence1= f"... kept values = {kept_uniq_unit_List}"
            print(sentence1) 
            
            dictnameSpdl=f'dict_Spindleprop_{Cortex}' #change dict according to the cortex 
            dictSpdl = globals()[dictnameSpdl]

            SpipropO=dictSpdl[i]
            SpipropM=SpipropO.copy()
            SWRpropO=dict_SWRprop[i]
            SWRpropM=SWRpropO.copy()
            SpipropM[["peak time", "start time", "end time"]] = SpipropM[["peak time", "start time", "end time"]]-(StartTime*1000)
            SWRpropM[["peak time", "start time", "end time"]] = SWRpropM[["peak time", "start time", "end time"]]-(StartTime*1000)        

            timeSpdl = range(int(durationSpdl*2*minian_freq))
            HalfSpdl = int(durationSpdl*minian_freq)
            
            timeSWR = range(int(durationSWR*2*minian_freq))
            HalfSWR = int(durationSWR*minian_freq)

            TimeStamps_miniscope=list(dict_StampsMiniscope[i]["Time Stamp (ms)"]) # + (StartTime*1000))

            SpipropTrunc = SpipropM[SpipropM["start time"]>0]
            SpipropTrunc = SpipropTrunc[SpipropTrunc["start time"]< (EndTime-StartTime)*1000]
            SWRpropTrunc = SWRpropM[SWRpropM["start time"]>0]
            SWRpropTrunc = SWRpropTrunc[SWRpropTrunc["start time"] < (EndTime-StartTime)*1000]
        
            nb_spindle = SpipropTrunc.shape[0]
            nb_swr = SWRpropTrunc.shape[0]

            for ii, unit in enumerate(units): # for each kept units (cause Cseries/Sseries only have kept units)

                lCseries = np.array(Cseries)[(unit)*(rec_dur + numbdropfr):(unit+1)*(rec_dur + numbdropfr)]
                lSseries = np.array(Sseries)[(unit)*(rec_dur + numbdropfr):(unit+1)*(rec_dur + numbdropfr)]
                ActivityCa_Spin = [] #For each unit  
                ActivityCa_Spin_Precoupled= [] #For each unit 
                ActivityCa_Spin_Postcoupled= [] #For each unit 
                ActivityCa_Spin_Uncoupled= [] #For each unit 
                ActivitySp_Spin = [] #For each unit  
                ActivitySp_Spin_Precoupled= [] #For each unit 
                ActivitySp_Spin_Postcoupled= [] #For each unit 
                ActivitySp_Spin_Uncoupled= [] #For each unit 
                startSpiList = list(pd.Series(SpipropTrunc["start time"]))
                endSpiList = list(pd.Series(SpipropTrunc["end time"]))

                for Pspin in range(nb_spindle): 
                    
                    # Get the calcium and spike trace associated with the spdl
                    startSpi=startSpiList[Pspin]
                    endSpi=endSpiList[Pspin]
                    Frame_Spindle_start = take_closest2(TimeStamps_miniscope, startSpi)
                    index = TimeStamps_miniscope.index(Frame_Spindle_start)
                    CaTrace = list(lCseries[index-HalfSpdl:index+HalfSpdl])
                    SpTraceO = list(lSseries[index-HalfSpdl:index+HalfSpdl]) 
                    peaks, _ = find_peaks(SpTraceO, height=np.std(SpTraceO))
                    SpTrace=np.zeros(len(SpTraceO))
                    SpTrace[peaks]=1

                    if len(CaTrace)<len(timeSpdl): 
                        print(" /!\ Spindle too close to the begining/end of the recording,", i, ", Spdl n°", Pspin, ", Start Spdl =", int(startSpi), "ms, Start Rec =", int(Frame_Spindle_start), 'ms') if ii==0 else None            
                    else:
                        ActivityCa_Spin.append(CaTrace)
                        ActivitySp_Spin.append(SpTraceO)

                        # Define if that spindle is coupled with a SWR or not
                        Spdl_statut=[]
                        startSWRList = list(pd.Series(SWRpropTrunc["start time"]))
                        if len(startSWRList)>0:
                            startClosest_SWR = take_closest2(startSWRList, startSpi)
                            distance = startClosest_SWR - startSpi
                            if (distance > (- before)) and (distance <  0):
                                Spdl_statut = ['PreCoupled']
                                cPreCoupled+=1 if ii==1 else 0
                                ActivityCa_Spin_Precoupled.append(CaTrace)
                                ActivitySp_Spin_Precoupled.append(SpTraceO)
                            elif (distance > (0)) and (distance <  after):
                                Spdl_statut = ['PostCoupled']
                                cPostCoupled+=1 if ii==1 else 0
                                ActivityCa_Spin_Postcoupled.append(CaTrace)
                                ActivitySp_Spin_Postcoupled.append(SpTraceO)
                            else:
                                Spdl_statut= ['UnCoupled']
                                cUnCoupled+=1 if ii==1 else 0
                                ActivityCa_Spin_Uncoupled.append(CaTrace)
                                ActivitySp_Spin_Uncoupled.append(SpTraceO)
                        else:
                            Spdl_statut= ['UnCoupled']
                            cUnCoupled+=1 if ii==1 else 0
                            ActivityCa_Spin_Uncoupled.append(CaTrace)
                            ActivitySp_Spin_Uncoupled.append(SpTraceO)

                        Spindles_GlobalResults.loc[counter, 'Mice'] = os.path.basename(folder_base)
                        Spindles_GlobalResults.loc[counter, 'Session'] = i 
                        Spindles_GlobalResults.loc[counter, 'Session_Time'] = None 
                        mapping['session'].columns.tolist()

                        indexMapp = np.where(B[i] == C_upd_unit_id[unit])[0]
                        Spindles_GlobalResults.loc[counter, 'Unique_Unit'] = indexMapp 
                        Spindles_GlobalResults.loc[counter, 'UnitNumber'] = unit 
                        Spindles_GlobalResults.loc[counter, 'UnitValue'] = C_upd_unit_id[unit] 
                        Spindles_GlobalResults.loc[counter, 'SpdlStatut'] = Spdl_statut
                        Spindles_GlobalResults.loc[counter, 'SpdlNumber'] = Pspin
                        Spindles_GlobalResults.loc[counter, 'SpdlDuration (ms)'] = endSpi- startSpi
                        IsTrue=is_between(startSWRList,startSpi, endSpi)
                        Spindles_GlobalResults.loc[counter, 'SWR inside Spdl'] = IsTrue
                        
                        if np.mean(CaTrace[:HalfSpdl],0) > np.mean(CaTrace[HalfSpdl:],0):
                            pref='Before'
                        elif np.mean(CaTrace[:HalfSpdl],0) < np.mean(CaTrace[HalfSpdl:],0):
                            pref='After' 
                        else:
                            pref='None'
                        Spindles_GlobalResults.loc[counter, 'CalciumActivityPreference'] = pref
                        Spindles_GlobalResults.loc[counter, 'CalciumActivityBefore'] = np.mean(CaTrace[:HalfSpdl],0)
                        Spindles_GlobalResults.loc[counter, 'CalciumActivityAfter'] = np.mean(CaTrace[HalfSpdl:],0)
                        Spindles_GlobalResults.loc[counter, 'AUC_calciumBefore'] = np.trapz(CaTrace[:HalfSpdl],np.arange(0,len(CaTrace[:HalfSpdl]),1))
                        Spindles_GlobalResults.loc[counter, 'AUC_calciumAfter'] = np.trapz(CaTrace[HalfSpdl:],np.arange(0,len(CaTrace[HalfSpdl:]),1))          
                        
                        if np.sum(SpTrace[:HalfSpdl],0) > np.sum(SpTrace[HalfSpdl:],0):
                            pref='Before'
                        elif np.sum(SpTrace[:HalfSpdl],0) < np.sum(SpTrace[HalfSpdl:],0):
                            pref='After' 
                        else:
                            pref='None'
                        Spindles_GlobalResults.loc[counter, 'SpikeActivityPreference'] = pref
                        Spindles_GlobalResults.loc[counter, 'SpikeActivityBefore'] = np.sum(SpTrace[:HalfSpdl],0)
                        Spindles_GlobalResults.loc[counter, 'SpikeActivityAfter'] = np.sum(SpTrace[HalfSpdl:],0)
                        counter+=1  

                # All Ca traces for each spindles per Unique unit (according to cross-registration)

                list_ActivityCa= ['ActivityCa_Spin', 'ActivityCa_Spin_Precoupled', 'ActivityCa_Spin_Postcoupled', 'ActivityCa_Spin_Uncoupled']
                list_dict_All_ActivityCa= ['dict_All_ActivityCa_spin', 'dict_All_ActivityCa_spin_Precoupled', 'dict_All_ActivityCa_spin_Postcoupled', 'dict_All_ActivityCa_spin_Uncoupled']
                for it, ActivityCaNames in enumerate(list_ActivityCa): # for each Spdl types
                    if len(indexMapp) > 0: #not empty --> cause some units are not in the cross registration..! Need to know why 

                        ActivityCa = locals()[ActivityCaNames]
                        dict_All_ActivityCa = locals()[list_dict_All_ActivityCa[it]]       
                        if len(ActivityCa)>0 :    
                            if np.shape(np.array(ActivityCa))[1] == int(norm_freq*durationSpdl*2):  #normalize traces to the same frequency rate         
                                ActivityCa= np.reshape(np.array(ActivityCa), (-1, len(np.array(ActivityCa)))) if np.ndim(ActivityCa) == 1 else np.array(ActivityCa)    
                                dict_All_ActivityCa[str(indexMapp)] = np.append(dict_All_ActivityCa[str(indexMapp)], np.array(ActivityCa), axis=0) if str(indexMapp) in dict_All_ActivityCa else np.array(ActivityCa)
                            else:
                                dataO = np.array(ActivityCa)
                                data= np.repeat(dataO, 2, axis=0) if dataO.shape[0] == 1 else dataO
                                x_mesh, y_mesh = np.meshgrid(np.arange(data.shape[1]), np.arange(data.shape[0]))
                                x_new_mesh, y_new_mesh = np.meshgrid(np.linspace(0, data.shape[1] - 1, int(norm_freq*durationSpdl*2)), np.linspace(0, data.shape[0] - 1, np.shape(data)[0]))
                                resampled_dataO = griddata((x_mesh.flatten(), y_mesh.flatten()), data.flatten(), (x_new_mesh, y_new_mesh), method='linear')
                                resampled_data= resampled_dataO[0,:] if dataO.shape[0] == 1 else resampled_dataO
                                resampled_data= np.reshape(resampled_data, (-1, len(resampled_data))) if np.ndim(resampled_data) == 1 else resampled_data
                                dict_All_ActivityCa[str(indexMapp)] = np.append(dict_All_ActivityCa[str(indexMapp)], np.array(resampled_data), axis=0) if str(indexMapp) in dict_All_ActivityCa else np.array(resampled_data)
                        sentence1bis=""
                    else: 
                        sentence1bis=f"/!\ Cell idx {unit} not in the cross registration" if it==1 else ""
                        print(sentence1bis) if it==1 else None

                # All Sp traces for each spindles per Unique unit (according to cross-registration)

                list_ActivitySp= ['ActivitySp_Spin', 'ActivitySp_Spin_Precoupled', 'ActivitySp_Spin_Postcoupled', 'ActivitySp_Spin_Uncoupled']
                list_dict_All_ActivitySp= ['dict_All_ActivitySp_spin', 'dict_All_ActivitySp_spin_Precoupled', 'dict_All_ActivitySp_spin_Postcoupled', 'dict_All_ActivitySp_spin_Uncoupled']
                for it, ActivitySpNames in enumerate(list_ActivitySp): # for each Spdl types
                    if len(indexMapp) > 0: #not empty --> cause some units are not in the cross registration..! Need to know why 

                        ActivitySp = locals()[ActivitySpNames]
                        dict_All_ActivitySp = locals()[list_dict_All_ActivitySp[it]]       
                        if len(ActivitySp)>0 :    
                            if np.shape(np.array(ActivitySp))[1] == int(norm_freq*durationSpdl*2):  #normalize traces to the same frequency rate         
                                ActivitySp= np.reshape(np.array(ActivitySp), (-1, len(np.array(ActivitySp)))) if np.ndim(ActivitySp) == 1 else np.array(ActivitySp)    
                                for cl in range(len(ActivitySp)):
                                    resampled_data_col=ActivitySp[cl]
                                    peaks, _ = find_peaks(resampled_data_col, height=np.std(resampled_data_col))
                                    SpTrace=np.zeros(len(resampled_data_col))
                                    SpTrace[peaks]=1
                                    ActivitySp[cl]=SpTrace                                
                                dict_All_ActivitySp[str(indexMapp)] = np.append(dict_All_ActivitySp[str(indexMapp)], np.array(ActivitySp), axis=0) if str(indexMapp) in dict_All_ActivitySp else np.array(ActivitySp)
                            else:
                                dataO = np.array(ActivitySp)
                                data= np.repeat(dataO, 2, axis=0) if dataO.shape[0] == 1 else dataO
                                x_mesh, y_mesh = np.meshgrid(np.arange(data.shape[1]), np.arange(data.shape[0]))
                                x_new_mesh, y_new_mesh = np.meshgrid(np.linspace(0, data.shape[1] - 1, int(norm_freq*durationSpdl*2)), np.linspace(0, data.shape[0] - 1, np.shape(data)[0]))
                                resampled_dataO = griddata((x_mesh.flatten(), y_mesh.flatten()), data.flatten(), (x_new_mesh, y_new_mesh), method='linear')
                                resampled_data= resampled_dataO[0,:] if dataO.shape[0] == 1 else resampled_dataO
                                resampled_data= np.reshape(resampled_data, (-1, len(resampled_data))) if np.ndim(resampled_data) == 1 else resampled_data
                                for cl in range(len(resampled_data)):
                                    resampled_data_col=resampled_data[cl]
                                    peaks, _ = find_peaks(resampled_data_col, height=np.std(resampled_data_col))
                                    SpTrace=np.zeros(len(resampled_data_col))
                                    SpTrace[peaks]=1
                                    resampled_data[cl]=SpTrace    
                                dict_All_ActivitySp[str(indexMapp)] = np.append(dict_All_ActivitySp[str(indexMapp)], np.array(resampled_data), axis=0) if str(indexMapp) in dict_All_ActivitySp else np.array(resampled_data)
                                
                ActivityCa_swr = [] #For each unit  
                ActivityCa_swr_Precoupled= [] #For each unit 
                ActivityCa_swr_Postcoupled= [] #For each unit 
                ActivityCa_swr_Uncoupled= [] #For each unit 

                ActivitySp_swr = [] #For each unit  
                ActivitySp_swr_Precoupled= [] #For each unit 
                ActivitySp_swr_Postcoupled= [] #For each unit 
                ActivitySp_swr_Uncoupled= [] #For each unit 

                startSwrList = list(pd.Series(SWRpropTrunc["start time"]))
                endSwrList = list(pd.Series(SWRpropTrunc["end time"]))
                for Pswr in range(nb_swr): 

                    # Get the calcium and spike trace associated with the swr

                    startSwr=startSwrList[Pswr]
                    endSwr=endSwrList[Pswr]
                    Frame_SWR_start = take_closest2(TimeStamps_miniscope, startSwr)
                    index = TimeStamps_miniscope.index(Frame_SWR_start)
                    CaTrace = list(lCseries[index-HalfSWR:index+HalfSWR])
                    SpTraceO = list(lSseries[index-HalfSWR:index+HalfSWR]) 
                    peaks, _ = find_peaks(SpTraceO, height=np.std(SpTraceO))
                    SpTrace=np.zeros(len(SpTraceO))
                    SpTrace[peaks]=1

                    if len(CaTrace)<len(timeSWR): 
                        print("/!\ SWR too close to the begining/end of the recording,", i, ", SWR n°", Pswr, ", Start SWR =",  int(startSwr), "ms, Start Rec =", int(Frame_SWR_start), 'ms') if ii==0 else None 
                    else:
                        ActivityCa_swr.append(CaTrace) 
                        ActivitySp_swr.append(SpTraceO) 

                        # Define if that SWR is coupled with a spdl or not

                        SWR_statut=[]
                        startSpiList = list(pd.Series(SpipropTrunc["start time"]))
                        endSpiList = list(pd.Series(SpipropTrunc["end time"]))
                        if len(startSpiList)>0:
                            startClosest_Spi = take_closest2(startSpiList, startSwr)# + StartTimeIndexSpi])
                            indexSpi = startSpiList.index(startClosest_Spi)
                            endClosest_Spi=endSpiList[indexSpi]
                            distance = startClosest_Spi - startSwr #  + StartTimeIndexSpi]  
                            IsTrue = 'False'             
                            if (distance > (- before)) and (distance <  0):
                                SWR_statut = ['Postcoupled']
                                cPostCoupledSWR+=1 if ii==1 else 0
                                ActivityCa_swr_Postcoupled.append(CaTrace)
                                ActivitySp_swr_Postcoupled.append(SpTraceO)
                                if startSwr<endClosest_Spi:
                                    IsTrue = 'True' #SWR inside the Spindle
                            elif (distance > (0)) and (distance <  after):
                                SWR_statut = ['Precoupled']
                                cPreCoupledSWR+=1 if ii==1 else 0
                                ActivityCa_swr_Precoupled.append(CaTrace)
                                ActivitySp_swr_Precoupled.append(SpTraceO)
                            else:
                                SWR_statut= ['UnCoupled']
                                cUnCoupledSWR+=1 if ii==1 else 0
                                ActivityCa_swr_Uncoupled.append(CaTrace)
                                ActivitySp_swr_Uncoupled.append(SpTraceO)
                        else: 
                            SWR_statut= ['UnCoupled']
                            cUnCoupledSWR+=1 if ii==1 else 0
                            ActivityCa_swr_Uncoupled.append(CaTrace)
                            ActivitySp_swr_Uncoupled.append(SpTraceO)

                        SWR_GlobalResults.loc[counter2, 'Mice'] = os.path.basename(folder_base)
                        SWR_GlobalResults.loc[counter2, 'Session'] = i 
                        SWR_GlobalResults.loc[counter2, 'Session_Time'] = None 
                        indexMapp = np.where(B[i] == C_upd_unit_id[unit])[0]
                        SWR_GlobalResults.loc[counter2, 'Unique_Unit'] = indexMapp 
                        SWR_GlobalResults.loc[counter2, 'UnitNumber'] = unit 
                        SWR_GlobalResults.loc[counter2, 'UnitValue'] = C_upd_unit_id[unit] 
                        SWR_GlobalResults.loc[counter2, 'SWRStatut'] = SWR_statut
                        SWR_GlobalResults.loc[counter2, 'SWRNumber'] = Pswr
                        SWR_GlobalResults.loc[counter2, 'SWRDuration (ms)'] = endSwr- startSwr
                        SWR_GlobalResults.loc[counter2, 'SWR inside Spdl'] = IsTrue
                        
                        if np.mean(CaTrace[:HalfSWR],0) > np.mean(CaTrace[HalfSWR:],0):
                            pref='Before'
                        elif np.mean(CaTrace[:HalfSWR],0) < np.mean(CaTrace[HalfSWR:],0):
                            pref='After' 
                        else:
                            pref='None'
                        SWR_GlobalResults.loc[counter2, 'CalciumActivityPreference'] = pref
                        SWR_GlobalResults.loc[counter2, 'CalciumActivityBefore'] = np.mean(CaTrace[:HalfSWR],0)
                        SWR_GlobalResults.loc[counter2, 'CalciumActivityAfter'] = np.mean(CaTrace[HalfSWR:],0)
                        SWR_GlobalResults.loc[counter2, 'AUC_calciumBefore'] = np.trapz(CaTrace[:HalfSWR],np.arange(0,len(CaTrace[:HalfSWR]),1))
                        SWR_GlobalResults.loc[counter2, 'AUC_calciumAfter'] = np.trapz(CaTrace[HalfSWR:],np.arange(0,len(CaTrace[HalfSWR:]),1))          

                        if np.sum(SpTrace[:HalfSWR],0) > np.sum(SpTrace[HalfSWR:],0):
                            pref='Before'
                        elif np.sum(SpTrace[:HalfSWR],0) < np.sum(SpTrace[HalfSWR:],0):
                            pref='After' 
                        else:
                            pref='None'
                        SWR_GlobalResults.loc[counter2, 'SpikeActivityPreference'] = pref
                        SWR_GlobalResults.loc[counter2, 'SpikeActivityBefore'] = np.sum(SpTrace[:HalfSWR],0)
                        SWR_GlobalResults.loc[counter2, 'SpikeActivityAfter'] = np.sum(SpTrace[HalfSWR:],0)
                        counter2+=1  

                # All Ca traces for each spindles per Unique unit (according to cross-registration)

                list_ActivityCa= ['ActivityCa_swr', 'ActivityCa_swr_Precoupled', 'ActivityCa_swr_Postcoupled', 'ActivityCa_swr_Uncoupled']
                list_dict_All_ActivityCa= ['dict_All_ActivityCa_swr', 'dict_All_ActivityCa_swr_Precoupled', 'dict_All_ActivityCa_swr_Postcoupled', 'dict_All_ActivityCa_swr_Uncoupled']
                for it, ActivityCaNames in enumerate(list_ActivityCa): 
                    if len(indexMapp) > 0: #not empty --> cause some units are not in the cross registration..! Need to know why 

                        ActivityCa = locals()[ActivityCaNames]
                        dict_All_ActivityCa = locals()[list_dict_All_ActivityCa[it]]                
                        if len(ActivityCa)>0 :  
                            if np.shape(np.array(ActivityCa))[1] == int(norm_freq*durationSWR*2):   #normalize traces to the same frequency rate    
                                ActivityCa= np.reshape(np.array(ActivityCa), (-1, len(np.array(ActivityCa)))) if np.ndim(ActivityCa) == 1 else np.array(ActivityCa)    
                                dict_All_ActivityCa[str(indexMapp)] = np.append(dict_All_ActivityCa[str(indexMapp)], np.array(ActivityCa), axis=0) if str(indexMapp) in dict_All_ActivityCa else np.array(ActivityCa)
                            else:
                                dataO = np.array(ActivityCa)
                                data= np.repeat(dataO, 2, axis=0) if dataO.shape[0] == 1 else dataO
                                x_mesh, y_mesh = np.meshgrid(np.arange(data.shape[1]), np.arange(data.shape[0]))
                                x_new_mesh, y_new_mesh = np.meshgrid(np.linspace(0, data.shape[1] - 1, int(norm_freq*durationSWR*2)), np.linspace(0, data.shape[0] - 1, np.shape(data)[0]))
                                resampled_dataO = griddata((x_mesh.flatten(), y_mesh.flatten()), data.flatten(), (x_new_mesh, y_new_mesh), method='linear')
                                resampled_data= resampled_dataO[0,:] if dataO.shape[0] == 1 else resampled_dataO
                                resampled_data= np.reshape(resampled_data, (-1, len(resampled_data))) if np.ndim(resampled_data) == 1 else resampled_data
                                dict_All_ActivityCa[str(indexMapp)] = np.append(dict_All_ActivityCa[str(indexMapp)], np.array(resampled_data), axis=0) if str(indexMapp) in dict_All_ActivityCa else np.array(resampled_data)
                
                # All Sp traces for each spindles per Unique unit (according to cross-registration)

                list_ActivitySp= ['ActivitySp_swr', 'ActivitySp_swr_Precoupled', 'ActivitySp_swr_Postcoupled', 'ActivitySp_swr_Uncoupled']
                list_dict_All_ActivitySp= ['dict_All_ActivitySp_swr', 'dict_All_ActivitySp_swr_Precoupled', 'dict_All_ActivitySp_swr_Postcoupled', 'dict_All_ActivitySp_swr_Uncoupled']
                for it, ActivitySpNames in enumerate(list_ActivitySp): 
                    if len(indexMapp) > 0: #not empty --> cause some units are not in the cross registration..! Need to know why 

                        ActivitySp = locals()[ActivitySpNames]
                        dict_All_ActivitySp = locals()[list_dict_All_ActivitySp[it]]                
                        if len(ActivitySp)>0 :  
                            if np.shape(np.array(ActivitySp))[1] == int(norm_freq*durationSWR*2):   
                                ActivitySp= np.reshape(np.array(ActivitySp), (-1, len(np.array(ActivitySp)))) if np.ndim(ActivitySp) == 1 else np.array(ActivitySp)    
                                for cl in range(len(ActivitySp)):
                                    resampled_data_col=ActivitySp[cl]
                                    peaks, _ = find_peaks(resampled_data_col, height=np.std(resampled_data_col))
                                    SpTrace=np.zeros(len(resampled_data_col))
                                    SpTrace[peaks]=1
                                    ActivitySp[cl]=SpTrace
                                dict_All_ActivitySp[str(indexMapp)] = np.append(dict_All_ActivitySp[str(indexMapp)], np.array(ActivitySp), axis=0) if str(indexMapp) in dict_All_ActivitySp else np.array(ActivitySp)
                            else: #normalize traces to the same frequency rate    
                                dataO = np.array(ActivitySp)
                                data= np.repeat(dataO, 2, axis=0) if dataO.shape[0] == 1 else dataO
                                x_mesh, y_mesh = np.meshgrid(np.arange(data.shape[1]), np.arange(data.shape[0]))
                                x_new_mesh, y_new_mesh = np.meshgrid(np.linspace(0, data.shape[1] - 1, int(norm_freq*durationSWR*2)), np.linspace(0, data.shape[0] - 1, np.shape(data)[0]))
                                resampled_dataO = griddata((x_mesh.flatten(), y_mesh.flatten()), data.flatten(), (x_new_mesh, y_new_mesh), method='linear')
                                resampled_data= resampled_dataO[0,:] if dataO.shape[0] == 1 else resampled_dataO
                                resampled_data= np.reshape(resampled_data, (-1, len(resampled_data))) if np.ndim(resampled_data) == 1 else resampled_data
                                for cl in range(len(resampled_data)):
                                    resampled_data_col=resampled_data[cl]
                                    peaks, _ = find_peaks(resampled_data_col, height=np.std(resampled_data_col))
                                    SpTrace=np.zeros(len(resampled_data_col))
                                    SpTrace[peaks]=1
                                    resampled_data[cl]=SpTrace
                                dict_All_ActivitySp[str(indexMapp)] = np.append(dict_All_ActivitySp[str(indexMapp)], np.array(resampled_data), axis=0) if str(indexMapp) in dict_All_ActivitySp else np.array(resampled_data)
            sentence2=f"... in {Cortex}: {nb_spindle} spindles ({cPreCoupled} Pre, {cPostCoupled} Post & {cUnCoupled} Uncoupled Spdl) and {nb_swr} SWR detected ({cPreCoupledSWR} Pre, {cPostCoupledSWR} Post & {cUnCoupledSWR} Uncoupled SWR)"
            print(sentence2)
        sentence3=f"Nb of unique units for {os.path.basename(folder_base)} = {len(dict_All_ActivityCa_spin)}"
        print(sentence3)

        # Save Spindles_GlobalResults

        mice=os.path.basename(folder_base) 
        filenameOut = folder_base / f'Spindles_{Cortex}_ABdetection_GlobalResultsAB_{mice}.xlsx'
        writer = pd.ExcelWriter(filenameOut)
        Spindles_GlobalResults.to_excel(writer)
        writer.close()
        
        # Do average Calcium results for Spindles 

        AVG_dict_All_ActivityCa_spin = {key: np.mean(matrix,0) for key, matrix in dict_All_ActivityCa_spin.items()}
        Array=list(AVG_dict_All_ActivityCa_spin.values())

        AVG_dict_All_ActivityCa_spin_Uncoupled = {key: np.mean(matrix,0) for key, matrix in dict_All_ActivityCa_spin_Uncoupled.items()}
        ArrayUn=list(AVG_dict_All_ActivityCa_spin_Uncoupled.values())

        AVG_dict_All_ActivityCa_spin_Precoupled = {key: np.mean(matrix,0) for key, matrix in dict_All_ActivityCa_spin_Precoupled.items()}
        ArrayPre=list(AVG_dict_All_ActivityCa_spin_Precoupled.values())

        AVG_dict_All_ActivityCa_spin_Postcoupled = {key: np.mean(matrix,0) for key, matrix in dict_All_ActivityCa_spin_Postcoupled.items()}
        ArrayPost=list(AVG_dict_All_ActivityCa_spin_Postcoupled.values())

        filenameOut = folder_base / f'Spindles_{Cortex}_ABdetection_CalciumAvgResultsAB_{mice}.xlsx'
        excel_writer = pd.ExcelWriter(filenameOut)

        Array=pd.DataFrame(Array)
        ArrayUn=pd.DataFrame(ArrayUn)
        ArrayPre=pd.DataFrame(ArrayPre)
        ArrayPost=pd.DataFrame(ArrayPost)

        Array.to_excel(excel_writer, sheet_name='All_Spindles', index=True, header=False)
        ArrayUn.to_excel(excel_writer, sheet_name='Uncoupled_Spindles', index=True, header=False)
        ArrayPre.to_excel(excel_writer, sheet_name='Precoupled_Spindles', index=True, header=False)
        ArrayPost.to_excel(excel_writer, sheet_name='Postcoupled_Spindles', index=True, header=False)

        excel_writer.close()

        # Do average Spike results for Spindles 

        AVG_dict_All_ActivitySp_spin = {key: np.mean(matrix,0) for key, matrix in dict_All_ActivitySp_spin.items()}
        Array=list(AVG_dict_All_ActivitySp_spin.values())

        AVG_dict_All_ActivitySp_spin_Uncoupled = {key: np.mean(matrix,0) for key, matrix in dict_All_ActivitySp_spin_Uncoupled.items()}
        ArrayUn=list(AVG_dict_All_ActivityCa_spin_Uncoupled.values())

        AVG_dict_All_ActivitySp_spin_Precoupled = {key: np.mean(matrix,0) for key, matrix in dict_All_ActivitySp_spin_Precoupled.items()}
        ArrayPre=list(AVG_dict_All_ActivitySp_spin_Precoupled.values())

        AVG_dict_All_ActivitySp_spin_Postcoupled = {key: np.mean(matrix,0) for key, matrix in dict_All_ActivitySp_spin_Postcoupled.items()}
        ArrayPost=list(AVG_dict_All_ActivitySp_spin_Postcoupled.values())

        filenameOut = folder_base / f'Spindles_{Cortex}_ABdetection_SpikeAvgResultsAB_{mice}.xlsx'
        excel_writer = pd.ExcelWriter(filenameOut)

        Array=pd.DataFrame(Array)
        ArrayUn=pd.DataFrame(ArrayUn)
        ArrayPre=pd.DataFrame(ArrayPre)
        ArrayPost=pd.DataFrame(ArrayPost)

        Array.to_excel(excel_writer, sheet_name='All_Spindles', index=True, header=False)
        ArrayUn.to_excel(excel_writer, sheet_name='Uncoupled_Spindles', index=True, header=False)
        ArrayPre.to_excel(excel_writer, sheet_name='Precoupled_Spindles', index=True, header=False)
        ArrayPost.to_excel(excel_writer, sheet_name='Postcoupled_Spindles', index=True, header=False)

        excel_writer.close()

               
        # Save SWR_GlobalResults

        mice=os.path.basename(folder_base) 
        filenameOut = folder_base / f'SWR_{Cortex}_ABdetection_GlobalResultsAB_{mice}.xlsx'
        writer = pd.ExcelWriter(filenameOut)
        SWR_GlobalResults.to_excel(writer)
        writer.close()

        # Do average Calcium results for SWR

        AVG_dict_All_ActivityCa_swr = {key: np.mean(matrix,0) for key, matrix in dict_All_ActivityCa_swr.items()}
        Array=list(AVG_dict_All_ActivityCa_swr.values())

        AVG_dict_All_ActivityCa_swr_Uncoupled = {key: np.mean(matrix,0) for key, matrix in dict_All_ActivityCa_swr_Uncoupled.items()}
        ArrayUn=list(AVG_dict_All_ActivityCa_swr_Uncoupled.values())

        AVG_dict_All_ActivityCa_swr_Precoupled = {key: np.mean(matrix,0) for key, matrix in dict_All_ActivityCa_swr_Precoupled.items()}
        ArrayPre=list(AVG_dict_All_ActivityCa_swr_Precoupled.values())

        AVG_dict_All_ActivityCa_swr_Postcoupled = {key: np.mean(matrix,0) for key, matrix in dict_All_ActivityCa_swr_Postcoupled.items()}
        ArrayPost=list(AVG_dict_All_ActivityCa_swr_Postcoupled.values())

        filenameOut = folder_base / f'SWR_{Cortex}_ABdetection_CalciumAvgResultsAB_{mice}.xlsx'
        excel_writer = pd.ExcelWriter(filenameOut)

        Array=pd.DataFrame(Array)
        ArrayUn=pd.DataFrame(ArrayUn)
        ArrayPre=pd.DataFrame(ArrayPre)
        ArrayPost=pd.DataFrame(ArrayPost)

        Array.to_excel(excel_writer, sheet_name='All_SWR', index=True, header=False)
        ArrayUn.to_excel(excel_writer, sheet_name='Uncoupled_SWR', index=True, header=False)
        ArrayPre.to_excel(excel_writer, sheet_name='Precoupled_SWR', index=True, header=False)
        ArrayPost.to_excel(excel_writer, sheet_name='Postcoupled_SWR', index=True, header=False)

        excel_writer.close()

        # Do average Spike results for SWR

        AVG_dict_All_ActivitySp_swr = {key: np.mean(matrix,0) for key, matrix in dict_All_ActivitySp_swr.items()}
        Array=list(AVG_dict_All_ActivitySp_swr.values())

        AVG_dict_All_ActivitySp_swr_Uncoupled = {key: np.mean(matrix,0) for key, matrix in dict_All_ActivitySp_swr_Uncoupled.items()}
        ArrayUn=list(AVG_dict_All_ActivitySp_swr_Uncoupled.values())

        AVG_dict_All_ActivitySp_swr_Precoupled = {key: np.mean(matrix,0) for key, matrix in dict_All_ActivitySp_swr_Precoupled.items()}
        ArrayPre=list(AVG_dict_All_ActivitySp_swr_Precoupled.values())

        AVG_dict_All_ActivitySp_swr_Postcoupled = {key: np.mean(matrix,0) for key, matrix in dict_All_ActivitySp_swr_Postcoupled.items()}
        ArrayPost=list(AVG_dict_All_ActivitySp_swr_Postcoupled.values())

        filenameOut = folder_base / f'SWR_{Cortex}_ABdetection_SpikeAvgResultsAB_{mice}.xlsx'
        excel_writer = pd.ExcelWriter(filenameOut)

        Array=pd.DataFrame(Array)
        ArrayUn=pd.DataFrame(ArrayUn)
        ArrayPre=pd.DataFrame(ArrayPre)
        ArrayPost=pd.DataFrame(ArrayPost)

        Array.to_excel(excel_writer, sheet_name='All_SWR', index=True, header=False)
        ArrayUn.to_excel(excel_writer, sheet_name='Uncoupled_SWR', index=True, header=False)
        ArrayPre.to_excel(excel_writer, sheet_name='Precoupled_SWR', index=True, header=False)
        ArrayPost.to_excel(excel_writer, sheet_name='Postcoupled_SWR', index=True, header=False)

        excel_writer.close() 


