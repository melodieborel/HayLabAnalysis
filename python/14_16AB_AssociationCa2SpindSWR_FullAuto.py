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
    dict_Spindleprop = {}
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
                dict_Spindleprop[subsession]  = pd.read_csv(Spindleproperties_PFC)
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
            dict_Spindleprop[session]  = pd.read_csv(Spindleproperties_PFC)
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
    for c in range(len(B)):
        print('unit n°', c)
        for sess in list(dict_Stamps.keys()):
            print('= unit', int(B[sess][c]), 'in', sess) if math.isnan (float(B[sess][c])) == False else None

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

    data = {}
    Struct = "PFC"
    before = 1000 # Max distance in ms between a SWR and a spindle to be considered as Precoupled
    after = 1000 # Max distance in ms between a spindle and a SWR to be considered as Postcoupled
    duration = 0.5 # number of sec before and after the Spdl onset taken into acount
    counter=0
    counter2=0

    norm_freq=10

    Spindles_GlobalResults= pd.DataFrame(data, columns=['Mice', 'Session','Session_Time','Unique_Unit','UnitNumber','UnitValue','SpdlStatut','SpdlNumber','SpdlDuration (ms)','SWR inside Spdl','CalciumActivityPreference', 'CalciumActivityBefore','CalciumActivityAfter','AUC_calciumBefore','AUC_calciumAfter','SpikeActivityPreference','SpikeActivityBefore','SpikeActivityAfter','AUC_spikeBefore', 'AUC_spikeAfter'])
    SWR_GlobalResults= pd.DataFrame(data, columns=['Mice', 'Session','Session_Time','Unique_Unit','UnitNumber','UnitValue','SWRStatut','SWRNumber','SWRDuration (ms)','SWR inside Spdl','CalciumActivityPreference', 'CalciumActivityBefore','CalciumActivityAfter','AUC_calciumBefore','AUC_calciumAfter','SpikeActivityPreference','SpikeActivityBefore','SpikeActivityAfter','AUC_spikeBefore', 'AUC_spikeAfter'])

    dict_All_ActivityCa_PFCspin={}
    dict_All_ActivityCa_PFCspin_Precoupled={}
    dict_All_ActivityCa_PFCspin_Postcoupled={}
    dict_All_ActivityCa_PFCspin_Uncoupled={}

    dict_All_ActivityCa_swr={}
    dict_All_ActivityCa_swr_Precoupled={}
    dict_All_ActivityCa_swr_Postcoupled={}
    dict_All_ActivityCa_swr_Uncoupled={}

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
        Supd = C.loc[:, :] #S.loc[:, :]
        
        if len(list_droppedframes) > 0:
            numbdropfr = len(list(item for item in range(numbdropfr) if list_droppedframes[item] < rec_dur))
        else:
            numbdropfr = 0
        EndTime = StartTime + ((rec_dur + numbdropfr)/minian_freq) # in seconds

        PFCspipropO=dict_Spindleprop[i]
        PFCspipropM=PFCspipropO.copy()
        SWRpropO=dict_SWRprop[i]
        SWRpropM=SWRpropO.copy()
        PFCspipropM[["peak time", "start time", "end time"]] = PFCspipropM[["peak time", "start time", "end time"]]-(StartTime*1000)
        SWRpropM[["peak time", "start time", "end time"]] = SWRpropM[["peak time", "start time", "end time"]]-(StartTime*1000)

        AA = C['unit_id']
        copyAA = list(AA.copy())
        unit_to_drop=dict_TodropFile[i]    
        for u in unit_to_drop:
            copyAA.remove(u)
        unit_to_keep = copyAA
        Cupd = Cupd.loc[unit_to_keep,:]
        Supd = Supd.loc[unit_to_keep,:]
        nb_unit = Cupd.shape[0]

        Cseries = Cupd.to_series()
        C_upd_unit_id = Cupd['unit_id'].values
        Sseries = Supd.to_series()

        sentence1= f"In {i} kept values = {C_upd_unit_id}"
        print(sentence1)

        time = range(int(duration*2*minian_freq))
        Half = int(duration*minian_freq)
        TimeStamps_miniscope=list(dict_StampsMiniscope[i]["Time Stamp (ms)"]) # + (StartTime*1000))


        PFCspipropTrunc = PFCspipropM[PFCspipropM["start time"]>0]
        PFCspipropTrunc = PFCspipropTrunc[PFCspipropTrunc["start time"]< TimeStamps_miniscope[-1]]
        SWRpropTrunc = SWRpropM[SWRpropM["start time"]>0]
        SWRpropTrunc = SWRpropTrunc[SWRpropTrunc["start time"] < TimeStamps_miniscope[-1]]
    
        units = range(nb_unit)
        nb_spindle = PFCspipropTrunc.shape[0]
        nb_swr = SWRpropTrunc.shape[0]
        for ii, unit in enumerate(units): # for each kept units (cause Cseries/Sseries only have kept units)
            lCseries = np.array(Cseries)[(unit)*rec_dur:(unit+1)*rec_dur]
            lSseries = np.array(Sseries)[(unit)*rec_dur:(unit+1)*rec_dur]
            ActivityCa_PFCspin = [] #For each unit  
            ActivityCa_PFCspin_Precoupled= [] #For each unit 
            ActivityCa_PFCspin_Postcoupled= [] #For each unit 
            ActivityCa_PFCspin_Uncoupled= [] #For each unit 
            startSpiList = list(pd.Series(PFCspipropTrunc["start time"]))
            endSpiList = list(pd.Series(PFCspipropTrunc["end time"]))
            for Pspin in range(nb_spindle): 

                # Get the calcium and spike trace associated with the spdl
                startSpi=startSpiList[Pspin]
                endSpi=endSpiList[Pspin]

                Frame_Spindle_start = take_closest2(TimeStamps_miniscope, startSpi)
                index = TimeStamps_miniscope.index(Frame_Spindle_start)
                CaTrace = list(lCseries[index-Half:index+Half])
                SpTrace = list(lSseries[index-Half:index+Half])            
                if len(CaTrace)<len(time): 
                    print("Spindle too close to the begining/end of the recording,", i, ", Spdl n°", Pspin, ", Start Spdl =", int(startSpi), "ms, Start Rec =", int(Frame_Spindle_start), 'ms')
                else:
                    ActivityCa_PFCspin.append(CaTrace)

                    # Define if that spindle is coupled with a SWR or not
                    Spdl_statut=[]
                    startSWRList = list(pd.Series(SWRpropTrunc["start time"]))
                    startClosest_SWR = take_closest2(startSWRList, startSpi)# + StartTimeIndexSpi])
                    distance = startClosest_SWR - startSpi#  + StartTimeIndexSpi]
                    if (distance > (- before)) and (distance <  0):
                        Spdl_statut = ['PreCoupled']
                        cPreCoupled+=1 if ii==1 else 0
                        ActivityCa_PFCspin_Precoupled.append(CaTrace)
                    elif (distance > (0)) and (distance <  after):
                        Spdl_statut = ['PostCoupled']
                        cPostCoupled+=1 if ii==1 else 0
                        ActivityCa_PFCspin_Postcoupled.append(CaTrace)
                    else:
                        Spdl_statut= ['UnCoupled']
                        cUnCoupled+=1 if ii==1 else 0
                        ActivityCa_PFCspin_Uncoupled.append(CaTrace)

                    Spindles_GlobalResults.loc[counter, 'Mice'] = os.path.basename(folder_base)
                    Spindles_GlobalResults.loc[counter, 'Session'] = i 
                    Spindles_GlobalResults.loc[counter, 'Session_Time'] = None 
                    indexMapp = np.where(B[i] == C_upd_unit_id[unit])[0]
                    Spindles_GlobalResults.loc[counter, 'Unique_Unit'] = indexMapp 
                    Spindles_GlobalResults.loc[counter, 'UnitNumber'] = unit 
                    Spindles_GlobalResults.loc[counter, 'UnitValue'] = C_upd_unit_id[unit] 
                    Spindles_GlobalResults.loc[counter, 'SpdlStatut'] = Spdl_statut
                    Spindles_GlobalResults.loc[counter, 'SpdlNumber'] = Pspin
                    Spindles_GlobalResults.loc[counter, 'SpdlDuration (ms)'] = endSpi- startSpi
                    IsTrue=is_between(startSWRList,startSpi, endSpi)
                    Spindles_GlobalResults.loc[counter, 'SWR inside Spdl'] = IsTrue
                    
                    if np.mean(CaTrace[:Half],0) > np.mean(CaTrace[Half:],0):
                        pref='Before'
                    elif np.mean(CaTrace[:Half],0) < np.mean(CaTrace[Half:],0):
                        pref='After' 
                    else:
                        pref='None'
                    Spindles_GlobalResults.loc[counter, 'CalciumActivityPreference'] = pref
                    Spindles_GlobalResults.loc[counter, 'CalciumActivityBefore'] = np.mean(CaTrace[:Half],0)
                    Spindles_GlobalResults.loc[counter, 'CalciumActivityAfter'] = np.mean(CaTrace[Half:],0)
                    Spindles_GlobalResults.loc[counter, 'AUC_calciumBefore'] = np.trapz(CaTrace[:Half],np.arange(0,len(CaTrace[:Half]),1))
                    Spindles_GlobalResults.loc[counter, 'AUC_calciumAfter'] = np.trapz(CaTrace[Half:],np.arange(0,len(CaTrace[Half:]),1))          

                    if np.mean(SpTrace[:Half],0) > np.mean(SpTrace[Half:],0):
                        pref='Before'
                    elif np.mean(SpTrace[:Half],0) < np.mean(SpTrace[Half:],0):
                        pref='After' 
                    else:
                        pref='None'
                    Spindles_GlobalResults.loc[counter, 'SpikeActivityPreference'] = pref
                    Spindles_GlobalResults.loc[counter, 'SpikeActivityBefore'] = np.mean(SpTrace[:Half],0)
                    Spindles_GlobalResults.loc[counter, 'SpikeActivityAfter'] = np.mean(SpTrace[Half:],0)
                    Spindles_GlobalResults.loc[counter, 'AUC_spikeBefore'] = np.trapz(SpTrace[:Half],np.arange(0,len(SpTrace[:Half]),1))
                    Spindles_GlobalResults.loc[counter, 'AUC_spikeAfter'] = np.trapz(SpTrace[Half:],np.arange(0,len(SpTrace[Half:]),1))
                    counter+=1  

            list_ActivityCa= ['ActivityCa_PFCspin', 'ActivityCa_PFCspin_Precoupled', 'ActivityCa_PFCspin_Postcoupled', 'ActivityCa_PFCspin_Uncoupled']
            list_dict_All_ActivityCa= ['dict_All_ActivityCa_PFCspin', 'dict_All_ActivityCa_PFCspin_Precoupled', 'dict_All_ActivityCa_PFCspin_Postcoupled', 'dict_All_ActivityCa_PFCspin_Uncoupled']
            for it, ActivityCaNames in enumerate(list_ActivityCa): # for each Spdl types
                if len(indexMapp) > 0: #not empty --> cause some units are not in the cross registration..! Need to know why 
                    # All traces for each spindles per Unique unit (according to cross-registration)
                    ActivityCa = locals()[ActivityCaNames]
                    dict_All_ActivityCa = locals()[list_dict_All_ActivityCa[it]]       
                    #print(indexMapp, '= cell nb in crossReg,', unit, '= cell index,', ActivityCaNames)
                    if len(ActivityCa)>0 :    
                        if np.shape(np.array(ActivityCa))[1] == norm_freq:  #normalize traces to the same frequency rate         
                            ActivityCa= np.reshape(np.array(ActivityCa), (-1, len(np.array(ActivityCa)))) if np.ndim(ActivityCa) == 1 else np.array(ActivityCa)    
                            dict_All_ActivityCa[str(indexMapp)] = np.append(dict_All_ActivityCa[str(indexMapp)], np.array(ActivityCa), axis=0) if str(indexMapp) in dict_All_ActivityCa else np.array(ActivityCa)
                        else:
                            dataO = np.array(ActivityCa)
                            data= np.repeat(dataO, 2, axis=0) if dataO.shape[0] == 1 else dataO
                            x_mesh, y_mesh = np.meshgrid(np.arange(data.shape[1]), np.arange(data.shape[0]))
                            x_new_mesh, y_new_mesh = np.meshgrid(np.linspace(0, data.shape[1] - 1, norm_freq), np.linspace(0, data.shape[0] - 1, np.shape(data)[0]))
                            resampled_dataO = griddata((x_mesh.flatten(), y_mesh.flatten()), data.flatten(), (x_new_mesh, y_new_mesh), method='linear')
                            resampled_data= resampled_dataO[0,:] if dataO.shape[0] == 1 else resampled_dataO
                            resampled_data= np.reshape(resampled_data, (-1, len(resampled_data))) if np.ndim(resampled_data) == 1 else resampled_data
                            dict_All_ActivityCa[str(indexMapp)] = np.append(dict_All_ActivityCa[str(indexMapp)], np.array(resampled_data), axis=0) if str(indexMapp) in dict_All_ActivityCa else np.array(resampled_data)
                    sentence1bis=""
                else: 
                    sentence1bis=f"Cell idx {unit} not in the cross registration!!!"
                    print(sentence1bis)
            
            ActivityCa_swr = [] #For each unit  
            ActivityCa_swr_Precoupled= [] #For each unit 
            ActivityCa_swr_Postcoupled= [] #For each unit 
            ActivityCa_swr_Uncoupled= [] #For each unit 

            startSwrList = list(pd.Series(SWRpropTrunc["start time"]))
            endSwrList = list(pd.Series(SWRpropTrunc["end time"]))
            for Pswr in range(nb_swr): 
                # Get the calcium and spike trace associated with the swr
                startSwr=startSwrList[Pswr]
                endSwr=endSwrList[Pswr]
                Frame_SWR_start = take_closest2(TimeStamps_miniscope, startSwr)
                index = TimeStamps_miniscope.index(Frame_SWR_start)
                CaTrace = list(lCseries[index-Half:index+Half])
                SpTrace = list(lSseries[index-Half:index+Half])     
                if len(CaTrace)<len(time): 
                    print("SWR too close to the begining/end of the recording,", i, ", SWR n°", Pswr, ", Start SWR =",  int(startSwr), "ms, Start Rec =", int(Frame_SWR_start), 'ms')
                else:
                    ActivityCa_swr.append(CaTrace) 

                    # Define if that SWR is coupled with a spdl or not
                    SWR_statut=[]
                    startSpiList = list(pd.Series(PFCspipropTrunc["start time"]))
                    endSpiList = list(pd.Series(PFCspipropTrunc["end time"]))
                    startClosest_Spi = take_closest2(startSpiList, startSwr)# + StartTimeIndexSpi])
                    indexSpi = startSpiList.index(startClosest_Spi)
                    endClosest_Spi=endSpiList[indexSpi]
                    distance = startClosest_Spi - startSwr #  + StartTimeIndexSpi]  
                    IsTrue = 'False'             
                    if (distance > (- before)) and (distance <  0):
                        SWR_statut = ['Postcoupled']
                        cPostCoupledSWR+=1 if ii==1 else 0
                        ActivityCa_swr_Postcoupled.append(CaTrace)
                        if startSwr<endClosest_Spi:
                            IsTrue = 'True' #SWR inside the Spindle
                    elif (distance > (0)) and (distance <  after):
                        SWR_statut = ['Precoupled']
                        cPreCoupledSWR+=1 if ii==1 else 0
                        ActivityCa_swr_Precoupled.append(CaTrace)
                    else:
                        SWR_statut= ['UnCoupled']
                        cUnCoupledSWR+=1 if ii==1 else 0
                        ActivityCa_swr_Uncoupled.append(CaTrace)

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
                    
                    if np.mean(CaTrace[:Half],0) > np.mean(CaTrace[Half:],0):
                        pref='Before'
                    elif np.mean(CaTrace[:Half],0) < np.mean(CaTrace[Half:],0):
                        pref='After' 
                    else:
                        pref='None'
                    SWR_GlobalResults.loc[counter2, 'CalciumActivityPreference'] = pref
                    SWR_GlobalResults.loc[counter2, 'CalciumActivityBefore'] = np.mean(CaTrace[:Half],0)
                    SWR_GlobalResults.loc[counter2, 'CalciumActivityAfter'] = np.mean(CaTrace[Half:],0)
                    SWR_GlobalResults.loc[counter2, 'AUC_calciumBefore'] = np.trapz(CaTrace[:Half],np.arange(0,len(CaTrace[:Half]),1))
                    SWR_GlobalResults.loc[counter2, 'AUC_calciumAfter'] = np.trapz(CaTrace[Half:],np.arange(0,len(CaTrace[Half:]),1))          

                    if np.mean(SpTrace[:Half],0) > np.mean(SpTrace[Half:],0):
                        pref='Before'
                    elif np.mean(SpTrace[:Half],0) < np.mean(SpTrace[Half:],0):
                        pref='After' 
                    else:
                        pref='None'
                    SWR_GlobalResults.loc[counter2, 'SpikeActivityPreference'] = pref
                    SWR_GlobalResults.loc[counter2, 'SpikeActivityBefore'] = np.mean(SpTrace[:Half],0)
                    SWR_GlobalResults.loc[counter2, 'SpikeActivityAfter'] = np.mean(SpTrace[Half:],0)
                    SWR_GlobalResults.loc[counter2, 'AUC_spikeBefore'] = np.trapz(SpTrace[:Half],np.arange(0,len(SpTrace[:Half]),1))
                    SWR_GlobalResults.loc[counter2, 'AUC_spikeAfter'] = np.trapz(SpTrace[Half:],np.arange(0,len(SpTrace[Half:]),1))
                    counter2+=1  

            list_ActivityCa= ['ActivityCa_swr', 'ActivityCa_swr_Precoupled', 'ActivityCa_swr_Postcoupled', 'ActivityCa_swr_Uncoupled']
            list_dict_All_ActivityCa= ['dict_All_ActivityCa_swr', 'dict_All_ActivityCa_swr_Precoupled', 'dict_All_ActivityCa_swr_Postcoupled', 'dict_All_ActivityCa_swr_Uncoupled']
            for it, ActivityCaNames in enumerate(list_ActivityCa): 
                if len(indexMapp) > 0: #not empty --> cause some units are not in the cross registration..! Need to know why 
                    # All traces for each spindles per Unique unit (according to cross-registration)
                    ActivityCa = locals()[ActivityCaNames]
                    dict_All_ActivityCa = locals()[list_dict_All_ActivityCa[it]]                
                    if len(ActivityCa)>0 :  
                        if np.shape(np.array(ActivityCa))[1] == norm_freq:   #normalize traces to the same frequency rate    
                            ActivityCa= np.reshape(np.array(ActivityCa), (-1, len(np.array(ActivityCa)))) if np.ndim(ActivityCa) == 1 else np.array(ActivityCa)    
                            dict_All_ActivityCa[str(indexMapp)] = np.append(dict_All_ActivityCa[str(indexMapp)], np.array(ActivityCa), axis=0) if str(indexMapp) in dict_All_ActivityCa else np.array(ActivityCa)
                        else:
                            dataO = np.array(ActivityCa)
                            data= np.repeat(dataO, 2, axis=0) if dataO.shape[0] == 1 else dataO
                            x_mesh, y_mesh = np.meshgrid(np.arange(data.shape[1]), np.arange(data.shape[0]))
                            x_new_mesh, y_new_mesh = np.meshgrid(np.linspace(0, data.shape[1] - 1, norm_freq), np.linspace(0, data.shape[0] - 1, np.shape(data)[0]))
                            resampled_dataO = griddata((x_mesh.flatten(), y_mesh.flatten()), data.flatten(), (x_new_mesh, y_new_mesh), method='linear')
                            resampled_data= resampled_dataO[0,:] if dataO.shape[0] == 1 else resampled_dataO
                            resampled_data= np.reshape(resampled_data, (-1, len(resampled_data))) if np.ndim(resampled_data) == 1 else resampled_data
                            dict_All_ActivityCa[str(indexMapp)] = np.append(dict_All_ActivityCa[str(indexMapp)], np.array(resampled_data), axis=0) if str(indexMapp) in dict_All_ActivityCa else np.array(resampled_data)
        
        sentence2=f"... with {nb_spindle} spindles ({cPreCoupled} Pre, {cPostCoupled} Post & {cUnCoupled} Uncoupled Spdl) and {nb_swr} SWR detected ({cPreCoupledSWR} Pre, {cPostCoupledSWR} Post & {cUnCoupledSWR} Uncoupled SWR)"
        print(sentence2)
    sentence3=f"Nb of unique units for {os.path.basename(folder_base)} = {len(dict_All_ActivityCa_PFCspin)}"
    print(sentence3)

    # Save to Excel

    mice=os.path.basename(folder_base) 
    filenameOut = folder_base / f'Spindles_ABdetection_GlobalResultsAB_{mice}.xlsx'
    writer = pd.ExcelWriter(filenameOut)
    Spindles_GlobalResults.to_excel(writer)
    writer.close()

    # Do average results for this mouse

    AVG_dict_All_ActivityCa_PFCspin = {key: np.mean(matrix,0) for key, matrix in dict_All_ActivityCa_PFCspin.items()}
    Array=list(AVG_dict_All_ActivityCa_PFCspin.values())

    AVG_dict_All_ActivityCa_PFCspin_Uncoupled = {key: np.mean(matrix,0) for key, matrix in dict_All_ActivityCa_PFCspin_Uncoupled.items()}
    ArrayUn=list(AVG_dict_All_ActivityCa_PFCspin_Uncoupled.values())

    AVG_dict_All_ActivityCa_PFCspin_Precoupled = {key: np.mean(matrix,0) for key, matrix in dict_All_ActivityCa_PFCspin_Precoupled.items()}
    ArrayPre=list(AVG_dict_All_ActivityCa_PFCspin_Precoupled.values())

    AVG_dict_All_ActivityCa_PFCspin_Postcoupled = {key: np.mean(matrix,0) for key, matrix in dict_All_ActivityCa_PFCspin_Postcoupled.items()}
    ArrayPost=list(AVG_dict_All_ActivityCa_PFCspin_Postcoupled.values())

    filenameOut = folder_base / f'Spindles_ABdetection_AverageResultsAB_{mice}.xlsx'
    excel_writer = pd.ExcelWriter(filenameOut)

    Array=pd.DataFrame(Array)
    ArrayUn=pd.DataFrame(ArrayUn)
    ArrayPre=pd.DataFrame(ArrayPre)
    ArrayPost=pd.DataFrame(ArrayPost)

    Array.to_excel(excel_writer, sheet_name='All_Spindles', index=True, header=False)
    ArrayUn.to_excel(excel_writer, sheet_name='Uncoupled_Spindles', index=True, header=False)
    ArrayPre.to_excel(excel_writer, sheet_name='Precoupled_Spindles', index=True, header=False)
    ArrayPost.to_excel(excel_writer, sheet_name='Postcoupled_Spindles', index=True, header=False)

    # Save the Excel file
    excel_writer.close()

    mice=os.path.basename(folder_base) 
    filenameOut = folder_base / f'SWR_ABdetection_GlobalResultsAB_{mice}.xlsx'
    writer = pd.ExcelWriter(filenameOut)
    Spindles_GlobalResults.to_excel(writer)
    writer.close()


    AVG_dict_All_ActivityCa_swr = {key: np.mean(matrix,0) for key, matrix in dict_All_ActivityCa_swr.items()}
    Array=list(AVG_dict_All_ActivityCa_swr.values())

    AVG_dict_All_ActivityCa_swr_Uncoupled = {key: np.mean(matrix,0) for key, matrix in dict_All_ActivityCa_swr_Uncoupled.items()}
    ArrayUn=list(AVG_dict_All_ActivityCa_swr_Uncoupled.values())

    AVG_dict_All_ActivityCa_swr_Precoupled = {key: np.mean(matrix,0) for key, matrix in dict_All_ActivityCa_swr_Precoupled.items()}
    ArrayPre=list(AVG_dict_All_ActivityCa_swr_Precoupled.values())

    AVG_dict_All_ActivityCa_swr_Postcoupled = {key: np.mean(matrix,0) for key, matrix in dict_All_ActivityCa_swr_Postcoupled.items()}
    ArrayPost=list(AVG_dict_All_ActivityCa_swr_Postcoupled.values())



    filenameOut = folder_base / f'SWR_ABdetection_AverageResultsAB_{mice}.xlsx'
    excel_writer = pd.ExcelWriter(filenameOut)

    Array=pd.DataFrame(Array)
    ArrayUn=pd.DataFrame(ArrayUn)
    ArrayPre=pd.DataFrame(ArrayPre)
    ArrayPost=pd.DataFrame(ArrayPost)

    Array.to_excel(excel_writer, sheet_name='All_SWR', index=True, header=False)
    ArrayUn.to_excel(excel_writer, sheet_name='Uncoupled_SWR', index=True, header=False)
    ArrayPre.to_excel(excel_writer, sheet_name='Precoupled_SWR', index=True, header=False)
    ArrayPost.to_excel(excel_writer, sheet_name='Postcoupled_SWR', index=True, header=False)

    # Save the Excel file
    excel_writer.close()


