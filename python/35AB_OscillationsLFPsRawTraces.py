# # Associate Ca2+ signal with spindles for each session & subsessions using crossregistration

#######################################################################################
                                # Load packages #
#######################################################################################

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
import shutil
from bisect import bisect_left
from ast import literal_eval

from itertools import groupby
from ephyviewer import mkQApp, MainViewer, TraceViewer
from IPython.display import display
from ipyfilechooser import FileChooser
from datetime import datetime

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

def take_closest(myList, myNumber):        
    #Assumes myList is sorted. Returns closest value to myNumber.
    #If two numbers are equally close, return the smallest number.        
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

def is_overlapping(starttime, endtime, starttimeList, endtimeList):
    IsTrue='False'
    for ind in range(len(starttimeList)):
        if starttime<=starttimeList[ind] and starttimeList[ind]<=endtime: # event n°2 begins after the start n°1               
            if (endtime-starttimeList[ind])>=int(0.5*(endtime-starttime)): # overlapp > to 50% of the duration of the event n°1
                IsTrue='True'
        elif starttime<=endtimeList[ind] and endtimeList[ind]<=endtime: # event n°2 ends before the end n°1 
            if (endtimeList[ind]-starttime)>=int(0.5*(endtime-starttime)): # overlapp > to 50% of the duration of the event n°1
                IsTrue='True'
    return IsTrue

def Convert(string):
    li = list(string.split(", "))
    li2 = len(li)
    return li2
    
#######################################################################################
                # Load sleep score and Ca2+ time series numpy arrays #
#######################################################################################

MiceList=['BlackLinesOK', 'BlueLinesOK', 'GreenDotsOK', 'GreenLinesOK', 'Purple', 'RedLinesOK','ThreeColDotsOK', 'ThreeBlueCrossesOK']

# Get the current date and time
FolderNameSave=str(datetime.now())
FolderNameSave = FolderNameSave.replace(" ", "_").replace(".", "_").replace(":", "_").replace("-", "_")
destination_folder= f"//10.69.168.1/crnldata/waking/audrey_hay/L1imaging/AnalysedMarch2023/Gaelle/Baseline_recording_ABmodified/AB_Analysis/RawLFPOscillations_PerMouse_{FolderNameSave}"
os.makedirs(destination_folder)
folder_to_save=Path(destination_folder)

# Copy the script file to the destination folder
source_script = "C:/Users/Manip2/SCRIPTS/Code python audrey/code python aurelie/HayLabAnalysis/python/35AB_OscillationsLFPsRawTraces.py"
destination_file_path = f"{destination_folder}/35AB_OscillationsLFPsRawTraces.txt"
shutil.copy(source_script, destination_file_path)

for micename in MiceList:

    dpath0 = "//10.69.168.1/crnldata/waking/audrey_hay/L1imaging/AnalysedMarch2023/Gaelle/Baseline_recording_ABmodified/"
    dpath=dpath0 + micename
    print(f"####################################################################################")
    print(f"################################### {micename} ####################################")
    print(f"####################################################################################")
    print(f"Path to the folder : {dpath}")
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

    subsessions = []
    dict_Calcium = {}
    dict_Spike = {}
    dict_SWRprop = {}
    dict_Spindleprop_PFC = {}
    dict_Spindleprop_S1 = {}
    dict_Stamps = {}
    dict_StampsMiniscope = {}
    dict_TodropFile = {}
    dict_DSprop_S1={}
    dict_DSprop_PFC={}

    dict_LFP_S1={}
    dict_LFP_PFC={}
    dict_LFP_CA1={}

    sessions = [folder.name for folder in folder_base.iterdir() if folder.is_dir() and "session" in folder.name]
    print(sessions)

    for session in sessions:
        folder_mini = folder_base / session / f'V4_Miniscope'
        nb_subsessions = sum(1 for p in folder_mini.iterdir() if p.is_dir() and p.name.startswith("session"))
        SWRproperties = folder_base /session / f'OpenEphys/SWRproperties_8sd_AB.xlsx'
        Spindleproperties_PFC = folder_base / session / f'OpenEphys/Spindlesproperties_PFC_7sd_AB.xlsx'
        Spindleproperties_S1 = folder_base / session / f'OpenEphys/Spindlesproperties_S1_7sd_AB.xlsx'
        DownStatesproperties_S1 = folder_base / session / f'OpenEphys/DownStatesproperties_S1_2Pro3Height_AB.xlsx'
        DownStatesproperties_PFC = folder_base / session / f'OpenEphys/DownStatesproperties_PFC_2Pro3Height_AB.xlsx'
        StampsFile = folder_base / session / f'SynchroFile.xlsx'
        StampsMiniscopeFile = folder_mini / f'timeStamps.csv'
        
        LFPFile = folder_base / session / f'OpenEphys/RawDataChannelExtractedDS.npy'
        All = np.load(LFPFile, mmap_mode= 'r')
        
        Channels = '//10.69.168.1/crnldata/waking/audrey_hay/L1imaging/AnalysedMarch2023/LFPChannels_perMice.xlsx' 

        allchannels = pd.read_excel(Channels)
        
        CA1ch1=int(allchannels[micename][2].split(',')[0])
        CA1ch2=int(allchannels[micename][2].split(',')[1])
        CA1  =  All[:, CA1ch1]-All[:, CA1ch2] 

        PFCch1=int(allchannels[micename][0].split(',')[0])
        PFCch2=int(allchannels[micename][0].split(',')[1])
        PFC  =  All[:, PFCch1]-All[:, PFCch2] 

        S1ch1=int(allchannels[micename][1].split(',')[0])
        S1ch2=int(allchannels[micename][1].split(',')[1])
        S1  =  All[:, S1ch1]-All[:, S1ch2] 


        if nb_subsessions!=0:
            for x in range(1, nb_subsessions+1):            
                subsession= session + str(x)
                subsessions.append(subsession)    
                minian_ds = open_minian(folder_mini / subsession / f'minian')      # OR minianAB
                dict_SWRprop[subsession]  = pd.read_excel(SWRproperties)
                dict_Spindleprop_PFC[subsession]  = pd.read_excel(Spindleproperties_PFC)
                dict_Spindleprop_S1[subsession]  = pd.read_excel(Spindleproperties_S1)            
                dict_DSprop_S1[subsession]  = pd.read_excel(DownStatesproperties_S1)            
                dict_DSprop_PFC[subsession]  = pd.read_excel(DownStatesproperties_PFC)            
                dict_Calcium[subsession] = minian_ds['C'] # calcium traces 
                dict_Spike[subsession] = minian_ds['S'] # estimated spikes
                dict_Stamps[subsession]  = pd.read_excel(StampsFile)
                dict_StampsMiniscope[subsession]  = pd.read_csv(StampsMiniscopeFile)
                dict_LFP_S1[subsession]  = S1
                dict_LFP_PFC[subsession]  = PFC
                dict_LFP_CA1[subsession]  = CA1
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
        else:
            minian_ds = open_minian(folder_mini / f'minian')            # OR minianAB
            dict_Calcium[session] = minian_ds['C'] # calcium traces 
            dict_Spike[session] = minian_ds['S'] # estimated spikes
            dict_SWRprop[session]  = pd.read_excel(SWRproperties)
            dict_Spindleprop_PFC[session]  = pd.read_excel(Spindleproperties_PFC)
            dict_Spindleprop_S1[session]  = pd.read_excel(Spindleproperties_S1)
            dict_DSprop_S1[session]  = pd.read_excel(DownStatesproperties_S1)            
            dict_DSprop_PFC[session]  = pd.read_excel(DownStatesproperties_PFC)            
            dict_Stamps[session]  = pd.read_excel(StampsFile)
            dict_StampsMiniscope[session]  = pd.read_csv(StampsMiniscopeFile)
            dict_LFP_S1[session]  = S1
            dict_LFP_PFC[session]  = PFC
            dict_LFP_CA1[session]  = CA1
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
    
    #######################################################################################
                             # Detect Global Spindles #
    #######################################################################################

    for session in sessions:        
        folder_mini = folder_base / session / f'V4_Miniscope'
        nb_subsessions = sum(1 for p in folder_mini.iterdir() if p.is_dir() and p.name.startswith("session"))
        if nb_subsessions!=0:
            for x in range(1, nb_subsessions+1):            
                subsession= session + str(x)
                listSpdlPFC= dict_Spindleprop_PFC[subsession]
                listSpdlS1= dict_Spindleprop_S1[subsession]
                listSpdlPFCstarts=listSpdlPFC["start time"]
                listSpdlPFCends=listSpdlPFC["end time"]
                listSpdlS1starts=listSpdlS1["start time"]
                listSpdlS1ends=listSpdlS1["end time"]
                for ss in range(len(listSpdlPFCstarts)): # for PFC 
                    startPFC=listSpdlPFCstarts[ss]
                    endPFC=listSpdlPFCends[ss]
                    Istrue=is_overlapping(startPFC, endPFC, listSpdlS1starts, listSpdlS1ends)
                    listSpdlPFC.loc[ss, 'GlobalSpindle'] =Istrue
                dict_Spindleprop_PFC[subsession]=listSpdlPFC
                if x==1 : #no need to overwritte for each subsession
                    filenameOut = folder_base / session / f'OpenEphys/Spindlesproperties_PFC_sd5bis_AB.xlsx'
                    print(filenameOut)
                    writer = pd.ExcelWriter(filenameOut)
                    dict_Spindleprop_PFC[subsession].to_excel(writer)
                    writer.close()
                for ss in range(len(listSpdlS1starts)): # for S1 
                    startS1=listSpdlS1starts[ss]
                    endS1=listSpdlS1ends[ss]
                    Istrue=is_overlapping(startS1, endS1, listSpdlPFCstarts, listSpdlPFCends)
                    listSpdlS1.loc[ss, 'GlobalSpindle'] =Istrue       
                dict_Spindleprop_S1[subsession]=listSpdlS1
                if x==1 : #no need to overwritte for each subsession
                    filenameOut = folder_base /session / f'OpenEphys/Spindlesproperties_S1_sd5bis_AB.xlsx'
                    print(filenameOut)
                    writer = pd.ExcelWriter(filenameOut)
                    dict_Spindleprop_S1[subsession].to_excel(writer)
                    writer.close()
        else: 
            listSpdlPFC= dict_Spindleprop_PFC[session]
            listSpdlS1= dict_Spindleprop_S1[session]
            listSpdlPFCstarts=listSpdlPFC["start time"]
            listSpdlPFCends=listSpdlPFC["end time"]
            listSpdlS1starts=listSpdlS1["start time"]
            listSpdlS1ends=listSpdlS1["end time"]
            for ss in range(len(listSpdlPFCstarts)): # for PFC 
                startPFC=listSpdlPFCstarts[ss]
                endPFC=listSpdlPFCends[ss]
                Istrue=is_overlapping(startPFC, endPFC, listSpdlS1starts, listSpdlS1ends)
                listSpdlPFC.loc[ss, 'GlobalSpindle'] =Istrue
            dict_Spindleprop_PFC[session]=listSpdlPFC
            filenameOut = folder_base / session / f'OpenEphys/Spindlesproperties_PFC_sd5bis_AB.xlsx'
            print(filenameOut)
            writer = pd.ExcelWriter(filenameOut)
            dict_Spindleprop_PFC[session].to_excel(writer)
            writer.close()
            for ss in range(len(listSpdlS1starts)): # for S1 
                startS1=listSpdlS1starts[ss]
                endS1=listSpdlS1ends[ss]
                Istrue=is_overlapping(startS1, endS1, listSpdlPFCstarts, listSpdlPFCends)
                listSpdlS1.loc[ss, 'GlobalSpindle'] =Istrue       
            dict_Spindleprop_S1[session]=listSpdlS1
            filenameOut = folder_base /session / f'OpenEphys/Spindlesproperties_S1_sd5bis_AB.xlsx'
            print(filenameOut)
            writer = pd.ExcelWriter(filenameOut)
            dict_Spindleprop_S1[session].to_excel(writer)
            writer.close()

    #######################################################################################
      # Distribute Ca2+ intensity & spikes to oscillations for each sessions/subsessions #
    #######################################################################################

    CortexList= ['PFC', 'S1']

    for Cortex in CortexList:

        data = {}        
        before = 1000 # Max distance in ms between a SWR and a spindle to be considered as Precoupled
        after = 1000 # Max distance in ms between a spindle and a SWR to be considered as Postcoupled
        durationSpdl = int(2*1000) # number of sec before and after the Spdl onset taken into acount
        durationSWR = int(0.5*1000) # number of sec before and after the SWR onset taken into acount
        durationDS = int(1*1000) # number of sec before and after the DownStates onset taken into acount

        norm_freq=20 # final miniscope frequency used for all recordings
        previousEndTime=0
        InitialStartTime=0

        Spdls=[] 
        SWRs=[]      
        DSs=[]      
        
        for session in list(dict_Stamps.keys()):  
            

            # Start time & freq miniscope

            StartTime = list(dict_Stamps[session][0])[0] # in seconds
            minian_freq=list(dict_Stamps[session][0])[2] # in Hz

            if minian_freq>=20: # should only remove 1 session                

                dictnameSpdl=f'dict_Spindleprop_{Cortex}' #change dict according to the cortex 
                dictnameDS=f'dict_DSprop_{Cortex}' #change dict according to the cortex 
                dictSpdl = globals()[dictnameSpdl]
                dictDS = globals()[dictnameDS]

                ### Spindles ###

                SpipropO=dictSpdl[session]
                nb_spindle = SpipropO.shape[0]

                startSpiList = list(pd.Series(SpipropO["start time"]))
                endSpiList = list(pd.Series(SpipropO["end time"]))
                
                for Pspin in range(nb_spindle):                         
                    startSpi=int(startSpiList[Pspin])
                    endSpi=endSpiList[Pspin]  
                    LFPctxn=f'dict_LFP_{Cortex}'
                    LFPctx=globals()[LFPctxn]
                    LFPctx=LFPctx[session] 
                    Spdls.append(LFPctx[startSpi-durationSpdl:startSpi+durationSpdl])  
                
                ### SWRs ###

                SWRpropO=dict_SWRprop[session]  
                nb_swr = SWRpropO.shape[0]
                
                startSwrList = list(pd.Series(SWRpropO["start time"]))
                endSwrList = list(pd.Series(SWRpropO["end time"]))

                for Pswr in range(nb_swr): 
                    startSwr=int(startSwrList[Pswr])
                    endSwr=endSwrList[Pswr]
                    LFPctx=dict_LFP_CA1[session]
                    SWRs.append(LFPctx[startSwr-durationSWR:startSwr+durationSWR])             

                ### DSs ###

                DSpropO=dictDS[session]
                nb_ds = DSpropO.shape[0]

                startdsList = list(pd.Series(DSpropO["start time"]))
                enddsList = list(pd.Series(DSpropO["end time"]))

                for Pds in range(nb_ds): 
                    startds=int(startdsList[Pds])                
                    endds=enddsList[Pds]
                    LFPctxn=f'dict_LFP_{Cortex}'
                    LFPctx=globals()[LFPctxn]
                    LFPctx=LFPctx[session] 
                    DSs.append(LFPctx[startds-durationDS:startds+durationDS])
 

        #######################################################################################
                                # Save Spindles analysis #
        #######################################################################################

        # Save the big summary table Spindles_GlobalResults

        mice=os.path.basename(folder_base) 
        filenameOut = folder_to_save / f'Spindles_{Cortex}_RawLFPs_{mice}.xlsx'
        writer = pd.ExcelWriter(filenameOut)
        Spdls=pd.DataFrame(Spdls)
        Spdls.to_excel(writer)
        writer.close()

        filenameOut = folder_to_save / f'SWR_{Cortex}_RawLFPs_{mice}.xlsx'
        writer = pd.ExcelWriter(filenameOut)
        SWRs=pd.DataFrame(SWRs)
        SWRs.to_excel(writer)
        writer.close()

        filenameOut = folder_to_save / f'DS_{Cortex}_RawLFPs_{mice}.xlsx'
        writer = pd.ExcelWriter(filenameOut)
        DSs=pd.DataFrame(DSs)
        DSs.to_excel(writer)
        writer.close()