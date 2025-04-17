# Stats on Ca2+ imaging with miniscope and Vigilance States

#######################################################################################
                            # Define Experiment type #
#######################################################################################

AnalysisID='_Correlations' #to identify this analysis from another
DrugExperiment=0 # 0 if Baseline, 1 if CGP, 2 if Baseline & CGP

saveexcel=1
Local=1

#choosed_folder='VigSt_2024-07-22_18_21_32_AB_FINAL' if DrugExperiment else 'VigSt_2024-07-22_17_16_28_AB_FINAL'
choosed_folder1='VigSt_2025-04-16_18_40_04_CellAssembly' # for Baseline Expe
choosed_folder2='VigSt_' # for CGP Expe

desired_order = ['AW','QW', 'NREM', 'IS', 'REM', 'undefined']   

#######################################################################################
                                # Load packages #
#######################################################################################

import statsmodels.api as sm
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
import pickle
import os
import statsmodels.api as sm
import statsmodels.formula.api as smf
from datetime import datetime
import shutil
from scipy.stats import ttest_ind
import statsmodels.api as sm
import re
import copy
from collections import defaultdict

import warnings
warnings.filterwarnings("ignore")

def divide_keys(data, startkey, everykey):
    for it in range(startkey, len(data), everykey):        
        key2 = list(data.keys())[it-1]
        key3 = list(data.keys())[it]
        d=data[key3]
        data[key3]=d.replace(0, np.nan)
        if startkey>1:
            key1 = list(data.keys())[it-2]
            data[key1] = data[key1] / data[key3]
        data[key2] = data[key2] / data[key3]
    keys_to_delete = list(data.keys())[startkey::everykey]
    for key in keys_to_delete:
        del data[key]
    return data   

########################################################################
        # SCRIPT 23AB_GrandAverages&Stats_for_VigilanceStates
########################################################################

# Specify the directory containing the Excel files
InitialDirectory1 = "//10.69.168.1/crnldata/waking/audrey_hay/L1imaging/Analysed2025_AB/_baseline_analysis" if Local else "/crnldata/waking/audrey_hay/L1imaging/Analysed2025_AB/_baseline_analysis" 
directory1= f'{InitialDirectory1}/{choosed_folder1}'
InitialDirectory2 ="//10.69.168.1/crnldata/waking/audrey_hay/L1imaging/Analysed2025_AB/_CGP_analysis" if Local else "/crnldata/waking/audrey_hay/L1imaging/Analysed2025_AB/_CGP_analysis"
directory2= f'{InitialDirectory2}/{choosed_folder2}'

# Get the current date and time
FolderNameSave=str(datetime.now())[:19]
FolderNameSave = FolderNameSave.replace(" ", "_").replace(".", "_").replace(":", "_")
destination_folder= f"//10.69.168.1/crnldata/waking/audrey_hay/L1imaging/Analysed2025_AB/_global_analysis/AVG_VigSt_{FolderNameSave}{AnalysisID}" if Local else f"/crnldata/waking/audrey_hay/L1imaging/Analysed2025_AB/_global_analysis/AVG_VigSt_{FolderNameSave}{AnalysisID}"
os.makedirs(destination_folder)
folder_to_save=Path(destination_folder)

# Copy the script file to the destination folder
source_script = "C:/Users/Manip2/SCRIPTS/CodePythonAudrey/CodePythonAurelie/HayLabAnalysis/python/_MINI&OE_3_average_correlations.py" if Local else "/python/_MINI&OE_3_average_correlations.py" 
destination_file_path = f"{destination_folder}/_MINI&OE_3_average_correlations.txt"
shutil.copy(source_script, destination_file_path)

directories= [directory1, directory2] if DrugExperiment else [directory1]

NrSubtypeList=['L1NDNF_mice','L2_3_mice']

for NrSubtype in NrSubtypeList:  
    dfs2_per_sheet = {}
    dfs3_per_sheet = {}
    dfs4_per_sheet = {}
    dfs5_per_sheet = {}
    dfs6_per_sheet = {}

    if NrSubtype=='L1NDNF_mice':
        MiceList=['BlackLines', 'BlueLines', 'GreenDots', 'GreenLines', 'RedLines']
    else:
        MiceList=['Purple', 'ThreeColDots', 'ThreeBlueCrosses']
    
    nametofind2='VigSt_CaCorr'
    nametofind3='VigSt_SpCorr'      
    nametofind4='TotCaCorr'
    nametofind5='TotSpCorr'
    nametofind6='StatesCaCorr'

    # Recursively traverse the directory structure
    for directory in directories:
        for root, _, files in os.walk(directory):
            for filename in files:
                # Check if the file is an pkl file and contains the specified name
                if filename.endswith('.pkl') and nametofind2 in filename: 
                    if any(name in filename for name in MiceList): 
                        # Construct the full path to the file
                        filepath = os.path.join(root, filename)
                        with open(filepath, 'rb') as pickle_file:
                            df = pickle.load(pickle_file)
                        for key, value in df.items():
                            if key in dfs2_per_sheet:
                                dfs2_per_sheet[key]=pd.concat([dfs2_per_sheet[key],value],axis=0)
                            else:
                                dfs2_per_sheet[key]=value
                        print(filename)
                if filename.endswith('.pkl') and nametofind3 in filename: 
                    if any(name in filename for name in MiceList): 
                        # Construct the full path to the file
                        filepath = os.path.join(root, filename)
                        with open(filepath, 'rb') as pickle_file:
                            df = pickle.load(pickle_file)
                        for key, value in df.items():
                            if key in dfs3_per_sheet:
                                dfs3_per_sheet[key]=pd.concat([dfs3_per_sheet[key],value],axis=0)
                            else:
                                dfs3_per_sheet[key]=value
                        print(filename)
                if filename.endswith('.pkl') and nametofind4 in filename: 
                    if any(name in filename for name in MiceList): 
                        # Construct the full path to the file
                        filepath = os.path.join(root, filename)
                        with open(filepath, 'rb') as pickle_file:
                            df = pickle.load(pickle_file)
                        for key, value in df.items():
                            if key in dfs4_per_sheet:
                                dfs4_per_sheet[key]=pd.concat([dfs4_per_sheet[key],value],axis=0)
                            else:
                                dfs4_per_sheet[key]=value
                        print(filename)
                if filename.endswith('.pkl') and nametofind5 in filename: 
                    if any(name in filename for name in MiceList): 
                        # Construct the full path to the file
                        filepath = os.path.join(root, filename)
                        with open(filepath, 'rb') as pickle_file:
                            df = pickle.load(pickle_file)
                        for key, value in df.items():
                            if key in dfs5_per_sheet:
                                dfs5_per_sheet[key]=pd.concat([dfs5_per_sheet[key],value],axis=0)
                            else:
                                dfs5_per_sheet[key]=value
                        print(filename)
                if filename.endswith('.pkl') and nametofind6 in filename: 
                    if any(name in filename for name in MiceList): 
                        # Construct the full path to the file
                        filepath = os.path.join(root, filename)
                        with open(filepath, 'rb') as pickle_file:
                            try : df = pickle.load(pickle_file)
                            except: pass
                        for key, value in df.items():
                            if key in dfs6_per_sheet:
                                dfs6_per_sheet[key]=pd.concat([dfs6_per_sheet[key],value],axis=0)
                            else:
                                dfs6_per_sheet[key]=value
                        print(filename)

    ######### Save the SubStates Ca correlation matrix   ########

    dfs6_per_sheet = {sheet_name: df.groupby(df.index).sum() for sheet_name, df in dfs6_per_sheet.items()} #cause was concatenated in the 0 axis
    dfs6_per_sheet=divide_keys(dfs6_per_sheet, 1, 2)
    for sheet_name, df in dfs6_per_sheet.items():
        df = df.sort_index(axis=1)
        df = df.sort_index(axis=0)
        dfs6_per_sheet[sheet_name]=df

    if saveexcel:
        file_path = f'{folder_to_save}/{NrSubtype}_SubSt_CaCorr.xlsx'
        with pd.ExcelWriter(file_path) as writer:        
            for sheet_name, df in dfs6_per_sheet.items():
                df.to_excel(writer, sheet_name=sheet_name, index=True, header=True)

    filenameOut = f'{folder_to_save}/{NrSubtype}_SubSt_CaCorr.pkl'
    with open(filenameOut, 'wb') as pickle_file:
        pickle.dump(dfs6_per_sheet, pickle_file)

    ######### Save the Ca correlation matrix   ########

    dfs2_per_sheet = {sheet_name: df.groupby(df.index).sum() for sheet_name, df in dfs2_per_sheet.items()} #cause was concatenated in the 0 axis
    dfs2_per_sheet=divide_keys(dfs2_per_sheet, 2, 3)
    for sheet_name, df in dfs2_per_sheet.items():
        df = df.sort_index(axis=1)
        df = df.sort_index(axis=0)
        dfs2_per_sheet[sheet_name]=df

    if saveexcel:
        file_path = f'{folder_to_save}/{NrSubtype}_VigSt_CaCorr.xlsx'
        with pd.ExcelWriter(file_path) as writer:        
            for sheet_name, df in dfs2_per_sheet.items():
                df.to_excel(writer, sheet_name=sheet_name, index=True, header=True)

    filenameOut = f'{folder_to_save}/{NrSubtype}_VigSt_CaCorr.pkl'
    with open(filenameOut, 'wb') as pickle_file:
        pickle.dump(dfs2_per_sheet, pickle_file)

    ######### Save the Sp correlation matrix  ########

    dfs3_per_sheet = {sheet_name: df.groupby(df.index).sum() for sheet_name, df in dfs3_per_sheet.items()}
    dfs3_per_sheet=divide_keys(dfs3_per_sheet, 2, 3)
    for sheet_name, df in dfs3_per_sheet.items():
        df = df.sort_index(axis=1)
        df = df.sort_index(axis=0)
        dfs3_per_sheet[sheet_name]=df

    if saveexcel:    
        file_path = f'{folder_to_save}/{NrSubtype}_VigSt_SpCorr.xlsx'
        with pd.ExcelWriter(file_path) as writer:        
            for sheet_name, df in dfs3_per_sheet.items():
                df.to_excel(writer, sheet_name=sheet_name, index=True, header=True)

    filenameOut = f'{folder_to_save}/{NrSubtype}_VigSt_SpCorr.pkl'
    with open(filenameOut, 'wb') as pickle_file:
        pickle.dump(dfs3_per_sheet, pickle_file)

    ####### Save the TOT Ca correlation matrix  ########

    dfs4_per_sheet = {sheet_name: df.groupby(df.index).sum() for sheet_name, df in dfs4_per_sheet.items()} #cause was concatenated in the 0 axis
    dfs4_per_sheet=divide_keys(dfs4_per_sheet, 2, 3)
    for sheet_name, df in dfs4_per_sheet.items():
        df = df.sort_index(axis=1)
        df = df.sort_index(axis=0)
        dfs4_per_sheet[sheet_name]=df

    if saveexcel:
        file_path = f'{folder_to_save}/{NrSubtype}_Tot_CaCorr.xlsx'      
        with pd.ExcelWriter(file_path) as writer:  
            for sheet_name, df in dfs4_per_sheet.items():
                df.to_excel(writer, sheet_name=sheet_name, index=True, header=True)

    filenameOut = f'{folder_to_save}/{NrSubtype}_Tot_CaCorr.pkl'
    with open(filenameOut, 'wb') as pickle_file:
        pickle.dump(dfs4_per_sheet, pickle_file)

    ######### Save the TOT Sp correlation matrix   ########

    dfs5_per_sheet = {sheet_name: df.groupby(df.index).sum() for sheet_name, df in dfs5_per_sheet.items()}
    dfs5_per_sheet=divide_keys(dfs5_per_sheet, 2, 3)
    for sheet_name, df in dfs5_per_sheet.items():
        df = df.sort_index(axis=1)
        df = df.sort_index(axis=0)
        dfs5_per_sheet[sheet_name]=df

    if saveexcel:   
        file_path = f'{folder_to_save}/{NrSubtype}_Tot_SpCorr.xlsx' 
        with pd.ExcelWriter(file_path) as writer:        
            for sheet_name, df in dfs5_per_sheet.items():
                df.to_excel(writer, sheet_name=sheet_name, index=True, header=True)

    filenameOut = f'{folder_to_save}/{NrSubtype}_Tot_SpCorr.pkl'
    with open(filenameOut, 'wb') as pickle_file:
        pickle.dump(dfs5_per_sheet, pickle_file)