# Stats on Ca2+ imaging with miniscope and Vigilance States

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


########################################################################
        # SCRIPT 23AB_GrandAverages&Stats_for_VigilanceStates
########################################################################

# Specify the directory containing the Excel files
directory = "//10.69.168.1/crnldata/waking/audrey_hay/L1imaging/AnalysedMarch2023/Gaelle/Baseline_recording_ABmodified/AB_analysis/Analysis_VigStates_2024-06-10_18_18_30_990530_noN2/"
# Get the current date and time
FolderNameSave=str(datetime.now())
FolderNameSave = FolderNameSave.replace(" ", "_").replace(".", "_").replace(":", "_")
destination_folder= f"//10.69.168.1/crnldata/waking/audrey_hay/L1imaging/AnalysedMarch2023/Gaelle/Baseline_recording_ABmodified/AB_analysis/Analysis_AVG_VigStates_{FolderNameSave}/"
os.makedirs(destination_folder)
folder_to_save=Path(destination_folder)

# Copy the script file to the destination folder
source_script = "C:/Users/Manip2/SCRIPTS/Code python audrey/code python aurelie/HayLabAnalysis/python/23_24AB_GrandAverage&Stats_for_VigilanceStates_FullAuto.py"
destination_file_path = f"{destination_folder}/23_24AB_GrandAverage&Stats_for_VigilanceStates_FullAuto.txt"
shutil.copy(source_script, destination_file_path)

NrSubtypeList=['All', 'L1']

for NrSubtype in NrSubtypeList:  

    # Initialize an empty list to store the dataframes
    dfs = []
    df= []

    if NrSubtype=='L1':
        MiceList=['BlackLinesOK', 'BlueLinesOK', 'GreenDotsOK', 'GreenLinesOK', 'RedLinesOK']
    else:
        MiceList=['Purple', 'ThreeColDotsOK', 'ThreeBlueCrossesOK']
    
    nametofind='VigilanceState_GlobalResults'

    # Recursively traverse the directory structure
    for root, _, files in os.walk(directory):
        for filename in files:
            # Check if the file is an Excel file and contains the specified name
            if filename.endswith('.xlsx') and nametofind in filename : 
                if any(name in filename for name in MiceList): 
                    # Construct the full path to the file
                    filepath = os.path.join(root, filename)
                    # Read the Excel file into a dataframe and append it to the list
                    df = pd.read_excel(filepath, index_col=0)      
                    dfs.append(df)

    # Concatenate all dataframes into a single dataframe
    combined_df = pd.concat(dfs, ignore_index=True)

    combined_df['Unique_Unit'] = combined_df['Unique_Unit'].astype(str)
    combined_df['UnitNumber'] = combined_df['UnitNumber'].astype(str)
    combined_df['UnitValue'] = combined_df['UnitValue'].astype(str)

    combined_df['Unit_ID'] = combined_df['Mice'] + combined_df['Unique_Unit']

    #Remove non defined Unique Units 
    combined_df = combined_df[combined_df['Unique_Unit'] != '[]']
    
    combined_df['NormalizedAUC_calcium'] = combined_df['AUC_calcium'] / combined_df['DurationSubstate']

    combined_df['Substate_ID'] = combined_df['Mice'] + combined_df['Session'] + combined_df['Substate'] + combined_df['SubstateNumber'].astype(str)
    combined_df['Session_ID'] = combined_df['Mice'] + combined_df['Session'].astype(str)

    # Save big dataset for stats

    filenameOut = f'{folder_to_save}/{NrSubtype}_VigilanceStates_GrandGlobalAB.xlsx'
    writer = pd.ExcelWriter(filenameOut)
    combined_df.to_excel(writer)
    writer.close()

########################################################################
            # SCRIPT 24AB_Load&Stats_for_VigilanceStates
########################################################################

def max_column_name(row):
    return row.idxmax()

def normalize_row(row):
    max_col = row.idxmax()  # Find the column with the maximum value
    max_val = row[max_col]  # Get the maximum value
    return row / max_val   # Normalize the row by dividing by the maximum value    

for NrSubtype in NrSubtypeList:

    analysisfile='VigilanceStates_GrandGlobalAB'
    directory="//10.69.168.1/crnldata/waking/audrey_hay/L1imaging/AnalysedMarch2023/Gaelle/Baseline_recording_ABmodified/AB_analysis/"
    combined_df = pd.read_excel(f'{folder_to_save}/{NrSubtype}_{analysisfile}.xlsx', index_col=0)
    
    ######################
    # CHOOSE OPTIONS
    ######################

    # MINIMUM VIG STATES DURATION #
    combined_df = combined_df[combined_df['DurationSubstate'] >= 20] 

    desired_order = ['Wake','NREM', 'REM']   
    #desired_order = ['Wake', 'N2', 'NREM', 'REM'] 

    # NO LOW FIRING RATE #
    #combined_df = combined_df[combined_df['Avg_SpikeActivityHz'] >= 0.05] 

    #######################
    # DETERMINATION COEFF #
    #######################

    CorrCoeff_calcium_perUnit = combined_df.pivot_table(index='Unit_ID', columns='Substate', values='Rsquared_CaCorrCoeff', aggfunc='mean')
    CorrCoeff_calcium_perUnit = CorrCoeff_calcium_perUnit[desired_order]

    # Save CorrCoeff_calcium_perUnit
    filenameOut = f'{folder_to_save}/{NrSubtype}_VigSt_Rsquared_CaCorrCoeff_perUnitAB.xlsx'
    writer = pd.ExcelWriter(filenameOut)
    CorrCoeff_calcium_perUnit.to_excel(writer)
    writer.close()

    CorrCoeff_spike_perUnit = combined_df.pivot_table(index='Unit_ID', columns='Substate', values='Rsquared_SpCorrCoeff', aggfunc='mean')
    CorrCoeff_spike_perUnit = CorrCoeff_spike_perUnit[desired_order]

    # Save CorrCoeff_spike_perUnit 
    filenameOut = f'{folder_to_save}/{NrSubtype}_VigSt_Rsquared_SpCorrCoeff_perUnitAB.xlsx'
    writer = pd.ExcelWriter(filenameOut)
    CorrCoeff_spike_perUnit.to_excel(writer)
    writer.close()

    #####################    
    # CORRELATION COEFF #
    #####################

    CorrCoeff_calcium_perUnit = combined_df.pivot_table(index='Unit_ID', columns='Substate', values='CaCorrCoeff', aggfunc='mean')
    CorrCoeff_calcium_perUnit = CorrCoeff_calcium_perUnit[desired_order]

    # Save CorrCoeff_calcium_perUnit
    filenameOut = f'{folder_to_save}/{NrSubtype}_VigSt_CaCorrCoeff_perUnitAB.xlsx'
    writer = pd.ExcelWriter(filenameOut)
    CorrCoeff_calcium_perUnit.to_excel(writer)
    writer.close()

    CorrCoeff_spike_perUnit = combined_df.pivot_table(index='Unit_ID', columns='Substate', values='SpCorrCoeff', aggfunc='mean')
    CorrCoeff_spike_perUnit = CorrCoeff_spike_perUnit[desired_order]

    # Save CorrCoeff_spike_perUnit 
    filenameOut = f'{folder_to_save}/{NrSubtype}_VigSt_SpCorrCoeff_perUnitAB.xlsx'
    writer = pd.ExcelWriter(filenameOut)
    CorrCoeff_spike_perUnit.to_excel(writer)
    writer.close()
    
    #####################
    # AUC CALCIUM #
    #####################
    
    resultNormalizedAUC_calcium_perUnit = combined_df.pivot_table(index='Unit_ID', columns='Substate', values='NormalizedAUC_calcium', aggfunc='mean')
    resultNormalizedAUC_calcium_perUnit = resultNormalizedAUC_calcium_perUnit[desired_order]

    max_values_per_row = resultNormalizedAUC_calcium_perUnit.max(axis=1)
    Diff_resultNormalizedAUC_calcium_perUnit = (resultNormalizedAUC_calcium_perUnit.sub(max_values_per_row, axis=0) +100)
    Diff_resultNormalizedAUC_calcium_perUnit['Activated_by'] = resultNormalizedAUC_calcium_perUnit.apply(max_column_name, axis=1)

    # Add a new column with the column name that has the maximum value
    resultNormalizedAUC_calcium_perUnit['Activated_by'] = resultNormalizedAUC_calcium_perUnit.apply(max_column_name, axis=1)

    # Save  df
    filenameOut = f'{folder_to_save}/{NrSubtype}_VigSt_AUC_perUnitAB_N.xlsx'
    writer = pd.ExcelWriter(filenameOut)
    resultNormalizedAUC_calcium_perUnit.to_excel(writer)
    writer.close()

    # Save  df
    filenameOut = f'{folder_to_save}/{NrSubtype}_DiffVigSt_AUC_perUnitAB_N.xlsx'
    writer = pd.ExcelWriter(filenameOut)
    Diff_resultNormalizedAUC_calcium_perUnit.to_excel(writer)
    writer.close()

    proportions = resultNormalizedAUC_calcium_perUnit['Activated_by'].value_counts(normalize=True)*100
    #print(proportions)

    resultNormalizedAUC_calcium_perMouse = combined_df.pivot_table(index='Mice', columns='Substate', values='NormalizedAUC_calcium', aggfunc='mean')
    resultNormalizedAUC_calcium_perMouse = resultNormalizedAUC_calcium_perMouse[desired_order]

    filenameOut = f'{folder_to_save}/{NrSubtype}_VigSt_AUC_perMouseAB.xlsx'
    writer = pd.ExcelWriter(filenameOut)
    resultNormalizedAUC_calcium_perMouse.to_excel(writer)
    writer.close()

    #####################
    # SPIKE ACTIVITY #
    #####################

    resultSpikeActivity_perUnit = combined_df.pivot_table(index='Unit_ID', columns='Substate', values='SpikeActivityHz', aggfunc='mean')    
    resultSpikeActivity_perUnit = resultSpikeActivity_perUnit[desired_order]

    max_values_per_row = resultSpikeActivity_perUnit.max(axis=1)
    Diff_resultSpikeActivity_perUnit = (resultSpikeActivity_perUnit.sub(max_values_per_row, axis=0) +100)
    Diff_resultSpikeActivity_perUnit['Activated_by'] = resultSpikeActivity_perUnit.apply(max_column_name, axis=1)

    # Add a new column with the column name that has the maximum value
    resultSpikeActivity_perUnit['Activated_by'] = resultSpikeActivity_perUnit.apply(max_column_name, axis=1)

    # Save  df
    filenameOut = f'{folder_to_save}/{NrSubtype}_VigSt_Spike_perUnitAB.xlsx'
    writer = pd.ExcelWriter(filenameOut)
    resultSpikeActivity_perUnit.to_excel(writer)
    writer.close()

    # Save  df
    filenameOut = f'{folder_to_save}/{NrSubtype}_DiffVigSt_Spike_perUnitAB_N.xlsx'
    writer = pd.ExcelWriter(filenameOut)
    Diff_resultSpikeActivity_perUnit.to_excel(writer)
    writer.close()

    resultSpikeActivity_perMouse = combined_df.pivot_table(index='Mice', columns='Substate', values='SpikeActivityHz', aggfunc='mean')
    resultSpikeActivity_perMouse = resultSpikeActivity_perMouse[desired_order]

    filenameOut = f'{folder_to_save}/{NrSubtype}_VigSt_Spike_perMouseAB.xlsx'
    writer = pd.ExcelWriter(filenameOut)
    resultSpikeActivity_perMouse.to_excel(writer)
    writer.close()

    proportions = resultSpikeActivity_perUnit['Activated_by'].value_counts(normalize=True)*100
    filenameOut = f'{folder_to_save}/{NrSubtype}_ActivityVigSt_proportions.xlsx'
    writer = pd.ExcelWriter(filenameOut)
    proportions.to_excel(writer)
    writer.close()

    #####################
    # ProportionVigStates
    #####################

    ProportionVigStates = combined_df.pivot_table(index='Session_ID', columns='Substate', values='DurationSubstate', aggfunc='sum')
    ProportionVigStates = ProportionVigStates[desired_order]
    
    ProportionVigStates2=ProportionVigStates.div(ProportionVigStates.sum(axis=1), axis=0)*100

    filenameOut = f'{folder_to_save}/{NrSubtype}_ProportionVigStatesAB.xlsx'
    writer = pd.ExcelWriter(filenameOut)
    ProportionVigStates2.to_excel(writer)
    writer.close()

    #####################
    # PREFERENCE #
    #####################
    #"""
    combined_dfNREM=combined_df.copy()
    combined_dfNREM['Substate'] = combined_dfNREM['Substate'].apply(lambda x: 'NREM' if x == 'NREM' else 'other')
    grouped = combined_dfNREM.groupby('Unit_ID')
    NREMprefUnits=[]
    NREMmatrix={}
    for key, group in grouped:
        nrem_values = group[group['Substate'] == 'NREM']['SpikeActivityHz'].tolist()
        other_values = group[group['Substate'] == 'other']['SpikeActivityHz'].tolist()
        NREMmatrix[key]=[nrem_values, other_values]
        if len(nrem_values) > 1 and len(other_values) > 1:  # Ensure there are enough values to perform the test
            t_stat, p_value = ttest_ind(nrem_values, other_values, equal_var=False)
            if p_value<0.2 and np.mean(nrem_values)>np.mean(other_values): 
                NREMprefUnits.append(key)
    
    combined_dfREM=combined_df.copy()
    combined_dfREM['Substate'] = combined_dfREM['Substate'].apply(lambda x: 'REM' if x == 'REM' else 'other')
    grouped = combined_dfREM.groupby('Unit_ID')
    REMprefUnits=[]
    REMmatrix={}
    for key, group in grouped:
        rem_values = group[group['Substate'] == 'REM']['SpikeActivityHz'].tolist()
        other_values = group[group['Substate'] == 'other']['SpikeActivityHz'].tolist()
        REMmatrix[key]=[rem_values, other_values]
        if len(rem_values) > 1 and len(other_values) > 1:  # Ensure there are enough values to perform the test
            t_stat, p_value = ttest_ind(rem_values, other_values, equal_var=False)
            if p_value<0.2  and np.mean(rem_values)>np.mean(other_values):  
                REMprefUnits.append(key)
 
    combined_dfWake=combined_df.copy()
    combined_dfWake['Substate'] = combined_dfWake['Substate'].apply(lambda x: 'Wake' if x == 'Wake' else 'other')
    grouped = combined_dfWake.groupby('Unit_ID')
    WakeprefUnits=[]
    Wmatrix={}
    for key, group in grouped:
        w_values = group[group['Substate'] == 'Wake']['SpikeActivityHz'].tolist()
        other_values = group[group['Substate'] == 'other']['SpikeActivityHz'].tolist()
        Wmatrix[key]=[w_values, other_values]
        if len(w_values) > 1 and len(other_values) > 1:  # Ensure there are enough values to perform the test
            t_stat, p_value = ttest_ind(w_values, other_values, equal_var=False)
            if p_value<0.2 and np.mean(w_values)>np.mean(other_values):  
                WakeprefUnits.append(key)

    # Save  df
    filenameOut = f'{folder_to_save}/{NrSubtype}_SignFiringPreference.xlsx'
    writer = pd.ExcelWriter(filenameOut)
    NREMprefUnits = pd.DataFrame(NREMprefUnits)
    REMprefUnits = pd.DataFrame(REMprefUnits)
    WakeprefUnits = pd.DataFrame(WakeprefUnits)
    NREMprefUnits.to_excel(writer, sheet_name='NREMpref', index=True, header=False) 
    REMprefUnits.to_excel(writer, sheet_name='REMpref', index=True, header=False) 
    WakeprefUnits.to_excel(writer, sheet_name='Wakepref', index=True, header=False) 
    writer.close()

    excel_file =  f'{folder_to_save}/{NrSubtype}_NREMmatrix.xlsx'
    with pd.ExcelWriter(excel_file) as writer:
        for key, df in NREMmatrix.items():
            df=pd.DataFrame(df)
            df.index = ['NREM', 'Other']
            df.to_excel(writer, sheet_name=key, header=False)

    excel_file =  f'{folder_to_save}/{NrSubtype}_REMmatrix.xlsx'
    with pd.ExcelWriter(excel_file) as writer:
        for key, df in REMmatrix.items():
            df=pd.DataFrame(df)
            df.index = ['REM', 'Other']
            df.to_excel(writer, sheet_name=key, header=False)

    excel_file =  f'{folder_to_save}/{NrSubtype}_Wakematrix.xlsx'
    with pd.ExcelWriter(excel_file) as writer:
        for key, df in Wmatrix.items():
            df=pd.DataFrame(df)
            df.index = ['Wake', 'Other']
            df.to_excel(writer, sheet_name=key, header=False)
    #"""
