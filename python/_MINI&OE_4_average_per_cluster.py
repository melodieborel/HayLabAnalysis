# Stats on Ca2+ imaging with miniscope and Vigilance States

#######################################################################################
                            # Define Experiment type #
#######################################################################################

AnalysisID='_CGP' #to identify this analysis from another
DrugExperiment=1 # 0 if Baseline, 1 if CGP, 2 if Baseline & CGP

saveexcel=0
Local=1

desired_order = ['AW','QW', 'NREM', 'IS', 'REM', 'undefined']   

#choosed_folder='VigSt_2024-07-22_18_21_32_AB_FINAL' if DrugExperiment else 'VigSt_2024-07-22_17_16_28_AB_FINAL'
choosed_folder1='VigSt_2025-05-03_10_01_32' # for Baseline Expe
choosed_folder2='VigSt_2025-05-03_12_01_21' # for CGP Expe
choosed_folder3='Corr_VigSt_2025-05-03_14_47_15' # for Correlations

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

def extract_micename(index_value):
    match = re.match(r"([a-zA-Z]+)", index_value)
    return match.group(1) if match else index_value

def max_column_name(row):
    return row.idxmax()
    

#######################################################################################
                            # Define Directories #
#######################################################################################

# Specify the directory containing the Excel files
InitialDirectory1 = "//10.69.168.1/crnldata/waking/audrey_hay/L1imaging/Analysed2025_AB/_baseline_analysis" if Local else "/crnldata/waking/audrey_hay/L1imaging/Analysed2025_AB/_baseline_analysis" 
directory1= f'{InitialDirectory1}/{choosed_folder1}'
InitialDirectory2 ="//10.69.168.1/crnldata/waking/audrey_hay/L1imaging/Analysed2025_AB/_CGP_analysis" if Local else "/crnldata/waking/audrey_hay/L1imaging/Analysed2025_AB/_CGP_analysis"
directory2= f'{InitialDirectory2}/{choosed_folder2}'

InitialDirectory3 ="//10.69.168.1/crnldata/waking/audrey_hay/L1imaging/Analysed2025_AB/_global_analysis" if Local else "/crnldata/waking/audrey_hay/L1imaging/Analysed2025_AB/_global_analysis"
directory3= f'{InitialDirectory3}/{choosed_folder3}'

# Get the current date and time
FolderNameSave=str(datetime.now())[:19]
FolderNameSave = FolderNameSave.replace(" ", "_").replace(".", "_").replace(":", "_")
destination_folder= f"//10.69.168.1/crnldata/waking/audrey_hay/L1imaging/Analysed2025_AB/_global_analysis/AVG_VigSt_{FolderNameSave}{AnalysisID}" if Local else f"/crnldata/waking/audrey_hay/L1imaging/Analysed2025_AB/_global_analysis/AVG_VigSt_{FolderNameSave}{AnalysisID}"
os.makedirs(destination_folder)
folder_to_save=Path(destination_folder)

# Copy the script file to the destination folder
source_script = "C:/Users/Manip2/SCRIPTS/CodePythonAudrey/CodePythonAurelie/HayLabAnalysis/python/_MINI&OE_4_average_per_cluster.py" if Local else "/python/_MINI&OE_4_average_per_cluster.py" 
destination_file_path = f"{destination_folder}/_MINI&OE_4_average_per_cluster.txt"
shutil.copy(source_script, destination_file_path)


#######################################################################################
                            # Do the average per cluster #
#######################################################################################

# Load Global Table 
with open(f'{directory1}/VigStates_Global_cluster.pkl', 'rb') as pickle_file:
    combined_dfO = pickle.load(pickle_file)

AllProportionVigStates=pd.DataFrame()
AllDurationVigStates=pd.DataFrame()
AllTotDurationVigStates=pd.DataFrame()

NrSubtypeList=['L1NDNF_mice','L2_3_mice']

for NrSubtype in NrSubtypeList:

    analysisfileCa='VigSt_CaCorr'
    with open(f'{directory3}/{NrSubtype}_{analysisfileCa}.pkl', 'rb') as pickle_file:
        combined_dfCa = pickle.load(pickle_file)
        del combined_dfCa['baseline_IS']
        del combined_dfCa['baseline_undefined']        
        del combined_dfCa['Z_baseline_IS']
        del combined_dfCa['Z_baseline_undefined']

    analysisfileSp='VigSt_SpCorr'
    with open(f'{directory3}/{NrSubtype}_{analysisfileSp}.pkl', 'rb') as pickle_file:
        combined_dfSp = pickle.load(pickle_file)   
        del combined_dfSp['baseline_IS']
        del combined_dfSp['baseline_undefined']        
        del combined_dfSp['Z_baseline_IS']
        del combined_dfSp['Z_baseline_undefined']

    analysisfileCa='Tot_CaCorr'
    with open(f'{directory3}/{NrSubtype}_{analysisfileCa}.pkl', 'rb') as pickle_file:
        combined_dfCaTot = pickle.load(pickle_file)

    analysisfileSp='Tot_SpCorr'
    with open(f'{directory3}/{NrSubtype}_{analysisfileSp}.pkl', 'rb') as pickle_file:
        combined_dfSpTot = pickle.load(pickle_file)
    
    analysisfileCa='SubSt_CaCorr'
    with open(f'{directory3}/{NrSubtype}_{analysisfileCa}.pkl', 'rb') as pickle_file:
        combined_dfCaSub = pickle.load(pickle_file)
        del combined_dfCaSub['baseline_IS']
        del combined_dfCaSub['baseline_undefined']
        del combined_dfCaSub['CGP_IS']
        del combined_dfCaSub['CGP_undefined']

    ######################
    # CHOOSE OPTIONS
    ######################

    # MINIMUM VIG STATES DURATION #
    #combined_df = combined_dfO[combined_dfO['DurationSubstate'] >= 0] #10
    
    Drugs= ['baseline', 'CGP'] if DrugExperiment else ['baseline']

    # NO LOW FIRING RATE #
    #combined_df = combined_df[combined_df['Avg_SpikeActivityHz'] >= 0.05] 
    
    #####################
    # LOAD CLUSTERS #
    #####################

    # /!/ The ones from Baseline (not CGP)

    combined_df = combined_dfO.copy()
    combined_df = combined_df[combined_df['NeuronType'] == NrSubtype]
    combined_df_Drug = combined_df.copy()
    combined_df_Drug = combined_df[combined_df['Drug'] == 'baseline']

    AllBaselineUnits = combined_df_Drug['Unit_ID'].unique()
    Cluster0units = combined_df_Drug[combined_df_Drug['ClusterHDBSCAN'] == 0]['Unit_ID'].unique()
    Cluster1units = combined_df_Drug[combined_df_Drug['ClusterHDBSCAN'] == 1]['Unit_ID'].unique()
    Cluster2units = combined_df_Drug[combined_df_Drug['ClusterHDBSCAN'] == 2]['Unit_ID'].unique()
    
    # Only keep units that appears in CGP & Baseline    
    if DrugExperiment>=1 :
        combined_df_CGP=combined_df.copy()
        combined_df_CGP = combined_df_CGP[combined_df_CGP['Drug'] == 'CGP'] 
        AllCGPUnits = combined_df_CGP['Unit_ID'].unique()

        AllBaselineUnits= np.intersect1d(AllBaselineUnits,AllCGPUnits)
        Cluster0units= np.intersect1d(Cluster0units,AllBaselineUnits)
        Cluster1units= np.intersect1d(Cluster1units,AllBaselineUnits)
        Cluster2units= np.intersect1d(Cluster2units,AllBaselineUnits)
        
    # Save the List of significant Unit more active in one vigilance state
    if NrSubtype=='L1NDNF_mice':
        os.makedirs(f'{folder_to_save}/baseline/')

    for Drug in Drugs:

        combined_df_DrugO=combined_df.copy() # no min vig states durations == for vig states stats
        combined_df_DrugO = combined_df_DrugO[combined_df_DrugO['Drug'] == Drug] 

        combined_df_Drug=combined_df.copy()
        combined_df_Drug = combined_df_Drug[combined_df_Drug['Drug'] == Drug] 
        AllUnits= combined_df_Drug['Unit_ID'].unique()

        folder_to_save2= f'{folder_to_save}/{Drug}/'
        if NrSubtype=='L1NDNF_mice' and Drug=='CGP':
            os.makedirs(folder_to_save2)

        List_SignFiringPreference=[Cluster0units, Cluster1units, Cluster2units, AllBaselineUnits, AllUnits] if DrugExperiment else [Cluster0units, Cluster1units, Cluster2units, AllUnits]
        SecondaryList=[Cluster1units, Cluster2units, Cluster0units, AllBaselineUnits, AllUnits] if DrugExperiment else [Cluster1units, Cluster2units, Cluster0units, AllUnits]
        List_Names=['Cluster0units', 'Cluster1units', 'Cluster2units', 'AllBaselineUnits', 'All']if DrugExperiment else ['Cluster0units', 'Cluster1units','Cluster2units', 'All']
        SecondaryList_Names=['Cluster1units', 'Cluster2units', 'Cluster0units', 'AllBaselineUnits', 'All' ] if DrugExperiment else ['Cluster1units', 'Cluster2units', 'Cluster0units', 'All']
        
        for listnb, listI  in enumerate(List_SignFiringPreference):

            if len(listI)>0:
            
                filtered_df = combined_df_Drug[combined_df_Drug['Unit_ID'].isin(listI)]
                List_name=List_Names[listnb]
                SecondaryList_name=SecondaryList_Names[listnb]
                listII=SecondaryList[listnb]

                filtered_df_AllDrug = combined_df[combined_df['Unit_ID'].isin(listI)]
                filenameOut = f'{folder_to_save}/GLM_{NrSubtype}_{List_name}_VigSt_Global.xlsx'
                writer = pd.ExcelWriter(filenameOut)
                filtered_df_AllDrug.to_excel(writer)
                writer.close()

                if NrSubtype=='L1NDNF_mice':
                    new_folder= f"{folder_to_save2}/{List_name}/"
                    os.makedirs(new_folder)

                if Drug==Drug: #'baseline':
                    
                    #####################################################
                    ## TOTAL Ca correlation with neuron from same population ##
                    #####################################################

                    # Keep only neurons from the list 
                    dfCaTot_filtered={}
                    for sheet_name, dfCa in combined_dfCaTot.items():
                        dfCa=pd.DataFrame(dfCa)
                        indices_to_keep_existing = [idx for idx in listI if idx in dfCa.index] #from first list
                        columns_to_keep_existing = [col for col in listI if col in dfCa.columns] #from second list
                        dfCaTot_filtered[sheet_name] = dfCa.loc[indices_to_keep_existing, columns_to_keep_existing]
                    
                    if DrugExperiment>=1:
                        # Keep only correlation pairs that occurs for each Drug condition
                        for sheet_name, df in dfCaTot_filtered.items(): #remove inactive/non recorded neurons
                            df = df[~(df.fillna(0) == 0).all(axis=1)]
                            df = df.loc[:, ~(df.fillna(0) == 0).all(axis=0)]
                            dfCaTot_filtered[sheet_name] = df    
                        first_key = list(dfCaTot_filtered.keys())[0]
                        common_columns = dfCaTot_filtered[first_key].columns
                        common_indices = dfCaTot_filtered[first_key].index
                        for df in dfCaTot_filtered.values():
                            common_columns = common_columns.intersection(df.columns)
                            common_indices = common_indices.intersection(df.index)
                        dfCaTot_Doublefiltered = {name: df.loc[common_indices, common_columns] for name, df in dfCaTot_filtered.items()}
                    else:
                        dfCaTot_Doublefiltered=dfCaTot_filtered

                    # Save Corr Pairs
                    file_path = f'{folder_to_save2}/{List_name}/{NrSubtype}_Tot_PairCorrCa_{List_name}.xlsx'
                    with pd.ExcelWriter(file_path) as writer:
                        for sheet_name, dfCa in dfCaTot_Doublefiltered.items():
                            if 'Z_' not in sheet_name:                            
                                dfCa = dfCa.sort_index(axis=1)
                                dfCa = dfCa.sort_index(axis=0)
                                dfCa.to_excel(writer, sheet_name=sheet_name, index=True, header=True)

                    # Flat correlations
                    SummaryMatrixCa= pd.DataFrame()
                    for sheet_name, df in dfCaTot_Doublefiltered.items():    
                        series_flattened = df.stack().reset_index()
                        series_flattened['index'] = series_flattened['level_0'] + '_' + series_flattened['level_1']
                        series_flattened = series_flattened.set_index('index')[0]
                        series_flattened_cleaned = series_flattened[series_flattened.index.str.split('_').str[0] != series_flattened.index.str.split('_').str[1]] #remove Neuron1 vs Neuron1
                        series_flattened_cleaned.name = sheet_name
                        SummaryMatrixCa = pd.concat([SummaryMatrixCa,series_flattened_cleaned], axis=1)
                    
                    SummaryMatrixCa_cleaned = SummaryMatrixCa.round(5) # to better detect duplicate                   
                    SummaryMatrixCa_cleaned = SummaryMatrixCa.drop_duplicates(subset=SummaryMatrixCa.columns[1:]) 
                    #SummaryMatrixCa_cleaned = SummaryMatrixCa_cleaned.drop(columns=[SummaryMatrixCa_cleaned.columns[0], SummaryMatrixCa_cleaned.columns[2], SummaryMatrixCa_cleaned.columns[4], SummaryMatrixCa_cleaned.columns[6], SummaryMatrixCa_cleaned.columns[8], SummaryMatrixCa_cleaned.columns[10]]) if DrugExperiment else SummaryMatrixCa_cleaned.drop(columns=[SummaryMatrixCa_cleaned.columns[0], SummaryMatrixCa_cleaned.columns[2], SummaryMatrixCa_cleaned.columns[4]])

                    filenameOut = f'{folder_to_save2}/{List_name}/{NrSubtype}_Tot_FlatPairCaCorr_{List_name}.xlsx'
                    SummaryMatrixCa_cleaned.to_excel(filenameOut, index=True, header=True)  
                
                    #####################################################
                    ## TOTAL Ca correlation with neuron from different population ##
                    #####################################################

                    # Keep only neurons from the list 
                    dfCaTot_filtered={}
                    for sheet_name, dfCa in combined_dfCaTot.items():
                        dfCa=pd.DataFrame(dfCa)
                        indices_to_keep_existing = [idx for idx in listI if idx in dfCa.index] #from first list
                        columns_to_keep_existing = [col for col in listII if col in dfCa.columns] #from second list
                        dfCaTot_filtered[sheet_name] = dfCa.loc[indices_to_keep_existing, columns_to_keep_existing]
                    
                    if DrugExperiment>=1:
                        # Keep only correlation pairs that occurs for each Drug condition
                        for sheet_name, df in dfCaTot_filtered.items(): #remove inactive/non recorded neurons
                            df = df[~(df.fillna(0) == 0).all(axis=1)]
                            df = df.loc[:, ~(df.fillna(0) == 0).all(axis=0)]
                            dfCaTot_filtered[sheet_name] = df    
                        first_key = list(dfCaTot_filtered.keys())[0]
                        common_columns = dfCaTot_filtered[first_key].columns
                        common_indices = dfCaTot_filtered[first_key].index
                        for df in dfCaTot_filtered.values():
                            common_columns = common_columns.intersection(df.columns)
                            common_indices = common_indices.intersection(df.index)
                        dfCaTot_Doublefiltered = {name: df.loc[common_indices, common_columns] for name, df in dfCaTot_filtered.items()}
                    else:
                        dfCaTot_Doublefiltered=dfCaTot_filtered

                    # Save Corr Pairs
                    file_path = f'{folder_to_save2}/{List_name}/{NrSubtype}_Tot_PairCorrCa_{SecondaryList_name}.xlsx'
                    with pd.ExcelWriter(file_path) as writer:
                        for sheet_name, dfCa in dfCaTot_Doublefiltered.items():
                            if 'Z_' not in sheet_name:
                                dfCa = dfCa.sort_index(axis=1)
                                dfCa = dfCa.sort_index(axis=0)
                                dfCa.to_excel(writer, sheet_name=sheet_name, index=True, header=True)

                    # Flat correlations
                    SummaryMatrixCa= pd.DataFrame()
                    for sheet_name, df in dfCaTot_Doublefiltered.items():    
                        series_flattened = df.stack().reset_index()
                        series_flattened['index'] = series_flattened['level_0'] + '_' + series_flattened['level_1']
                        series_flattened = series_flattened.set_index('index')[0]
                        series_flattened_cleaned = series_flattened[series_flattened.index.str.split('_').str[0] != series_flattened.index.str.split('_').str[1]] #remove Neuron1 vs Neuron1
                        series_flattened_cleaned.name = sheet_name
                        SummaryMatrixCa = pd.concat([SummaryMatrixCa,series_flattened_cleaned], axis=1)
                    
                    SummaryMatrixCa_cleaned = SummaryMatrixCa.round(5) # to better detect duplicate                   
                    SummaryMatrixCa_cleaned = SummaryMatrixCa.drop_duplicates(subset=SummaryMatrixCa.columns[1:]) 
                    #SummaryMatrixCa_cleaned = SummaryMatrixCa_cleaned.drop(columns=[SummaryMatrixCa_cleaned.columns[0], SummaryMatrixCa_cleaned.columns[2], SummaryMatrixCa_cleaned.columns[4], SummaryMatrixCa_cleaned.columns[6], SummaryMatrixCa_cleaned.columns[8], SummaryMatrixCa_cleaned.columns[10]]) if DrugExperiment else SummaryMatrixCa_cleaned.drop(columns=[SummaryMatrixCa_cleaned.columns[0], SummaryMatrixCa_cleaned.columns[2], SummaryMatrixCa_cleaned.columns[4]])

                    filenameOut = f'{folder_to_save2}/{List_name}/{NrSubtype}_Tot_FlatPairCaCorr_{SecondaryList_name}.xlsx'
                    SummaryMatrixCa_cleaned.to_excel(filenameOut, index=True, header=True)  
                

                    #####################################################
                    ## Vig St Ca correlation with neuron from same population ##
                    #####################################################

                    # Keep only neurons from the list 
                    dfCa_filtered={}
                    for sheet_name, dfCa in combined_dfCa.items():
                        dfCa=pd.DataFrame(dfCa)
                        indices_to_keep_existing = [idx for idx in listI if idx in dfCa.index] #from first list
                        columns_to_keep_existing = [col for col in listI if col in dfCa.columns] #from second list
                        dfCa_filtered[sheet_name] = dfCa.loc[indices_to_keep_existing, columns_to_keep_existing]

                    if DrugExperiment==1: 
                        # Keep only correlation pairs that occurs for each Drug condition                
                        dfCa_filtered2 = copy.deepcopy(dfCa_filtered)
                        keys_to_keep =['baseline_AW', 'Z_baseline_AW', 'CGP_AW', 'Z_CGP_AW']
                        dfCa_filtered2 = {k: dfCa_filtered2[k] for k in keys_to_keep if k in dfCa_filtered2}
                        for sheet_name, df in dfCa_filtered2.items(): #remove inactive/non recorded neurons
                            df = df[~(df.fillna(0) == 0).all(axis=1)]
                            df = df.loc[:, ~(df.fillna(0) == 0).all(axis=0)]
                            dfCa_filtered2[sheet_name] = df    
                        first_key = list(dfCa_filtered2.keys())[0]
                        common_columns = dfCa_filtered2[first_key].columns
                        common_indices = dfCa_filtered2[first_key].index
                        for df in dfCa_filtered2.values():
                            common_columns = common_columns.intersection(df.columns)
                            common_indices = common_indices.intersection(df.index)
                        dfCa_DoubleFiltered = {name: df.loc[common_indices, common_columns] for name, df in dfCa_filtered2.items()}

                        dfCa_filtered2 = copy.deepcopy(dfCa_filtered)
                        keys_to_keep =['baseline_QW', 'Z_baseline_QW', 'CGP_QW', 'Z_CGP_QW']
                        dfCa_filtered2 = {k: dfCa_filtered2[k] for k in keys_to_keep if k in dfCa_filtered2}
                        for sheet_name, df in dfCa_filtered2.items(): #remove inactive/non recorded neurons
                            df = df[~(df.fillna(0) == 0).all(axis=1)]
                            df = df.loc[:, ~(df.fillna(0) == 0).all(axis=0)]
                            dfCa_filtered2[sheet_name] = df    
                        first_key = list(dfCa_filtered2.keys())[0]
                        common_columns = dfCa_filtered2[first_key].columns
                        common_indices = dfCa_filtered2[first_key].index
                        for df in dfCa_filtered2.values():
                            common_columns = common_columns.intersection(df.columns)
                            common_indices = common_indices.intersection(df.index)
                        dfCa_DoubleFiltered2 = {name: df.loc[common_indices, common_columns] for name, df in dfCa_filtered2.items()}
                        
                        dfCa_filtered2 = copy.deepcopy(dfCa_filtered)
                        keys_to_keep =['baseline_NREM', 'Z_baseline_NREM', 'CGP_NREM', 'Z_CGP_NREM']
                        dfCa_filtered2 = {k: dfCa_filtered2[k] for k in keys_to_keep if k in dfCa_filtered2}
                        for sheet_name, df in dfCa_filtered2.items(): #remove inactive/non recorded neurons
                            df = df[~(df.fillna(0) == 0).all(axis=1)]
                            df = df.loc[:, ~(df.fillna(0) == 0).all(axis=0)]
                            dfCa_filtered2[sheet_name] = df    
                        first_key = list(dfCa_filtered2.keys())[0]
                        common_columns = dfCa_filtered2[first_key].columns
                        common_indices = dfCa_filtered2[first_key].index
                        for df in dfCa_filtered2.values():
                            common_columns = common_columns.intersection(df.columns)
                            common_indices = common_indices.intersection(df.index)
                        dfCa_DoubleFiltered2 = {name: df.loc[common_indices, common_columns] for name, df in dfCa_filtered2.items()}

                        dfCa_filtered2 = copy.deepcopy(dfCa_filtered)
                        keys_to_keep =['baseline_REM', 'Z_baseline_REM', 'CGP_REM', 'Z_CGP_REM']
                        dfCa_filtered2 = {k: dfCa_filtered2[k] for k in keys_to_keep if k in dfCa_filtered2}
                        for sheet_name, df in dfCa_filtered2.items(): #remove inactive/non recorded neurons
                            df = df[~(df.fillna(0) == 0).all(axis=1)]
                            df = df.loc[:, ~(df.fillna(0) == 0).all(axis=0)]
                            dfCa_filtered2[sheet_name] = df    
                        first_key = list(dfCa_filtered2.keys())[0]
                        common_columns = dfCa_filtered2[first_key].columns
                        common_indices = dfCa_filtered2[first_key].index
                        for df in dfCa_filtered2.values():
                            common_columns = common_columns.intersection(df.columns)
                            common_indices = common_indices.intersection(df.index)
                        dfCa_DoubleFiltered3 = {name: df.loc[common_indices, common_columns] for name, df in dfCa_filtered2.items()}
            
                        dfCa_DoubleFiltered.update(dfCa_DoubleFiltered2)    # modifies z with keys and values of y
                        dfCa_DoubleFiltered.update(dfCa_DoubleFiltered3)    # modifies z with keys and values of y

                    elif DrugExperiment==0:
                        # Keep only correlation pairs that occurs for each Vig States  
                        dfCa_filtered2 = copy.deepcopy(dfCa_filtered)
                        keys_to_keep =['baseline_AW', 'Z_baseline_AW','baseline_QW', 'Z_baseline_QW', 'baseline_NREM', 'Z_baseline_NREM', 'baseline_REM', 'Z_baseline_REM']
                        dfCa_filtered2 = {k: dfCa_filtered2[k] for k in keys_to_keep if k in dfCa_filtered2}
                        for sheet_name, df in dfCa_filtered2.items(): #remove inactive/non recorded neurons
                            df = df[~(df.fillna(0) == 0).all(axis=1)]
                            df = df.loc[:, ~(df.fillna(0) == 0).all(axis=0)]
                            dfCa_filtered2[sheet_name] = df    
                        first_key = list(dfCa_filtered2.keys())[0]
                        common_columns = dfCa_filtered2[first_key].columns
                        common_indices = dfCa_filtered2[first_key].index
                        for df in dfCa_filtered2.values():
                            common_columns = common_columns.intersection(df.columns)
                            common_indices = common_indices.intersection(df.index)
                        dfCa_DoubleFiltered = {name: df.loc[common_indices, common_columns] for name, df in dfCa_filtered2.items()}
                                            
                    elif DrugExperiment ==2: 
                        # Keep only correlation pairs that occurs for each Vig States  
                        dfCa_filtered2 = copy.deepcopy(dfCa_filtered)
                        keys_to_keep =['CGP_AW', 'Z_CGP_AW','CGP_QW', 'Z_CGP_QW', 'CGP_NREM', 'Z_CGP_NREM', 'CGP_REM', 'Z_CGP_REM']
                        dfCa_filtered2 = {k: dfCa_filtered2[k] for k in keys_to_keep if k in dfCa_filtered2}
                        for sheet_name, df in dfCa_filtered2.items(): #remove inactive/non recorded neurons
                            df = df[~(df.fillna(0) == 0).all(axis=1)]
                            df = df.loc[:, ~(df.fillna(0) == 0).all(axis=0)]
                            dfCa_filtered2[sheet_name] = df    
                        first_key = list(dfCa_filtered2.keys())[0]
                        common_columns = dfCa_filtered2[first_key].columns
                        common_indices = dfCa_filtered2[first_key].index
                        for df in dfCa_filtered2.values():
                            common_columns = common_columns.intersection(df.columns)
                            common_indices = common_indices.intersection(df.index)
                        dfCa_DoubleFiltered = {name: df.loc[common_indices, common_columns] for name, df in dfCa_filtered2.items()}
                        
                        dfCa_filtered2 = copy.deepcopy(dfCa_filtered)
                        keys_to_keep =['baseline_AW', 'Z_baseline_AW','baseline_QW', 'Z_baseline_QW', 'baseline_NREM', 'Z_baseline_NREM', 'baseline_REM', 'Z_baseline_REM']
                        dfCa_filtered2 = {k: dfCa_filtered2[k] for k in keys_to_keep if k in dfCa_filtered2}       
                        for sheet_name, df in dfCa_filtered2.items(): #remove inactive/non recorded neurons
                            df = df[~(df.fillna(0) == 0).all(axis=1)]
                            df = df.loc[:, ~(df.fillna(0) == 0).all(axis=0)]
                            dfCa_filtered2[sheet_name] = df    
                        first_key = list(dfCa_filtered2.keys())[0]
                        common_columns = dfCa_filtered2[first_key].columns
                        common_indices = dfCa_filtered2[first_key].index
                        for df in dfCa_filtered2.values():
                            common_columns = common_columns.intersection(df.columns)
                            common_indices = common_indices.intersection(df.index)
                        dfCa_DoubleFiltered2 = {name: df.loc[common_indices, common_columns] for name, df in dfCa_filtered2.items()}

                        dfCa_DoubleFiltered.update(dfCa_DoubleFiltered2)    # modifies z with keys and values of y
                    
                    # Save Corr Pairs
                    file_path = f'{folder_to_save2}/{List_name}/{NrSubtype}_VigSt_PairCorrCa_{List_name}.xlsx'
                    with pd.ExcelWriter(file_path) as writer:
                        for sheet_name, dfCa in dfCa_DoubleFiltered.items():
                            if 'Z_' not in sheet_name:                            
                                dfCa = dfCa.sort_index(axis=1)
                                dfCa = dfCa.sort_index(axis=0)
                                dfCa.to_excel(writer, sheet_name=sheet_name, index=True, header=True)

                    # Flat correlations
                    SummaryMatrixCa= pd.DataFrame()
                    for sheet_name, df in dfCa_DoubleFiltered.items():    
                        series_flattened = df.stack().reset_index()
                        series_flattened['index'] = series_flattened['level_0'] + '_' + series_flattened['level_1']
                        series_flattened = series_flattened.set_index('index')[0]
                        series_flattened_cleaned = series_flattened[series_flattened.index.str.split('_').str[0] != series_flattened.index.str.split('_').str[1]] #remove Neuron1 vs Neuron1
                        series_flattened_cleaned.name = sheet_name
                        SummaryMatrixCa = pd.concat([SummaryMatrixCa,series_flattened_cleaned], axis=1)
                    
                    SummaryMatrixCa_cleaned = SummaryMatrixCa.round(5) # to better detect duplicate                   
                    SummaryMatrixCa_cleaned = SummaryMatrixCa.drop_duplicates(subset=SummaryMatrixCa.columns[1:]) 
                    #SummaryMatrixCa_cleaned = SummaryMatrixCa_cleaned.drop(columns=[SummaryMatrixCa_cleaned.columns[0], SummaryMatrixCa_cleaned.columns[2], SummaryMatrixCa_cleaned.columns[4], SummaryMatrixCa_cleaned.columns[6], SummaryMatrixCa_cleaned.columns[8], SummaryMatrixCa_cleaned.columns[10]]) if DrugExperiment else SummaryMatrixCa_cleaned.drop(columns=[SummaryMatrixCa_cleaned.columns[0], SummaryMatrixCa_cleaned.columns[2], SummaryMatrixCa_cleaned.columns[4]])
                    """
                    df_reset = SummaryMatrixCa_cleaned.reset_index()       
                    if len(df_reset)>0:
                        melted_df = pd.melt(df_reset, id_vars=['index'], var_name='VigilanceSt', value_name='CorrCoeff')
                        split_columns = melted_df['VigilanceSt'].str.split('_', expand=True)
                        split_columns.columns = ['Transformation','Drug','Substate']
                        melted_df = pd.concat([melted_df, split_columns], axis=1)
                        extracted_micename = [extract_micename(idx) for idx in melted_df['index']]
                        melted_df['Mice']=extracted_micename            
                    else: 
                        melted_df = pd.DataFrame() 

                    filenameOut = f'{folder_to_save2}/{List_name}/GLM_{NrSubtype}_FlatPairCaCorr_{List_name}.xlsx'
                    melted_df.to_excel(filenameOut, index=True, header=True)
                    """
                    
                    filenameOut = f'{folder_to_save2}/{List_name}/{NrSubtype}_VigSt_FlatPairCaCorr_{List_name}.xlsx'
                    SummaryMatrixCa_cleaned.to_excel(filenameOut, index=True, header=True) 

                    ###########################################################
                    ## Vig St Ca correlation with neurons from different population ##
                    ###########################################################

                    # Keep only neurons from the list 
                    dfCa_filtered={}
                    for sheet_name, dfCa in combined_dfCa.items():
                        dfCa=pd.DataFrame(dfCa)
                        indices_to_keep_existing = [idx for idx in listI if idx in dfCa.index] #from first list
                        columns_to_keep_existing = [col for col in listII if col in dfCa.columns] #from second list
                        dfCa_filtered[sheet_name] = dfCa.loc[indices_to_keep_existing, columns_to_keep_existing]
                        
                    if DrugExperiment==1: 
                        # Keep only correlation pairs that occurs for each Drug condition                     
                        dfCa_filtered2 = copy.deepcopy(dfCa_filtered)
                        keys_to_keep =['baseline_AW', 'Z_baseline_AW', 'CGP_AW', 'Z_CGP_AW']
                        dfCa_filtered2 = {k: dfCa_filtered2[k] for k in keys_to_keep if k in dfCa_filtered2}
                        for sheet_name, df in dfCa_filtered2.items(): #remove inactive/non recorded neurons
                            df = df[~(df.fillna(0) == 0).all(axis=1)]
                            df = df.loc[:, ~(df.fillna(0) == 0).all(axis=0)]
                            dfCa_filtered2[sheet_name] = df    
                        first_key = list(dfCa_filtered2.keys())[0]
                        common_columns = dfCa_filtered2[first_key].columns
                        common_indices = dfCa_filtered2[first_key].index
                        for df in dfCa_filtered2.values():
                            common_columns = common_columns.intersection(df.columns)
                            common_indices = common_indices.intersection(df.index)
                        dfCa_DoubleFiltered = {name: df.loc[common_indices, common_columns] for name, df in dfCa_filtered2.items()}

                        dfCa_filtered2 = copy.deepcopy(dfCa_filtered)
                        keys_to_keep =['baseline_QW', 'Z_baseline_QW', 'CGP_QW', 'Z_CGP_QW']
                        dfCa_filtered2 = {k: dfCa_filtered2[k] for k in keys_to_keep if k in dfCa_filtered2}
                        for sheet_name, df in dfCa_filtered2.items(): #remove inactive/non recorded neurons
                            df = df[~(df.fillna(0) == 0).all(axis=1)]
                            df = df.loc[:, ~(df.fillna(0) == 0).all(axis=0)]
                            dfCa_filtered2[sheet_name] = df    
                        first_key = list(dfCa_filtered2.keys())[0]
                        common_columns = dfCa_filtered2[first_key].columns
                        common_indices = dfCa_filtered2[first_key].index
                        for df in dfCa_filtered2.values():
                            common_columns = common_columns.intersection(df.columns)
                            common_indices = common_indices.intersection(df.index)
                        dfCa_DoubleFiltered2 = {name: df.loc[common_indices, common_columns] for name, df in dfCa_filtered2.items()}

                        dfCa_filtered2 = copy.deepcopy(dfCa_filtered)
                        keys_to_keep =['baseline_NREM', 'Z_baseline_NREM', 'CGP_NREM', 'Z_CGP_NREM']
                        dfCa_filtered2 = {k: dfCa_filtered2[k] for k in keys_to_keep if k in dfCa_filtered2}
                        for sheet_name, df in dfCa_filtered2.items(): #remove inactive/non recorded neurons
                            df = df[~(df.fillna(0) == 0).all(axis=1)]
                            df = df.loc[:, ~(df.fillna(0) == 0).all(axis=0)]
                            dfCa_filtered2[sheet_name] = df    
                        first_key = list(dfCa_filtered2.keys())[0]
                        common_columns = dfCa_filtered2[first_key].columns
                        common_indices = dfCa_filtered2[first_key].index
                        for df in dfCa_filtered2.values():
                            common_columns = common_columns.intersection(df.columns)
                            common_indices = common_indices.intersection(df.index)
                        dfCa_DoubleFiltered2 = {name: df.loc[common_indices, common_columns] for name, df in dfCa_filtered2.items()}
                        
                        dfCa_filtered2 = copy.deepcopy(dfCa_filtered)
                        keys_to_keep =['baseline_REM', 'Z_baseline_REM', 'CGP_REM', 'Z_CGP_REM']
                        dfCa_filtered2 = {k: dfCa_filtered2[k] for k in keys_to_keep if k in dfCa_filtered2}
                        for sheet_name, df in dfCa_filtered2.items(): #remove inactive/non recorded neurons
                            df = df[~(df.fillna(0) == 0).all(axis=1)]
                            df = df.loc[:, ~(df.fillna(0) == 0).all(axis=0)]
                            dfCa_filtered2[sheet_name] = df    
                        first_key = list(dfCa_filtered2.keys())[0]
                        common_columns = dfCa_filtered2[first_key].columns
                        common_indices = dfCa_filtered2[first_key].index
                        for df in dfCa_filtered2.values():
                            common_columns = common_columns.intersection(df.columns)
                            common_indices = common_indices.intersection(df.index)
                        dfCa_DoubleFiltered3 = {name: df.loc[common_indices, common_columns] for name, df in dfCa_filtered2.items()}
            
                        dfCa_DoubleFiltered.update(dfCa_DoubleFiltered2)    # modifies z with keys and values of y
                        dfCa_DoubleFiltered.update(dfCa_DoubleFiltered3)    # modifies z with keys and values of y

                    elif DrugExperiment==0:
                        # Keep only correlation pairs that occurs for each Vig States  
                        dfCa_filtered2 = copy.deepcopy(dfCa_filtered)
                        for key in ['CGP_AW', 'Z_CGP_AW','CGP_QW', 'Z_CGP_QW', 'CGP_NREM', 'Z_CGP_NREM', 'CGP_REM', 'Z_CGP_REM']:
                            if key in dfCa_filtered2:
                                del dfCa_filtered2[key]
                        for sheet_name, df in dfCa_filtered2.items(): #remove inactive/non recorded neurons
                            df = df[~(df.fillna(0) == 0).all(axis=1)]
                            df = df.loc[:, ~(df.fillna(0) == 0).all(axis=0)]
                            dfCa_filtered2[sheet_name] = df    
                        first_key = list(dfCa_filtered2.keys())[0]
                        common_columns = dfCa_filtered2[first_key].columns
                        common_indices = dfCa_filtered2[first_key].index
                        for df in dfCa_filtered2.values():
                            common_columns = common_columns.intersection(df.columns)
                            common_indices = common_indices.intersection(df.index)
                        dfCa_DoubleFiltered = {name: df.loc[common_indices, common_columns] for name, df in dfCa_filtered2.items()}
                    
                    elif DrugExperiment ==2: 
                        # Keep only correlation pairs that occurs for each Vig States  
                        dfCa_filtered2 = copy.deepcopy(dfCa_filtered)
                        keys_to_keep =['CGP_AW', 'Z_CGP_AW','CGP_QW', 'Z_CGP_QW', 'CGP_NREM', 'Z_CGP_NREM', 'CGP_REM', 'Z_CGP_REM']
                        dfCa_filtered2 = {k: dfCa_filtered2[k] for k in keys_to_keep if k in dfCa_filtered2}
                        for sheet_name, df in dfCa_filtered2.items(): #remove inactive/non recorded neurons
                            df = df[~(df.fillna(0) == 0).all(axis=1)]
                            df = df.loc[:, ~(df.fillna(0) == 0).all(axis=0)]
                            dfCa_filtered2[sheet_name] = df    
                        first_key = list(dfCa_filtered2.keys())[0]
                        common_columns = dfCa_filtered2[first_key].columns
                        common_indices = dfCa_filtered2[first_key].index
                        for df in dfCa_filtered2.values():
                            common_columns = common_columns.intersection(df.columns)
                            common_indices = common_indices.intersection(df.index)
                        dfCa_DoubleFiltered = {name: df.loc[common_indices, common_columns] for name, df in dfCa_filtered2.items()}
                        
                        dfCa_filtered2 = copy.deepcopy(dfCa_filtered)
                        keys_to_keep =['baseline_AW', 'Z_baseline_AW','baseline_QW', 'Z_baseline_QW', 'baseline_NREM', 'Z_baseline_NREM', 'baseline_REM', 'Z_baseline_REM']
                        dfCa_filtered2 = {k: dfCa_filtered2[k] for k in keys_to_keep if k in dfCa_filtered2}              
                        for sheet_name, df in dfCa_filtered2.items(): #remove inactive/non recorded neurons
                            df = df[~(df.fillna(0) == 0).all(axis=1)]
                            df = df.loc[:, ~(df.fillna(0) == 0).all(axis=0)]
                            dfCa_filtered2[sheet_name] = df    
                        first_key = list(dfCa_filtered2.keys())[0]
                        common_columns = dfCa_filtered2[first_key].columns
                        common_indices = dfCa_filtered2[first_key].index
                        for df in dfCa_filtered2.values():
                            common_columns = common_columns.intersection(df.columns)
                            common_indices = common_indices.intersection(df.index)
                        dfCa_DoubleFiltered2 = {name: df.loc[common_indices, common_columns] for name, df in dfCa_filtered2.items()}

                        dfCa_DoubleFiltered.update(dfCa_DoubleFiltered2)    # modifies z with keys and values of y
                    
                    # Save Corr Pairs
                    file_path = f'{folder_to_save2}/{List_name}/{NrSubtype}_VigSt_PairCorrCa_{SecondaryList_name}.xlsx'
                    with pd.ExcelWriter(file_path) as writer:
                        for sheet_name, dfCa in dfCa_DoubleFiltered.items():
                            if 'Z_' not in sheet_name:
                                dfCa = dfCa.sort_index(axis=1)
                                dfCa = dfCa.sort_index(axis=0)
                                dfCa.to_excel(writer, sheet_name=sheet_name, index=True, header=True)
                    
                    # Flat correlations
                    SummaryMatrixCa= pd.DataFrame()
                    for sheet_name, df in dfCa_DoubleFiltered.items():    
                        series_flattened = df.stack().reset_index()
                        series_flattened['index'] = series_flattened['level_0'] + '_' + series_flattened['level_1']
                        series_flattened = series_flattened.set_index('index')[0]
                        series_flattened_cleaned = series_flattened[series_flattened.index.str.split('_').str[0] != series_flattened.index.str.split('_').str[1]] #remove Neuron1 vs Neuron1
                        series_flattened_cleaned.name = sheet_name
                        SummaryMatrixCa = pd.concat([SummaryMatrixCa,series_flattened_cleaned], axis=1)
                    
                    SummaryMatrixCa_cleaned = SummaryMatrixCa.round(5) # to better detect duplicate                   
                    SummaryMatrixCa_cleaned = SummaryMatrixCa.drop_duplicates(subset=SummaryMatrixCa.columns[1:]) 
                    #SummaryMatrixCa_cleaned = SummaryMatrixCa_cleaned.drop(columns=[SummaryMatrixCa_cleaned.columns[0], SummaryMatrixCa_cleaned.columns[2], SummaryMatrixCa_cleaned.columns[4], SummaryMatrixCa_cleaned.columns[6], SummaryMatrixCa_cleaned.columns[8], SummaryMatrixCa_cleaned.columns[10]]) if DrugExperiment else SummaryMatrixCa_cleaned.drop(columns=[SummaryMatrixCa_cleaned.columns[0], SummaryMatrixCa_cleaned.columns[2], SummaryMatrixCa_cleaned.columns[4]])

                    """
                    df_reset = SummaryMatrixCa_cleaned.reset_index()       
                    if len(df_reset)>0:
                        melted_df = pd.melt(df_reset, id_vars=['index'], var_name='VigilanceSt', value_name='CorrCoeff')
                        split_columns = melted_df['VigilanceSt'].str.split('_', expand=True)
                        split_columns.columns = ['Transformation','Drug','Substate']
                        melted_df = pd.concat([melted_df, split_columns], axis=1)
                        extracted_micename = [extract_micename(idx) for idx in melted_df['index']]
                        melted_df['Mice']=extracted_micename            
                    else: 
                        melted_df = pd.DataFrame()

                    filenameOut = f'{folder_to_save2}/{List_name}/GLM_{NrSubtype}_FlatPairCaCorr_{SecondaryList_name}.xlsx'
                    melted_df.to_excel(filenameOut, index=True, header=True)  
                    """
                    filenameOut = f'{folder_to_save2}/{List_name}/{NrSubtype}_VigSt_FlatPairCaCorr_{SecondaryList_name}.xlsx'
                    SummaryMatrixCa_cleaned.to_excel(filenameOut, index=True, header=True)   
                    
            if len(filtered_df)>0:

                filenameOut = f'{folder_to_save2}/{List_name}/GLM_{NrSubtype}_VigSt_Global.xlsx'
                filtered_df.to_excel(filenameOut)
                
                #####################    
                # POPULATION COUPLING #
                #####################

                CaPopCoupling_perUnit = filtered_df.pivot_table(index='Unit_ID', columns='Session_ID', values='TotZ_CaPopCoupling', aggfunc='mean')
                CaPopCoupling_perUnit['AllSession']=CaPopCoupling_perUnit.mean(axis=1)  
                # Save CaPopCoupling_perUnit
                filenameOutCaPopCoupling = f'{folder_to_save2}/{List_name}/{NrSubtype}_Tot_ZCaPopCoupling.xlsx'
                writerCaPopCoupling = pd.ExcelWriter(filenameOutCaPopCoupling)
                CaPopCoupling_perUnit.to_excel(writerCaPopCoupling)
                writerCaPopCoupling.close()

                #####################
                # AUC CALCIUM #
                #####################

                resultNormalizedAUC_calcium_perUnit = filtered_df.pivot_table(index='Unit_ID', columns='Substate', values='NormalizedAUC_calcium', aggfunc='mean')   
                try : resultNormalizedAUC_calcium_perUnit = resultNormalizedAUC_calcium_perUnit[desired_order]
                except: pass
                resultNormalizedAUC_calcium_perUnit['Activated_by'] = resultNormalizedAUC_calcium_perUnit.apply(max_column_name, axis=1)
                proportions = resultNormalizedAUC_calcium_perUnit['Activated_by'].value_counts(normalize=True)*100
                
                filenameOut = f'{folder_to_save2}/{List_name}/{NrSubtype}_VigStAct_nAUC_proportions.xlsx'
                proportions.to_excel(filenameOut)

                filenameOutAUC = f'{folder_to_save2}/{List_name}/{NrSubtype}_VigSt_nAUC.xlsx'
                resultNormalizedAUC_calcium_perUnit.to_excel(filenameOutAUC)

                if Drug == 'baseline':
                    BaselineResultNormalizedAUC=resultNormalizedAUC_calcium_perUnit
                    
                resultNormalizedAUC_calcium_perMouse = filtered_df.pivot_table(index='Mice', columns='Substate', values='NormalizedAUC_calcium', aggfunc='mean')
                try: resultNormalizedAUC_calcium_perMouse = resultNormalizedAUC_calcium_perMouse[desired_order]
                except: pass
                filenameOutAUCM = f'{folder_to_save2}/{List_name}/{NrSubtype}_VigSt_AUC_perMouse.xlsx'
                resultNormalizedAUC_calcium_perMouse.to_excel(filenameOutAUCM)

                #####################
                # DECONV ACTIVITY #
                #####################
                
                resultSpikeActivity_perUnit = filtered_df.pivot_table(index='Unit_ID', columns='Substate', values='DeconvSpikeMeanActivity', aggfunc='mean')    
                try: resultSpikeActivity_perUnit = resultSpikeActivity_perUnit[desired_order]
                except: pass
                resultSpikeActivity_perUnit['Activated_by'] = resultSpikeActivity_perUnit.apply(max_column_name, axis=1)
                proportions = resultSpikeActivity_perUnit['Activated_by'].value_counts(normalize=True)*100
                filenameOut = f'{folder_to_save2}/{List_name}/{NrSubtype}_VigStAct_Deconv_proportions.xlsx'
                proportions.to_excel(filenameOut)

                # Save  df
                filenameOut = f'{folder_to_save2}/{List_name}/{NrSubtype}_VigSt_DeconvSpike.xlsx'
                resultSpikeActivity_perUnit.to_excel(filenameOut)

                resultSpikeActivity_perMouse = filtered_df.pivot_table(index='Mice', columns='Substate', values='DeconvSpikeMeanActivity', aggfunc='mean')
                try: resultSpikeActivity_perMouse = resultSpikeActivity_perMouse[desired_order]
                except: pass
                filenameOut = f'{folder_to_save2}/{List_name}/{NrSubtype}_VigSt_DeconvSpike_perMouse.xlsx'
                writer = pd.ExcelWriter(filenameOut)
                resultSpikeActivity_perMouse.to_excel(writer)
                writer.close()
            
                
                #####################
                # Spike ACTIVITY #
                #####################
                
                resultSpikeActivity_perUnit = filtered_df.pivot_table(index='Unit_ID', columns='Substate', values='SpikeActivityHz', aggfunc='mean')    
                try: resultSpikeActivity_perUnit = resultSpikeActivity_perUnit[desired_order]
                except: pass
                resultSpikeActivity_perUnit['Activated_by'] = resultSpikeActivity_perUnit.apply(max_column_name, axis=1)
                proportions = resultSpikeActivity_perUnit['Activated_by'].value_counts(normalize=True)*100
                filenameOut = f'{folder_to_save2}/{List_name}/{NrSubtype}_VigStAct_SpikeHz_proportions.xlsx'
                proportions.to_excel(filenameOut)

                # Save  df
                filenameOut = f'{folder_to_save2}/{List_name}/{NrSubtype}_VigSt_SpikeActivityHz.xlsx'
                resultSpikeActivity_perUnit.to_excel(filenameOut)

                resultSpikeActivity_perMouse = filtered_df.pivot_table(index='Mice', columns='Substate', values='SpikeActivityHz', aggfunc='mean')
                try: resultSpikeActivity_perMouse = resultSpikeActivity_perMouse[desired_order]
                except: pass
                filenameOut = f'{folder_to_save2}/{List_name}/{NrSubtype}_VigSt_SpikeActivityHz_perMouse.xlsx'
                writer = pd.ExcelWriter(filenameOut)
                resultSpikeActivity_perMouse.to_excel(writer)
                writer.close()

            
                #####################    
                # POPULATION COUPLING #
                #####################

                CaPopCoupling_perUnit = filtered_df.pivot_table(index='Unit_ID', columns='Substate', values='Z_CaPopCoupling', aggfunc='mean')
                try: CaPopCoupling_perUnit = CaPopCoupling_perUnit[desired_order]
                except: pass
                # Save CaPopCoupling_perUnit
                filenameOutCaPopCoupling = f'{folder_to_save2}/{List_name}/{NrSubtype}_VigSt_Z_CaPopCoupling.xlsx'
                writerCaPopCoupling = pd.ExcelWriter(filenameOutCaPopCoupling)
                CaPopCoupling_perUnit.to_excel(writerCaPopCoupling)
                writerCaPopCoupling.close()
                
                CaPopCoupling_perUnit = filtered_df.pivot_table(index='Unit_ID', columns=[filtered_df['Substate'], filtered_df['Session_ID']], values='Z_CaPopCoupling', aggfunc='mean')   
                try : CaPopCoupling_perUnit = CaPopCoupling_perUnit[desired_order]
                except: pass

                filenameOutAUC = f'{folder_to_save2}/{List_name}/{NrSubtype}_VarAcrossVigSt_CaPopCoup.xlsx'
                CaPopCoupling_perUnit.to_excel(filenameOutAUC)

        if DrugExperiment:
            resultNormalizedAUC_calcium_perUnit = resultNormalizedAUC_calcium_perUnit.rename(columns={'Wake': 'CGP Wake', 'NREM': 'CGP NREM', 'REM': 'CGP REM', 'Activated_by':'CGP Activated_by', 'RatioNREM_REM':'CGP RatioNREM_REM'})
            mergeRes=pd.concat([BaselineResultNormalizedAUC,resultNormalizedAUC_calcium_perUnit], axis=1)
            filenameOut = f'{folder_to_save}/{NrSubtype}_SelectivityIndex.xlsx'
            mergeRes.to_excel(filenameOut)

#######################
# Propreties VigStates
#######################

filenameOut = f'{folder_to_save}/VigStPropreties.xlsx'
writer = pd.ExcelWriter(filenameOut)

combined_df2 = combined_dfO.drop_duplicates(subset='Substate_ID', keep='first')

DurationVigStates = combined_df2.pivot_table(index='Mice', columns=[combined_df2['Drug'], combined_df2['Substate']], values='DurationSubstate', aggfunc='mean', fill_value=None)
try: DurationVigStates = DurationVigStates[desired_order]
except: pass        
AllDurationVigStates=pd.concat([AllDurationVigStates, DurationVigStates], axis=0)

TotDurationVigStates = combined_df2.pivot_table(index='Mice', columns=[combined_df2['Drug'], combined_df2['Substate']], values='DurationSubstate', aggfunc='sum', fill_value=None)
try: TotDurationVigStates = TotDurationVigStates[desired_order]
except: pass
AllTotDurationVigStates=pd.concat([AllTotDurationVigStates, TotDurationVigStates], axis=0)

AllDurationVigStates.to_excel(writer, sheet_name=f'MeanEpisodeDurations')
AllTotDurationVigStates.to_excel(writer, sheet_name=f'TotalEpisodeDurations')
writer.close()