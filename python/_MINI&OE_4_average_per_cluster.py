

# Stats on Ca2+ imaging with miniscope and Vigilance States

#######################################################################################
                            # Define Experiment type #
#######################################################################################

AnalysisID='_CellAssembly' #to identify this analysis from another
DrugExperiment=0 # 0 if Baseline, 1 if CGP, 2 if Baseline & CGP

saveexcel=0
Local=1

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

def extract_micename(index_value):
    match = re.match(r"([a-zA-Z]+)", index_value)
    return match.group(1) if match else index_value

def max_column_name(row):
    return row.idxmax()
    
# Define a function to replace values based on conditions
def replace_values(arr, idx ):
    result = []
    for value in arr:
        if value < -0.5:
            result.append(f'REM{idx}')
        elif -0.5 <= value <= 0.5:
            result.append(f'NonSel{idx}')
        elif 0.5 < value:
            result.append(f'NREM{idx}')
        else:
            result.append(value)  # Keep value if it's outside the defined ranges
    return result

def normalize_row(row):
    max_col = row.idxmax()  # Find the column with the maximum value
    max_val = row[max_col]  # Get the maximum value
    return row / max_val   # Normalize the row by dividing by the maximum value

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

#######################################################################################
                            # Define Directories #
#######################################################################################

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
source_script = "C:/Users/Manip2/SCRIPTS/CodePythonAudrey/CodePythonAurelie/HayLabAnalysis/python/_MINI&OE_4_average_per_cluster.py" if Local else "/python/_MINI&OE_4_average_per_cluster.py" 
destination_file_path = f"{destination_folder}/_MINI&OE_4_average_per_cluster.txt"
shutil.copy(source_script, destination_file_path)


#######################################################################################
                            # Do the average per cluster #
#######################################################################################

# Load Global Table 
with open(f'{directory1}/VigStates_Global.pkl', 'rb') as pickle_file:
    combined_df = pickle.load(pickle_file)


AllProportionVigStates=pd.DataFrame()
AllDurationVigStates=pd.DataFrame()
AllTotDurationVigStates=pd.DataFrame()

NrSubtypeList=['L1NDNF_mice','L2_3_mice']

for NrSubtype in NrSubtypeList:

    analysisfileCa='VigSt_CaCorr'
    with open(f'{folder_to_save}/{NrSubtype}_{analysisfileCa}.pkl', 'rb') as pickle_file:
        combined_dfCa = pickle.load(pickle_file)

    analysisfileSp='VigSt_SpCorr'
    with open(f'{folder_to_save}/{NrSubtype}_{analysisfileSp}.pkl', 'rb') as pickle_file:
        combined_dfSp = pickle.load(pickle_file)   
        
    analysisfileCa='Tot_CaCorr'
    with open(f'{folder_to_save}/{NrSubtype}_{analysisfileCa}.pkl', 'rb') as pickle_file:
        combined_dfCaTot = pickle.load(pickle_file)

    analysisfileSp='Tot_SpCorr'
    with open(f'{folder_to_save}/{NrSubtype}_{analysisfileSp}.pkl', 'rb') as pickle_file:
        combined_dfSpTot = pickle.load(pickle_file)
    
    analysisfileCa='SubSt_CaCorr'
    with open(f'{folder_to_save}/{NrSubtype}_{analysisfileCa}.pkl', 'rb') as pickle_file:
        combined_dfCaSub = pickle.load(pickle_file)


    ######################
    # CHOOSE OPTIONS
    ######################

    # MINIMUM VIG STATES DURATION #
    combined_df = combined_dfO[combined_dfO['DurationSubstate'] >= 0] #10
    
    Drugs= ['baseline', 'CGP'] if DrugExperiment else ['baseline']

    # NO LOW FIRING RATE #
    #combined_df = combined_df[combined_df['Avg_SpikeActivityHz'] >= 0.05] 
    
    #####################
    # LOAD CLUSTERS #
    #####################

    # /!/ The ones from Baseline (not CGP)

    combined_df_Drug=combined_df.copy()
    combined_df_Drug = combined_df_Drug[combined_df_Drug['Drug'] == 'baseline']
    AllBaselineUnits = combined_df_Drug['Unit_ID'].unique()
  
    # Load computed clusters
    dfclusterpath= "//10.69.168.1/crnldata/waking/audrey_hay/L1imaging/Analysed2025_AB/_global_analysis/ClusterAnalysis_dfL1.csv" if NrSubtype=='L1' else "//10.69.168.1/crnldata/waking/audrey_hay/L1imaging/Analysed2025_AB/_global_analysis/ClusterAnalysis_dfL23.csv"
    dfcluster=pd.read_csv(dfclusterpath, index_col=0)
    Cluster0untis = dfcluster[dfcluster['Cluster'] ==0].index
    Cluster1untis = dfcluster[dfcluster['Cluster'] ==1].index
    Cluster2untis = dfcluster[dfcluster['Cluster'] ==2].index
    
    # Only keep units that appears in CGP & Baseline
    
    if DrugExperiment>=1 :
        combined_df_CGP=combined_df.copy()
        combined_df_CGP = combined_df_CGP[combined_df_CGP['Drug'] == 'CGP'] 
        AllCGPUnits = combined_df_CGP['Unit_ID'].unique()

        AllBaselineUnits= np.intersect1d(AllBaselineUnits,AllCGPUnits)
        Cluster0untis= np.intersect1d(Cluster0untis,AllBaselineUnits)
        Cluster1untis= np.intersect1d(Cluster1untis,AllBaselineUnits)
        Cluster2untis= np.intersect1d(Cluster2untis,AllBaselineUnits)
        
    # Save the List of significant Unit more active in one vigilance state
    if NrSubtype=='L1NDNF_mice':
        os.makedirs(f'{folder_to_save}/Baseline/')

    filenameOut = f'{folder_to_save}/Baseline/{NrSubtype}_ActivityPreference.xlsx'
    writer = pd.ExcelWriter(filenameOut)    


    AllBaselineUnitsDF = pd.DataFrame(AllBaselineUnits)
    Cluster0untisDF= pd.DataFrame(Cluster0untis)
    Cluster1untisDF= pd.DataFrame(Cluster1untis)
    Cluster2untisDF= pd.DataFrame(Cluster2untis)
    AllBaselineUnitsDF.to_excel(writer, sheet_name='AllBaselineUnits', index=True, header=False) 
    Cluster0untisDF.to_excel(writer, sheet_name='Cluster0untis', index=True, header=False) 
    Cluster1untisDF.to_excel(writer, sheet_name='Cluster1untis', index=True, header=False) 
    Cluster2untisDF.to_excel(writer, sheet_name='Cluster2untis', index=True, header=False) 
    
    writer.close()

    for Drug in Drugs:

        combined_df_DrugO=combined_dfO.copy() # no min vig states durations == for vig states stats
        combined_df_DrugO = combined_df_DrugO[combined_df_DrugO['Drug'] == Drug] 

        combined_df_Drug=combined_df.copy()
        combined_df_Drug = combined_df_Drug[combined_df_Drug['Drug'] == Drug] 
        AllUnits= combined_df_Drug['Unit_ID'].unique()

        folder_to_save2= f'{folder_to_save}/{Drug}/'
        if NrSubtype=='L1NDNF_mice' and Drug=='CGP':
            os.makedirs(folder_to_save2)


        List_SignFiringPreference=[Cluster0untis, Cluster1untis, Cluster2untis, AllBaselineUnits, AllUnits] if DrugExperiment else [Cluster0untis, Cluster1untis, Cluster2untis, AllUnits]
        SecondaryList=[Cluster1untis, Cluster2untis, Cluster0untis, AllBaselineUnits, AllUnits] if DrugExperiment else [Cluster1untis, Cluster2untis, Cluster0untis, AllUnits]
        List_Names=['Cluster0untis', 'Cluster1untis', 'Cluster2untis', 'AllBaselineUnits', 'All']if DrugExperiment else ['Cluster0untis', 'Cluster1untis','Cluster2untis', 'All']
        SecondaryList_Names=['Cluster1untis', 'Cluster2untis', 'Cluster0untis', 'AllBaselineUnits', 'All' ] if DrugExperiment else ['Cluster1untis', 'Cluster2untis', 'Cluster0untis', 'All']
        
        for listnb, listI  in enumerate(List_SignFiringPreference):
            
            filtered_df = combined_df_Drug[combined_df_Drug['Unit_ID'].isin(listI)]
            List_name=List_Names[listnb]
            SecondaryList_name=SecondaryList_Names[listnb]
            listII=SecondaryList[listnb]

            filtered_df_AllDrug = combined_df[combined_df['Unit_ID'].isin(listI)]
            filenameOut = f'{folder_to_save}/GLM_{NrSubtype}_{List_name}_VigSt_Global.xlsx'
            writer = pd.ExcelWriter(filenameOut)
            filtered_df_AllDrug.to_excel(writer)
            writer.close()

            if NrSubtype=='L1':
                new_folder= f"{folder_to_save2}/{List_name}/"
                os.makedirs(new_folder)

            if Drug==Drug: #'Baseline':
                
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
                    # Keep only correlation pairs that occurs for the 2 Drugs
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
                    # Keep only correlation pairs that occurs for the 2 Drugs
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
                    # Keep only correlation pairs that occurs for each Drug                 
                    dfCa_filtered2 = copy.deepcopy(dfCa_filtered)
                    for key in ['Baseline_NREM', 'Z_Baseline_NREM', 'Baseline_REM', 'Z_Baseline_REM','CGP_NREM', 'Z_CGP_NREM', 'CGP_REM', 'Z_CGP_REM']:
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

                    dfCa_filtered2 = copy.deepcopy(dfCa_filtered)
                    for key in ['Baseline_Wake', 'Z_Baseline_Wake', 'Baseline_REM', 'Z_Baseline_REM','CGP_Wake', 'Z_CGP_Wake', 'CGP_REM', 'Z_CGP_REM']:
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
                    dfCa_DoubleFiltered2 = {name: df.loc[common_indices, common_columns] for name, df in dfCa_filtered2.items()}
                    
                    dfCa_filtered2 = copy.deepcopy(dfCa_filtered)
                    for key in ['Baseline_Wake', 'Z_Baseline_Wake', 'Baseline_NREM', 'Z_Baseline_NREM','CGP_Wake', 'Z_CGP_Wake', 'CGP_NREM', 'Z_CGP_NREM']:
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
                    dfCa_DoubleFiltered3 = {name: df.loc[common_indices, common_columns] for name, df in dfCa_filtered2.items()}
        
                    dfCa_DoubleFiltered.update(dfCa_DoubleFiltered2)    # modifies z with keys and values of y
                    dfCa_DoubleFiltered.update(dfCa_DoubleFiltered3)    # modifies z with keys and values of y

                elif DrugExperiment==0:
                    # Keep only correlation pairs that occurs for each Vig States  
                    dfCa_filtered2 = copy.deepcopy(dfCa_filtered)
                    for key in ['CGP_Wake', 'Z_CGP_Wake', 'CGP_NREM', 'Z_CGP_NREM', 'CGP_REM', 'Z_CGP_REM']:
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
                    for key in ['Baseline_Wake', 'Z_Baseline_Wake', 'Baseline_NREM', 'Z_Baseline_NREM', 'Baseline_REM', 'Z_Baseline_REM']:
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
                    
                    dfCa_filtered2 = copy.deepcopy(dfCa_filtered)
                    for key in ['CGP_Wake', 'Z_CGP_Wake', 'CGP_NREM', 'Z_CGP_NREM', 'CGP_REM', 'Z_CGP_REM']:
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
                    # Keep only correlation pairs that occurs for each Drug                      
                    dfCa_filtered2 = copy.deepcopy(dfCa_filtered)
                    for key in ['Baseline_NREM', 'Z_Baseline_NREM', 'Baseline_REM', 'Z_Baseline_REM','CGP_NREM', 'Z_CGP_NREM', 'CGP_REM', 'Z_CGP_REM']:
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

                    dfCa_filtered2 = copy.deepcopy(dfCa_filtered)
                    for key in ['Baseline_Wake', 'Z_Baseline_Wake', 'Baseline_REM', 'Z_Baseline_REM','CGP_Wake', 'Z_CGP_Wake', 'CGP_REM', 'Z_CGP_REM']:
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
                    dfCa_DoubleFiltered2 = {name: df.loc[common_indices, common_columns] for name, df in dfCa_filtered2.items()}
                    
                    dfCa_filtered2 = copy.deepcopy(dfCa_filtered)
                    for key in ['Baseline_Wake', 'Z_Baseline_Wake', 'Baseline_NREM', 'Z_Baseline_NREM','CGP_Wake', 'Z_CGP_Wake', 'CGP_NREM', 'Z_CGP_NREM']:
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
                    dfCa_DoubleFiltered3 = {name: df.loc[common_indices, common_columns] for name, df in dfCa_filtered2.items()}
        
                    dfCa_DoubleFiltered.update(dfCa_DoubleFiltered2)    # modifies z with keys and values of y
                    dfCa_DoubleFiltered.update(dfCa_DoubleFiltered3)    # modifies z with keys and values of y

                elif DrugExperiment==0:
                    # Keep only correlation pairs that occurs for each Vig States  
                    dfCa_filtered2 = copy.deepcopy(dfCa_filtered)
                    for key in ['CGP_Wake', 'Z_CGP_Wake', 'CGP_NREM', 'Z_CGP_NREM', 'CGP_REM', 'Z_CGP_REM']:
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
                    for key in ['Baseline_Wake', 'Z_Baseline_Wake', 'Baseline_NREM', 'Z_Baseline_NREM', 'Baseline_REM', 'Z_Baseline_REM']:
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
                    
                    dfCa_filtered2 = copy.deepcopy(dfCa_filtered)
                    for key in ['CGP_Wake', 'Z_CGP_Wake', 'CGP_NREM', 'Z_CGP_NREM', 'CGP_REM', 'Z_CGP_REM']:
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

                CaPopCoupling_perUnit = filtered_df.pivot_table(index='Unit_ID', columns='Session', values='TotZ_CaPopCoupling', aggfunc='mean')
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
                resultNormalizedAUC_calcium_perUnit['RatioNREM_REM'] =discrimination_index(resultNormalizedAUC_calcium_perUnit)

                filenameOutAUC = f'{folder_to_save2}/{List_name}/{NrSubtype}_VigSt_nAUC.xlsx'
                resultNormalizedAUC_calcium_perUnit.to_excel(filenameOutAUC)

                if Drug == 'Baseline':
                    BaselineResultNormalizedAUC=resultNormalizedAUC_calcium_perUnit

                proportions = resultNormalizedAUC_calcium_perUnit['Activated_by'].value_counts(normalize=True)*100

                resultNormalizedAUC_calcium_perMouse = filtered_df.pivot_table(index='Mice', columns='Substate', values='NormalizedAUC_calcium', aggfunc='mean')
                try: resultNormalizedAUC_calcium_perMouse = resultNormalizedAUC_calcium_perMouse[desired_order]
                except: pass
                filenameOutAUCM = f'{folder_to_save2}/{List_name}/{NrSubtype}_VigSt_AUC_perMouse.xlsx'
                resultNormalizedAUC_calcium_perMouse.to_excel(filenameOutAUCM)

                #####################
                # SI evolution #
                #####################
                """
                combined_df_Drug['TrSession'] = combined_df_Drug['Session'].str[:8]
                resultNormalizedAUC_calcium_perUnit = combined_df_Drug.pivot_table(index='Unit_ID', columns=[combined_df_Drug['Substate'],combined_df_Drug['TrSession']], values='NormalizedAUC_calcium', aggfunc='mean')   
                A=(resultNormalizedAUC_calcium_perUnit['NREM']-resultNormalizedAUC_calcium_perUnit['REM'])/(resultNormalizedAUC_calcium_perUnit['NREM']+resultNormalizedAUC_calcium_perUnit['REM'])
                A = A[A.notna().sum(axis=1)>1]
                A['first_non_nan'] = A.bfill(axis=1).iloc[:, 0]
                A['last_non_nan'] = A.ffill(axis=1).iloc[:, -2]

                # Apply the function
                A['first_non_nan_Names'] = replace_values(A['first_non_nan'], '_1st')
                A['last_non_nan_Names'] = replace_values(A['last_non_nan'], '_last')

                dftest=A['first_non_nan_Names'] + A['first_non_nan_Names']

                
                filenameOutSI = f'{folder_to_save2}/{List_name}/{NrSubtype}_SIevolution.xlsx'
                A.to_excel(filenameOutSI)
                """
                #####################
                # DECONV ACTIVITY #
                #####################
                
                resultSpikeActivity_perUnit = filtered_df.pivot_table(index='Unit_ID', columns='Substate', values='DeconvSpikeMeanActivity', aggfunc='mean')    
                try: resultSpikeActivity_perUnit = resultSpikeActivity_perUnit[desired_order]
                except: pass
                resultSpikeActivity_perUnit['Activated_by'] = resultSpikeActivity_perUnit.apply(max_column_name, axis=1)

                # Save  df
                filenameOut = f'{folder_to_save2}/{List_name}/{NrSubtype}_VigSt_DeconvSpike.xlsx'
                writer = pd.ExcelWriter(filenameOut)
                resultSpikeActivity_perUnit.to_excel(writer)
                writer.close()

                resultSpikeActivity_perMouse = filtered_df.pivot_table(index='Mice', columns='Substate', values='DeconvSpikeMeanActivity', aggfunc='mean')
                try: resultSpikeActivity_perMouse = resultSpikeActivity_perMouse[desired_order]
                except: pass
                filenameOut = f'{folder_to_save2}/{List_name}/{NrSubtype}_VigSt_DeconvSpike_perMouse.xlsx'
                writer = pd.ExcelWriter(filenameOut)
                resultSpikeActivity_perMouse.to_excel(writer)
                writer.close()

                
                proportions = resultSpikeActivity_perUnit['Activated_by'].value_counts(normalize=True)*100
                filenameOut = f'{folder_to_save2}/{List_name}/{NrSubtype}_ActivityVigSt_DeconvSpike_proportions.xlsx'
                writer = pd.ExcelWriter(filenameOut)
                proportions.to_excel(writer)
                writer.close()
                
                #####################
                # Spike ACTIVITY #
                #####################
                
                resultSpikeActivity_perUnit = filtered_df.pivot_table(index='Unit_ID', columns='Substate', values='SpikeActivityHz', aggfunc='mean')    
                try: resultSpikeActivity_perUnit = resultSpikeActivity_perUnit[desired_order]
                except: pass
                resultSpikeActivity_perUnit['Activated_by'] = resultSpikeActivity_perUnit.apply(max_column_name, axis=1)

                # Save  df
                filenameOut = f'{folder_to_save2}/{List_name}/{NrSubtype}_VigSt_SpikeActivityHz.xlsx'
                writer = pd.ExcelWriter(filenameOut)
                resultSpikeActivity_perUnit.to_excel(writer)
                writer.close()

                resultSpikeActivity_perMouse = filtered_df.pivot_table(index='Mice', columns='Substate', values='SpikeActivityHz', aggfunc='mean')
                try: resultSpikeActivity_perMouse = resultSpikeActivity_perMouse[desired_order]
                except: pass
                filenameOut = f'{folder_to_save2}/{List_name}/{NrSubtype}_VigSt_SpikeActivityHz_perMouse.xlsx'
                writer = pd.ExcelWriter(filenameOut)
                resultSpikeActivity_perMouse.to_excel(writer)
                writer.close()

                
                proportions = resultSpikeActivity_perUnit['Activated_by'].value_counts(normalize=True)*100
                filenameOut = f'{folder_to_save2}/{List_name}/{NrSubtype}_ActivityVigSt_SpikeActivityHz_proportions.xlsx'
                writer = pd.ExcelWriter(filenameOut)
                proportions.to_excel(writer)
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
                
                
                CaPopCoupling_perUnit = filtered_df.pivot_table(index='Unit_ID', columns=[filtered_df['Substate'], filtered_df['Session']], values='Z_CaPopCoupling', aggfunc='mean')   
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