# Stats on Ca2+ imaging with miniscope and Vigilance States

#######################################################################################
                            # Define Experiment type #
#######################################################################################

AnalysisID='_AB_TEST' #to identify this analysis from another
DrugExperiment=1 

saveexcel=0
Local=1

#choosed_folder='VigSt_2024-07-22_18_21_32_AB_FINAL' if DrugExperiment else 'VigSt_2024-07-22_17_16_28_AB_FINAL'
choosed_folder1='VigSt_2024-08-07_15_22_29_AB_ALL' # for Baseline Expe
choosed_folder2='VigSt_2024-08-07_14_50_31_AB_ALL' # for CGP Expe

desired_order = ['Wake','NREM', 'REM']   
#desired_order = ['Wake', 'N2', 'NREM', 'REM'] 

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
    
def discrimination_index(df):
    #Index = (AUCa-AUCb)/(AUCa+AUCb)
    if 'REM' in df.columns and 'NREM' in df.columns :
        disc_index =(df['NREM']-df['REM'])/(df['NREM']+df['REM'])
        nan_rows = df['NREM'].isna() | df['REM'].isna()
        disc_index[nan_rows] = np.nan
    else:
        disc_index=np.full(len(df), np.nan)
    #min=df[['NREM', 'REM']].min(axis=1).abs()+1
    #NREMadj=df['NREM']#+min
    #REMadj=df['REM']#+min
    #disc_index=NREMadj/REMadj
    return disc_index

def normalize_row(row):
    max_col = row.idxmax()  # Find the column with the maximum value
    max_val = row[max_col]  # Get the maximum value
    return row / max_val   # Normalize the row by dividing by the maximum value

def divide_keys(data):
    for it in range(2, len(data), 3):        
        key1 = list(data.keys())[it-2]
        key2 = list(data.keys())[it-1]
        key3 = list(data.keys())[it]
        d=data[key3]
        data[key3]=d.replace(0, np.nan)
        data[key1] = data[key1] / data[key3]
        data[key2] = data[key2] / data[key3]
    keys_to_delete = list(data.keys())[2::3]
    for key in keys_to_delete:
        del data[key]
    return data    

########################################################################
        # SCRIPT 23AB_GrandAverages&Stats_for_VigilanceStates
########################################################################

# Specify the directory containing the Excel files
InitialDirectory1 = "//10.69.168.1/crnldata/waking/audrey_hay/L1imaging/AnalysedMarch2023/Gaelle/Baseline_recording_ABmodified/AB_Analysis" if Local else "/crnldata/waking/audrey_hay/L1imaging/AnalysedMarch2023/Gaelle/Baseline_recording_ABmodified/AB_Analysis" 
directory1= f'{InitialDirectory1}/{choosed_folder1}'
InitialDirectory2 ="//10.69.168.1/crnldata/waking/audrey_hay/L1imaging/AnalysedMarch2023/Gaelle/CGP/AB_Analysis" if Local else "/crnldata/waking///L1imaging/AnalysedMarch2023/Gaelle/CGP/AB_Analysis"
directory2= f'{InitialDirectory2}/{choosed_folder2}'

# Get the current date and time
FolderNameSave=str(datetime.now())[:19]
FolderNameSave = FolderNameSave.replace(" ", "_").replace(".", "_").replace(":", "_")
destination_folder= f"//10.69.168.1/crnldata/waking/audrey_hay/L1imaging/AnalysedMarch2023/Gaelle/AB_GlobalAnalysis/AVG_VigSt_{FolderNameSave}{AnalysisID}" if Local else f"/crnldata/waking/audrey_hay/L1imaging/AnalysedMarch2023/Gaelle/AB_GlobalAnalysis/AVG_VigSt_{FolderNameSave}{AnalysisID}"
os.makedirs(destination_folder)
folder_to_save=Path(destination_folder)

# Copy the script file to the destination folder
source_script = "C:/Users/Manip2/SCRIPTS/CodePythonAudrey/CodePythonAurelie/HayLabAnalysis/python/23_24_GrandAverage&Stats_for_VigilanceStates_FullAuto.py" if Local else "/HayLabAnalysis/python/23_24.py" 
destination_file_path = f"{destination_folder}/23_24_GrandAverage&Stats_for_VigilanceStates_FullAuto.txt"
shutil.copy(source_script, destination_file_path)

NrSubtypeList=['L1','L2&3']

for NrSubtype in NrSubtypeList:  

    # Initialize an empty df/dict to store the dataframes
    dfs = []
    dfs2_per_sheet = {}
    dfs3_per_sheet = {}

    if NrSubtype=='L1':
        MiceList=['BlackLinesOK', 'BlueLinesOK', 'GreenDotsOK', 'GreenLinesOK', 'RedLinesOK']
    else:
        MiceList=['Purple', 'ThreeColDotsOK', 'ThreeBlueCrossesOK']
    
    nametofind='Global'
    nametofind2='CaCorr'
    nametofind3='SpCorr'

    # Recursively traverse the directory structure
    for directory in [directory1, directory2]:
        for root, _, files in os.walk(directory):
            for filename in files:
                # Check if the file is an Excel file and contains the specified name
                if filename.endswith('.pkl') and nametofind in filename : 
                    if any(name in filename for name in MiceList): 
                        # Construct the full path to the file
                        filepath = os.path.join(root, filename)
                        # Read the file and append it to the list
                        with open(filepath, 'rb') as pickle_file:
                            df = pickle.load(pickle_file)
                        dfs.append(df)
                        print(filename)
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

    # Concatenate all dataframes into a single dataframe
    combined_df = pd.concat(dfs, ignore_index=True)

    # Remove non defined Unique Units 
    combined_df = combined_df[combined_df['Unique_Unit'] != '[]']
    combined_df = combined_df.dropna(subset=['Unique_Unit'])

    combined_df['Unique_Unit'] = combined_df['Unique_Unit'].astype(int).astype(str)
    combined_df['UnitNumber'] = combined_df['UnitNumber'].astype(str)
    combined_df['UnitValue'] = combined_df['UnitValue'].astype(str)

    combined_df['Unit_ID'] = combined_df['Mice'] + combined_df['Unique_Unit']
    
    combined_df['NormalizedAUC_calcium'] = combined_df['AUC_calcium'] / combined_df['DurationSubstate']

    combined_df['Substate_ID'] = combined_df['Mice'] + combined_df['Session'] + combined_df['Substate'] + combined_df['SubstateNumber'].astype(str)
    combined_df['Session_ID'] = combined_df['Mice'] + combined_df['Session'].astype(str)

    # Save big dataset for stats

    filenameOut = f'{folder_to_save}/{NrSubtype}_VigSt_Global.xlsx'
    writer = pd.ExcelWriter(filenameOut)
    combined_df.to_excel(writer)
    writer.close()

    # Save the Ca correlation matrix  

    file_path = f'{folder_to_save}/{NrSubtype}_VigSt_CaCorr.xlsx'
    dfs2_per_sheet = {sheet_name: df.groupby(df.index).sum() for sheet_name, df in dfs2_per_sheet.items()} #cause was concatenated in the 0 axis
    dfs2_per_sheet=divide_keys(dfs2_per_sheet)
    for sheet_name, df in dfs2_per_sheet.items():
        df = df.sort_index(axis=1)
        df = df.sort_index(axis=0)
        dfs2_per_sheet[sheet_name]=df

    if saveexcel:
        with pd.ExcelWriter(file_path) as writer:        
            for sheet_name, df in dfs2_per_sheet.items():
                df.to_excel(writer, sheet_name=sheet_name, index=True, header=True)

    filenameOut = f'{folder_to_save}/{NrSubtype}_VigSt_CaCorr.pkl'
    with open(filenameOut, 'wb') as pickle_file:
        pickle.dump(dfs2_per_sheet, pickle_file)

    # Save the Sp correlation matrix  

    file_path = f'{folder_to_save}/{NrSubtype}_VigSt_SpCorr.xlsx'
    dfs3_per_sheet = {sheet_name: df.groupby(df.index).sum() for sheet_name, df in dfs3_per_sheet.items()}
    dfs3_per_sheet=divide_keys(dfs3_per_sheet)
    for sheet_name, df in dfs3_per_sheet.items():
        df = df.sort_index(axis=1)
        df = df.sort_index(axis=0)
        dfs3_per_sheet[sheet_name]=df

    if saveexcel:    
        with pd.ExcelWriter(file_path) as writer:        
            for sheet_name, df in dfs3_per_sheet.items():
                df.to_excel(writer, sheet_name=sheet_name, index=True, header=True)

    filenameOut = f'{folder_to_save}/{NrSubtype}_VigSt_SpCorr.pkl'
    with open(filenameOut, 'wb') as pickle_file:
        pickle.dump(dfs3_per_sheet, pickle_file)


########################################################################
            # SCRIPT 24AB_Load&Stats_for_VigilanceStates
########################################################################

AllProportionVigStates=pd.DataFrame()
AllDurationVigStates=pd.DataFrame()
AllTotDurationVigStates=pd.DataFrame()

for NrSubtype in NrSubtypeList:

    analysisfile='VigSt_Global'
    combined_dfO = pd.read_excel(f'{folder_to_save}/{NrSubtype}_{analysisfile}.xlsx', index_col=0)
    
    analysisfileCa='VigSt_CaCorr'
    with open(f'{folder_to_save}/{NrSubtype}_{analysisfileCa}.pkl', 'rb') as pickle_file:
        combined_dfCa = pickle.load(pickle_file)

    analysisfileSp='VigSt_SpCorr'
    with open(f'{folder_to_save}/{NrSubtype}_{analysisfileSp}.pkl', 'rb') as pickle_file:
        combined_dfSp = pickle.load(pickle_file)

    ######################
    # CHOOSE OPTIONS
    ######################

    # MINIMUM VIG STATES DURATION #
    combined_df = combined_dfO[combined_dfO['DurationSubstate'] >= 20] 
    
    Drugs= ['Baseline', 'CGP'] if DrugExperiment else ['Baseline']

    # NO LOW FIRING RATE #
    #combined_df = combined_df[combined_df['Avg_SpikeActivityHz'] >= 0.05] 
    
    #####################
    # PREFERENCE #
    #####################
    
    # /!/ The ones from Baseline (not CGP)

    combined_df_Drug=combined_df.copy()
    combined_df_Drug = combined_df_Drug[combined_df_Drug['Drug'] == 'Baseline'] if DrugExperiment else combined_df_Drug

    # the highest mean = the preference OR discrimination factor 
    AllBaselineUnits = combined_df_Drug['Unit_ID'].unique()
        
    AresultActivity_perUnit = combined_df_Drug.pivot_table(index='Unit_ID', columns='Substate', values='NormalizedAUC_calcium', aggfunc='mean')   
    try : AresultActivity_perUnit = AresultActivity_perUnit[desired_order]
    except: pass
    AresultActivity_perUnit['Activated_by'] = AresultActivity_perUnit.apply(max_column_name, axis=1)
    AresultActivity_perUnit['RatioNREM_REM'] =discrimination_index(AresultActivity_perUnit)
    lower_threshold = -0.5 # 5 times more active in REM  #AresultActivity_perUnit['DiscriminationIndex'].quantile(0.20)
    upper_threshold = .5 # 5 times more active in NREM #AresultActivity_perUnit['DiscriminationIndex'].quantile(0.80)
    REMspeunits = AresultActivity_perUnit[AresultActivity_perUnit['RatioNREM_REM'] <= lower_threshold].index
    NREMspeunits = AresultActivity_perUnit[AresultActivity_perUnit['RatioNREM_REM'] >= upper_threshold].index
    NotSpeunits = AresultActivity_perUnit[(AresultActivity_perUnit['RatioNREM_REM'] > lower_threshold) & (AresultActivity_perUnit['RatioNREM_REM'] < upper_threshold)].index
    
    # Only keep units that appears in CGP & Baseline
    if DrugExperiment:
        combined_df_CGP=combined_df.copy()
        combined_df_CGP = combined_df_CGP[combined_df_CGP['Drug'] == 'CGP'] 
        AllCGPUnits = combined_df_CGP['Unit_ID'].unique()
        AllBaselineUnits= np.intersect1d(AllCGPUnits, AllBaselineUnits)
        REMspeunits=np.intersect1d(REMspeunits, AllCGPUnits)
        NREMspeunits=np.intersect1d(NREMspeunits, AllCGPUnits)
        NotSpeunits=np.intersect1d(NotSpeunits, AllCGPUnits)
       
    # Save the List of significant Unit more active in one vigilance state
    if NrSubtype=='L1':
        os.makedirs(f'{folder_to_save}/Baseline/')
    filenameOut = f'{folder_to_save}/Baseline/{NrSubtype}_ActivityPreference.xlsx'
    writer = pd.ExcelWriter(filenameOut)    
    AllBaselineUnitsDF = pd.DataFrame(AllBaselineUnits)
    REMspeunitsDF = pd.DataFrame(REMspeunits)
    NREMspeunitsDF = pd.DataFrame(NREMspeunits)
    NotSpeunitsDF = pd.DataFrame(NotSpeunits)
    AllBaselineUnitsDF.to_excel(writer, sheet_name='AllBaselineUnits', index=True, header=False) 
    REMspeunitsDF.to_excel(writer, sheet_name='REMspe', index=True, header=False) 
    NREMspeunitsDF.to_excel(writer, sheet_name='NREMspe', index=True, header=False) 
    NotSpeunitsDF.to_excel(writer, sheet_name='NotSpe', index=True, header=False) 

    writer.close()

    for Drug in Drugs:

        combined_df_DrugO=combined_dfO.copy()
        combined_df_DrugO = combined_df_DrugO[combined_df_DrugO['Drug'] == Drug]  if DrugExperiment else combined_df_DrugO
        combined_df_Drug=combined_df.copy()
        combined_df_Drug = combined_df_Drug[combined_df_Drug['Drug'] == Drug] if DrugExperiment else combined_df_Drug

        folder_to_save2= f'{folder_to_save}/{Drug}/'
        if NrSubtype=='L1' and Drug=='CGP':
            os.makedirs(folder_to_save2)

        if DrugExperiment: 
            combined_df_CGP=combined_df.copy()
            combined_df_CGP = combined_df_CGP[combined_df_CGP['Drug'] == Drug]
            AllUnits= combined_df_CGP['Unit_ID'].unique()

        List_SignFiringPreference=[NREMspeunits, REMspeunits, NotSpeunits, AllBaselineUnits, AllUnits] if DrugExperiment else [NREMspeunits, REMspeunits, NotSpeunits, AllBaselineUnits] #NREMprefUnits, REMprefUnits, WakeprefUnits] 
        SecondaryList=[REMspeunits, NotSpeunits, NREMspeunits, AllBaselineUnits, AllUnits] if DrugExperiment else [REMspeunits, NotSpeunits, NREMspeunits, AllBaselineUnits] #NREMprefUnits, REMprefUnits, WakeprefUnits] 
        List_Names=['NREMspe','REMspe','NotSpe', 'AllBaseline', 'All'] if DrugExperiment else ['NREMspe','REMspe','NotSpe', 'All'] #'NREMpref', 'REMpref', 'Wakepref' ] 
        SecondaryList_Names=['REMspe','NotSpe','NREMspe', 'AllBaseline', 'All'] if DrugExperiment else ['REMspe','NotSpe','NREMspe', 'All'] #'NREMpref', 'REMpref', 'Wakepref' ] 

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

            if Drug=='Baseline':

                ## Ca correlation with neuron from same population ##
                #####################################################

                # Keep only neurons from the list 
                dfCa_filtered={}
                for sheet_name, dfCa in combined_dfCa.items():
                    dfCa=pd.DataFrame(dfCa)
                    indices_to_keep_existing = [idx for idx in listI if idx in dfCa.index] #from first list
                    columns_to_keep_existing = [col for col in listI if col in dfCa.columns] #from second list
                    dfCa_filtered[sheet_name] = dfCa.loc[indices_to_keep_existing, columns_to_keep_existing]

                if DrugExperiment: 
                    
                    # Keep only correlation pairs that occurs for the 3 vigilances states / the 2 Drugs
                    """
                    for sheet_name, df in dfCa_filtered.items(): #remove inactive/non recorded neurons
                        df = df[~(df.fillna(0) == 0).all(axis=1)]
                        df = df.loc[:, ~(df.fillna(0) == 0).all(axis=0)]
                        dfCa_filtered[sheet_name] = df    
                    first_key = list(dfCa_filtered.keys())[0]
                    common_columns = dfCa_filtered[first_key].columns
                    common_indices = dfCa_filtered[first_key].index
                    for df in dfCa_filtered.values():
                        common_columns = common_columns.intersection(df.columns)
                        common_indices = common_indices.intersection(df.index)
                    dfCa_DoubleFiltered = {name: df.loc[common_indices, common_columns] for name, df in dfCa_filtered.items()}
                    """
                    # Per Drug
                    """
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

                    dfCa_DoubleFiltered = dfCa_DoubleFiltered.copy()   # start with keys and values of x
                    dfCa_DoubleFiltered.update(dfCa_DoubleFiltered2)    # modifies z with keys and values of y
                    """

                    # Per vigilance state
                    
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
                    
                    # All correlations
                    """
                    dfCa_filtered2 = copy.deepcopy(dfCa_filtered)
                    """

                else:
                    # Keep only correlation pairs that occurs for the 3 vigilances states / the 2 Drugs
                    for sheet_name, df in dfCa_filtered.items(): #remove inactive/non recorded neurons
                        df = df[~(df.fillna(0) == 0).all(axis=1)]
                        df = df.loc[:, ~(df.fillna(0) == 0).all(axis=0)]
                        dfCa_filtered[sheet_name] = df    
                    first_key = list(dfCa_filtered.keys())[0]
                    common_columns = dfCa_filtered[first_key].columns
                    common_indices = dfCa_filtered[first_key].index
                    for df in dfCa_filtered.values():
                        common_columns = common_columns.intersection(df.columns)
                        common_indices = common_indices.intersection(df.index)
                    dfCa_DoubleFiltered = {name: df.loc[common_indices, common_columns] for name, df in dfCa_filtered.items()}

                file_path = f'{folder_to_save2}/{List_name}/{NrSubtype}_VigSt_PairCorrCa_{List_name}.xlsx'
                with pd.ExcelWriter(file_path) as writer:
                    for sheet_name, dfCa in dfCa_DoubleFiltered.items():
                        if 'Z_' not in sheet_name:
                            dfCa.to_excel(writer, sheet_name=sheet_name, index=True, header=True)

                SummaryMatrixCa= pd.DataFrame()
                for sheet_name, df in dfCa_DoubleFiltered.items():    
                    series_flattened = df.stack().reset_index()
                    series_flattened['combined_index'] = series_flattened['level_0'] + '_' + series_flattened['level_1']
                    series_flattened = series_flattened.set_index('combined_index')[0]
                    series_flattened_cleaned = series_flattened[series_flattened.index.str.split('_').str[0] != series_flattened.index.str.split('_').str[1]] #remove Neuron1 vs Neuron1
                    SummaryMatrixCa[sheet_name] = series_flattened_cleaned
                
                SummaryMatrixCa = SummaryMatrixCa.round(5) # to better detect duplicate                   
                SummaryMatrixCa_cleaned = SummaryMatrixCa.drop_duplicates(subset=SummaryMatrixCa.columns[1:]) 
                SummaryMatrixCa_cleaned = SummaryMatrixCa_cleaned.drop(columns=[SummaryMatrixCa_cleaned.columns[0], SummaryMatrixCa_cleaned.columns[2], SummaryMatrixCa_cleaned.columns[4], SummaryMatrixCa_cleaned.columns[6], SummaryMatrixCa_cleaned.columns[8], SummaryMatrixCa_cleaned.columns[10]]) if DrugExperiment else SummaryMatrixCa_cleaned.drop(columns=[SummaryMatrixCa_cleaned.columns[0], SummaryMatrixCa_cleaned.columns[2], SummaryMatrixCa_cleaned.columns[4]])

                df_reset = SummaryMatrixCa_cleaned.reset_index()       
                if len(df_reset)>0:
                    melted_df = pd.melt(df_reset, id_vars=['combined_index'], var_name='VigilanceSt', value_name='CorrCoeff')
                    split_columns = melted_df['VigilanceSt'].str.split('_', expand=True)
                    split_columns.columns = ['Transformation','Drug','Substate']
                    melted_df = pd.concat([melted_df, split_columns], axis=1)
                    extracted_micename = [extract_micename(idx) for idx in melted_df['combined_index']]
                    melted_df['Mice']=extracted_micename            
                else: 
                    melted_df = pd.DataFrame()

                filenameOut = f'{folder_to_save2}/{List_name}/{NrSubtype}_VigSt_FlatPairCaCorr_{List_name}.xlsx'
                SummaryMatrixCa_cleaned.to_excel(filenameOut, index=True, header=True)  

                filenameOut = f'{folder_to_save2}/{List_name}/GLM_{NrSubtype}_FlatPairCaCorr_{List_name}.xlsx'
                melted_df.to_excel(filenameOut, index=True, header=True)
                
                ## Ca correlation with neurons from different population ##
                ###########################################################

                # Keep only neurons from the list 
                dfCa_filtered={}
                for sheet_name, dfCa in combined_dfCa.items():
                    dfCa=pd.DataFrame(dfCa)
                    indices_to_keep_existing = [idx for idx in listI if idx in dfCa.index] #from first list
                    columns_to_keep_existing = [col for col in listII if col in dfCa.columns] #from second list
                    dfCa_filtered[sheet_name] = dfCa.loc[indices_to_keep_existing, columns_to_keep_existing]
                    
                if DrugExperiment: 
                    # Keep only correlation pairs that occurs for the 3 vigilances states / the 2 Drugs
                    """
                    for sheet_name, df in dfCa_filtered.items(): #remove inactive/non recorded neurons
                        df = df[~(df.fillna(0) == 0).all(axis=1)]
                        df = df.loc[:, ~(df.fillna(0) == 0).all(axis=0)]
                        dfCa_filtered[sheet_name] = df    
                    first_key = list(dfCa_filtered.keys())[0]
                    common_columns = dfCa_filtered[first_key].columns
                    common_indices = dfCa_filtered[first_key].index
                    for df in dfCa_filtered.values():
                        common_columns = common_columns.intersection(df.columns)
                        common_indices = common_indices.intersection(df.index)
                    dfCa_DoubleFiltered = {name: df.loc[common_indices, common_columns] for name, df in dfCa_filtered.items()}
                    """
                    # Per Drug
                    """
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

                    dfCa_DoubleFiltered = dfCa_DoubleFiltered.copy()   # start with keys and values of x
                    dfCa_DoubleFiltered.update(dfCa_DoubleFiltered2)    # modifies z with keys and values of y
                    """

                    # Per vigilance state
                    
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
                    
                    # All correlations
                    """
                    dfCa_filtered2 = copy.deepcopy(dfCa_filtered)
                    """

                else:
                    # Keep only correlation pairs that occurs for the 3 vigilances states / the 2 Drugs
                    for sheet_name, df in dfCa_filtered.items(): #remove inactive/non recorded neurons
                        df = df[~(df.fillna(0) == 0).all(axis=1)]
                        df = df.loc[:, ~(df.fillna(0) == 0).all(axis=0)]
                        dfCa_filtered[sheet_name] = df    
                    first_key = list(dfCa_filtered.keys())[0]
                    common_columns = dfCa_filtered[first_key].columns
                    common_indices = dfCa_filtered[first_key].index
                    for df in dfCa_filtered.values():
                        common_columns = common_columns.intersection(df.columns)
                        common_indices = common_indices.intersection(df.index)
                    dfCa_DoubleFiltered = {name: df.loc[common_indices, common_columns] for name, df in dfCa_filtered.items()}

                file_path = f'{folder_to_save2}/{List_name}/{NrSubtype}_VigSt_PairCorrCa_{SecondaryList_name}.xlsx'
                with pd.ExcelWriter(file_path) as writer:
                    for sheet_name, dfCa in dfCa_DoubleFiltered.items():
                        if 'Z_' not in sheet_name:
                            dfCa.to_excel(writer, sheet_name=sheet_name, index=True, header=True)

                SummaryMatrixCa= pd.DataFrame()
                for sheet_name, df in dfCa_DoubleFiltered.items():    
                    series_flattened = df.stack().reset_index()
                    series_flattened['combined_index'] = series_flattened['level_0'] + '_' + series_flattened['level_1']
                    series_flattened = series_flattened.set_index('combined_index')[0]
                    series_flattened_cleaned = series_flattened[series_flattened.index.str.split('_').str[0] != series_flattened.index.str.split('_').str[1]] #remove Neuron1 vs Neuron1
                    SummaryMatrixCa[sheet_name] = series_flattened_cleaned
                
                SummaryMatrixCa = SummaryMatrixCa.round(5) # to better detect duplicate                   
                SummaryMatrixCa_cleaned = SummaryMatrixCa.drop_duplicates(subset=SummaryMatrixCa.columns[1:]) 
                SummaryMatrixCa_cleaned = SummaryMatrixCa_cleaned.drop(columns=[SummaryMatrixCa_cleaned.columns[0], SummaryMatrixCa_cleaned.columns[2], SummaryMatrixCa_cleaned.columns[4], SummaryMatrixCa_cleaned.columns[6], SummaryMatrixCa_cleaned.columns[8], SummaryMatrixCa_cleaned.columns[10]]) if DrugExperiment else SummaryMatrixCa_cleaned.drop(columns=[SummaryMatrixCa_cleaned.columns[0], SummaryMatrixCa_cleaned.columns[2], SummaryMatrixCa_cleaned.columns[4]])

                df_reset = SummaryMatrixCa_cleaned.reset_index()       
                if len(df_reset)>0:
                    melted_df = pd.melt(df_reset, id_vars=['combined_index'], var_name='VigilanceSt', value_name='CorrCoeff')
                    split_columns = melted_df['VigilanceSt'].str.split('_', expand=True)
                    split_columns.columns = ['Transformation','Drug','Substate']
                    melted_df = pd.concat([melted_df, split_columns], axis=1)
                    extracted_micename = [extract_micename(idx) for idx in melted_df['combined_index']]
                    melted_df['Mice']=extracted_micename            
                else: 
                    melted_df = pd.DataFrame()

                filenameOut = f'{folder_to_save2}/{List_name}/{NrSubtype}_VigSt_FlatPairCaCorr_{SecondaryList_name}.xlsx'
                SummaryMatrixCa_cleaned.to_excel(filenameOut, index=True, header=True)  

                filenameOut = f'{folder_to_save2}/{List_name}/GLM_{NrSubtype}_FlatPairCaCorr_{SecondaryList_name}.xlsx'
                melted_df.to_excel(filenameOut, index=True, header=True)   
                
            if len(filtered_df)>0:

                filenameOut = f'{folder_to_save2}/{List_name}/GLM_{NrSubtype}_VigSt_Global.xlsx'
                writer = pd.ExcelWriter(filenameOut)
                filtered_df.to_excel(writer)
                writer.close()

                #####################    
                # CORRELATION COEFF #
                #####################
                """
                CaPopCoupling_perUnit = filtered_df.pivot_table(index='Unit_ID', columns='Substate', values='CaPopCoupling', aggfunc='mean')
                try: CaPopCoupling_perUnit = CaPopCoupling_perUnit[desired_order]
                except: pass
                # Save CaPopCoupling_perUnit
                filenameOut = f'{folder_to_save2}/{List_name}/{NrSubtype}_VigSt_CaPopCoupling.xlsx'
                writer = pd.ExcelWriter(filenameOut)
                CaPopCoupling_perUnit.to_excel(writer)
                writer.close()
                """
                CaPopCoupling_perUnit = filtered_df.pivot_table(index='Unit_ID', columns='Substate', values='Z_CaPopCoupling', aggfunc='mean')
                try: CaPopCoupling_perUnit = CaPopCoupling_perUnit[desired_order]
                except: pass
                # Save CaPopCoupling_perUnit
                filenameOut = f'{folder_to_save2}/{List_name}/{NrSubtype}_VigSt_Z_CaPopCoupling.xlsx'
                writer = pd.ExcelWriter(filenameOut)
                CaPopCoupling_perUnit.to_excel(writer)
                writer.close()
                
                #####################
                # AUC CALCIUM #
                #####################

                resultNormalizedAUC_calcium_perUnit = filtered_df.pivot_table(index='Unit_ID', columns='Substate', values='NormalizedAUC_calcium', aggfunc='mean')   
                try : resultNormalizedAUC_calcium_perUnit = resultNormalizedAUC_calcium_perUnit[desired_order]
                except: pass
                resultNormalizedAUC_calcium_perUnit['Activated_by'] = resultNormalizedAUC_calcium_perUnit.apply(max_column_name, axis=1)
                resultNormalizedAUC_calcium_perUnit['RatioNREM_REM'] =discrimination_index(resultNormalizedAUC_calcium_perUnit)

                filenameOut = f'{folder_to_save2}/{List_name}/{NrSubtype}_VigSt_nAUC.xlsx'
                resultNormalizedAUC_calcium_perUnit.to_excel(filenameOut)

                if Drug == 'Baseline':
                    BaselineResultNormalizedAUC=resultNormalizedAUC_calcium_perUnit

                proportions = resultNormalizedAUC_calcium_perUnit['Activated_by'].value_counts(normalize=True)*100

                resultNormalizedAUC_calcium_perMouse = filtered_df.pivot_table(index='Mice', columns='Substate', values='NormalizedAUC_calcium', aggfunc='mean')
                try: resultNormalizedAUC_calcium_perMouse = resultNormalizedAUC_calcium_perMouse[desired_order]
                except: pass
                filenameOut = f'{folder_to_save2}/{List_name}/{NrSubtype}_VigSt_AUC_perMouse.xlsx'
                resultNormalizedAUC_calcium_perMouse.to_excel(filenameOut)

                #####################
                # DECONV ACTIVITY #
                #####################
                """
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
                """

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

    ProportionVigStates = combined_df2.pivot_table(index='Session_ID', columns=[combined_df2['Drug'], combined_df2['Substate']], values='DurationSubstate', aggfunc='sum', fill_value=None)
    try: ProportionVigStates = ProportionVigStates[desired_order]
    except: pass
    ProportionVigStates=ProportionVigStates.div(ProportionVigStates.sum(axis=1), axis=0)*100
    AllProportionVigStates=pd.concat([AllProportionVigStates, ProportionVigStates], axis=0)

    DurationVigStates = combined_df2.pivot_table(index='Session_ID', columns=[combined_df2['Drug'], combined_df2['Substate']], values='DurationSubstate', aggfunc='mean', fill_value=None)
    try: DurationVigStates = DurationVigStates[desired_order]
    except: pass        
    AllDurationVigStates=pd.concat([AllDurationVigStates, DurationVigStates], axis=0)

    TotDurationVigStates = combined_df2.pivot_table(index='Session_ID', columns=[combined_df2['Drug'], combined_df2['Substate']], values='DurationSubstate', aggfunc='sum', fill_value=None)
    try: TotDurationVigStates = TotDurationVigStates[desired_order]
    except: pass
    AllTotDurationVigStates=pd.concat([AllTotDurationVigStates, TotDurationVigStates], axis=0)

AllProportionVigStates.to_excel(writer, sheet_name=f'ProportionVigStates')        
AllDurationVigStates.to_excel(writer, sheet_name=f'MeanEpisodeDurations')
AllTotDurationVigStates.to_excel(writer, sheet_name=f'TotalEpisodeDurations')
writer.close()

