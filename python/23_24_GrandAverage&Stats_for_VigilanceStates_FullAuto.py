# Stats on Ca2+ imaging with miniscope and Vigilance States

#######################################################################################
                            # Define Experiment type #
#######################################################################################

#DrugExperiment=1 if CGP Experiment // DrugExperiment=0 if Baseline Experiment
DrugExperiment=1

AnalysisID='_AB_FINAL' #to identify this analysis from another

choosed_folder='VigSt_2024-07-22_18_21_32_AB_FINAL'

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

########################################################################
        # SCRIPT 23AB_GrandAverages&Stats_for_VigilanceStates
########################################################################

# Specify the directory containing the Excel files
InitialDirectory = "//10.69.168.1/crnldata/waking/audrey_hay/L1imaging/AnalysedMarch2023/Gaelle/CGP/AB_Analysis" if DrugExperiment else f"//10.69.168.1/crnldata/waking/audrey_hay/L1imaging/AnalysedMarch2023/Gaelle/Baseline_recording_ABmodified/AB_Analysis"
directory= f'{InitialDirectory}/{choosed_folder}'

# Get the current date and time
FolderNameSave=str(datetime.now())[:19]
FolderNameSave = FolderNameSave.replace(" ", "_").replace(".", "_").replace(":", "_")
destination_folder= f"//10.69.168.1/crnldata/waking/audrey_hay/L1imaging/AnalysedMarch2023/Gaelle/CGP/AB_Analysis/AVG_VigSt_{FolderNameSave}{AnalysisID}" if DrugExperiment else f"//10.69.168.1/crnldata/waking/audrey_hay/L1imaging/AnalysedMarch2023/Gaelle/Baseline_recording_ABmodified/AB_Analysis/AVG_VigSt_{FolderNameSave}{AnalysisID}"
os.makedirs(destination_folder)
folder_to_save=Path(destination_folder)

# Copy the script file to the destination folder
source_script = "C:/Users/Manip2/SCRIPTS/CodePythonAudrey/CodePythonAurelie/HayLabAnalysis/python/23_24_GrandAverage&Stats_for_VigilanceStates_FullAuto.py"
destination_file_path = f"{destination_folder}/23_24_GrandAverage&Stats_for_VigilanceStates_FullAuto.txt"
shutil.copy(source_script, destination_file_path)

NrSubtypeList=['All', 'L1']

for NrSubtype in NrSubtypeList:  

    # Initialize an empty list to store the dataframes
    dfs = []
    df= []
    df2=[]
    dfs2_per_sheet = {}
    dfs3 = []
    df3=[]
    dfs3_per_sheet = {}

    if NrSubtype=='L1':
        MiceList=['BlackLinesOK', 'BlueLinesOK', 'GreenDotsOK', 'GreenLinesOK', 'RedLinesOK']
    else:
        MiceList=['Purple', 'ThreeColDotsOK', 'ThreeBlueCrossesOK']
    
    nametofind='Global'
    nametofind2='CaCorr'
    nametofind3='SpCorr'

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
            if filename.endswith('.xlsx') and nametofind2 in filename: 
                if any(name in filename for name in MiceList): 
                    # Construct the full path to the file
                    filepath = os.path.join(root, filename)
                    # Read the Excel file into a dataframe and append it to the list
                    excel_data = pd.read_excel(filepath, sheet_name=None, index_col=0)           
                    for sheet_name, df2 in excel_data.items():
                        if len(df2)>0:
                            if sheet_name in dfs2_per_sheet:                                       
                                updated_matrix = pd.concat([dfs2_per_sheet[sheet_name], df2], axis=0)                    
                                dfs2_per_sheet[sheet_name] = updated_matrix    
                            else:
                                dfs2_per_sheet[sheet_name] = df2 
                    # Keep only correlation that occurs for the 3 vigilances states / the 2 Drugs
                    first_key = next(iter(dfs2_per_sheet))
                    common_columns = dfs2_per_sheet[first_key].columns
                    for df in dfs2_per_sheet.values():
                        common_columns = common_columns.intersection(df.columns)

                    common_indices = dfs2_per_sheet[first_key].index
                    for df in dfs2_per_sheet.values():
                        common_indices = common_indices.intersection(df.index)

                    filtered_df_dict2 = {name: df.loc[common_indices, common_columns] for name, df in dfs2_per_sheet.items()}
            if filename.endswith('.xlsx') and nametofind3 in filename: 
                if any(name in filename for name in MiceList): 
                    # Construct the full path to the file
                    filepath = os.path.join(root, filename)
                    # Read the Excel file into a dataframe and append it to the list
                    excel_data = pd.read_excel(filepath, sheet_name=None, index_col=0)           
                    for sheet_name, df3 in excel_data.items():
                        if len(df3)>0:
                            if sheet_name in dfs3_per_sheet:   
                                updated_matrix = pd.concat((dfs3_per_sheet[sheet_name],df3), axis=0)                
                                dfs3_per_sheet[sheet_name] = updated_matrix                    
                            else:                    
                                dfs3_per_sheet[sheet_name] = df3 #one average trace per unique unit, len(df3)==nb unit recorded for that mouse
                    # Keep only correlation that occurs for the 3 vigilances states / the 2 Drugs
                    first_key = next(iter(dfs3_per_sheet))
                    common_columns = dfs3_per_sheet[first_key].columns
                    for df in dfs3_per_sheet.values():
                        common_columns = common_columns.intersection(df.columns)

                    common_indices = dfs3_per_sheet[first_key].index
                    for df in dfs3_per_sheet.values():
                        common_indices = common_indices.intersection(df.index)

                    filtered_df_dict3 = {name: df.loc[common_indices, common_columns] for name, df in dfs3_per_sheet.items()}

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
    with pd.ExcelWriter(file_path) as writer:
        for sheet_name, df in dfs2_per_sheet.items():
            df.to_excel(writer, sheet_name=sheet_name, index=True, header=True)

    # Save the Sp correlation matrix  

    file_path = f'{folder_to_save}/{NrSubtype}_VigSt_SpCorr.xlsx'
    with pd.ExcelWriter(file_path) as writer:
        for sheet_name, df in dfs3_per_sheet.items():
            df.to_excel(writer, sheet_name=sheet_name, index=True, header=True)

########################################################################
            # SCRIPT 24AB_Load&Stats_for_VigilanceStates
########################################################################

for NrSubtype in NrSubtypeList:

    analysisfile='VigSt_Global'
    combined_dfO = pd.read_excel(f'{folder_to_save}/{NrSubtype}_{analysisfile}.xlsx', index_col=0)
    
    combined_dfCa = {}
    analysisfileCa='VigSt_CaCorr'
    excel_data = pd.read_excel(f'{folder_to_save}/{NrSubtype}_{analysisfileCa}.xlsx', sheet_name=None, index_col=0)           
    for sheet_name, df in excel_data.items():
        if len(df)>0:
            combined_dfCa[sheet_name] = df  

    combined_dfSp = {}
    analysisfileSp='VigSt_SpCorr'
    excel_data = pd.read_excel(f'{folder_to_save}/{NrSubtype}_{analysisfileSp}.xlsx', sheet_name=None, index_col=0)           
    for sheet_name, df in excel_data.items():
        if len(df)>0:
            combined_dfSp[sheet_name] = df  

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
    
    # /!\ The ones from Baseline (not CGP)

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
    
    """
    index_lists = {}
    for unitindex in AresultActivity_perUnit['Activated_by'].unique():
        indices = AresultActivity_perUnit.index[AresultActivity_perUnit['Activated_by'] == unitindex].tolist()
        index_lists[unitindex] = indices    
    NREMprefUnits=index_lists['NREM']
    try: REMprefUnits=index_lists['REM'] 
    except: REMprefUnits=[]
    WakeprefUnits=index_lists['Wake']
    """
    
    # Save the List of significant Unit more active in one vigilance state
    if NrSubtype=='All':
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
    """
    NREMprefUnitsDF = pd.DataFrame(NREMprefUnits)
    REMprefUnitsDF = pd.DataFrame(REMprefUnits)
    WakeprefUnitsDF = pd.DataFrame(WakeprefUnits)
    NREMprefUnitsDF.to_excel(writer, sheet_name='NREMpref', index=True, header=False) 
    REMprefUnitsDF.to_excel(writer, sheet_name='REMpref', index=True, header=False) 
    WakeprefUnitsDF.to_excel(writer, sheet_name='Wakepref', index=True, header=False)     
    """
    writer.close()

    for Drug in Drugs:

        combined_df_DrugO=combined_dfO.copy()
        combined_df_DrugO = combined_df_DrugO[combined_df_DrugO['Drug'] == Drug] if DrugExperiment else combined_df_DrugO
        combined_df_Drug=combined_df.copy()
        combined_df_Drug = combined_df_Drug[combined_df_Drug['Drug'] == Drug] if DrugExperiment else combined_df_Drug

        folder_to_save2= f'{folder_to_save}/{Drug}/'
        if NrSubtype=='All' and Drug=='CGP':
            os.makedirs(folder_to_save2)

        if DrugExperiment: 
            combined_df_CGP=combined_df.copy()
            combined_df_CGP = combined_df_CGP[combined_df_CGP['Drug'] == Drug]
            AllUnits= combined_df_CGP['Unit_ID'].unique()

        List_SignFiringPreference=[NREMspeunits, REMspeunits, NotSpeunits, AllBaselineUnits, AllUnits] if DrugExperiment else [NREMspeunits, REMspeunits, NotSpeunits, AllBaselineUnits] #NREMprefUnits, REMprefUnits, WakeprefUnits] 
        List_Names=['NREMspeUnits','REMspeUnits','NotSpeUnits', 'AllBaselineUnits', 'AllUnits'] if DrugExperiment else ['NREMspeUnits','REMspeUnits','NotSpeUnits', 'AllUnits'] #'NREMprefUnits', 'REMprefUnits', 'WakeprefUnits' ] 

        for listnb, listI  in enumerate(List_SignFiringPreference):
            
            filtered_df = combined_df_Drug[combined_df_Drug['Unit_ID'].isin(listI)]
            List_name=List_Names[listnb]

            filtered_df_AllDrug = combined_df[combined_df['Unit_ID'].isin(listI)]
            filenameOut = f'{folder_to_save}/GLM_{NrSubtype}_{List_name}_VigSt_Global.xlsx'
            writer = pd.ExcelWriter(filenameOut)
            filtered_df_AllDrug.to_excel(writer)
            writer.close()

            if NrSubtype=='All':
                new_folder= f"{folder_to_save2}/{List_name}/"
                os.makedirs(new_folder)

            if Drug=='Baseline':

                ## Ca correlation

                # Keep only neurons from the list 
                dfCa_filtered={}
                for sheet_name, dfCa in combined_dfCa.items():
                    dfCa=pd.DataFrame(dfCa)
                    indices_to_keep_existing = [idx for idx in listI if idx in dfCa.index]
                    columns_to_keep_existing = [col for col in listI if col in dfCa.columns]
                    dfCa_filtered[sheet_name] = dfCa.loc[indices_to_keep_existing, columns_to_keep_existing]

                # Keep only correlation that occurs for the 3 vigilances states / the 2 Drugs
                first_key = next(iter(dfCa_filtered))
                common_columns = dfCa_filtered[first_key].columns
                for df in dfCa_filtered.values():
                    common_columns = common_columns.intersection(df.columns)
                common_indices = dfCa_filtered[first_key].index
                for df in dfCa_filtered.values():
                    common_indices = common_indices.intersection(df.index)
                dfCa_DoubleFiltered = {name: df.loc[common_indices, common_columns] for name, df in dfCa_filtered.items()}

                file_path = f'{folder_to_save2}/{List_name}/{NrSubtype}_VigSt_Pairwise_CaCorr.xlsx'
                with pd.ExcelWriter(file_path) as writer:
                    for sheet_name, dfCa in dfCa_DoubleFiltered.items():
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
                                
                df_reset = SummaryMatrixCa_cleaned.reset_index()       
                columns_to_keep = df_reset.columns[[0, 2, 4, 6, 8, 10, 12]] if DrugExperiment else df_reset.columns[[0, 2, 4, 6]] # Columns 2, 4, and 6 (0-based index)
                df_reset = df_reset[columns_to_keep]
                if len(df_reset)>0:
                    melted_df = pd.melt(df_reset, id_vars=['combined_index'], var_name='VigilanceSt', value_name='CorrCoeff')
                    split_columns = melted_df['VigilanceSt'].str.split('_', expand=True)
                    split_columns.columns = ['Transformation','Drug','Substate']
                    melted_df = pd.concat([melted_df, split_columns], axis=1)
                    extracted_micename = [extract_micename(idx) for idx in melted_df['combined_index']]
                    melted_df['Mice']=extracted_micename            
                else: 
                    melted_df = pd.DataFrame()

                filenameOut = f'{folder_to_save2}/{List_name}/{NrSubtype}_VigSt_FlatPairwise_CaCorr.xlsx'
                SummaryMatrixCa_cleaned.to_excel(filenameOut, index=True, header=True)  

                filenameOut = f'{folder_to_save2}/{List_name}/GLM_{NrSubtype}_FlatPairwise_CaCorr.xlsx'
                melted_df.to_excel(filenameOut, index=True, header=True)   
                
                ## Sp correlation

                # Keep only neurons from the list 
                dfSp_filtered={}
                for sheet_name, dfSp in combined_dfSp.items():
                    dfSp=pd.DataFrame(dfSp)
                    indices_to_keep_existing = [idx for idx in listI if idx in dfSp.index]
                    columns_to_keep_existing = [col for col in listI if col in dfSp.columns]
                    dfSp_filtered[sheet_name] = dfSp.loc[indices_to_keep_existing, columns_to_keep_existing]

                # Keep only correlation that occurs for the 3 vigilances states / the 2 Drugs
                first_key = next(iter(dfSp_filtered))
                common_columns = dfSp_filtered[first_key].columns
                for df in dfSp_filtered.values():
                    common_columns = common_columns.intersection(df.columns)
                common_indices = dfSp_filtered[first_key].index
                for df in dfSp_filtered.values():
                    common_indices = common_indices.intersection(df.index)
                dfSp_DoubleFiltered = {name: df.loc[common_indices, common_columns] for name, df in dfSp_filtered.items()}

                file_path = f'{folder_to_save2}/{List_name}/{NrSubtype}_VigSt_Pairwise_SpCorr.xlsx'
                with pd.ExcelWriter(file_path) as writer:
                    for sheet_name, dfSp in dfSp_DoubleFiltered.items():
                        dfSp.to_excel(writer, sheet_name=sheet_name, index=True, header=True)

                SummaryMatrixSp= pd.DataFrame()
                for sheet_name, df in dfSp_DoubleFiltered.items():    
                    series_flattened = df.stack().reset_index()
                    series_flattened['combined_index'] = series_flattened['level_0'] + '_' + series_flattened['level_1']
                    series_flattened = series_flattened.set_index('combined_index')[0]
                    series_flattened_cleaned = series_flattened[series_flattened.index.str.split('_').str[0] != series_flattened.index.str.split('_').str[1]] #remove Neuron1 vs Neuron1
                    SummaryMatrixSp[sheet_name] = series_flattened_cleaned

                SummaryMatrixSp = SummaryMatrixSp.round(5) # to better detect duplicate                   
                SummaryMatrixSp_cleaned = SummaryMatrixSp.drop_duplicates(subset=SummaryMatrixSp.columns[1:])
        
                df_reset = SummaryMatrixSp_cleaned.reset_index()       
                columns_to_keep = df_reset.columns[[0, 2, 4, 6, 8, 10, 12]] if DrugExperiment else df_reset.columns[[0, 2, 4, 6]] # Columns 2, 4, and 6 (0-based index)
                df_reset = df_reset[columns_to_keep]
                if len(df_reset)>0:
                    melted_df = pd.melt(df_reset, id_vars=['combined_index'], var_name='VigilanceSt', value_name='CorrCoeff')
                    split_columns = melted_df['VigilanceSt'].str.split('_', expand=True)
                    split_columns.columns = ['Transformation','Drug','Substate']
                    extracted_micename = [extract_micename(idx) for idx in melted_df['combined_index']]
                    melted_df['Mice']=extracted_micename            
                else: 
                    melted_df = pd.DataFrame()
                    
                filenameOut = f'{folder_to_save2}/{List_name}/{NrSubtype}_VigSt_FlatPairwise_SpCorr.xlsx'
                SummaryMatrixSp_cleaned.to_excel(filenameOut, index=True, header=True)   
                filenameOutGLM = f'{folder_to_save2}/{List_name}/GLM_{NrSubtype}_FlatPairwise_SpCorr.xlsx'
                melted_df.to_excel(filenameOutGLM, index=True, header=True)   

            if len(filtered_df)>0:

                filenameOut = f'{folder_to_save2}/{List_name}/GLM_{NrSubtype}_VigSt_Global.xlsx'
                writer = pd.ExcelWriter(filenameOut)
                filtered_df.to_excel(writer)
                writer.close()

                #######################
                # DETERMINATION COEFF #
                #######################

                CorrCoeff_calcium_perUnit = filtered_df.pivot_table(index='Unit_ID', columns='Substate', values='Rsquared_CaCorrCoeff', aggfunc='mean')
                try: CorrCoeff_calcium_perUnit = CorrCoeff_calcium_perUnit[desired_order]
                except: pass

                # Save CorrCoeff_calcium_perUnit
                filenameOut = f'{folder_to_save2}/{List_name}/{NrSubtype}_VigSt_Rsq_CaCorrCoeff.xlsx'
                writer = pd.ExcelWriter(filenameOut)
                CorrCoeff_calcium_perUnit.to_excel(writer)
                writer.close()

                CorrCoeff_spike_perUnit = filtered_df.pivot_table(index='Unit_ID', columns='Substate', values='Rsquared_SpCorrCoeff', aggfunc='mean')
                try: CorrCoeff_spike_perUnit = CorrCoeff_spike_perUnit[desired_order]
                except: pass
                # Save CorrCoeff_spike_perUnit 
                filenameOut = f'{folder_to_save2}/{List_name}/{NrSubtype}_VigSt_Rsq_SpCorrCoeff.xlsx'
                writer = pd.ExcelWriter(filenameOut)
                CorrCoeff_spike_perUnit.to_excel(writer)
                writer.close()

                #####################    
                # CORRELATION COEFF #
                #####################

                CorrCoeff_calcium_perUnit = filtered_df.pivot_table(index='Unit_ID', columns='Substate', values='CaCorrCoeff', aggfunc='mean')
                try: CorrCoeff_calcium_perUnit = CorrCoeff_calcium_perUnit[desired_order]
                except: pass
                # Save CorrCoeff_calcium_perUnit
                filenameOut = f'{folder_to_save2}/{List_name}/{NrSubtype}_VigSt_CaCorrCoeff.xlsx'
                writer = pd.ExcelWriter(filenameOut)
                CorrCoeff_calcium_perUnit.to_excel(writer)
                writer.close()

                CorrCoeff_spike_perUnit = filtered_df.pivot_table(index='Unit_ID', columns='Substate', values='SpCorrCoeff', aggfunc='mean')
                try: CorrCoeff_spike_perUnit = CorrCoeff_spike_perUnit[desired_order]
                except: pass
                # Save CorrCoeff_spike_perUnit 
                filenameOut = f'{folder_to_save2}/{List_name}/{NrSubtype}_VigSt_SpCorrCoeff.xlsx'
                writer = pd.ExcelWriter(filenameOut)
                CorrCoeff_spike_perUnit.to_excel(writer)
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

                """
                SpPopCoupling_perUnit = filtered_df.pivot_table(index='Unit_ID', columns='Substate', values='SpPopCoupling', aggfunc='mean')
                try: SpPopCoupling_perUnit = SpPopCoupling_perUnit[desired_order]
                except: pass
                # Save SpPopCoupling_perUnit 
                filenameOut = f'{folder_to_save2}/{List_name}/{NrSubtype}_VigSt_SpPopCoupling.xlsx'
                writer = pd.ExcelWriter(filenameOut)
                SpPopCoupling_perUnit.to_excel(writer)
                writer.close()
                """

                SpPopCoupling_perUnit = filtered_df.pivot_table(index='Unit_ID', columns='Substate', values='Z_SpPopCoupling', aggfunc='mean')
                try: SpPopCoupling_perUnit = SpPopCoupling_perUnit[desired_order]
                except: pass
                # Save SpPopCoupling_perUnit 
                filenameOut = f'{folder_to_save2}/{List_name}/{NrSubtype}_VigSt_Z_SpPopCoupling.xlsx'
                writer = pd.ExcelWriter(filenameOut)
                SpPopCoupling_perUnit.to_excel(writer)
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
                # SPIKE ACTIVITY #
                #####################

                resultSpikeActivity_perUnit = filtered_df.pivot_table(index='Unit_ID', columns='Substate', values='SpikeActivityHz', aggfunc='mean')    
                try: resultSpikeActivity_perUnit = resultSpikeActivity_perUnit[desired_order]
                except: pass
                resultSpikeActivity_perUnit['Activated_by'] = resultSpikeActivity_perUnit.apply(max_column_name, axis=1)

                # Save  df
                filenameOut = f'{folder_to_save2}/{List_name}/{NrSubtype}_VigSt_Spike.xlsx'
                writer = pd.ExcelWriter(filenameOut)
                resultSpikeActivity_perUnit.to_excel(writer)
                writer.close()

                resultSpikeActivity_perMouse = filtered_df.pivot_table(index='Mice', columns='Substate', values='SpikeActivityHz', aggfunc='mean')
                try: resultSpikeActivity_perMouse = resultSpikeActivity_perMouse[desired_order]
                except: pass
                filenameOut = f'{folder_to_save2}/{List_name}/{NrSubtype}_VigSt_Spike_perMouse.xlsx'
                writer = pd.ExcelWriter(filenameOut)
                resultSpikeActivity_perMouse.to_excel(writer)
                writer.close()
                """
                proportions = resultSpikeActivity_perUnit['Activated_by'].value_counts(normalize=True)*100
                filenameOut = f'{folder_to_save2}/{List_name}/{NrSubtype}_ActivityVigSt_proportions.xlsx'
                writer = pd.ExcelWriter(filenameOut)
                proportions.to_excel(writer)
                writer.close()
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

                """
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

        filenameOut = f'{folder_to_save}/{Drug}_{NrSubtype}_VigStPropreties.xlsx'
        writer = pd.ExcelWriter(filenameOut)

        combined_df_Drug2 = combined_df_DrugO.drop_duplicates(subset='Substate_ID', keep='first')

        ProportionVigStates = combined_df_Drug2.pivot_table(index='Session_ID', columns='Substate', values='DurationSubstate', aggfunc='sum')
        try: ProportionVigStates = ProportionVigStates[desired_order]
        except: pass
        ProportionVigStates2=ProportionVigStates.div(ProportionVigStates.sum(axis=1), axis=0)*100
        ProportionVigStates2.to_excel(writer, sheet_name='ProportionVigStates')

        DurationVigStates = combined_df_Drug2.pivot_table(index='Session_ID', columns='Substate', values='DurationSubstate', aggfunc='mean')
        DurationVigStates = combined_df_Drug2.pivot_table(index='Substate_ID', columns='Substate', values='DurationSubstate', fill_value=None)
        
        try: DurationVigStates = DurationVigStates[desired_order]
        except: pass
        DurationVigStates.to_excel(writer, sheet_name='EpisodeDurations')

        writer.close()
     