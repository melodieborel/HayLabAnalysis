# Stats on Ca2+ imaging with miniscope and Vigilance States

#######################################################################################
                            # Define Experiment type #
#######################################################################################

#DrugExperiment=1 if CGP Experiment // DrugExperiment=0 if Baseline Experiment
DrugExperiment=1 

suffix={} #to identify this analysis from another

choosed_folder='Analysis_VigStates_2024-07-09_14_52_13_353475_AB'

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

import warnings
warnings.filterwarnings("ignore")


########################################################################
        # SCRIPT 23AB_GrandAverages&Stats_for_VigilanceStates
########################################################################

# Specify the directory containing the Excel files
InitialDirectory = "//10.69.168.1/crnldata/waking/audrey_hay/L1imaging/AnalysedMarch2023/Gaelle/CGP/AB_analysis" if DrugExperiment else f"//10.69.168.1/crnldata/waking/audrey_hay/L1imaging/AnalysedMarch2023/Gaelle/Baseline_recording_ABmodified/AB_Analysis"
directory= f'{InitialDirectory}/{choosed_folder}'

# Get the current date and time
FolderNameSave=str(datetime.now())
FolderNameSave = FolderNameSave.replace(" ", "_").replace(".", "_").replace(":", "_")
destination_folder= f"//10.69.168.1/crnldata/waking/audrey_hay/L1imaging/AnalysedMarch2023/Gaelle/CGP/AB_Analysis/Analysis_AVG_VigStates_{FolderNameSave}{suffix}" if DrugExperiment else f"//10.69.168.1/crnldata/waking/audrey_hay/L1imaging/AnalysedMarch2023/Gaelle/Baseline_recording_ABmodified/AB_Analysis/Analysis_AVG_VigStates_{FolderNameSave}{suffix}"
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
    
    nametofind='VigilanceState_GlobalResults'
    nametofind2='VigilanceState_CaCorrelation'
    nametofind3='VigilanceState_SpCorrelation'

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
                                dfs2_per_sheet[sheet_name] = df2  #one average trace per unique unit, len(df2)==nb unit recorded for that mouse
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

    filenameOut = f'{folder_to_save}/{NrSubtype}_VigilanceStates_GrandGlobalAB.xlsx'
    writer = pd.ExcelWriter(filenameOut)
    combined_df.to_excel(writer)
    writer.close()

    # Save the Ca correlation matrix  

    file_path = f'{folder_to_save}/{NrSubtype}_VigilanceStates_CaCorrelationAB.xlsx'
    with pd.ExcelWriter(file_path) as writer:
        for sheet_name, df in dfs2_per_sheet.items():
            df.to_excel(writer, sheet_name=sheet_name, index=True, header=True)

    # Save the Sp correlation matrix  

    file_path = f'{folder_to_save}/{NrSubtype}_VigilanceStates_SpCorrelationAB.xlsx'
    with pd.ExcelWriter(file_path) as writer:
        for sheet_name, df in dfs3_per_sheet.items():
            df.to_excel(writer, sheet_name=sheet_name, index=True, header=True)

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
    combined_df = pd.read_excel(f'{folder_to_save}/{NrSubtype}_{analysisfile}.xlsx', index_col=0)
    
    combined_dfCa = {}
    analysisfileCa='VigilanceStates_CaCorrelationAB'
    excel_data = pd.read_excel(f'{folder_to_save}/{NrSubtype}_{analysisfileCa}.xlsx', sheet_name=None, index_col=0)           
    for sheet_name, df in excel_data.items():
        if len(df)>0:
            combined_dfCa[sheet_name] = df  

    combined_dfSp = {}
    analysisfileSp='VigilanceStates_SpCorrelationAB'
    excel_data = pd.read_excel(f'{folder_to_save}/{NrSubtype}_{analysisfileSp}.xlsx', sheet_name=None, index_col=0)           
    for sheet_name, df in excel_data.items():
        if len(df)>0:
            combined_dfSp[sheet_name] = df  

    ######################
    # CHOOSE OPTIONS
    ######################

    # MINIMUM VIG STATES DURATION #
    combined_df = combined_df[combined_df['DurationSubstate'] >= 20] 
    
    Drugs= ['Baseline', 'CGP'] if DrugExperiment else ['Baseline']

    # NO LOW FIRING RATE #
    #combined_df = combined_df[combined_df['Avg_SpikeActivityHz'] >= 0.05] 
    
    #####################
    # GLM #
    #####################         
    print(NrSubtype)

    mixedlm_model = sm.MixedLM.from_formula("SpikeActivityHz ~ Substate", groups='Unit_ID', data=combined_df)
    result = mixedlm_model.fit()
    print(result.summary())      
    
    mixedlm_model3 = sm.MixedLM.from_formula("NormalizedAUC_calcium ~ Substate", groups='Unit_ID', data=combined_df)
    result3 = mixedlm_model3.fit()
    print(result3.summary())      
    try: 
        mixedlm_model4 = sm.MixedLM.from_formula("DeconvSpikeMeanActivity ~ Substate", groups='Unit_ID', data=combined_df)
        result4 = mixedlm_model4.fit()
        print(result4.summary())   
    except: 
        print('Error on DeconvSpikeMeanActivity') 

    #####################
    # PREFERENCE #
    #####################
    #"""
    
    for Drug in Drugs: 
        
        folder_to_save2= f'{folder_to_save}/{Drug}/'
        if NrSubtype=='All':
            os.makedirs(folder_to_save2)

        combined_df_Drug=combined_df.copy()

        combined_df_Drug = combined_df_Drug[combined_df_Drug['Drug'] == Drug] if DrugExperiment else combined_df_Drug

        AllUnits = combined_df_Drug['Unit_ID'].unique()

        combined_df_DrugNREM=combined_df_Drug.copy()
        combined_df_DrugNREM['Substate'] = combined_df_DrugNREM['Substate'].apply(lambda x: 'NREM' if x == 'NREM' else 'other')
        grouped = combined_df_DrugNREM.groupby('Unit_ID')
        NREMprefUnits=[]
        NREMmatrix={}
        for key, group in grouped:
            nrem_values = group[group['Substate'] == 'NREM']['DeconvSpikeMeanActivity'].tolist()
            other_values = group[group['Substate'] == 'other']['DeconvSpikeMeanActivity'].tolist()
            NREMmatrix[key]=[nrem_values, other_values]
            if len(nrem_values) > 1 and len(other_values) > 1:  # Ensure there are enough values to perform the test
                t_stat, p_value = ttest_ind(nrem_values, other_values, equal_var=False)
                if p_value<0.2 and np.mean(nrem_values)>np.mean(other_values): 
                    NREMprefUnits.append(key)
        
        combined_df_DrugREM=combined_df_Drug.copy()
        combined_df_DrugREM['Substate'] = combined_df_DrugREM['Substate'].apply(lambda x: 'REM' if x == 'REM' else 'other')
        grouped = combined_df_DrugREM.groupby('Unit_ID')
        REMprefUnits=[]
        REMmatrix={}
        for key, group in grouped:
            rem_values = group[group['Substate'] == 'REM']['DeconvSpikeMeanActivity'].tolist()
            other_values = group[group['Substate'] == 'other']['DeconvSpikeMeanActivity'].tolist()
            REMmatrix[key]=[rem_values, other_values]
            if len(rem_values) > 1 and len(other_values) > 1:  # Ensure there are enough values to perform the test
                t_stat, p_value = ttest_ind(rem_values, other_values, equal_var=False)
                if p_value<0.2  and np.mean(rem_values)>np.mean(other_values):  
                    REMprefUnits.append(key)
    
        combined_df_DrugWake=combined_df_Drug.copy()
        combined_df_DrugWake['Substate'] = combined_df_DrugWake['Substate'].apply(lambda x: 'Wake' if x == 'Wake' else 'other')
        grouped = combined_df_DrugWake.groupby('Unit_ID')
        WakeprefUnits=[]
        Wmatrix={}
        for key, group in grouped:
            w_values = group[group['Substate'] == 'Wake']['DeconvSpikeMeanActivity'].tolist()
            other_values = group[group['Substate'] == 'other']['DeconvSpikeMeanActivity'].tolist()
            Wmatrix[key]=[w_values, other_values]
            if len(w_values) > 1 and len(other_values) > 1:  # Ensure there are enough values to perform the test
                t_stat, p_value = ttest_ind(w_values, other_values, equal_var=False)
                if p_value<0.2 and np.mean(w_values)>np.mean(other_values):  
                    WakeprefUnits.append(key)

        # Save the List of significant Unit more active in one vigilance state
        filenameOut = f'{folder_to_save2}/{NrSubtype}_SignDeconvPreference.xlsx'
        writer = pd.ExcelWriter(filenameOut)
        AllUnits = pd.DataFrame(AllUnits)
        NREMprefUnits = pd.DataFrame(NREMprefUnits)
        REMprefUnits = pd.DataFrame(REMprefUnits)
        WakeprefUnits = pd.DataFrame(WakeprefUnits)
        NREMprefUnits.to_excel(writer, sheet_name='NREMpref', index=True, header=False) 
        REMprefUnits.to_excel(writer, sheet_name='REMpref', index=True, header=False) 
        WakeprefUnits.to_excel(writer, sheet_name='Wakepref', index=True, header=False) 
        AllUnits.to_excel(writer, sheet_name='AllUnits', index=True, header=False) 
        writer.close()

        # Save the Matrix with activity values for one state compare to all the other per Unit
        excel_file =  f'{folder_to_save2}/{NrSubtype}_NREMmatrix.xlsx'
        with pd.ExcelWriter(excel_file) as writer:
            for key, df in NREMmatrix.items():
                df=pd.DataFrame(df)
                df.index = ['NREM', 'Other']
                df.to_excel(writer, sheet_name=key, header=False)

        excel_file =  f'{folder_to_save2}/{NrSubtype}_REMmatrix.xlsx'
        with pd.ExcelWriter(excel_file) as writer:
            for key, df in REMmatrix.items():
                df=pd.DataFrame(df)
                df.index = ['REM', 'Other']
                df.to_excel(writer, sheet_name=key, header=False)

        excel_file =  f'{folder_to_save2}/{NrSubtype}_Wakematrix.xlsx'
        with pd.ExcelWriter(excel_file) as writer:
            for key, df in Wmatrix.items():
                df=pd.DataFrame(df)
                df.index = ['Wake', 'Other']
                df.to_excel(writer, sheet_name=key, header=False)

        #"""

        # Just the highest mean = the preference

        AresultSpikeActivity_perUnit = combined_df_Drug.pivot_table(index='Unit_ID', columns='Substate', values='DeconvSpikeMeanActivity', aggfunc='mean')   
        try : AresultSpikeActivity_perUnit = AresultSpikeActivity_perUnit[desired_order]
        except: pass
        AresultSpikeActivity_perUnit['Activated_by'] = AresultSpikeActivity_perUnit.apply(max_column_name, axis=1)

        index_lists = {}
        for unitindex in AresultSpikeActivity_perUnit['Activated_by'].unique():
            indices = AresultSpikeActivity_perUnit.index[AresultSpikeActivity_perUnit['Activated_by'] == unitindex].tolist()
            index_lists[unitindex] = indices
        
        NREMprefUnits=index_lists['NREM']
        try: REMprefUnits=index_lists['REM'] 
        except: REMprefUnits=[]
        WakeprefUnits=index_lists['Wake']
        AllUnits = combined_df_Drug['Unit_ID'].unique()

        # Save the List of significant Unit more active in one vigilance state
        filenameOut = f'{folder_to_save2}/{NrSubtype}_AverageDeconvPreference.xlsx'
        writer = pd.ExcelWriter(filenameOut)
        AllUnitsDF = pd.DataFrame(AllUnits)
        NREMprefUnitsDF = pd.DataFrame(NREMprefUnits)
        REMprefUnitsDF = pd.DataFrame(REMprefUnits)
        WakeprefUnitsDF = pd.DataFrame(WakeprefUnits)
        NREMprefUnitsDF.to_excel(writer, sheet_name='NREMpref', index=True, header=False) 
        REMprefUnitsDF.to_excel(writer, sheet_name='REMpref', index=True, header=False) 
        WakeprefUnitsDF.to_excel(writer, sheet_name='Wakepref', index=True, header=False) 
        AllUnitsDF.to_excel(writer, sheet_name='AllUnits', index=True, header=False) 
        writer.close()


        List_SignFiringPreference=[NREMprefUnits, REMprefUnits, WakeprefUnits, AllUnits]
        List_Names=['NREMprefUnits', 'REMprefUnits', 'WakeprefUnits', 'AllUnits']    

        for listnb, listI  in enumerate(List_SignFiringPreference):
            
            #list=listI[0] #convert df to list if with pvalue
            filtered_df = combined_df_Drug[combined_df_Drug['Unit_ID'].isin(listI)]
            List_name=List_Names[listnb]

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

                file_path = f'{folder_to_save2}/{List_name}/{NrSubtype}_VigSt_Pairwise_CaCorrelationAB.xlsx'
                with pd.ExcelWriter(file_path) as writer:
                    for sheet_name, dfCa in dfCa_DoubleFiltered.items():
                        dfCa.to_excel(writer, sheet_name=sheet_name, index=True, header=True)

                SummaryMatrixCa= pd.DataFrame()
                for sheet_name, df in dfCa_DoubleFiltered.items():    
                    series_flattened = df.stack().reset_index()
                    series_flattened['combined_index'] = series_flattened['level_0'] + '_' + series_flattened['level_1']
                    series_flattened = series_flattened.set_index('combined_index')[0]
                    SummaryMatrixCa[sheet_name] = series_flattened

                SummaryMatrixCa_cleaned = SummaryMatrixCa[~(SummaryMatrixCa == 1).all(axis=1)]
                Sum = SummaryMatrixCa_cleaned.sum(axis=1)
                mask=(Sum == 3)
                SummaryMatrixCa_cleaned_filtered = SummaryMatrixCa_cleaned[~mask]
                SummaryMatrixCa_cleaned_filtered = SummaryMatrixCa_cleaned_filtered.drop_duplicates() #cause Neuron1_Neuron2 & Neuron2_Neuron1

                filenameOut = f'{folder_to_save2}/{List_name}/{NrSubtype}_VigSt_FlattenPairwise_CaCorrelationAB.xlsx'
                SummaryMatrixCa_cleaned_filtered.to_excel(filenameOut, index=True, header=True)

                # Sp correlation

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

                file_path = f'{folder_to_save2}/{List_name}/{NrSubtype}_VigSt_Pairwise_SpCorrelationAB.xlsx'
                with pd.ExcelWriter(file_path) as writer:
                    for sheet_name, dfSp in dfSp_DoubleFiltered.items():
                        dfSp.to_excel(writer, sheet_name=sheet_name, index=True, header=True)

                SummaryMatrixSp= pd.DataFrame()
                for sheet_name, df in dfSp_DoubleFiltered.items():    
                    series_flattened = df.stack().reset_index()
                    series_flattened['combined_index'] = series_flattened['level_0'] + '_' + series_flattened['level_1']
                    series_flattened = series_flattened.set_index('combined_index')[0]
                    SummaryMatrixSp[sheet_name] = series_flattened

                SummaryMatrixSp_cleaned = SummaryMatrixSp[~(SummaryMatrixSp == 1).all(axis=1)]
                Sum = SummaryMatrixSp_cleaned.sum(axis=1)
                mask=(Sum == 3)
                SummaryMatrixSp_cleaned_filtered = SummaryMatrixSp_cleaned[~mask]
                SummaryMatrixSp_cleaned_filtered = SummaryMatrixSp_cleaned_filtered.drop_duplicates() #cause Neuron1_Neuron2 & Neuron2_Neuron1

                filenameOut = f'{folder_to_save2}/{List_name}/{NrSubtype}_VigSt_FlattenPairwise_SpCorrelationAB.xlsx'
                SummaryMatrixSp_cleaned_filtered.to_excel(filenameOut, index=True, header=True)


            if len(filtered_df)>0:
                #######################
                # DETERMINATION COEFF #
                #######################

                CorrCoeff_calcium_perUnit = filtered_df.pivot_table(index='Unit_ID', columns='Substate', values='Rsquared_CaCorrCoeff', aggfunc='mean')
                try: CorrCoeff_calcium_perUnit = CorrCoeff_calcium_perUnit[desired_order]
                except: pass

                # Save CorrCoeff_calcium_perUnit
                filenameOut = f'{folder_to_save2}/{List_name}/{NrSubtype}_VigSt_Rsquared_CaCorrCoeff_perUnitAB.xlsx'
                writer = pd.ExcelWriter(filenameOut)
                CorrCoeff_calcium_perUnit.to_excel(writer)
                writer.close()

                CorrCoeff_spike_perUnit = filtered_df.pivot_table(index='Unit_ID', columns='Substate', values='Rsquared_SpCorrCoeff', aggfunc='mean')
                try: CorrCoeff_spike_perUnit = CorrCoeff_spike_perUnit[desired_order]
                except: pass
                # Save CorrCoeff_spike_perUnit 
                filenameOut = f'{folder_to_save2}/{List_name}/{NrSubtype}_VigSt_Rsquared_SpCorrCoeff_perUnitAB.xlsx'
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
                filenameOut = f'{folder_to_save2}/{List_name}/{NrSubtype}_VigSt_CaCorrCoeff_perUnitAB.xlsx'
                writer = pd.ExcelWriter(filenameOut)
                CorrCoeff_calcium_perUnit.to_excel(writer)
                writer.close()

                CorrCoeff_spike_perUnit = filtered_df.pivot_table(index='Unit_ID', columns='Substate', values='SpCorrCoeff', aggfunc='mean')
                try: CorrCoeff_spike_perUnit = CorrCoeff_spike_perUnit[desired_order]
                except: pass
                # Save CorrCoeff_spike_perUnit 
                filenameOut = f'{folder_to_save2}/{List_name}/{NrSubtype}_VigSt_SpCorrCoeff_perUnitAB.xlsx'
                writer = pd.ExcelWriter(filenameOut)
                CorrCoeff_spike_perUnit.to_excel(writer)
                writer.close()
                
                #####################
                # AUC CALCIUM #
                #####################
                
                resultNormalizedAUC_calcium_perUnit = filtered_df.pivot_table(index='Unit_ID', columns='Substate', values='NormalizedAUC_calcium', aggfunc='mean')
                try: resultNormalizedAUC_calcium_perUnit = resultNormalizedAUC_calcium_perUnit[desired_order]
                except: pass
                max_values_per_row = resultNormalizedAUC_calcium_perUnit.max(axis=1)
                Diff_resultNormalizedAUC_calcium_perUnit = (resultNormalizedAUC_calcium_perUnit.sub(max_values_per_row, axis=0) +100)
                Diff_resultNormalizedAUC_calcium_perUnit['Activated_by'] = resultNormalizedAUC_calcium_perUnit.apply(max_column_name, axis=1)

                # Add a new column with the column name that has the maximum value
                resultNormalizedAUC_calcium_perUnit['Activated_by'] = resultNormalizedAUC_calcium_perUnit.apply(max_column_name, axis=1)

                # Save  df
                filenameOut = f'{folder_to_save2}/{List_name}/{NrSubtype}_VigSt_AUC_perUnitAB_N.xlsx'
                writer = pd.ExcelWriter(filenameOut)
                resultNormalizedAUC_calcium_perUnit.to_excel(writer)
                writer.close()

                # Save  df
                filenameOut = f'{folder_to_save2}/{List_name}/{NrSubtype}_DiffVigSt_AUC_perUnitAB_N.xlsx'
                writer = pd.ExcelWriter(filenameOut)
                Diff_resultNormalizedAUC_calcium_perUnit.to_excel(writer)
                writer.close()

                proportions = resultNormalizedAUC_calcium_perUnit['Activated_by'].value_counts(normalize=True)*100
                #print(proportions)

                resultNormalizedAUC_calcium_perMouse = filtered_df.pivot_table(index='Mice', columns='Substate', values='NormalizedAUC_calcium', aggfunc='mean')
                try: resultNormalizedAUC_calcium_perMouse = resultNormalizedAUC_calcium_perMouse[desired_order]
                except: pass
                filenameOut = f'{folder_to_save2}/{List_name}/{NrSubtype}_VigSt_AUC_perMouseAB.xlsx'
                writer = pd.ExcelWriter(filenameOut)
                resultNormalizedAUC_calcium_perMouse.to_excel(writer)
                writer.close()

                #####################
                # SPIKE ACTIVITY #
                #####################

                resultSpikeActivity_perUnit = filtered_df.pivot_table(index='Unit_ID', columns='Substate', values='SpikeActivityHz', aggfunc='mean')    
                try: resultSpikeActivity_perUnit = resultSpikeActivity_perUnit[desired_order]
                except: pass
                max_values_per_row = resultSpikeActivity_perUnit.max(axis=1)
                Diff_resultSpikeActivity_perUnit = (resultSpikeActivity_perUnit.sub(max_values_per_row, axis=0) +100)
                Diff_resultSpikeActivity_perUnit['Activated_by'] = resultSpikeActivity_perUnit.apply(max_column_name, axis=1)

                # Add a new column with the column name that has the maximum value
                resultSpikeActivity_perUnit['Activated_by'] = resultSpikeActivity_perUnit.apply(max_column_name, axis=1)

                # Save  df
                filenameOut = f'{folder_to_save2}/{List_name}/{NrSubtype}_VigSt_Spike_perUnitAB.xlsx'
                writer = pd.ExcelWriter(filenameOut)
                resultSpikeActivity_perUnit.to_excel(writer)
                writer.close()

                # Save  df
                filenameOut = f'{folder_to_save2}/{List_name}/{NrSubtype}_DiffVigSt_Spike_perUnitAB_N.xlsx'
                writer = pd.ExcelWriter(filenameOut)
                Diff_resultSpikeActivity_perUnit.to_excel(writer)
                writer.close()

                resultSpikeActivity_perMouse = filtered_df.pivot_table(index='Mice', columns='Substate', values='SpikeActivityHz', aggfunc='mean')
                try: resultSpikeActivity_perMouse = resultSpikeActivity_perMouse[desired_order]
                except: pass
                filenameOut = f'{folder_to_save2}/{List_name}/{NrSubtype}_VigSt_Spike_perMouseAB.xlsx'
                writer = pd.ExcelWriter(filenameOut)
                resultSpikeActivity_perMouse.to_excel(writer)
                writer.close()

                proportions = resultSpikeActivity_perUnit['Activated_by'].value_counts(normalize=True)*100
                filenameOut = f'{folder_to_save2}/{List_name}/{NrSubtype}_ActivityVigSt_proportions.xlsx'
                writer = pd.ExcelWriter(filenameOut)
                proportions.to_excel(writer)
                writer.close()

                #####################
                # DECONV ACTIVITY #
                #####################

                resultSpikeActivity_perUnit = filtered_df.pivot_table(index='Unit_ID', columns='Substate', values='DeconvSpikeMeanActivity', aggfunc='mean')    
                try: resultSpikeActivity_perUnit = resultSpikeActivity_perUnit[desired_order]
                except: pass
                max_values_per_row = resultSpikeActivity_perUnit.max(axis=1)
                Diff_resultSpikeActivity_perUnit = (resultSpikeActivity_perUnit.sub(max_values_per_row, axis=0) +100)
                Diff_resultSpikeActivity_perUnit['Activated_by'] = resultSpikeActivity_perUnit.apply(max_column_name, axis=1)

                # Add a new column with the column name that has the maximum value
                resultSpikeActivity_perUnit['Activated_by'] = resultSpikeActivity_perUnit.apply(max_column_name, axis=1)

                # Save  df
                filenameOut = f'{folder_to_save2}/{List_name}/{NrSubtype}_VigSt_DeconvSpike_perUnitAB.xlsx'
                writer = pd.ExcelWriter(filenameOut)
                resultSpikeActivity_perUnit.to_excel(writer)
                writer.close()

                # Save  df
                filenameOut = f'{folder_to_save2}/{List_name}/{NrSubtype}_DiffVigSt_DeconvSpike_perUnitAB_N.xlsx'
                writer = pd.ExcelWriter(filenameOut)
                Diff_resultSpikeActivity_perUnit.to_excel(writer)
                writer.close()

                resultSpikeActivity_perMouse = filtered_df.pivot_table(index='Mice', columns='Substate', values='DeconvSpikeMeanActivity', aggfunc='mean')
                try: resultSpikeActivity_perMouse = resultSpikeActivity_perMouse[desired_order]
                except: pass
                filenameOut = f'{folder_to_save2}/{List_name}/{NrSubtype}_VigSt_DeconvSpike_perMouseAB.xlsx'
                writer = pd.ExcelWriter(filenameOut)
                resultSpikeActivity_perMouse.to_excel(writer)
                writer.close()

                proportions = resultSpikeActivity_perUnit['Activated_by'].value_counts(normalize=True)*100
                filenameOut = f'{folder_to_save2}/{List_name}/{NrSubtype}_ActivityVigSt_DeconvSpike_proportions.xlsx'
                writer = pd.ExcelWriter(filenameOut)
                proportions.to_excel(writer)
                writer.close()

        #######################
        # Propreties VigStates
        #######################

        filenameOut = f'{folder_to_save}/{NrSubtype}_PropretiesVigStatesAB.xlsx'
        writer = pd.ExcelWriter(filenameOut)

        ProportionVigStates = filtered_df.pivot_table(index='Session_ID', columns='Substate', values='DurationSubstate', aggfunc='sum')
        try: ProportionVigStates = ProportionVigStates[desired_order]
        except: pass
        ProportionVigStates2=ProportionVigStates.div(ProportionVigStates.sum(axis=1), axis=0)*100
        ProportionVigStates2.to_excel(writer, sheet_name='ProportionVigStates')

        DurationVigStates = combined_df.pivot_table(index='Session_ID', columns='Substate', values='DurationSubstate', aggfunc='mean')
        try: DurationVigStates = DurationVigStates[desired_order]
        except: pass
        DurationVigStates.to_excel(writer, sheet_name='EpisodeDurations')

        writer.close()
     