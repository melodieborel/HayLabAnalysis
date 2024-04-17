
# Stats on Ca2+ imaging with miniscope and Osc

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
from scipy.stats import zscore


########################################################################
        # SCRIPT 27AB_GrandAverages&Stats_for_Osc
########################################################################

# Specify the directory containing the Excel files
directory = "//10.69.168.1/crnldata/waking/audrey_hay/L1imaging/AnalysedMarch2023/Gaelle/Baseline_recording/"

NrSubtypeList=['All', 'L1']
CortexList=['PFC', 'S1']

OscList=['Spdl', 'SWR']
OscillationList=['Spindles', 'SWR']

for o, Osc in enumerate(OscList): 
    
    print(Osc, 'oscillations analysis...')

    for NrSubtype in NrSubtypeList: 
        
        print('... for', NrSubtype, 'neurons...')

        for Cortex in CortexList:
            
            print('... in the', Cortex)

            # Initialize an empty list to store the dataframes
            dfs = []
            df=[]
            dfs2 = []
            df2=[]
            dfs2_per_sheet = {}
            dfs3 = []
            df3=[]
            dfs3_per_sheet = {}

            if NrSubtype=='L1':
                MiceList=['BlackLinesOK', 'BlueLinesOK', 'GreenDotsOK', 'GreenLinesOK', 'RedLinesOK']
            else:
                MiceList=['Purple', 'ThreeColDotsOK', 'ThreeBlueCrossesOK']

            nametofind=f'{OscillationList[o]}_{Cortex}_ABdetection_GlobalResultsAB'
            nametofind2=f'{OscillationList[o]}_{Cortex}_ABdetection_CalciumAvgResultsAB'
            nametofind3=f'{OscillationList[o]}_{Cortex}_ABdetection_SpikeAvgResultsAB'

            # Recursively traverse the directory structure
            for root, _, files in os.walk(directory):
                for filename in files:
                    # Check if the file is an Excel file and contains the specified name
                    if filename.endswith('.xlsx') and nametofind in filename:
                        if any(name in root for name in MiceList):  
                            # Construct the full path to the file
                            filepath = os.path.join(root, filename)
                            # Read the Excel file into a dataframe and append it to the list
                            df = pd.read_excel(filepath, index_col=0)
                            dfs.append(df)
                    if filename.endswith('.xlsx') and nametofind2 in filename: 
                        if any(name in root for name in MiceList): 
                            # Construct the full path to the file
                            filepath = os.path.join(root, filename)
                            # Read the Excel file into a dataframe and append it to the list
                            excel_data = pd.read_excel(filepath, sheet_name=None, index_col=0)            
                            for sheet_name, df2 in excel_data.items():
                                if sheet_name in dfs2_per_sheet:   
                                    if len(df2)>0:
                                        updated_matrix = np.concatenate([dfs2_per_sheet[sheet_name], df2], axis=0)                    
                                        dfs2_per_sheet[sheet_name] = updated_matrix               
                                else:        
                                    dfs2_per_sheet[sheet_name] = df2 #one average trace per unique unit, len(df2)==nb unit recorded for that mouse

                    if filename.endswith('.xlsx') and nametofind3 in filename: 
                        if any(name in root for name in MiceList): 
                            # Construct the full path to the file
                            filepath = os.path.join(root, filename)
                            # Read the Excel file into a dataframe and append it to the list
                            excel_data = pd.read_excel(filepath, sheet_name=None, index_col=0, header=None)           
                            for sheet_name, df3 in excel_data.items():
                                if sheet_name in dfs3_per_sheet:   
                                    updated_matrix = pd.concat((dfs3_per_sheet[sheet_name],df3), ignore_index=True, axis=0)                
                                    dfs3_per_sheet[sheet_name] = updated_matrix                    
                                else:                    
                                    dfs3_per_sheet[sheet_name] = df3 #one average trace per unique unit, len(df3)==nb unit recorded for that mouse

            # Concatenate all dataframes into a single dataframe
            combined_df = pd.concat(dfs, ignore_index=True)

            combined_df['Unique_Unit'] = combined_df['Unique_Unit'].astype(str)
            combined_df['UnitNumber'] = combined_df['UnitNumber'].astype(str)
            combined_df['UnitValue'] = combined_df['UnitValue'].astype(str)

            combined_df['Unit_ID'] = combined_df['Mice'] + combined_df['Unique_Unit']

            unique_count = combined_df['Unit_ID'].nunique()
            print(unique_count, f'{NrSubtype} neurons recorded') 

            # Remove non defined Unique Units 
            combined_df = combined_df[combined_df['Unique_Unit'] != '[]']
            unique_count = combined_df['Unit_ID'].nunique()
            print(unique_count, f'{NrSubtype} neurons in the cross-registration') 
            
            combined_df[f'{Osc}_ID'] = combined_df['Mice'] + combined_df['Session'] + combined_df[f'{Osc}Number'].astype(str)
            
            unique_count = combined_df[f'{Osc}_ID'].nunique()
            print(unique_count, f'{Osc} recorded in total in the {Cortex}')

            # GLOBAL

            filenameOut = f'{directory}{NrSubtype}_{Osc}_{Cortex}_ABdetection_GrandGlobalAB.xlsx'
            writer = pd.ExcelWriter(filenameOut)
            combined_df.to_excel(writer)
            writer.close()

            # CALCIUM traces

            filenameOut = f'{directory}{NrSubtype}_{Osc}_{Cortex}_ABdetection_CalciumGrandAverageAB.xlsx'
            excel_writer = pd.ExcelWriter(filenameOut)

            Array=pd.DataFrame(dfs2_per_sheet[f'All_{OscillationList[o]}'])
            ArrayUn=pd.DataFrame(dfs2_per_sheet[f'Uncoupled_{OscillationList[o]}'])
            ArrayPre=pd.DataFrame(dfs2_per_sheet[f'Precoupled_{OscillationList[o]}'])
            ArrayPost=pd.DataFrame(dfs2_per_sheet[f'Postcoupled_{OscillationList[o]}'])

            Array.to_excel(excel_writer, sheet_name=f'All_{OscillationList[o]}', index=True, header=False)
            ArrayUn.to_excel(excel_writer, sheet_name=f'Uncoupled_{OscillationList[o]}', index=True, header=False)
            ArrayPre.to_excel(excel_writer, sheet_name=f'Precoupled_{OscillationList[o]}', index=True, header=False)
            ArrayPost.to_excel(excel_writer, sheet_name=f'Postcoupled_{OscillationList[o]}', index=True, header=False)

            excel_writer.close()

            # CALCIUM traces Z-score normalization

            filenameOut = f'{directory}{NrSubtype}_{Osc}_{Cortex}_ABdetection_Zscored_CalciumGrandAverageAB.xlsx'
            excel_writer = pd.ExcelWriter(filenameOut)

            Array=Array.T
            Array = Array.apply(zscore)
            ArrayUn=ArrayUn.apply(zscore)
            ArrayPre=ArrayPre.apply(zscore)
            ArrayPost=ArrayPost.apply(zscore)

            Array.to_excel(excel_writer, sheet_name=f'All_{OscillationList[o]}', index=True, header=False)
            ArrayUn.to_excel(excel_writer, sheet_name=f'Uncoupled_{OscillationList[o]}', index=True, header=False)
            ArrayPre.to_excel(excel_writer, sheet_name=f'Precoupled_{OscillationList[o]}', index=True, header=False)
            ArrayPost.to_excel(excel_writer, sheet_name=f'Postcoupled_{OscillationList[o]}', index=True, header=False)

            excel_writer.close()

            # SPIKE

            filenameOut = f'{directory}{NrSubtype}_{Osc}_{Cortex}_ABdetection_SpikeGrandAverageAB.xlsx'
            excel_writer = pd.ExcelWriter(filenameOut)

            Array=pd.DataFrame(dfs3_per_sheet[f'All_{OscillationList[o]}'])
            ArrayUn=pd.DataFrame(dfs3_per_sheet[f'Uncoupled_{OscillationList[o]}'])
            ArrayPre=pd.DataFrame(dfs3_per_sheet[f'Precoupled_{OscillationList[o]}'])
            ArrayPost=pd.DataFrame(dfs3_per_sheet[f'Postcoupled_{OscillationList[o]}'])

            Array.to_excel(excel_writer, sheet_name=f'All_{OscillationList[o]}', index=True, header=False)
            ArrayUn.to_excel(excel_writer, sheet_name=f'Uncoupled_{OscillationList[o]}', index=True, header=False)
            ArrayPre.to_excel(excel_writer, sheet_name=f'Precoupled_{OscillationList[o]}', index=True, header=False)
            ArrayPost.to_excel(excel_writer, sheet_name=f'Postcoupled_{OscillationList[o]}', index=True, header=False)

            excel_writer.close()

        
    ########################################################################
            # SCRIPT 28AB_GrandAverages&Stats_for_Osc
    ########################################################################


    for NrSubtype in NrSubtypeList:  

        for Cortex in CortexList:

            # Load the Excel file into a DataFrame

            analysisfile=f'{Osc}_{Cortex}_ABdetection_GrandGlobalAB'
            combined_df = pd.read_excel(f'{directory}{NrSubtype}_{analysisfile}.xlsx', index_col=0)            

            analysisfile2=f'{Osc}_{Cortex}_ABdetection_CalciumGrandAverageAB'
            excel_file = f'{directory}{NrSubtype}_{analysisfile2}.xlsx'        
            excel_data = pd.read_excel(excel_file, sheet_name=None, index_col=0, header=None) 
            dict_Osc_average = {}
            for sheet_name, sheet_data in excel_data.items():
                dict_Osc_average[sheet_name] = sheet_data

            analysisfile3=f'{Osc}_{Cortex}_ABdetection_SpikeGrandAverageAB'
            excel_file = f'{directory}{NrSubtype}_{analysisfile3}.xlsx'
            excel_data = pd.read_excel(excel_file, sheet_name=None, index_col=0, header=None) 
            dict_Osc_Spikeaverage = {}
            for sheet_name, sheet_data in excel_data.items():
                dict_Osc_Spikeaverage[sheet_name] = sheet_data

            def clean_string(s):
                return s.replace("[", "").replace("]", "").replace("'", "")

            # Apply the function to the specific column
            combined_df[f'{Osc}Statut'] = combined_df[f'{Osc}Statut'].apply(clean_string)

            Oscdurmean = combined_df.groupby(f'{Osc}_ID')[f'{Osc}Duration (ms)'].mean()
            filenameOut = f'{directory}{NrSubtype}_{Osc}durations_AB.xlsx'
            writer = pd.ExcelWriter(filenameOut)
            Oscdurmean.to_excel(writer)
            writer.close()
            print(Cortex, Osc, 'mean duration = ', np.mean(Oscdurmean), 'ms for', NrSubtype, 'neurons')

            # Create df for All Osc
            AUC_calciumBefore_perUnit = combined_df.groupby('Unit_ID')['AUC_calciumBefore'].mean()
            AUC_calciumAfter_perUnit = combined_df.groupby('Unit_ID')['AUC_calciumAfter'].mean()
            allAUCresult = pd.concat([AUC_calciumBefore_perUnit,AUC_calciumAfter_perUnit], axis=1)

            # Create dfs for each Osc Statut
            OscStatut_AUC_calciumBefore_perUnit = combined_df.pivot_table(index='Unit_ID', columns=f'{Osc}Statut', values='AUC_calciumBefore', aggfunc='mean')
            OscStatutNames=OscStatut_AUC_calciumBefore_perUnit.columns.tolist()
            new_columns = [f'Before_{col}' for col in OscStatut_AUC_calciumBefore_perUnit.columns]
            OscStatut_AUC_calciumBefore_perUnit.columns = new_columns
            OscStatut_AUC_calciumAfter_perUnit = combined_df.pivot_table(index='Unit_ID', columns=f'{Osc}Statut', values='AUC_calciumAfter', aggfunc='mean')
            new_columns = [f'After_{col}' for col in OscStatut_AUC_calciumAfter_perUnit.columns]
            OscStatut_AUC_calciumAfter_perUnit.columns = new_columns
            # Combine each Osc Statut df into one
            AUC_interleaved_cols = [col for pair in zip(OscStatut_AUC_calciumBefore_perUnit.columns, OscStatut_AUC_calciumAfter_perUnit.columns) for col in pair]
            AUCresult = pd.concat([OscStatut_AUC_calciumBefore_perUnit[col] if col in OscStatut_AUC_calciumBefore_perUnit.columns else OscStatut_AUC_calciumAfter_perUnit[col] for col in AUC_interleaved_cols], axis=1)

            # Create a dict containing each Osc Statut and All Osc

            AUCdata_dict = {}
            c=0
            for i in range(len(OscStatutNames)):
                columns = AUCresult.columns[c:c+2]
                c=c+2
                AUCdata_dict[OscStatutNames[i]] = AUCresult[columns].values.tolist()
            AUCdata_dict[f'All{Osc}']=allAUCresult.values.tolist()


            AUCdata_dict_ActiveOnly={}
            for key in AUCdata_dict:
                # Filter out rows with non-zero values
                filtered_rows = [row for row in AUCdata_dict[key] if any(val != 0 for val in row) and np.any(~np.isnan(row))]
                # Update the key with filtered rows
                AUCdata_dict_ActiveOnly[key] = filtered_rows


            filenameOut = f'{directory}{NrSubtype}_{Osc}_{Cortex}_AUC_perUnitAB.xlsx'
            with pd.ExcelWriter(filenameOut) as writer:
                # Iterate over each key-value pair in the dictionary
                for sheet_name, data_list in AUCdata_dict.items():
                    # Convert the list to a DataFrame
                    df = pd.DataFrame(data_list, columns=['Before', 'After'])
                    # Write the DataFrame to a separate sheet
                    df.to_excel(writer, sheet_name=sheet_name, index=False)

            filenameOut = f'{directory}{NrSubtype}_{Osc}_{Cortex}_AUC_perUnit_Active OnlyAB.xlsx'
            with pd.ExcelWriter(filenameOut) as writer:
                # Iterate over each key-value pair in the dictionary
                for sheet_name, data_list in AUCdata_dict_ActiveOnly.items():
                    # Convert the list to a DataFrame
                    df = pd.DataFrame(data_list, columns=['Before', 'After'])
                    # Write the DataFrame to a separate sheet
                    df.to_excel(writer, sheet_name=sheet_name, index=False)

            SpikeActivityPreference_perUnits = combined_df.pivot_table(index='Unit_ID', columns='SpikeActivityPreference', aggfunc='size', fill_value=0)
            SpikeActivityPreference_perUnits['%PrefBefore']=(SpikeActivityPreference_perUnits['Before'])/(SpikeActivityPreference_perUnits['Before']+SpikeActivityPreference_perUnits['After'])*100
            SpikeActivityPreference_perUnits['%PrefAfter']=(SpikeActivityPreference_perUnits['After'])/(SpikeActivityPreference_perUnits['Before']+SpikeActivityPreference_perUnits['After'])*100

            CalciumActivityPreference_perUnits = combined_df.pivot_table(index='Unit_ID', columns='CalciumActivityPreference', aggfunc='size', fill_value=0)
            CalciumActivityPreference_perUnits['%PrefBefore']=(CalciumActivityPreference_perUnits['Before'])/(CalciumActivityPreference_perUnits['Before']+CalciumActivityPreference_perUnits['After'])*100
            CalciumActivityPreference_perUnits['%PrefAfter']=(CalciumActivityPreference_perUnits['After'])/(CalciumActivityPreference_perUnits['Before']+CalciumActivityPreference_perUnits['After'])*100

            # Create df for All {Osc}
            SpikeActivityBefore_perUnit = combined_df.groupby('Unit_ID')['SpikeActivityBefore'].mean()
            SpikeActivityAfter_perUnit = combined_df.groupby('Unit_ID')['SpikeActivityAfter'].mean()
            allSPresult = pd.concat([SpikeActivityBefore_perUnit,SpikeActivityAfter_perUnit], axis=1)

            # Create dfs for each {Osc} Statut
            OscStatut_SpikeActivityBefore_perUnit = combined_df.pivot_table(index='Unit_ID', columns=f'{Osc}Statut', values='SpikeActivityBefore', aggfunc='mean')
            OscStatutNames=OscStatut_SpikeActivityBefore_perUnit.columns.tolist()
            new_columns = [f'Before_{col}' for col in OscStatut_SpikeActivityBefore_perUnit.columns]
            OscStatut_SpikeActivityBefore_perUnit.columns = new_columns
            OscStatut_SpikeActivityAfter_perUnit = combined_df.pivot_table(index='Unit_ID', columns=f'{Osc}Statut', values='SpikeActivityAfter', aggfunc='mean')
            new_columns = [f'After_{col}' for col in OscStatut_SpikeActivityAfter_perUnit.columns]
            OscStatut_SpikeActivityAfter_perUnit.columns = new_columns
            # Combine each Osc Statut df into one
            SP_interleaved_cols = [col for pair in zip(OscStatut_SpikeActivityBefore_perUnit.columns, OscStatut_SpikeActivityAfter_perUnit.columns) for col in pair]
            SPresult = pd.concat([OscStatut_SpikeActivityBefore_perUnit[col] if col in OscStatut_SpikeActivityBefore_perUnit.columns else OscStatut_SpikeActivityAfter_perUnit[col] for col in SP_interleaved_cols], axis=1)

            # Create a dict containing each Osc Statut and All Osc
            SPdata_dict = {}
            c=0
            for i in range(len(OscStatutNames)):
                columns = SPresult.columns[c:c+2]
                c=c+2
                SPdata_dict[OscStatutNames[i]] = SPresult[columns].values.tolist()
            SPdata_dict[f'All{Osc}']=allSPresult.values.tolist()


            SPdata_dict_ActiveOnly={}
            for key in SPdata_dict:
                # Filter out rows with non-zero values
                filtered_rows = [row for row in SPdata_dict[key] if any(val != 0 for val in row) and np.any(~np.isnan(row))]
                # Update the key with filtered rows
                SPdata_dict_ActiveOnly[key] = filtered_rows


            filenameOut = f'{directory}{NrSubtype}_{Osc}_{Cortex}_SpikeAct_perUnitAB.xlsx'
            with pd.ExcelWriter(filenameOut) as writer:
                # Iterate over each key-value pair in the dictionary
                for sheet_name, data_list in SPdata_dict.items():
                    # Convert the list to a DataFrame
                    df = pd.DataFrame(data_list, columns=['Before', 'After'])
                    # Write the DataFrame to a separate sheet
                    df.to_excel(writer, sheet_name=sheet_name, index=False)

            filenameOut = f'{directory}{NrSubtype}_{Osc}_{Cortex}_SpikeAct_perUnit_Active OnlyAB.xlsx'
            with pd.ExcelWriter(filenameOut) as writer:
                # Iterate over each key-value pair in the dictionary
                for sheet_name, data_list in SPdata_dict_ActiveOnly.items():
                    # Convert the list to a DataFrame
                    df = pd.DataFrame(data_list, columns=['Before', 'After'])
                    # Write the DataFrame to a separate sheet
                    df.to_excel(writer, sheet_name=sheet_name, index=False)
