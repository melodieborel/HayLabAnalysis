

# Stats on Ca2+ imaging with miniscope and Spindles

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


########################################################################
        # SCRIPT 25AB_GrandAverages&Stats_for_Spindles
########################################################################


# Specify the directory containing the Excel files
directory = "//10.69.168.1/crnldata/waking/audrey_hay/L1imaging/AnalysedMarch2023/Gaelle/Baseline_recording/"

NrSubtypeList=['All', 'L1']
CortexList=['PFC', 'S1']

for NrSubtype in NrSubtypeList:  

    for Cortex in CortexList:

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

        nametofind=f'Spindles_{Cortex}_ABdetection_GlobalResultsAB'
        nametofind2=f'Spindles_{Cortex}_ABdetection_CalciumAvgResultsAB'
        nametofind3=f'Spindles_{Cortex}_ABdetection_SpikeAvgResultsAB'

        # Recursively traverse the directory structure
        for root, _, files in os.walk(directory):
            for filename in files:
                # Check if the file is an Excel file and contains the specified name
                if filename.endswith('.xlsx') and nametofind in filename: 
                    # Construct the full path to the file
                    filepath = os.path.join(root, filename)
                    # Read the Excel file into a dataframe and append it to the list
                    df = pd.read_excel(filepath, index_col=0)
                    dfs.append(df)
                if filename.endswith('.xlsx') and nametofind2 in filename: 
                    # Construct the full path to the file
                    filepath = os.path.join(root, filename)
                    # Read the Excel file into a dataframe and append it to the list
                    excel_data = pd.read_excel(filepath, sheet_name=None, index_col=0)            
                    for sheet_name, df2 in excel_data.items():
                        if sheet_name in dfs2_per_sheet:   
                            updated_matrix = np.concatenate([dfs2_per_sheet[sheet_name], df2], axis=0)                    
                            dfs2_per_sheet[sheet_name] = updated_matrix                    
                        else:                    
                            dfs2_per_sheet[sheet_name] = df2 #one average trace per unique unit, len(df2)==nb unit recorded for that mouse
                if filename.endswith('.xlsx') and nametofind3 in filename: 
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

        unique_count = combined_df['Unique_Unit'].nunique()
        print(unique_count, 'units total') # weird that not the same nb as the len of average df2


        #Remove non defined Unique Units 
        combined_df = combined_df[combined_df['Unique_Unit'] != '[]']
        combined_df['Unique_Unit'] = combined_df['Unique_Unit'].astype(str)
        combined_df['UnitNumber'] = combined_df['UnitNumber'].astype(str)
        combined_df['UnitValue'] = combined_df['UnitValue'].astype(str)

        combined_df['Unit_ID'] = combined_df['Mice'] + combined_df['Unique_Unit']
        combined_df['Spdl_ID'] = combined_df['Mice'] + combined_df['Session'] +combined_df['SpdlNumber'].astype(str)

        combined_df

        unique_count = combined_df['Unique_Unit'].nunique()
        print(unique_count, 'units total') # weird that not the same nb as the len of average df2
        unique_count = combined_df['Spdl_ID'].nunique()
        print(unique_count, 'spdl total')

        # Save big dataset for stats

        filenameOut = f'{directory}{NrSubtype}_Spindles_{Cortex}_ABdetection_GrandGlobalAB.xlsx'
        writer = pd.ExcelWriter(filenameOut)
        combined_df.to_excel(writer)
        writer.close()

        filenameOut = f'{directory}{NrSubtype}_Spindles_{Cortex}_ABdetection_CalciumGrandAverageAB.xlsx'
        excel_writer = pd.ExcelWriter(filenameOut)

        Array=pd.DataFrame(dfs2_per_sheet['All_Spindles'])
        ArrayUn=pd.DataFrame(dfs2_per_sheet['Uncoupled_Spindles'])
        ArrayPre=pd.DataFrame(dfs2_per_sheet['Precoupled_Spindles'])
        ArrayPost=pd.DataFrame(dfs2_per_sheet['Postcoupled_Spindles'])

        Array.to_excel(excel_writer, sheet_name='All_Spindles', index=True, header=False)
        ArrayUn.to_excel(excel_writer, sheet_name='Uncoupled_Spindles', index=True, header=False)
        ArrayPre.to_excel(excel_writer, sheet_name='Precoupled_Spindles', index=True, header=False)
        ArrayPost.to_excel(excel_writer, sheet_name='Postcoupled_Spindles', index=True, header=False)

        # Save the Excel file
        excel_writer.close()

        plt.close()
        test=np.array(dfs2_per_sheet['All_Spindles']) #one average trace per unique unit 
        plt.plot(test.T)
        plt.show()

        filenameOut = f'{directory}{NrSubtype}_Spindles_{Cortex}_ABdetection_SpikeGrandAverageAB.xlsx'
        excel_writer = pd.ExcelWriter(filenameOut)

        Array=pd.DataFrame(dfs3_per_sheet['All_Spindles'])
        ArrayUn=pd.DataFrame(dfs3_per_sheet['Uncoupled_Spindles'])
        ArrayPre=pd.DataFrame(dfs3_per_sheet['Precoupled_Spindles'])
        ArrayPost=pd.DataFrame(dfs3_per_sheet['Postcoupled_Spindles'])

        Array.to_excel(excel_writer, sheet_name='All_Spindles', index=True, header=False)
        ArrayUn.to_excel(excel_writer, sheet_name='Uncoupled_Spindles', index=True, header=False)
        ArrayPre.to_excel(excel_writer, sheet_name='Precoupled_Spindles', index=True, header=False)
        ArrayPost.to_excel(excel_writer, sheet_name='Postcoupled_Spindles', index=True, header=False)

        # Save the Excel file
        excel_writer.close()

########################################################################
        # SCRIPT 26AB_GrandAverages&Stats_for_Spindles
########################################################################

for NrSubtype in NrSubtypeList:  

    for Cortex in CortexList:
        
        # Load the Excel file into a DataFrame

        analysisfile=f'Spindles_{Cortex}_ABdetection_GrandGlobalAB'
        combined_df = pd.read_excel(f'{directory}{NrSubtype}_{analysisfile}.xlsx', index_col=0)            

        analysisfile2=f'Spindles_{Cortex}_ABdetection_CalciumGrandAverageAB'
        excel_file = f'{directory}{NrSubtype}_{analysisfile2}.xlsx'
        excel_data = pd.read_excel(excel_file, sheet_name=None, index_col=0, header=None) 
        dict_Spindles_average = {}
        for sheet_name, sheet_data in excel_data.items():
            dict_Spindles_average[sheet_name] = sheet_data


        analysisfile3=f'Spindles_{Cortex}_ABdetection_SpikeGrandAverageAB'
        excel_file = f'{directory}{NrSubtype}_{analysisfile3}.xlsx'
        excel_data = pd.read_excel(excel_file, sheet_name=None, index_col=0, header=None) 
        dict_Spindles_Spikeaverage = {}
        for sheet_name, sheet_data in excel_data.items():
            dict_Spindles_Spikeaverage[sheet_name] = sheet_data

        def clean_string(s):
            return s.replace("[", "").replace("]", "").replace("'", "")

        # Apply the function to the specific column
        combined_df['SpdlStatut'] = combined_df['SpdlStatut'].apply(clean_string)

        Spdldurmean = combined_df.groupby('Spdl_ID')['SpdlDuration (ms)'].mean()
        filenameOut = f'{directory}{NrSubtype}_SpdlPFCdurations_AB.xlsx'
        writer = pd.ExcelWriter(filenameOut)
        Spdldurmean.to_excel(writer)
        writer.close()
        print(np.mean(Spdldurmean))

        # Create df for All Spdl
        AUC_calciumBefore_perUnit = combined_df.groupby('Unit_ID')['AUC_calciumBefore'].mean()
        AUC_calciumAfter_perUnit = combined_df.groupby('Unit_ID')['AUC_calciumAfter'].mean()
        allAUCresult = pd.concat([AUC_calciumBefore_perUnit,AUC_calciumAfter_perUnit], axis=1)

        # Create dfs for each Spdl Statut
        SpdlStatut_AUC_calciumBefore_perUnit = combined_df.pivot_table(index='Unit_ID', columns='SpdlStatut', values='AUC_calciumBefore', aggfunc='mean')
        SpdlStatutNames=SpdlStatut_AUC_calciumBefore_perUnit.columns.tolist()
        new_columns = [f'Before_{col}' for col in SpdlStatut_AUC_calciumBefore_perUnit.columns]
        SpdlStatut_AUC_calciumBefore_perUnit.columns = new_columns
        SpdlStatut_AUC_calciumAfter_perUnit = combined_df.pivot_table(index='Unit_ID', columns='SpdlStatut', values='AUC_calciumAfter', aggfunc='mean')
        new_columns = [f'After_{col}' for col in SpdlStatut_AUC_calciumAfter_perUnit.columns]
        SpdlStatut_AUC_calciumAfter_perUnit.columns = new_columns
        # Combine each Spdl Statut df into one
        AUC_interleaved_cols = [col for pair in zip(SpdlStatut_AUC_calciumBefore_perUnit.columns, SpdlStatut_AUC_calciumAfter_perUnit.columns) for col in pair]
        AUCresult = pd.concat([SpdlStatut_AUC_calciumBefore_perUnit[col] if col in SpdlStatut_AUC_calciumBefore_perUnit.columns else SpdlStatut_AUC_calciumAfter_perUnit[col] for col in AUC_interleaved_cols], axis=1)

        # Create a dict containing each Spdl Statut and All Spdl

        AUCdata_dict = {}
        c=0
        for i in range(len(SpdlStatutNames)):
            columns = AUCresult.columns[c:c+2]
            c=c+2
            AUCdata_dict[SpdlStatutNames[i]] = AUCresult[columns].values.tolist()
        AUCdata_dict['AllSpdl']=allAUCresult.values.tolist()


        AUCdata_dict_ActiveOnly={}
        for key in AUCdata_dict:
            # Filter out rows with non-zero values
            filtered_rows = [row for row in AUCdata_dict[key] if any(val != 0 for val in row) and np.any(~np.isnan(row))]
            # Update the key with filtered rows
            AUCdata_dict_ActiveOnly[key] = filtered_rows


        filenameOut = f'{directory}{NrSubtype}_Spdl_{Cortex}_AUC_perUnitAB.xlsx'
        with pd.ExcelWriter(filenameOut) as writer:
            # Iterate over each key-value pair in the dictionary
            for sheet_name, data_list in AUCdata_dict.items():
                # Convert the list to a DataFrame
                df = pd.DataFrame(data_list, columns=['Before', 'After'])
                # Write the DataFrame to a separate sheet
                df.to_excel(writer, sheet_name=sheet_name, index=False)

        filenameOut = f'{directory}{NrSubtype}_Spdl_{Cortex}_AUC_perUnit_Active_OnlyAB.xlsx'
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

        # Create df for All Spdl
        SpikeActivityBefore_perUnit = combined_df.groupby('Unit_ID')['SpikeActivityBefore'].mean()
        SpikeActivityAfter_perUnit = combined_df.groupby('Unit_ID')['SpikeActivityAfter'].mean()
        allSPresult = pd.concat([SpikeActivityBefore_perUnit,SpikeActivityAfter_perUnit], axis=1)

        # Create dfs for each Spdl Statut
        SpdlStatut_SpikeActivityBefore_perUnit = combined_df.pivot_table(index='Unit_ID', columns='SpdlStatut', values='SpikeActivityBefore', aggfunc='mean')
        SpdlStatutNames=SpdlStatut_SpikeActivityBefore_perUnit.columns.tolist()
        new_columns = [f'Before_{col}' for col in SpdlStatut_SpikeActivityBefore_perUnit.columns]
        SpdlStatut_SpikeActivityBefore_perUnit.columns = new_columns
        SpdlStatut_SpikeActivityAfter_perUnit = combined_df.pivot_table(index='Unit_ID', columns='SpdlStatut', values='SpikeActivityAfter', aggfunc='mean')
        new_columns = [f'After_{col}' for col in SpdlStatut_SpikeActivityAfter_perUnit.columns]
        SpdlStatut_SpikeActivityAfter_perUnit.columns = new_columns
        # Combine each Spdl Statut df into one
        SP_interleaved_cols = [col for pair in zip(SpdlStatut_SpikeActivityBefore_perUnit.columns, SpdlStatut_SpikeActivityAfter_perUnit.columns) for col in pair]
        SPresult = pd.concat([SpdlStatut_SpikeActivityBefore_perUnit[col] if col in SpdlStatut_SpikeActivityBefore_perUnit.columns else SpdlStatut_SpikeActivityAfter_perUnit[col] for col in SP_interleaved_cols], axis=1)


        # Create a dict containing each Spdl Statut and All Spdl
        SPdata_dict = {}
        c=0
        for i in range(len(SpdlStatutNames)):
            columns = SPresult.columns[c:c+2]
            c=c+2
            SPdata_dict[SpdlStatutNames[i]] = SPresult[columns].values.tolist()
        SPdata_dict['AllSpdl']=allSPresult.values.tolist()

        SPdata_dict_ActiveOnly={}
        for key in SPdata_dict:
            # Filter out rows with non-zero values
            filtered_rows = [row for row in SPdata_dict[key] if any(val != 0 for val in row) and np.any(~np.isnan(row))]
            # Update the key with filtered rows
            SPdata_dict_ActiveOnly[key] = filtered_rows


        filenameOut = f'{directory}{NrSubtype}_Spdl_{Cortex}_SpikeAct_perUnitAB.xlsx'
        with pd.ExcelWriter(filenameOut) as writer:
            # Iterate over each key-value pair in the dictionary
            for sheet_name, data_list in SPdata_dict.items():
                # Convert the list to a DataFrame
                df = pd.DataFrame(data_list, columns=['Before', 'After'])
                # Write the DataFrame to a separate sheet
                df.to_excel(writer, sheet_name=sheet_name, index=False)

        filenameOut = f'{directory}{NrSubtype}_Spdl_{Cortex}_SpikeAct_perUnit_Active_OnlyAB.xlsx'
        with pd.ExcelWriter(filenameOut) as writer:
            # Iterate over each key-value pair in the dictionary
            for sheet_name, data_list in SPdata_dict_ActiveOnly.items():
                # Convert the list to a DataFrame
                df = pd.DataFrame(data_list, columns=['Before', 'After'])
                # Write the DataFrame to a separate sheet
                df.to_excel(writer, sheet_name=sheet_name, index=False)
