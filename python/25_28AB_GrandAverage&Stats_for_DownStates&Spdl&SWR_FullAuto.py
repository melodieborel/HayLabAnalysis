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
from datetime import datetime
import shutil

########################################################################
        # SCRIPT 27AB_GrandAverages&Stats_for_Osc
########################################################################

# Specify the directory containing the Excel files
directory = "//10.69.168.1/crnldata/waking/audrey_hay/L1imaging/AnalysedMarch2023/Gaelle/Baseline_recording_ABmodified/AB_Analysis/OscillationsAnalysis_PerMouse_2024_06_13_16_33_03_913022"

# Get the current date and time
FolderNameSave=str(datetime.now())
FolderNameSave = FolderNameSave.replace(" ", "_").replace(".", "_").replace(":", "_")
destination_folder= f"//10.69.168.1/crnldata/waking/audrey_hay/L1imaging/AnalysedMarch2023/Gaelle/Baseline_recording_ABmodified/AB_Analysis/GlobalOscillationsAnalysis_{FolderNameSave}/"
os.makedirs(destination_folder)
folder_to_save=Path(destination_folder)

# Copy the script file to the destination folder
source_script = "C:/Users/Manip2/SCRIPTS/Code python audrey/code python aurelie/HayLabAnalysis/python/25_28AB_GrandAverage&Stats_for_DownStates&Spdl&SWR_FullAuto.py"
destination_file_path = f"{destination_folder}/25_28AB_GrandAverage&Stats_for_DownStates&Spdl&SWR_FullAuto.txt"
shutil.copy(source_script, destination_file_path)

NrSubtypeList=['All', 'L1']
CortexList=['PFC', 'S1']

OscList=['Spdl', 'SWR', 'DS']
OscillationList=['Spindles', 'SWR', 'DS']

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
            dfs4 = []
            df4=[]
            dfs4_per_sheet = {}
            combined_df=[]

            if NrSubtype=='L1':
                MiceList=['BlackLinesOK', 'BlueLinesOK', 'GreenDotsOK', 'GreenLinesOK', 'RedLinesOK']
            else:
                MiceList=['Purple', 'ThreeColDotsOK', 'ThreeBlueCrossesOK']

            nametofind=f'{OscillationList[o]}_{Cortex}_ABdetection_GlobalResultsAB'
            nametofind2=f'{OscillationList[o]}_{Cortex}_ABdetection_CalciumAvgResultsAB'
            nametofind3=f'{OscillationList[o]}_{Cortex}_ABdetection_SpikeAvgResultsAB'
            nametofind4=f'{OscillationList[o]}_{Cortex}_ABdetection_SpikeSumResultsAB'

            # Recursively traverse the directory structure
            for root, _, files in os.walk(directory):
                for filename in files:
                    # Check if the file is an Excel file and contains the specified name
                    if filename.endswith('.xlsx') and nametofind in filename:
                        if any(name in filename for name in MiceList):  
                            # Construct the full path to the file
                            filepath = os.path.join(root, filename)
                            # Read the Excel file into a dataframe and append it to the list
                            df = pd.read_excel(filepath, index_col=0)
                            dfs.append(df)
                            print(filename)
                    if filename.endswith('.xlsx') and nametofind2 in filename: 
                        if any(name in filename for name in MiceList): 
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
                            print(filename)
                    if filename.endswith('.xlsx') and nametofind3 in filename: 
                        if any(name in filename for name in MiceList): 
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
                            print(filename)
                    if filename.endswith('.xlsx') and nametofind4 in filename: 
                        if any(name in filename for name in MiceList): 
                            # Construct the full path to the file
                            filepath = os.path.join(root, filename)
                            # Read the Excel file into a dataframe and append it to the list
                            excel_data = pd.read_excel(filepath, sheet_name=None, index_col=0, header=None)           
                            for sheet_name, df4 in excel_data.items():
                                if sheet_name in dfs4_per_sheet:   
                                    updated_matrix = pd.concat((dfs4_per_sheet[sheet_name],df4), ignore_index=True, axis=0)                
                                    dfs4_per_sheet[sheet_name] = updated_matrix                    
                                else:                    
                                    dfs4_per_sheet[sheet_name] = df4 #one average trace per unique unit, len(df4)==nb unit recorded for that mouse
                            print(filename)

            # GLOBAL

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

            filenameOut = f'{folder_to_save}/{NrSubtype}_{Osc}_{Cortex}_ABdetection_GrandGlobalAB.xlsx'
            writer = pd.ExcelWriter(filenameOut)
            combined_df.to_excel(writer)
            writer.close()

            # CALCIUM traces dfs2_per_sheet

            filenameOut = f'{folder_to_save}/{NrSubtype}_{Osc}_{Cortex}_ABdetection_CalciumGrandAverageAB.xlsx'
            excel_writer = pd.ExcelWriter(filenameOut)

            Array=[]
            ArrayUn=[]
            ArrayPre=[]
            ArrayPost=[]
            
            Array=pd.DataFrame(dfs2_per_sheet[f'All_{OscillationList[o]}'])
            ArrayUn=pd.DataFrame(dfs2_per_sheet[f'Uncoupled_{OscillationList[o]}'])
            ArrayPre=pd.DataFrame(dfs2_per_sheet[f'Precoupled_{OscillationList[o]}'])
            ArrayPost=pd.DataFrame(dfs2_per_sheet[f'Postcoupled_{OscillationList[o]}'])
            
            Array.to_excel(excel_writer, sheet_name=f'All_{OscillationList[o]}', index=True, header=False)
            ArrayUn.to_excel(excel_writer, sheet_name=f'Uncoupled_{OscillationList[o]}', index=True, header=False)
            ArrayPre.to_excel(excel_writer, sheet_name=f'Precoupled_{OscillationList[o]}', index=True, header=False)
            ArrayPost.to_excel(excel_writer, sheet_name=f'Postcoupled_{OscillationList[o]}', index=True, header=False)

            if Osc=='Spdl':
                ArrayGlobal=[]
                ArrayLocal=[]
                ArrayGlobal=pd.DataFrame(dfs2_per_sheet[f'Global_Spindles'])
                ArrayLocal=pd.DataFrame(dfs2_per_sheet[f'Local_Spindles'])
                ArrayGlobal.to_excel(excel_writer, sheet_name=f'Global_Spindles', index=True, header=False)
                ArrayLocal.to_excel(excel_writer, sheet_name=f'Local_Spindles', index=True, header=False)

            excel_writer.close()

            # CALCIUM traces Normalization dfs2_per_sheet

            filenameOut = f'{folder_to_save}/{NrSubtype}_{Osc}_{Cortex}_ABdetection_Normalized_CalciumGrandAverageAB.xlsx'
            excel_writer = pd.ExcelWriter(filenameOut)

            row_sums = Array.sum(axis=1)
            Array = Array.div(row_sums, axis=0)
            row_sums = ArrayUn.sum(axis=1)
            ArrayUn = ArrayUn.div(row_sums, axis=0)
            row_sums = ArrayPre.sum(axis=1)
            ArrayPre = ArrayPre.div(row_sums, axis=0)
            row_sums = ArrayPost.sum(axis=1)
            ArrayPost = ArrayPost.div(row_sums, axis=0)

            Array.to_excel(excel_writer, sheet_name=f'All_{OscillationList[o]}', index=True, header=False)
            ArrayUn.to_excel(excel_writer, sheet_name=f'Uncoupled_{OscillationList[o]}', index=True, header=False)
            ArrayPre.to_excel(excel_writer, sheet_name=f'Precoupled_{OscillationList[o]}', index=True, header=False)
            ArrayPost.to_excel(excel_writer, sheet_name=f'Postcoupled_{OscillationList[o]}', index=True, header=False)

            if Osc=='Spdl':
                row_sums = ArrayGlobal.sum(axis=1)
                ArrayGlobal = ArrayGlobal.div(row_sums, axis=0)
                row_sums = ArrayLocal.sum(axis=1)
                ArrayLocal = ArrayLocal.div(row_sums, axis=0)
                ArrayGlobal.to_excel(excel_writer, sheet_name=f'Global_Spindles', index=True, header=False)
                ArrayLocal.to_excel(excel_writer, sheet_name=f'Local_Spindles', index=True, header=False)

            excel_writer.close()

            # SPIKE AVG dfs3_per_sheet

            filenameOut = f'{folder_to_save}/{NrSubtype}_{Osc}_{Cortex}_ABdetection_SpikeGrandAverageAB.xlsx'
            excel_writer = pd.ExcelWriter(filenameOut)

            Array=[]
            ArrayUn=[]
            ArrayPre=[]
            ArrayPost=[]

            Array=pd.DataFrame(dfs3_per_sheet[f'All_{OscillationList[o]}'])
            ArrayUn=pd.DataFrame(dfs3_per_sheet[f'Uncoupled_{OscillationList[o]}'])
            ArrayPre=pd.DataFrame(dfs3_per_sheet[f'Precoupled_{OscillationList[o]}'])
            ArrayPost=pd.DataFrame(dfs3_per_sheet[f'Postcoupled_{OscillationList[o]}'])

            Array.to_excel(excel_writer, sheet_name=f'All_{OscillationList[o]}', index=True, header=False)
            ArrayUn.to_excel(excel_writer, sheet_name=f'Uncoupled_{OscillationList[o]}', index=True, header=False)
            ArrayPre.to_excel(excel_writer, sheet_name=f'Precoupled_{OscillationList[o]}', index=True, header=False)
            ArrayPost.to_excel(excel_writer, sheet_name=f'Postcoupled_{OscillationList[o]}', index=True, header=False)
            
            if Osc=='Spdl':
                ArrayGlobal=[]
                ArrayLocal=[]
                ArrayGlobal=pd.DataFrame(dfs3_per_sheet[f'Global_Spindles'])
                ArrayLocal=pd.DataFrame(dfs3_per_sheet[f'Local_Spindles'])
                ArrayGlobal.to_excel(excel_writer, sheet_name=f'Global_Spindles', index=True, header=False)
                ArrayLocal.to_excel(excel_writer, sheet_name=f'Local_Spindles', index=True, header=False)

            excel_writer.close()
            
            # SPIKE SUM dfs4_per_sheet

            filenameOut = f'{folder_to_save}/{NrSubtype}_{Osc}_{Cortex}_ABdetection_SpikeGrandSumAB.xlsx'
            excel_writer = pd.ExcelWriter(filenameOut)

            Array=[]
            ArrayUn=[]
            ArrayPre=[]
            ArrayPost=[]

            Array=pd.DataFrame(dfs4_per_sheet[f'All_{OscillationList[o]}'])
            ArrayUn=pd.DataFrame(dfs4_per_sheet[f'Uncoupled_{OscillationList[o]}'])
            ArrayPre=pd.DataFrame(dfs4_per_sheet[f'Precoupled_{OscillationList[o]}'])
            ArrayPost=pd.DataFrame(dfs4_per_sheet[f'Postcoupled_{OscillationList[o]}'])

            Array.to_excel(excel_writer, sheet_name=f'All_{OscillationList[o]}', index=True, header=False)
            ArrayUn.to_excel(excel_writer, sheet_name=f'Uncoupled_{OscillationList[o]}', index=True, header=False)
            ArrayPre.to_excel(excel_writer, sheet_name=f'Precoupled_{OscillationList[o]}', index=True, header=False)
            ArrayPost.to_excel(excel_writer, sheet_name=f'Postcoupled_{OscillationList[o]}', index=True, header=False)

            if Osc=='Spdl':
                ArrayGlobal=[]
                ArrayLocal=[]
                ArrayGlobal=pd.DataFrame(dfs4_per_sheet[f'Global_Spindles'])
                ArrayLocal=pd.DataFrame(dfs4_per_sheet[f'Local_Spindles'])
                ArrayGlobal.to_excel(excel_writer, sheet_name=f'Global_Spindles', index=True, header=False)
                ArrayLocal.to_excel(excel_writer, sheet_name=f'Local_Spindles', index=True, header=False)

            excel_writer.close()

            # Spike SUM traces Normalization dfs4_per_sheet

            filenameOut = f'{folder_to_save}/{NrSubtype}_{Osc}_{Cortex}_ABdetection_Normalized_SpikeGrandSumAB.xlsx'
            excel_writer = pd.ExcelWriter(filenameOut)

            row_sums = Array.sum(axis=1)
            Array = Array.div(row_sums, axis=0)
            row_sums = ArrayUn.sum(axis=1)
            ArrayUn = ArrayUn.div(row_sums, axis=0)
            row_sums = ArrayPre.sum(axis=1)
            ArrayPre = ArrayPre.div(row_sums, axis=0)
            row_sums = ArrayPost.sum(axis=1)
            ArrayPost = ArrayPost.div(row_sums, axis=0)

            Array.to_excel(excel_writer, sheet_name=f'All_{OscillationList[o]}', index=True, header=False)
            ArrayUn.to_excel(excel_writer, sheet_name=f'Uncoupled_{OscillationList[o]}', index=True, header=False)
            ArrayPre.to_excel(excel_writer, sheet_name=f'Precoupled_{OscillationList[o]}', index=True, header=False)
            ArrayPost.to_excel(excel_writer, sheet_name=f'Postcoupled_{OscillationList[o]}', index=True, header=False)

            if Osc=='Spdl':
                row_sums = ArrayGlobal.sum(axis=1)
                ArrayGlobal = ArrayGlobal.div(row_sums, axis=0)
                row_sums = ArrayLocal.sum(axis=1)
                ArrayLocal = ArrayLocal.div(row_sums, axis=0)
                ArrayGlobal.to_excel(excel_writer, sheet_name=f'Global_Spindles', index=True, header=False)
                ArrayLocal.to_excel(excel_writer, sheet_name=f'Local_Spindles', index=True, header=False)

            excel_writer.close()