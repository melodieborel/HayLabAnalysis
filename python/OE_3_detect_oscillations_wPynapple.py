#######################################################################################
                            # Define Directory #
#######################################################################################

dir = "//10.69.168.1/crnldata/waking/audrey_hay/L1imaging/Analysed2025_AB/"

#######################################################################################
                                # Load Packages #
#######################################################################################

import os
import numpy as np
import pandas as pd
from pathlib import Path
import os
import numpy as np
import pynapple as nap
import pandas as pd
from pathlib import Path
import os
import warnings
from scipy.stats import zscore
warnings.filterwarnings("ignore")

#######################################################################################
                                # Define functions #
#######################################################################################
def restriction_parameter(All_Spindle):
    nb_spindle = All_Spindle.shape[0]
    listtodrop = []
    for tt in range(nb_spindle-1):
        # merge spdl that starts within a spdl
        if(All_Spindle['end time'][tt]>All_Spindle['start time'][tt + 1]):
            if(All_Spindle['Duration'][tt]<All_Spindle['Duration'][tt + 1]):
                if(All_Spindle['start time'][tt]<All_Spindle['start time'][tt + 1]):
                    All_Spindle['start time'][tt+1] = All_Spindle['start time'][tt]
                    listtodrop.append(tt)
                else:
                    listtodrop.append(tt)
            if(All_Spindle['Duration'][tt]>All_Spindle['Duration'][tt + 1]):
                if(All_Spindle['end time'][tt]<All_Spindle['end time'][tt + 1]):
                    All_Spindle['end time'][tt] = All_Spindle['end time'][tt + 1]
                    listtodrop.append(tt+1)
                else:
                    listtodrop.append(tt+1)
    """
    for tt in range(nb_spindle-1):
        # merge spdls that are 200ms apart
        if((All_Spindle['start time'][tt + 1] - All_Spindle['end time'][tt])<200):
            if((All_Spindle['Duration'][tt])<All_Spindle['Duration'][tt + 1]): #first spdl longer so remove/merge the second one
                All_Spindle['start time'][tt + 1] = min(All_Spindle['start time'][tt], All_Spindle['start time'][tt+1])
                All_Spindle['end time'][tt+ 1] = max(All_Spindle['end time'][tt + 1], All_Spindle['end time'][tt])
                listtodrop.append(tt)
            if((All_Spindle['Duration'][tt+1])<All_Spindle['Duration'][tt]): #second spdl longer so remove/merge the first one
                All_Spindle['start time'][tt] = min(All_Spindle['start time'][tt], All_Spindle['start time'][tt+1])
                All_Spindle['end time'][tt] = max(All_Spindle['end time'][tt + 1], All_Spindle['end time'][tt])
                listtodrop.append(tt+1)
    """
    for tt in range(nb_spindle):
        #Update duration because of the merging
        All_Spindle['Duration'][tt]=All_Spindle['end time'][tt]-All_Spindle['start time'][tt]

    for tt in range(nb_spindle): #All_Spindle.index:
        #Remove Spdl that last less than 500ms
        if (All_Spindle['Duration'][tt]<500):
            listtodrop.append(tt)        
    
    All_Spindle = All_Spindle.drop(listtodrop) 
    All_Spindle = All_Spindle.reset_index(drop=True)
    return All_Spindle


#######################################################################################
                                # Load Signals #
#######################################################################################

#Load LFP coordinates 
Channels = f'{os.getcwd()}\python\_LFP_coordinates_of_all_mice.csv'
all_LFPcoordinates = pd.read_csv(Channels, index_col=0)

for dpath in Path(dir).glob('**/DataFrame_rawdataDS.pkl'):
    

    folder_base = Path(dpath).parent
    print(folder_base)        
    
    isSWRfile=list(Path(folder_base).glob('**/SWRproperties.csv'))
    isSpdlfile=list(Path(folder_base).glob('**/Spindleproperties_S1.csv'))

    if  len(isSpdlfile)==0 : #len(isSWRfile)==0 and
                
        #Load signals

        LFPfile = Path(f'{folder_base}\DataFrame_rawdataDS.pkl')
        LFPs_df = pd.read_pickle(LFPfile)
        samplerate = 1000 
        numchannel = LFPs_df.shape[1]
        rec_ch_list = LFPs_df.columns.values
        # Load LFPs timestamps 
        for file_pathTS in folder_base.parent.parent.glob('**/continuous/*/timeStampsDS.npy'):
            print('LFPs timestamps file found')
            LFPtimestamps = np.load(file_pathTS)  
        print(round(LFPs_df.shape[0]/samplerate/60), 'min of recording')

        # Identify mouse & choose threshold for detection
        mouse = []
        pos_mice = []
        for mouse_name in all_LFPcoordinates.index:
            if mouse_name in LFPfile.__str__():
                mouse.append(mouse_name)
                pos_mice.append(LFPfile.__str__().find(mouse_name)) 
        mouse = [x for _, x in sorted(zip(pos_mice, mouse))] # sort mouse in the same order as they appear in the path
        mouse=mouse[0]   

        # Identify electrodes & create differential LFPs
        all_LFPcoordinates = all_LFPcoordinates.astype(str)
        for region in all_LFPcoordinates.loc[mouse].index:
            locals()[region] = []
            locals()[f'{region}_ch'] = []
        ID=0
        rec_ch_list_mouse = [value for value in rec_ch_list if 0+(ID*32) <= value <= 31+(ID*32)]
        for rec_ch in rec_ch_list_mouse:
            for idx, LFPcoord_str in enumerate(all_LFPcoordinates.loc[mouse]):
                region = all_LFPcoordinates.loc[mouse].index[idx]
                if LFPcoord_str != 'nan':
                    LFPcoord = LFPcoord_str.split('_')[:2] # only take into account the 2 first of electrode of that region 
                    num_ch = np.where(str(rec_ch-(ID*32)) == np.array(LFPcoord))[0]
                    if len(num_ch) > 0:
                        region = all_LFPcoordinates.loc[mouse].index[idx]
                        LFP = locals()[region]
                        LFP = LFP-np.array(LFPs_df[(rec_ch)]) if len(LFP) > 0 else np.array(LFPs_df[(rec_ch)])
                        locals()[region] = LFP
                        locals()[f'{region}_ch'].append(rec_ch)
                        break
                    continue    
        for region in all_LFPcoordinates.loc[mouse].index:
            LFP = locals()[region]
            LFP_ch = locals()[f'{region}_ch']


        
        # Load Sleep Scoring    
        ScoringFile= folder_base / f'Sleep_Scoring_6Stages_5sEpoch.csv'
        SleepScored = pd.read_csv(ScoringFile)
        SleepScored['label']= SleepScored['label'].str.extract(r'(\d+)', expand=False)
        SleepScoredTS=np.array(SleepScored['label'])
        scale_factor=samplerate/0.2  #cause scoring was done in 5 seconds bin, ie 0.2 Hz              
        SleepScoredTS_upscaled0 = np.repeat(SleepScoredTS, scale_factor, axis=0)
        SleepScoredTS_upscaled0=SleepScoredTS_upscaled0.astype(int)
        SleepScoredTS_upscaled=SleepScoredTS_upscaled0.copy()

        SleepScoredTS_upscaled[SleepScoredTS_upscaled==2]=1 # Transform QW in AW
        SleepScoredTS_upscaled[SleepScoredTS_upscaled==4]=1 # Transform REM in AW
        SleepScoredTS_upscaled[SleepScoredTS_upscaled==5]=1 # Transform IS in AW
        SleepScoredTS_upscaled[SleepScoredTS_upscaled==6]=1 # Transform undefined in AW

        if len(CA1) > len(SleepScoredTS_upscaled):
            SleepScoredTS_upscaled = np.pad(SleepScoredTS_upscaled, (0, (len(CA1) - len(SleepScoredTS_upscaled))), constant_values=SleepScoredTS_upscaled[-1])
        elif len(CA1) < len(SleepScoredTS_upscaled):
            SleepScoredTS_upscaled=SleepScoredTS_upscaled[:len(CA1)]

        notNREMBool = SleepScoredTS_upscaled==1 


        CA1=zscore(np.array(CA1))#-np.mean(np.array(CA1)))#/np.mean(np.array(CA1))
        S1=zscore(np.array(S1))#-np.mean(np.array(S1)))#/np.mean(np.array(S1))
        PFC=zscore(np.array(PFC))#-np.mean(np.array(PFC)))#/np.mean(np.array(PFC))

        CA1[notNREMBool] = 0
        S1[notNREMBool] = 0
        PFC[notNREMBool] = 0
        
        ###########################################
                # SWR in CA1: 120-200 Hz #
        ###########################################

        lfp = nap.Tsd(t=np.arange(len(CA1))/samplerate, d=CA1)
        freqs = np.geomspace(3, 200, 100)
        mwt_RUN = nap.compute_wavelet_transform(lfp, fs=samplerate, freqs=freqs, gaussian_width=4, window_length=1) # gw=4, wl=1 only very pretty SWRs
        ripple_freq_index = np.logical_and(freqs > 120, freqs < 200)
        ripple_power = np.mean(np.abs(mwt_RUN[:, ripple_freq_index]), 1)
        smoothed_ripple_power = ripple_power.smooth(0.01) #in seconds 0.005 
        threshold_ripple_power = smoothed_ripple_power.threshold(.15) #.25 only very pretty SWRs // .1 too much for TBC
        rip_ep = threshold_ripple_power.time_support
        rip_ep['dur_ms']=np.round((rip_ep['end']-rip_ep['start'])*1000)
        rip_ep=rip_ep[rip_ep['dur_ms']>10] #5ms

        emp=[]
        All_SWR = pd.DataFrame(emp, columns = ['start time', 'end time', 'Duration'])
        All_SWR['start time']=rip_ep['start']*1000
        All_SWR['end time']=rip_ep['end']*1000
        All_SWR['Duration']=rip_ep['dur_ms']
        All_SWR['toKeep']= 'True'

        # Store the results in All_SWR_prop pd dataframe and save as pkl/csv for post processing.
        filename = folder_base / f'SWR_detection.csv'
        All_SWR.to_csv(filename)

        print(len(All_SWR), 'SWR detected in CA1')

        ########################################
                    # Spdl: 9-18 Hz #
        ########################################
        
        #####################################
            ##         PFC         ##

        lfp = nap.Tsd(t=np.arange(len(PFC))/samplerate, d=PFC)
        freqs = np.geomspace(3, 200, 100)
        mwt_RUN = nap.compute_wavelet_transform(lfp, fs=samplerate, freqs=freqs, gaussian_width=3, window_length=1.5)
        spdl_freq_index = np.logical_and(freqs > 11, freqs < 17)
        spdl_power = np.mean(np.abs(mwt_RUN[:, spdl_freq_index]), 1)
        smoothed_spdl_power = spdl_power.smooth(0.1) #in seconds 0.005 
        threshold_spdl_power = smoothed_spdl_power.threshold(.4)
        spdl_ep = threshold_spdl_power.time_support
        spdl_ep['dur_ms']=np.round((spdl_ep['end']-spdl_ep['start'])*1000)
        #spdl_ep=spdl_ep[spdl_ep['dur_ms']>500]

        emp=[]
        All_SpindlePFC = pd.DataFrame(emp, columns = ['start time', 'end time', 'Duration'])
        All_SpindlePFC['start time']=spdl_ep['start']*1000
        All_SpindlePFC['end time']=spdl_ep['end']*1000
        All_SpindlePFC['Duration']=spdl_ep['dur_ms']
        All_SpindlePFC['toKeep']= 'True'

        filename = folder_base / f'SpindlesPFC_detection.csv'
        All_SpindlePFC.to_csv(filename)

        print(len(All_SpindlePFC), 'Spdl detected in PFC')

        #####################################
                ##         S1         ##

        lfp = nap.Tsd(t=np.arange(len(S1))/samplerate, d=S1)
        freqs = np.geomspace(3, 200, 100)
        mwt_RUN = nap.compute_wavelet_transform(lfp, fs=samplerate, freqs=freqs, gaussian_width=3, window_length=1.5)
        spdl_freq_index = np.logical_and(freqs > 11, freqs < 17)
        spdl_power = np.mean(np.abs(mwt_RUN[:, spdl_freq_index]), 1)
        smoothed_spdl_power = spdl_power.smooth(0.1) #in seconds 0.005 
        threshold_spdl_power = smoothed_spdl_power.threshold(.4)
        spdl_ep = threshold_spdl_power.time_support
        spdl_ep['dur_ms']=np.round((spdl_ep['end']-spdl_ep['start'])*1000)
        #spdl_ep=spdl_ep[spdl_ep['dur_ms']>500]

        emp=[]
        All_SpindleS1 = pd.DataFrame(emp, columns = ['start time', 'end time', 'Duration'])
        All_SpindleS1['start time']=spdl_ep['start']*1000
        All_SpindleS1['end time']=spdl_ep['end']*1000
        All_SpindleS1['Duration']=spdl_ep['dur_ms']
        All_SpindleS1['toKeep']= 'True'

        filename = folder_base / f'SpindlesS1_detection.csv'
        All_SpindleS1.to_csv(filename)
        
        print(len(All_SpindleS1), 'Spdl detected in S1')

    else: 
        print('SWR & Spdl already detected')

    ########################################
                # Merge Spdl #
    ########################################
    
    All_SpindleS1 = pd.read_csv(Path(f'{folder_base}\SpindlesS1_detection.csv'))
    All_SpindlePFC = pd.read_csv(Path(f'{folder_base}\SpindlesPFC_detection.csv'))

    All_SpindlePFC['CTX']='PFC'
    All_SpindlePFC['StartingLoc']='PFC'
    All_SpindlePFC = All_SpindlePFC.reset_index(drop=True)
    NewPFClist=restriction_parameter(All_SpindlePFC)

    All_SpindleS1['CTX']='S1'
    All_SpindleS1['StartingLoc']='S1'
    All_SpindleS1 = All_SpindleS1.reset_index(drop=True)
    NewS1list=restriction_parameter(All_SpindleS1)

    #Spdllist=pd.concat([All_SpindlePFC,All_SpindleS1],ignore_index=False)  
    Spdllist=pd.concat([NewPFClist,NewS1list],ignore_index=False)  

    Spdllist['LocalGlobal']='Local'
    Spdllist['DistanceClosestSpdl']=np.inf

    Spdllist = Spdllist.sort_values(by='start time')
    Spdllist = Spdllist.reset_index(drop=True)

    liststarts=Spdllist["start time"]
    listends=Spdllist["end time"]
    NewSpdllist=Spdllist.copy()
    NewSpdllist['toKeep'] = NewSpdllist['toKeep'].astype(str)

    for spdl1 in NewSpdllist.index : # range(len(Spdllist)) :
        start1=liststarts[spdl1]
        end1=listends[spdl1]
        otherunit_range = [x for x in NewSpdllist.index  if x != spdl1]
        #print(f'spdl1 n°{spdl1}')
        if NewSpdllist.loc[spdl1,'toKeep'] != 'False':
            for spdl2 in otherunit_range:
                if NewSpdllist.loc[spdl2,'toKeep'] != 'False':
                    start2=liststarts[spdl2]
                    end2=listends[spdl2]
                    if NewSpdllist.loc[spdl1,'CTX']!= NewSpdllist.loc[spdl2,'CTX']: # Needs to be from 2 differents brain areas

                        if start2>start1:
                            NewSpdllist.loc[spdl1, 'DistanceClosestSpdl'] = start1-start2 if start2-start1<NewSpdllist.loc[spdl1, 'DistanceClosestSpdl'] else NewSpdllist.loc[spdl1, 'DistanceClosestSpdl']
                            NewSpdllist.loc[spdl2, 'DistanceClosestSpdl'] = start2-start1 if start2-start1<NewSpdllist.loc[spdl2, 'DistanceClosestSpdl'] else NewSpdllist.loc[spdl2, 'DistanceClosestSpdl']
                        else:
                            NewSpdllist.loc[spdl1, 'DistanceClosestSpdl'] = start1-start2 if start1-start2<NewSpdllist.loc[spdl1, 'DistanceClosestSpdl'] else NewSpdllist.loc[spdl1, 'DistanceClosestSpdl']
                            NewSpdllist.loc[spdl2, 'DistanceClosestSpdl'] = start2-start1 if start1-start2<NewSpdllist.loc[spdl2, 'DistanceClosestSpdl'] else NewSpdllist.loc[spdl2, 'DistanceClosestSpdl']

                        if start1<=start2 and start2<=end1: # event n°2 begins after the start n°1               
                            if (end1-start2)>=int(0.5*(end1-start1)): # overlapp > to 50% of the duration of the event n°1                                
                                NewSpdllist.loc[spdl1, 'LocalGlobal']='Global'
                                NewSpdllist.loc[spdl1, 'StartingLoc']=NewSpdllist.loc[spdl1,'StartingLoc']
                                NewSpdllist.loc[spdl1, 'CTX']='S1PFC'
                                NewSpdllist.loc[spdl1, 'start time']=int((start1 + start2)/2)
                                NewSpdllist.loc[spdl1, 'end time']=max(end1, end2)
                                NewSpdllist.loc[spdl1, 'Duration']=max(end1, end2)-int((start1 + start2)/2)
                                NewSpdllist.loc[spdl2, 'toKeep']='False'
                                #print(f'Cdt n°1: Global, keep spdl n°{spdl1} {start1} instead of spdl n°{spdl2} {start2}')
                                break
                        elif start1<=end2 and end2<=end1: # event n°2 ends before the end n°1 
                            if (end2-start1)>=int(0.5*(end1-start1)): # overlapp > to 50% of the duration of the event n°1
                                NewSpdllist.loc[spdl2, 'LocalGlobal']='Global'
                                NewSpdllist.loc[spdl2, 'StartingLoc']=NewSpdllist.loc[spdl2,'StartingLoc']
                                NewSpdllist.loc[spdl2, 'CTX']='S1PFC'
                                NewSpdllist.loc[spdl2, 'start time']=int((start1 + start2)/2)
                                NewSpdllist.loc[spdl2, 'end time']=max(end1, end2)
                                NewSpdllist.loc[spdl2, 'Duration']=max(end1, end2)-int((start1 + start2)/2)
                                NewSpdllist.loc[spdl1, 'toKeep']='False'
                                #print(f'Cdt n°2: Global, keep spdl n°{spdl2} {start2}, instead of spdl n°{spdl1} {start1}')
                                break
    NewSpdllist = NewSpdllist[NewSpdllist['toKeep'].isin(['True', 'VRAI'])]
    NewSpdllist = NewSpdllist.sort_values(by='start time')
    NewSpdllist = NewSpdllist.reset_index(drop=True)
    NewSpdllist['DistanceClosestSpdl'] = NewSpdllist['DistanceClosestSpdl'] *-1 # to have in positive the spdl that arrives after the onset and vice versa
    NewSpdllist=restriction_parameter(NewSpdllist)

    filenameOutput = folder_base / f'SpindlesS1&PFC_detection.csv' 
    NewSpdllist.to_csv(filenameOutput, sep= ',')