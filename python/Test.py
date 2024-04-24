
# ## Load LFP and packages

from scipy import signal
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, Cursor
from scipy import fftpack
import pandas as pd
from pathlib import Path
import os
from IPython.display import display
from ipyfilechooser import FileChooser
from datetime import datetime
import shutil
from scipy.signal import find_peaks
from scipy.signal import chirp, find_peaks, peak_widths

# Perform analysis for each mouse

MiceList=['BlackLinesOK']
dpath0 = "//10.69.168.1/crnldata/waking/audrey_hay/L1imaging/AnalysedMarch2023/Gaelle/Baseline_recording_ABmodified/"

for micename in MiceList:
    
    dpath=Path(dpath0 + micename)
    # Load sleep score and Ca2+ time series numpy arrays
    nb_sessions = sum(1 for p in dpath.iterdir() if p.is_dir() and p.name.startswith("session"))    
    sessions = [folder.name for folder in dpath.iterdir() if folder.is_dir() and "session" in folder.name]
    print('Sessions :', sessions)
    print('dpath :', dpath)
    for session in sessions:  
        print('session :', session)
        folder_base = Path(dpath) / session / f'OpenEphys/'
        print('folder_base :', folder_base)
