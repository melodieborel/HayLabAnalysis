import spikeinterface.full as si
import numpy as np
import pickle
import os
import shutil
import time

start_time = time.time()

duration_extract = 2 #min

with open('/crnldata/waking/audrey_hay/NPX/NPXprobe.pkl', 'rb') as outp: 
    probe = pickle.load(outp)
print(probe)


probe.set_device_channel_indices(np.arange(384))

file = "/mnt/data/ahay/NP_spikes_2024-07-22T17_55_16.raw"

raw_rec = si.read_binary(file, dtype='uint16', num_channels=384, sampling_frequency=30_000.)
raw_rec = raw_rec.set_probe(probe)
raw_rec = raw_rec.frame_slice(0, 30_000 * 60 * duration_extract)
print(raw_rec)



dirpath = os.path.join(os.getcwd(), 'kilosort4_output')
print(dirpath)
if os.path.exists(dirpath) and os.path.isdir(dirpath):
    print("exists")
    shutil.rmtree(dirpath)
    print(f'{dirpath} was deleted')

sorting = si.run_sorter('kilosort4', raw_rec, verbose=True)

print("job done in " + str(time.time()-start_time) + " seconds")
