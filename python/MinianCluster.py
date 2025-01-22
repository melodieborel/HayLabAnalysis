##################################
        # SETTING-UP #
##################################

suffix = ''
intSave=False

import itertools as itt
import os
import sys
import numpy as np
import xarray as xr
from dask.distributed import Client, LocalCluster


minian_path = os.path.join(os.path.abspath('..'),'minian')
print("The folder used for minian procedures is : {}".format(minian_path))

sys.path.append(minian_path)
from minian.cnmf import (
    compute_AtC,
    compute_trace,
    get_noise_fft,
    smooth_sig,
    unit_merge,
    update_spatial,
    update_temporal,
    update_background,
)
from minian.initialization import (
    gmm_refine,
    initA,
    initC,
    intensity_refine,
    ks_refine,
    pnr_refine,
    seeds_init,
    seeds_merge,
)
from minian.motion_correction import apply_transform, estimate_motion
from minian.preprocessing import denoise, remove_background
from minian.utilities import (
    TaskAnnotation,
    get_optimal_chk,
    load_videos,
    open_minian,
    save_minian,
)

#dpath='crnldata/forgetting/Aurelie/Test'
dpath='/mnt/data/minianAB/'


minian_ds_path = os.path.join(dpath, f'minian{suffix}')
intpath = os.path.join(dpath, f'minian_intermediate{suffix}')
print(minian_ds_path)


##################################
        # PARAMETERS #
##################################
# Initial parameters
subset = dict(frame=slice(0, None))
subset_mc = None
output_size = 100
n_workers = int(os.getenv("MINIAN_NWORKERS", 4))
param_save_minian = {
    "dpath": minian_ds_path,
    "meta_dict": dict(session=-1, animal=-2),
    "overwrite": True,
}

# Pre-processing Parameters#
param_load_videos = {
    "pattern": "[0-9]+\.avi$",
    "dtype": np.uint8,
    "downsample": dict(frame=1, height=1, width=1),
    "downsample_strategy": "subset",
}
param_denoise ={"method": "median", "ksize": 7} #{"method": "median", "ksize": 5} #Default minian = {"method": "median", "ksize": 7}
param_background_removal = {"method": "tophat", "wnd": 15}

# Motion Correction Parameters#
subset_mc = None
param_estimate_motion = {"dim": "frame"}

# Initialization Parameters#
param_seeds_init = {
    "wnd_size": 1000, # 100, #Default minian = 1000
    "method": "rolling",
    "stp_size": 50, #50, #Default minian =500
    "max_wnd": 15, #20,#generally 10 updated here to 20 to account for L1 wide dendritic trees #Default minian =15
    "diff_thres": 3, #3
}
param_pnr_refine = {"noise_freq": 0.06, "thres": 1}
param_ks_refine = {"sig": 0.05}
param_seeds_merge = {"thres_dist": 10, "thres_corr": 0.8, "noise_freq": 0.06}
param_initialize = {"thres_corr": 0.8, "wnd": 10, "noise_freq": 0.06} 
param_init_merge = {"thres_corr": 0.8}

# CNMF Parameters# 0.025 for threecolordots
param_get_noise = {"noise_range": (0.06, 0.5)}
param_first_spatial = {
    "dl_wnd": 15, #15, #Default minian = 10
    "sparse_penal": 0.012, #0.012, #Default minian =0.01
    "size_thres": (25, None),
}
param_first_temporal = {
    "noise_freq": 0.06,
    "sparse_penal": 1,
    "p": 1,
    "add_lag": 20,
    "jac_thres": 0.2,
}
param_first_merge = {"thres_corr": 0.8}
param_second_spatial = {
    "dl_wnd": 10,
    "sparse_penal": 0.01, #0.005, #Default minian =0.01
    "size_thres": (25, None),
}
param_second_temporal = {
    "noise_freq": 0.06,
    "sparse_penal": 1,
    "p": 1,
    "add_lag": 20,
    "jac_thres": 0.4,
}

##################################
            # CLUSTER #
##################################
"""
cluster = LocalCluster(
    n_workers=n_workers,
    memory_limit="8GB", #4
    resources={"MEM": 1},
    threads_per_worker=4,
    dashboard_address=":8780" #port 8787 already used by jupyter
)

annt_plugin = TaskAnnotation()
cluster.scheduler.add_plugin(annt_plugin)
client = Client(cluster)
print(cluster)
print(client)
"""
##################################
        # PRE- PROCESSING #
##################################

# PRE- PROCESSING
varr = load_videos(dpath, **param_load_videos)
chk, _ = get_optimal_chk(varr, dtype=np.float64)
if intSave: 
    varr = save_minian(
        varr.chunk({"frame": chk["frame"], "height": -1, "width": -1}).rename("varr"),
        intpath,
        overwrite=True,
    )
varr_ref = varr.sel(height=slice(0, 600), width=slice(0, 600))

# CLEAN UP - GLOW REMOVAL 
varr_min = varr_ref.min("frame").compute()
varr_ref = varr_ref - varr_min

# DENOISE & BACKGROUND REMOVAL
varr_ref = denoise(varr_ref, **param_denoise)
varr_ref = remove_background(varr_ref, **param_background_removal)
if intSave: 
    varr_ref = save_minian(varr_ref.rename("varr_ref"), dpath=intpath, overwrite=True)

##################################
        # MOTION CORRECTION #
##################################

motion = estimate_motion(varr_ref.sel(subset_mc), **param_estimate_motion)
if intSave: 
    motion = save_minian(
        motion.rename("motion").chunk({"frame": chk["frame"]}), **param_save_minian
    )
Y = apply_transform(varr_ref, motion, fill=0)

Y_fm_chk =Y.astype(float)
if intSave: 
    Y_fm_chk = save_minian(Y.astype(float).rename("Y_fm_chk"), intpath, overwrite=True)
    Y_hw_chk = save_minian(
        Y_fm_chk.rename("Y_hw_chk"),  
        intpath,
        overwrite=True,
        chunks={"frame": -1, "height": chk["height"], "width": chk["width"]},
    )
vid_arr = xr.concat([varr_ref, Y_fm_chk], "width").chunk({"width": -1})
max_proj = save_minian(
    Y_fm_chk.max("frame").rename("max_proj"), **param_save_minian
).compute()

##################################
        # INTIALISATION #
##################################

# GENERATE SEEDS 
seeds = seeds_init(Y_fm_chk, **param_seeds_init)
seeds, pnr, gmm = pnr_refine(Y_hw_chk, seeds, **param_pnr_refine)
seeds = ks_refine(Y_hw_chk, seeds, **param_ks_refine)

# MERGE SEEDS
seeds_final = seeds[seeds["mask_ks"] & seeds["mask_pnr"]].reset_index(drop=True)
seeds_final = seeds_merge(Y_hw_chk, max_proj, seeds_final, **param_seeds_merge)
print("{} units found".format(seeds_final["mask_mrg"].count()))
print("{} units merged found".format(len(seeds_final[seeds_final["mask_mrg"]])))


# INITIAISE SPATIAL MATRIX
A_init = initA(Y_hw_chk, seeds_final[seeds_final["mask_mrg"]], **param_initialize)
if intSave: 
    A_init = save_minian(A_init.rename("A_init"), intpath, overwrite=True)

# INITIALISE TEMPORAL MATRIX
C_init = initC(Y_fm_chk, A_init)
if intSave: 
    C_init = save_minian(
        C_init.rename("C_init"), intpath, overwrite=True, chunks={"unit_id": 1, "frame": -1}
    )

# MERGE UNITS
A_merged, C_merged = unit_merge(A_init, C_init, **param_init_merge)
if intSave: 
    A_merged = save_minian(A_merged.rename("A"), intpath, overwrite=True)
    C_merged = save_minian(C_merged.rename("C"), intpath, overwrite=True)
    C_chk_merged = save_minian(
        C_merged.rename("C_chk"),
        intpath,
        overwrite=True,
        chunks={"unit_id": -1, "frame": chk["frame"]},
    )

# INITIALIZE BACKGROUND TERMS
b_init, f_init = update_background(Y_fm_chk, A_merged, C_chk_merged)
if intSave: 
    f_init = save_minian(f_init.rename("f"), intpath, overwrite=True)
    b_init = save_minian(b_init.rename("b"), intpath, overwrite=True)

##################################
            # CNMF #
##################################

# ESTIMATE SPATIAL NOISE
sn_spatial = get_noise_fft(Y_hw_chk, **param_get_noise)
if intSave: 
    sn_spatial = save_minian(sn_spatial.rename("sn_spatial"), intpath, overwrite=True)

# FIRST SPATIAL UPDATE
A_firstS, mask_firstS, norm_fac_firstS = update_spatial(
    Y_hw_chk, A_merged, C_merged, sn_spatial, **param_first_spatial
)
C_firstS = save_minian(
    (C_merged.sel(unit_id=mask_firstS) * norm_fac_firstS).rename("C_new"), intpath, overwrite=True
)
C_chk_firstS = save_minian(
    (C_chk_merged.sel(unit_id=mask_firstS) * norm_fac_firstS).rename("C_chk_new"), intpath, overwrite=True
)
b_firstS, f_firstS = update_background(Y_fm_chk, A_firstS, C_chk_firstS)

if intSave: 
    A_firstS = save_minian(
        A_firstS.rename("A"),
        intpath,
        overwrite=True,
        chunks={"unit_id": 1, "height": -1, "width": -1},
    )
    b_firstS = save_minian(b_firstS.rename("b"), intpath, overwrite=True)
    f_firstS = save_minian(
        f_firstS.chunk({"frame": chk["frame"]}).rename("f"), intpath, overwrite=True
    )
    C_firstS = save_minian(C_firstS.rename("C"), intpath, overwrite=True)
    C_chk_firstS = save_minian(C_chk_firstS.rename("C_chk"), intpath, overwrite=True)

# FIRST TEMPORAL UPDATE
YrA_firstT = save_minian(
    compute_trace(Y_fm_chk, A_firstS, b_firstS, C_chk_firstS, f_firstS).rename("YrA"),
    intpath,
    overwrite=True,
    chunks={"unit_id": 1, "frame": -1},
)
C_firstT, S_firstT, b0_firstT, c0_firstT, g_firstT, mask_firstT = update_temporal(
    A_firstS, C_firstS, YrA=YrA_firstT, **param_first_temporal
)
if intSave: 
    C_firstT = save_minian(
        C_firstT.rename("C").chunk({"unit_id": 1, "frame": -1}), intpath, overwrite=True
    )
    C_chk_firstT = save_minian(
        C_firstT.rename("C_chk"),
        intpath,
        overwrite=True,
        chunks={"unit_id": -1, "frame": chk["frame"]},
    )
    S_firstT = save_minian(
        S_firstT.rename("S").chunk({"unit_id": 1, "frame": -1}), intpath, overwrite=True
    )
    b0_firstT = save_minian(
        b0_firstT.rename("b0").chunk({"unit_id": 1, "frame": -1}), intpath, overwrite=True
    )
    c0_firstT = save_minian(
        c0_firstT.rename("c0").chunk({"unit_id": 1, "frame": -1}), intpath, overwrite=True
    )
A_firstT = A_firstS.sel(unit_id=C_firstT.coords["unit_id"].values)

# MERGE UNITS
A_mrg, C_mrg, [sig_mrg] = unit_merge(A_firstT, C_firstT, [C_firstT + b0_firstT + c0_firstT], **param_first_merge)

if intSave: 
    A_mrg = save_minian(A_mrg.rename("A_mrg"), intpath, overwrite=True)
    C_mrg = save_minian(C_mrg.rename("C_mrg"), intpath, overwrite=True)

C_chk_mrg = save_minian(
    C_mrg.rename("C_mrg_chk"),
    intpath,
    overwrite=True,
    chunks={"unit_id": -1, "frame": chk["frame"]},
)
sig_mrg = save_minian(sig_mrg.rename("sig_mrg"), intpath, overwrite=True)

# SECOND SPATIAL UPDATE
A_secS, mask_secS, norm_fac_secS = update_spatial(
    Y_hw_chk, A_mrg, C_mrg, sn_spatial, **param_second_spatial
)
C_secS = save_minian(
    (C_mrg.sel(unit_id=mask_secS) * norm_fac_secS).rename("C_new"), intpath, overwrite=True
)
C_chk_secS = save_minian(
    (C_chk_mrg.sel(unit_id=mask_secS) * norm_fac_secS).rename("C_chk_new"), intpath, overwrite=True
)

b_secS, f_secS = update_background(Y_fm_chk, A_secS, C_chk_secS)

if intSave: 
    A_secS = save_minian(
        A_secS.rename("A"),
        intpath,
        overwrite=True,
        chunks={"unit_id": 1, "height": -1, "width": -1},
    )
    b_secS = save_minian(b_secS.rename("b"), intpath, overwrite=True)
    f_secS = save_minian(
        f_secS.chunk({"frame": chk["frame"]}).rename("f"), intpath, overwrite=True
    )
    C_secS = save_minian(C_secS.rename("C"), intpath, overwrite=True)
    C_chk_secS = save_minian(C_chk_secS.rename("C_chk"), intpath, overwrite=True)


# SECOND TEMPORAL UPDATE
YrA_secT = save_minian(
    compute_trace(Y_fm_chk, A_secS, b_secS, C_chk_secS, f_secS).rename("YrA"),
    intpath,
    overwrite=True,
    chunks={"unit_id": 1, "frame": -1},
)
C_secT, S_secT, b0_secT, c0_secT, g_secT, mask_secT = update_temporal(
    A_secS, C_secS, YrA=YrA_secT, **param_second_temporal
)
if intSave: 
    C_secT = save_minian(
        C_secT.rename("C").chunk({"unit_id": 1, "frame": -1}), intpath, overwrite=True
    )
    C_chk_secT = save_minian(
        C_secT.rename("C_chk"),
        intpath,
        overwrite=True,
        chunks={"unit_id": -1, "frame": chk["frame"]},
    )
    S_secT = save_minian(
        S_secT.rename("S").chunk({"unit_id": 1, "frame": -1}), intpath, overwrite=True
    )
    b0_secT = save_minian(
        b0_secT.rename("b0").chunk({"unit_id": 1, "frame": -1}), intpath, overwrite=True
    )
    c0_secT = save_minian(
        c0_secT.rename("c0").chunk({"unit_id": 1, "frame": -1}), intpath, overwrite=True
    )
    A_secT = A_secS.sel(unit_id=C_secT.coords["unit_id"].values)


# SAVE FINAL RESULTS
A = save_minian(A_secT.rename("A"), **param_save_minian)
C = save_minian(C_secT.rename("C"), **param_save_minian)
S = save_minian(S_secT.rename("S"), **param_save_minian)
c0 = save_minian(c0_secT.rename("c0"), **param_save_minian)
b0 = save_minian(b0_secT.rename("b0"), **param_save_minian)
b = save_minian(b_init.rename("b"), **param_save_minian)
f = save_minian(f_init.rename("f"), **param_save_minian)

#DELETE INTERMEDIATE PATH
os.rmdir(intpath)

"""
# Close cluster
client.close()
cluster.close()
"""