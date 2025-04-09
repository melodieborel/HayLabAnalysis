##################################
        # SETTING-UP #
##################################

suffix = ''

import itertools as itt
import os
import sys

import holoviews as hv
import numpy as np
import xarray as xr
from dask.distributed import Client, LocalCluster
import shutil
import time
from dask import config

st = time.time()

##################################
        # PARAMETERS #
##################################
# Set up Initial Basic Parameters#
minian_path = "/home/clementine.robein/minian"

dpath='/mnt/data/ClemR_minian/'
mouse_name = [f for f in os.listdir(dpath) if os.path.isdir(os.path.join(dpath, f))]
mouse_name = mouse_name[0]
path_mouse = os.path.join(dpath, mouse_name)

session_name = [f for f in os.listdir(path_mouse) if os.path.isdir(os.path.join(path_mouse, f))]
session_name = session_name[0]
dpath = os.path.join(path_mouse, session_name)

minian_ds_path = os.path.join(dpath, f'minian{suffix}')
intpath = os.path.join(dpath, f'minian_intermediate{suffix}')
subset = dict(frame=slice(0, None))
subset_mc = None
interactive = False
output_size = 100
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
param_denoise = {"method": "median", "ksize": 5} #Default minian ={"method": "median", "ksize": 7}
param_background_removal = {"method": "tophat", "wnd": 15}

# Motion Correction Parameters#
subset_mc = None
param_estimate_motion = {"dim": "frame"}

# Initialization Parameters#
param_seeds_init = {
    "wnd_size": 100, # 100, #Default minian = 1000
    "method": "rolling",
    "stp_size": 50, #50, #Default minian =500
    "max_wnd": 20, #20,#generally 10 updated here to 20 to account for L1 wide dendritic trees #Default minian = 15 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    "diff_thres": 3, #3
}
param_pnr_refine = {"noise_freq": 0.06, "thres": 1}
param_ks_refine = {"sig": 0.05}
param_seeds_merge = {"thres_dist": 10, "thres_corr": 0.8, "noise_freq": 0.06}
param_initialize = {"thres_corr": 0.8, "wnd": 10, "noise_freq": 0.06}
param_init_merge = {"thres_corr": 0.8}

# CNMF Parameters#
param_get_noise = {"noise_range": (0.06, 0.5)}
param_first_spatial = {
    "dl_wnd": 20, #15, #Default minian = 10 #the window size of the morphological dilation operation
    "sparse_penal": 0.0015, #0.012, #Default minian =0.01 # the bigger, the smaller the ROI !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    "size_thres": (25, None), # range of area (number of non-zero pixels) of the spatial footprints that will be accepted
}
param_first_temporal = {
    "noise_freq": 0.06,
    "sparse_penal": 1,
    "p": 1,
    "add_lag": 20,
    "jac_thres": 0.2,
}
param_first_merge = {"thres_corr": 0.8}
param_second_spatial = param_first_spatial
param_second_temporal = param_first_temporal
""""
param_second_spatial = {
    "dl_wnd": 20,
    "sparse_penal": 0.002,
    "size_thres": (25, None),
}"

param_second_temporal = {
    "noise_freq": 0.06,
    "sparse_penal": 1,
    "p": 1,
    "add_lag": 20,
    "jac_thres": 0.4,
}
"""
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MINIAN_INTERMEDIATE"] = intpath

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
from minian.visualization import (
    CNMFViewer,
    VArrayViewer,
    generate_videos,
    visualize_gmm_fit,
    visualize_motion,
    visualize_preprocess,
    visualize_seeds,
    visualize_spatial_update,
    visualize_temporal_update,
    write_video,
)

##################################
            # CLUSTER #
##################################

if __name__ == "__main__": # needed if dask client runned into a .py script    

    dpath = os.path.abspath(dpath)
    
    cluster = LocalCluster(
        n_workers=int(os.getenv("MINIAN_NWORKERS", 40)), # /!\ max 40 or 64 CPUs per node in remote machine # /!\ 8 total cores in local machine 
        memory_limit="6GB", #per worker, /!\ max 95 or 256 GB per node in remote machine # /!\ 32GB total RAM in local machine 
        resources={"MEM": 1}, #set to 1 before
        threads_per_worker=2,
        dashboard_address=None,
        #processes=False, # to avoid distributed.nanny - WARNING - Restarting worker ?
    )

    config.set({'interface': 'lo'}) 
    annt_plugin = TaskAnnotation()
    cluster.scheduler.add_plugin(annt_plugin)
    client = Client(cluster) 

    ##################################
            # PRE- PROCESSING #
    ##################################
    # PRE- PROCESSING
    varr = load_videos(dpath, **param_load_videos)
    chk, _ = get_optimal_chk(varr, dtype=float)

    varr = save_minian(
        varr.chunk({"frame": chk["frame"], "height": -1, "width": -1}).rename("varr"),
        intpath,
        overwrite=True,
    )

    # SELECT A SUBSET OF THE VIDEO IF NEEDED
    #varr_ref = varr.sel(height=slice(120, 600))
    varr_ref = varr.sel(subset)

    # CLEAN UP - GLOW REMOVAL 
    varr_min = varr_ref.min("frame").compute()
    varr_ref = varr_ref - varr_min

    # DENOISE & BACKGROUND REMOVAL
    varr_ref = denoise(varr_ref, **param_denoise)

    varr_ref = remove_background(varr_ref, **param_background_removal)

    varr_ref = save_minian(varr_ref.rename("varr_ref"), dpath=intpath, overwrite=True)

    ##################################
            # MOTION CORRECTION #
    ##################################
    motion = estimate_motion(varr_ref.sel(subset_mc), **param_estimate_motion)

    motion = save_minian(
        motion.rename("motion").chunk({"frame": chk["frame"]}), **param_save_minian
    )

    Y = apply_transform(varr_ref, motion, fill=0)

    Y_fm_chk = save_minian(Y.astype(float).rename("Y_fm_chk"), intpath, overwrite=True)
    Y_hw_chk = save_minian(
        Y_fm_chk.rename("Y_hw_chk"),
        intpath,
        overwrite=True,
        chunks={"frame": -1, "height": chk["height"], "width": chk["width"]},
    )

    vid_arr = xr.concat([varr_ref, Y_fm_chk], "width").chunk({"width": -1})
    
    #write_video(vid_arr, "minian_mc.mp4", dpath)

    max_proj = save_minian(
        Y_fm_chk.max("frame").rename("max_proj"), **param_save_minian
    ).compute()

    ##################################
            # INTIALISATION #
    ##################################

    # GENERATE SEEDS 
    seeds = seeds_init(Y_fm_chk, **param_seeds_init)

    seeds.head()

    seeds, pnr, gmm = pnr_refine(Y_hw_chk, seeds, **param_pnr_refine)

    seeds = ks_refine(Y_hw_chk, seeds, **param_ks_refine)

    # MERGE SEEDS
    seeds_final = seeds[seeds["mask_ks"] & seeds["mask_pnr"]].reset_index(drop=True)
    seeds_final = seeds_merge(Y_hw_chk, max_proj, seeds_final, **param_seeds_merge)

    # INITIALISE SPATIAL MATRIX
    A_init = initA(Y_hw_chk, seeds_final[seeds_final["mask_mrg"]], **param_initialize)
    A_init = save_minian(A_init.rename("A_init"), intpath, overwrite=True)

    # INITIALISE TEMPORAL MATRIX
    C_init = initC(Y_fm_chk, A_init)
    C_init = save_minian(
        C_init.rename("C_init"), intpath, overwrite=True, chunks={"unit_id": 1, "frame": -1}
    )

    # MERGE UNITS
    A, C = unit_merge(A_init, C_init, **param_init_merge)
    A = save_minian(A.rename("A"), intpath, overwrite=True)
    C = save_minian(C.rename("C"), intpath, overwrite=True)
    C_chk = save_minian(
        C.rename("C_chk"),
        intpath,
        overwrite=True,
        chunks={"unit_id": -1, "frame": chk["frame"]},
    )

    # INITIALIZE BACKGROUND TERMS
    b, f = update_background(Y_fm_chk, A, C_chk)
    f = save_minian(f.rename("f"), intpath, overwrite=True)
    b = save_minian(b.rename("b"), intpath, overwrite=True)

    ##################################
                # CNMF #
    ##################################

    # ESTIMATE SPATIAL NOISE
    sn_spatial = get_noise_fft(Y_hw_chk, **param_get_noise)
    sn_spatial = save_minian(sn_spatial.rename("sn_spatial"), intpath, overwrite=True)
    
    # FIRST SPATIAL UPDATE
    A_new, mask, norm_fac = update_spatial(
        Y_hw_chk, A, C, sn_spatial, **param_first_spatial
    )
    C_new = save_minian(
        (C.sel(unit_id=mask) * norm_fac).rename("C_new"), intpath, overwrite=True
    )
    C_chk_new = save_minian(
        (C_chk.sel(unit_id=mask) * norm_fac).rename("C_chk_new"), intpath, overwrite=True
    )

    b_new, f_new = update_background(Y_fm_chk, A_new, C_chk_new)

    A = save_minian(
        A_new.rename("A"),
        intpath,
        overwrite=True,
        chunks={"unit_id": 1, "height": -1, "width": -1},
    )
    b = save_minian(b_new.rename("b"), intpath, overwrite=True)
    f = save_minian(
        f_new.chunk({"frame": chk["frame"]}).rename("f"), intpath, overwrite=True
    )
    C = save_minian(C_new.rename("C"), intpath, overwrite=True)
    C_chk = save_minian(C_chk_new.rename("C_chk"), intpath, overwrite=True)
    
    # FIRST TEMPORAL UPDATE
    YrA = save_minian(
        compute_trace(Y_fm_chk, A, b, C_chk, f).rename("YrA"),
        intpath,
        overwrite=True,
        chunks={"unit_id": 1, "frame": -1},
    )

    C_new, S_new, b0_new, c0_new, g, mask = update_temporal(
        A, C, YrA=YrA, **param_first_temporal
    )

    C = save_minian(
        C_new.rename("C").chunk({"unit_id": 1, "frame": -1}), intpath, overwrite=True
    )
    C_chk = save_minian(
        C.rename("C_chk"),
        intpath,
        overwrite=True,
        chunks={"unit_id": -1, "frame": chk["frame"]},
    )
    S = save_minian(
        S_new.rename("S").chunk({"unit_id": 1, "frame": -1}), intpath, overwrite=True
    )
    b0 = save_minian(
        b0_new.rename("b0").chunk({"unit_id": 1, "frame": -1}), intpath, overwrite=True
    )
    c0 = save_minian(
        c0_new.rename("c0").chunk({"unit_id": 1, "frame": -1}), intpath, overwrite=True
    )
    A = A.sel(unit_id=C.coords["unit_id"].values)

    # MERGE UNITS
    A_mrg, C_mrg, [sig_mrg] = unit_merge(A, C, [C + b0 + c0], **param_first_merge)


    A = save_minian(A_mrg.rename("A_mrg"), intpath, overwrite=True)
    C = save_minian(C_mrg.rename("C_mrg"), intpath, overwrite=True)
    C_chk = save_minian(
        C.rename("C_mrg_chk"),
        intpath,
        overwrite=True,
        chunks={"unit_id": -1, "frame": chk["frame"]},
    )
    sig = save_minian(sig_mrg.rename("sig_mrg"), intpath, overwrite=True)
    
    # SECOND SPATIAL UPDATE
    A_new, mask, norm_fac = update_spatial(
        Y_hw_chk, A, C, sn_spatial, **param_second_spatial
    )
    C_new = save_minian(
        (C.sel(unit_id=mask) * norm_fac).rename("C_new"), intpath, overwrite=True
    )
    C_chk_new = save_minian(
        (C_chk.sel(unit_id=mask) * norm_fac).rename("C_chk_new"), intpath, overwrite=True
    )

    b_new, f_new = update_background(Y_fm_chk, A_new, C_chk_new)

    A = save_minian(
        A_new.rename("A"),
        intpath,
        overwrite=True,
        chunks={"unit_id": 1, "height": -1, "width": -1},
    )
    b = save_minian(b_new.rename("b"), intpath, overwrite=True)
    f = save_minian(
        f_new.chunk({"frame": chk["frame"]}).rename("f"), intpath, overwrite=True
    )
    C = save_minian(C_new.rename("C"), intpath, overwrite=True)
    C_chk = save_minian(C_chk_new.rename("C_chk"), intpath, overwrite=True)
    
    # SECOND TEMPORAL UPDATE

    YrA = save_minian(
        compute_trace(Y_fm_chk, A, b, C_chk, f).rename("YrA"),
        intpath,
        overwrite=True,
        chunks={"unit_id": 1, "frame": -1},
    )

    C_new, S_new, b0_new, c0_new, g, mask = update_temporal(
        A, C, YrA=YrA, **param_second_temporal
    )
    C = save_minian(
        C_new.rename("C").chunk({"unit_id": 1, "frame": -1}), intpath, overwrite=True
    )
    C_chk = save_minian(
        C.rename("C_chk"),
        intpath,
        overwrite=True,
        chunks={"unit_id": -1, "frame": chk["frame"]},
    )
    S = save_minian(
        S_new.rename("S").chunk({"unit_id": 1, "frame": -1}), intpath, overwrite=True
    )
    b0 = save_minian(
        b0_new.rename("b0").chunk({"unit_id": 1, "frame": -1}), intpath, overwrite=True
    )
    c0 = save_minian(
        c0_new.rename("c0").chunk({"unit_id": 1, "frame": -1}), intpath, overwrite=True
    )
    A = A.sel(unit_id=C.coords["unit_id"].values)
   
    # SAVE FINAL RESULTS

    A = save_minian(A.rename("A"), **param_save_minian)
    C = save_minian(C.rename("C"), **param_save_minian)
    S = save_minian(S.rename("S"), **param_save_minian)
    c0 = save_minian(c0.rename("c0"), **param_save_minian)
    b0 = save_minian(b0.rename("b0"), **param_save_minian)
    b = save_minian(b.rename("b"), **param_save_minian)
    f = save_minian(f.rename("f"), **param_save_minian)
    
    #generate_videos(varr.sel(subset), Y_fm_chk, A=A, C=C_chk, vpath=dpath)

    client.close()
    cluster.close()
    
    # Copy the script file to the destination folder
    source_script = "/mnt/data/home/aurelie.brecier/HayLabAnalysis/python/MINI_1_run_minian_pipeline.py"
    destination_file_path = f"{minian_ds_path}/MINI_1_run_minian_pipeline.txt"
    shutil.copy(source_script, destination_file_path)

    elapsed_time = time.time() - st
    print('Execution time:', time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))