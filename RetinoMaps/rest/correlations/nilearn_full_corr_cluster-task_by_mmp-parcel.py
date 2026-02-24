#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Feb 6, 2025

---------------------------------------------------
Written by Marco Bedini (marco.bedini@univ-amu.fr)
---------------------------------------------------

Compute full correlations between clusters (seeds) and parcels (targets)
These full correlations are computed with Nilearn as a sanity check wrt to
Results obtained using Workbench
"""

import os
import sys
import numpy as np
import pandas as pd
from pathlib import Path
from nilearn.connectome import ConnectivityMeasure

USER = os.environ["USER"]

# ============================================================
# Paths
# ============================================================
USER = os.environ["USER"]

# Main folders
main_data = "/scratch/mszinte/data/RetinoMaps/derivatives/pp_data"
seed_folder = main_data

# Output folders
full_output_folder = "/scratch/mszinte/data/RetinoMaps/derivatives/pp_data/group/91k/rest/partial_corr"
os.makedirs(full_output_folder, exist_ok=True)

# General utils
main_codes = f"/home/{USER}/meso_H/projects"
gutils_path = os.path.join(main_codes, "pRF_analysis/analysis_code/utils")
sys.path.append(gutils_path)
from surface_utils import load_surface
from cifti_utils import from_91k_to_32k

# Custom utils
utils_path = os.path.join(main_codes, "pRF_analysis/RetinoMaps")
sys.path.append(utils_path)

# Personal imports
sys.path.append("{}/../../../../utils".format(os.getcwd()))
from settings_utils import load_settings

# Load settings
main_codes = os.path.abspath(os.path.join(os.getcwd(), "../../../../"))
settings_path = os.path.join(main_codes, utils_path, "settings.yml")
prf_settings_path = os.path.join(main_codes, utils_path, "prf-analysis.yml")
settings = load_settings([settings_path, prf_settings_path])
analysis_info = settings[0]
subjects = analysis_info['subjects']

# ============================================================
# ROIs
# ============================================================

# Load seed clusters and parcel to cluster assignments
clusters = analysis_info['rois-drawn']
seed_to_parcels = analysis_info['rois-group-mmp']

# Have mPCS as the first cluster instead of V1
clusters.reverse()

parcels = []
for cl in clusters:
    parcels.extend(seed_to_parcels[cl])

seed_to_number = {s: i+1 for i,s in enumerate(clusters)}

#%% Initialize storage

all_subject_full_matrices = []
all_subject_full_matrices_fisherz = []
all_subject_parcel_full = []
all_subject_parcel_full_fisherz = []

#%% Subject loop

for subject in subjects:
    print(f"\nProcessing {subject}")

    timeseries_fn = f'{main_data}/{subject}/91k/rest/timeseries/{subject}_ses-01_task-rest_space-fsLR_den-91k_desc-denoised_bold.dtseries.nii'
    ts_img, ts_data_raw = load_surface(timeseries_fn)
    res = from_91k_to_32k(ts_img, ts_data_raw, return_concat_hemis=True, return_32k_mask=True)
    ts_data = res['data_concat']
    n_time = ts_data.shape[0]

    # --- clusters ---
    cluster_ts_list = []
    cluster_names_used = []
    for roi in clusters:
        lh, rh = [load_surface(f'{seed_folder}/{subject}/91k/rest/seed/{subject}_91k_intertask_Sac_Pur_vision-pursuit-saccade_{hemi}_{roi}.shape.gii')[1]
                  for hemi in ('lh','rh')]
        mask = np.hstack((lh,rh)).ravel()
        if np.any(mask):
            ts = ts_data[:, mask > 0]
            ts_mean = ts.mean(axis=1) if ts.ndim > 1 else ts
            cluster_ts_list.append(ts_mean)
            cluster_names_used.append(roi)

    cluster_ts = np.column_stack(cluster_ts_list) if cluster_ts_list else np.empty((n_time,0))

    # --- parcels ---
    parcel_ts_list = []
    parcel_names_used = []
    for parcel in parcels:
        lh, rh = [load_surface(f'{main_codes}/pRF_analysis/RetinoMaps/rest/mmp1_clusters/parcels/{hemi}_{parcel}_ROI.shape.gii')[1]
                  for hemi in ('L','R')]
        mask = np.hstack((lh,rh)).ravel()
        if np.any(mask):
            ts = ts_data[:, mask > 0]
            ts_mean = ts.mean(axis=1) if ts.ndim > 1 else ts
            parcel_ts_list.append(ts_mean)
            parcel_names_used.append(parcel)

    parcel_ts = np.column_stack(parcel_ts_list) if parcel_ts_list else np.empty((n_time,0))

    print(f"{subject}: using {len(parcel_names_used)}/{len(parcels)} parcels")

    # --- parcel × parcel ---
    if parcel_ts.shape[1] > 1:
        parcel_corr = ConnectivityMeasure(kind='correlation').fit_transform([parcel_ts])[0]
        np.fill_diagonal(parcel_corr, np.nan)
        parcel_corr_fisherz = np.arctanh(parcel_corr)
    else:
        parcel_corr = np.full((parcel_ts.shape[1],parcel_ts.shape[1]),np.nan)
        parcel_corr_fisherz = parcel_corr.copy()

    n_clusters_present = cluster_ts.shape[1]
    n_parcels_present = parcel_ts.shape[1]

    if n_clusters_present == 0 or n_parcels_present == 0:
        full_matrix = np.empty((n_clusters_present,n_parcels_present))
        full_matrix_fisherz = full_matrix.copy()
    else:
        combined = np.hstack([cluster_ts,parcel_ts])
        corr_all = ConnectivityMeasure(kind='correlation').fit_transform([combined])[0]
        full_matrix = corr_all[:n_clusters_present, n_clusters_present:]
        full_matrix_fisherz = np.arctanh(full_matrix)

    # --- fill global matrices ---
    full_filled = np.full((len(clusters),len(parcels)),np.nan)
    full_filled_fz = full_filled.copy()

    parcel_full_filled = np.full((len(parcels),len(parcels)),np.nan)
    parcel_full_filled_fz = parcel_full_filled.copy()

    for i,p_i in enumerate(parcel_names_used):
        gi = parcels.index(p_i)
        for j,p_j in enumerate(parcel_names_used):
            gj = parcels.index(p_j)
            parcel_full_filled[gi,gj] = parcel_corr[i,j]
            parcel_full_filled_fz[gi,gj] = parcel_corr_fisherz[i,j]

    for i_cl,cl in enumerate(cluster_names_used):
        gr = clusters.index(cl)
        for j_pa,pa in enumerate(parcel_names_used):
            gc = parcels.index(pa)
            full_filled[gr,gc] = full_matrix[i_cl,j_pa]
            full_filled_fz[gr,gc] = full_matrix_fisherz[i_cl,j_pa]

    # --- subject output dirs ---
    sub_full = f'{main_data}/{subject}/91k/rest/corr/full_corr/nilearn_full_corr'
    os.makedirs(sub_full,exist_ok=True)

    # =========================
    # Save NUMPY arrays
    # =========================

    np.save(os.path.join(sub_full, "cluster_by_mmp-parcel_full.npy"), full_filled)
    np.save(os.path.join(sub_full, "cluster_by_mmp-parcel_full_fisherz.npy"), full_filled_fz)
    np.save(os.path.join(sub_full, "parcel_by_mmp-parcel_full.npy"), parcel_full_filled)
    np.save(os.path.join(sub_full, "parcel_by_mmp-parcel_full_fisherz.npy"), parcel_full_filled_fz)

    # =========================
    # Save as labeled DataFrames (CSV)
    # =========================

    df_full = pd.DataFrame(full_filled, index=clusters, columns=parcels)
    df_full_fz = pd.DataFrame(full_filled_fz, index=clusters, columns=parcels)

    df_parcel_full = pd.DataFrame(parcel_full_filled, index=parcels, columns=parcels)
    df_parcel_full_fz = pd.DataFrame(parcel_full_filled_fz, index=parcels, columns=parcels)

    df_full.to_csv(os.path.join(sub_full, "cluster_by_mmp-parcel_full.csv"))
    df_full_fz.to_csv(os.path.join(sub_full, "cluster_by_mmp-parcel_full_fisherz.csv"))
    df_parcel_full.to_csv(os.path.join(sub_full, "parcel_by_mmp-parcel_full.csv"))
    df_parcel_full_fz.to_csv(os.path.join(sub_full, "parcel_by_mmp-parcel_full_fisherz.csv"))

    # =========================
    # Append for group stats
    # =========================

    all_subject_full_matrices.append(full_filled)
    all_subject_full_matrices_fisherz.append(full_filled_fz)
    all_subject_parcel_full.append(parcel_full_filled)
    all_subject_parcel_full_fisherz.append(parcel_full_filled_fz)


#%% Group stats

all_subject_full_matrices = np.stack(all_subject_full_matrices, axis=0)
all_subject_full_matrices_fisherz = np.stack(all_subject_full_matrices_fisherz, axis=0)

all_subject_parcel_full = np.stack(all_subject_parcel_full, axis=0)
all_subject_parcel_full_fisherz = np.stack(all_subject_parcel_full_fisherz, axis=0)

# Cluster × Parcel
mean_full = np.nanmean(all_subject_full_matrices, axis=0)
mean_full_fz = np.nanmean(all_subject_full_matrices_fisherz, axis=0)

median_full = np.nanmedian(all_subject_full_matrices, axis=0)
median_full_fz = np.nanmedian(all_subject_full_matrices_fisherz, axis=0)

# Parcel × Parcel
mean_parcel_full = np.nanmean(all_subject_parcel_full, axis=0)
mean_parcel_full_fisherz = np.nanmean(all_subject_parcel_full_fisherz, axis=0)

median_parcel_full = np.nanmedian(all_subject_parcel_full, axis=0)
median_parcel_full_fisherz = np.nanmedian(all_subject_parcel_full_fisherz, axis=0)

# DataFrames
df_mean_full = pd.DataFrame(mean_full, index=clusters, columns=parcels)
df_median_full = pd.DataFrame(median_full, index=clusters, columns=parcels)

df_mean_parcel_full = pd.DataFrame(mean_parcel_full, index=parcels, columns=parcels)
df_median_parcel_full = pd.DataFrame(median_parcel_full, index=parcels, columns=parcels)

# Save CSVs
df_mean_full.to_csv(os.path.join(full_output_folder, "group_mean_cluster_by_mmp-parcel_full.csv"))
df_median_full.to_csv(os.path.join(full_output_folder, "group_median_cluster_by_mmp-parcel_full.csv"))

df_mean_parcel_full.to_csv(os.path.join(full_output_folder, "group_mean_parcel_by_mmp-parcel_full.csv"))
df_median_parcel_full.to_csv(os.path.join(full_output_folder, "group_median_parcel_by_mmp-parcel_full.csv"))

# Save NPZ
np.savez_compressed(
    os.path.join(full_output_folder, "group_full_corr.npz"),
    mean_cluster_parcel_full=mean_full,
    mean_cluster_parcel_full_fisherz=mean_full_fz,
    median_cluster_parcel_full=median_full,
    median_cluster_parcel_full_fisherz=median_full_fz,
    mean_parcel_full=mean_parcel_full,
    mean_parcel_full_fisherz=mean_parcel_full_fisherz,
    median_parcel_full=median_parcel_full,
    median_parcel_full_fisherz=median_parcel_full_fisherz,
    subjects=np.array(subjects),
    clusters=np.array(clusters),
    parcels=np.array(parcels),
)

print("Done.")
