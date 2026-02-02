#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 29 12:14:17 2025

---------------------------------------------------
Written by Marco Bedini (marco.bedini@univ-amu.fr)
---------------------------------------------------

Compute full and partial correlations between clusters (seeds) and parcels (targets),
using Nilearn for full correlations and computing partials per-seed while excluding
seed-local parcels from the conditioning set (those entries set to NaN).
"""

import os
import sys
import json
import numpy as np
from nilearn.connectome import ConnectivityMeasure

USER = os.environ["USER"]

# Main folders
main_data = "/scratch/mszinte/data/RetinoMaps/derivatives/pp_data"
seed_folder = main_data
output_folder = "/scratch/mszinte/data/RetinoMaps/derivatives/pp_data/group/91k/rest/nilearn_partial_corr"

# Custom utils
main_codes = f"/home/{USER}/projects"
utils_path = os.path.join(main_codes, "pRF_analysis/analysis_code/utils")
sys.path.append(utils_path)
from surface_utils import load_surface
from cifti_utils import from_91k_to_32k

# Settings
settings_filepath = f"/home/{USER}/projects/pRF_analysis/RetinoMaps/settings.json"
with open(settings_filepath, "r") as f:
    settings = json.load(f)

subjects = settings["subjects"]
task_name = settings["task_names"][0]

# Seeds (clusters) and flattened parcels
clusters = settings['rois']
parcels = [p for group in settings['rois_groups'] for p in group]

# Map: cluster -> list of parcel names to exclude from partial-conditioning
exclude_parcels_per_cluster = {
    'mPCS': ['SCEF', 'p32pr', '24dv'],
    'sPCS': ['FEF', 'i6-8', '6a', '6d', '6mp', '6ma'],
    'iPCS': ['PEF', 'IFJp', '6v', '6r', 'IFJa', '55b'],
    'sIPS': ['VIP', 'LIPv', 'LIPd', 'IP2', '7PC', 'AIP', '7AL', '7Am', '7Pm'],
    'iIPS': ['IP0', 'IPS1', 'V7', 'MIP', 'IP1', 'V6A', '7PL'],
    'hMT+': ['V8', 'PIT', 'PH', 'FFC', 'VMV1', 'VMV2', 'VMV3', 'VVC'],
    'VO': ['V4t', 'MST', 'MT', 'FST'],
    'LO': ['LO1', 'LO2', 'LO3'],
    'V3AB': ['V3CD', 'V3A', 'V3B'],
    'V3': ['V3', 'V4'],
    'V2': ['V2'],
    'V1': ['V1']
}

#%% Initialize storage
all_subject_full_matrices = []
all_subject_full_matrices_fisherz = []
all_subject_parcel_full = []
all_subject_parcel_full_fisherz = []
all_subject_partial_matrices = []
all_subject_partial_matrices_fisherz = []
all_subject_parcel_names = []

# Loop over subjects
for subject in subjects:
    print(f"\nProcessing {subject}")
    
    # Load timeseries
    timeseries_fn = f'{main_data}/{subject}/91k/rest/timeseries/{subject}_ses-01_task-{task_name}_space-fsLR_den-91k_desc-denoised_bold.dtseries.nii'
    ts_img, ts_data_raw = load_surface(timeseries_fn)
    res = from_91k_to_32k(ts_img, ts_data_raw, return_concat_hemis=True, return_32k_mask=True)
    ts_data = res['data_concat']  # (timepoints x vertices)
    n_time = ts_data.shape[0]

    # Load cluster timeseries (seeds)
    cluster_ts_list = []
    cluster_names_used = []
    for roi in clusters:
        lh, rh = [load_surface(f'{seed_folder}/{subject}/91k/rest/seed/{subject}_91k_intertask_Sac_Pur_vision-pursuit-saccade_{hemi}_{roi}.shape.gii')[1] for hemi in ('lh','rh')]
        mask = np.hstack((lh, rh)).ravel()
        if np.any(mask):
            ts = ts_data[:, mask > 0]
            # guard for single-vertex selection => ensure 2D
            ts_mean = ts.mean(axis=1) if ts.ndim > 1 else ts
            cluster_ts_list.append(ts_mean)
            cluster_names_used.append(roi)
    # column-stack or empty array
    cluster_ts = np.column_stack(cluster_ts_list) if cluster_ts_list else np.empty((n_time, 0))

    # Load parcel timeseries (targets) - we grab this from the MMP ATLAS
    parcel_ts_list = []
    parcel_names_used = []
    for parcel in parcels:
        lh, rh = [load_surface(f'{main_codes}/pRF_analysis/RetinoMaps/rest/mmp1_clusters/parcels/{hemi}_{parcel}_ROI.shape.gii')[1] for hemi in ('L','R')]
        mask = np.hstack((lh, rh)).ravel()
        if np.any(mask):
            ts = ts_data[:, mask > 0]
            ts_mean = ts.mean(axis=1) if ts.ndim > 1 else ts
            parcel_ts_list.append(ts_mean)
            parcel_names_used.append(parcel)
    parcel_ts = np.column_stack(parcel_ts_list) if parcel_ts_list else np.empty((n_time, 0))

    print(f"{subject}: using {len(parcel_names_used)}/{len(parcels)} parcels")
    
    # Parcel × parcel full correlation (symmetric, no exclusions; useful as a sanity check wrt Workbench results)
    if parcel_ts.shape[1] > 1:
        parcel_conn = ConnectivityMeasure(kind='correlation')
        parcel_corr = parcel_conn.fit_transform([parcel_ts])[0]
        parcel_corr_fisherz = np.arctanh(parcel_corr)
    else:
        parcel_corr = np.full((parcel_ts.shape[1], parcel_ts.shape[1]), np.nan)
        parcel_corr_fisherz = np.full((parcel_ts.shape[1], parcel_ts.shape[1]), np.nan)

    # If either side empty, create empty matrices and continue
    n_clusters_present = cluster_ts.shape[1]
    n_parcels_present = parcel_ts.shape[1]
    if n_clusters_present == 0 or n_parcels_present == 0:
        full_matrix = np.empty((n_clusters_present, n_parcels_present))
        full_matrix_fisherz = np.empty((n_clusters_present, n_parcels_present))
        partial_matrix = np.empty((n_clusters_present, n_parcels_present))
        partial_matrix_fisherz = np.empty((n_clusters_present, n_parcels_present))
    else:
        # Full correlation via Nilearn
        combined_for_full = np.hstack([cluster_ts, parcel_ts])  # time x (n_clusters + n_parcels)
        full_conn = ConnectivityMeasure(kind='correlation')
        corr_all = full_conn.fit_transform([combined_for_full])[0]  # (nvars, nvars)
        full_matrix = corr_all[:n_clusters_present, n_clusters_present:(n_clusters_present + n_parcels_present)]
        full_matrix_fisherz = np.arctanh(full_matrix)

        # Partial correlations computed per-seed while excluding within seed parcels
        partial_matrix = np.full((n_clusters_present, n_parcels_present), np.nan)
        partial_matrix_fisherz = np.full((n_clusters_present, n_parcels_present), np.nan)

        # Precompute index mapping: parcel name -> column index in parcel_ts / parcel_names_used
        parcel_name_to_idx = {name: idx for idx, name in enumerate(parcel_names_used)}

        # For each seed (row) compute partials using combined_ts that excludes within seed parcels
        for i_cl, cl_name in enumerate(cluster_names_used):
            # Determine which parcel names to exclude for this cluster (intersection with present parcels)
            exclude_list = exclude_parcels_per_cluster.get(cl_name, [])
            exclude_set = set(exclude_list)
            # Indices of parcels that are allowed (not excluded)
            included_parcel_indices = [j for j, pname in enumerate(parcel_names_used) if pname not in exclude_set]

            if len(included_parcel_indices) == 0:
                # nothing to compute for this seed (all parcels excluded)
                continue

            # Build combined timeseries containing: all clusters (same order) + included parcels only
            combined_seed = np.hstack([cluster_ts, parcel_ts[:, included_parcel_indices]])  # time x (n_clusters + n_included)
            
            # Compute partial correlation on this combined set
            partial_conn = ConnectivityMeasure(kind='partial correlation')
            tmp = partial_conn.fit_transform([combined_seed])[0]

            # Extract row corresponding to the seed (seed is among the first n_clusters_present columns)
            seed_col = i_cl
            n_included = len(included_parcel_indices)
            # The parcel block starts at column index n_clusters_present and spans n_included columns
            partial_values_included = tmp[seed_col, n_clusters_present:(n_clusters_present + n_included)]
            partial_values_included_fisherz = np.arctanh(partial_values_included)

            # Place values back into the full-sized partial_matrix at correct parcel positions
            for k, pj in enumerate(included_parcel_indices):
                partial_matrix[i_cl, pj] = partial_values_included[k]
                partial_matrix_fisherz[i_cl, pj] = partial_values_included_fisherz[k]
            # excluded parcel columns remain NaN for this seed (as desired)

    # Fill global matrix with NaNs for missing parcels and clusters
    # Note: NaNs here reflect parcels that are systematically excluded
    # from partial correlation conditioning for a given seed.
    # Group averages are therefore computed over the valid subset only.

    full_filled = np.full((len(clusters), len(parcels)), np.nan)
    full_filled_fisherz = np.full((len(clusters), len(parcels)), np.nan)
    partial_filled = np.full((len(clusters), len(parcels)), np.nan)
    partial_filled_fisherz = np.full((len(clusters), len(parcels)), np.nan)
    
    parcel_full_filled = np.full((len(parcels), len(parcels)), np.nan)
    parcel_full_filled_fisherz = np.full((len(parcels), len(parcels)), np.nan)

    for i, p_i in enumerate(parcel_names_used):
        gi = parcels.index(p_i)
        for j, p_j in enumerate(parcel_names_used):
            gj = parcels.index(p_j)
            parcel_full_filled[gi, gj] = parcel_corr[i, j]
            parcel_full_filled_fisherz[gi, gj] = parcel_corr_fisherz[i, j]
            
    all_subject_parcel_full.append(parcel_full_filled)
    all_subject_parcel_full_fisherz.append(parcel_full_filled_fisherz)
 
    # Per-subject output folder
    sub_out = f'{main_data}/{subject}/91k/rest/corr/partial_corr'
    os.makedirs(sub_out, exist_ok=True)

    # Map subject-local cluster names to global cluster indices
    for i_cl, cl in enumerate(cluster_names_used):
        global_r = clusters.index(cl)
        for j_pa, pa in enumerate(parcel_names_used):
            global_c = parcels.index(pa)
            # Guard indexing (full_matrix and partial_matrix are subject-local shapes)
            if full_matrix.size:
                full_filled[global_r, global_c] = full_matrix[i_cl, j_pa]
                full_filled_fisherz[global_r, global_c] = full_matrix_fisherz[i_cl, j_pa]
            if partial_matrix.size:
                # partial_matrix has NaN for excluded parcels by construction
                partial_filled[global_r, global_c] = partial_matrix[i_cl, j_pa]
                partial_filled_fisherz[global_r, global_c] = partial_matrix_fisherz[i_cl, j_pa]

    # Save per subject outputs
    np.save(os.path.join(sub_out, "cluster_by_mmp-parcel_full.npy"), full_filled)
    np.save(os.path.join(sub_out, "cluster_by_mmp-parcel_full_fisherz.npy"), full_filled_fisherz)
    np.save(os.path.join(sub_out, "cluster_by_mmp-parcel_partial.npy"), partial_filled)
    np.save(os.path.join(sub_out, "cluster_by_mmp-parcel_partial_fisherz.npy"), partial_filled_fisherz)
    np.save(os.path.join(sub_out, "parcel_by_mmp-parcel_full.npy"), parcel_full_filled)
    np.save(os.path.join(sub_out, "parcel_by_mmp-parcel_full_fisherz.npy"), parcel_full_filled_fisherz)

    # Append for group stats
    all_subject_full_matrices.append(full_filled)
    all_subject_full_matrices_fisherz.append(full_filled_fisherz)
    all_subject_partial_matrices.append(partial_filled)
    all_subject_partial_matrices_fisherz.append(partial_filled_fisherz)

#%% Group stats (n_subjects, n_clusters, n_parcels)
all_subject_full_matrices = np.stack(all_subject_full_matrices, axis=0)
all_subject_full_matrices_fisherz = np.stack(all_subject_full_matrices_fisherz, axis=0)
all_subject_partial_matrices = np.stack(all_subject_partial_matrices, axis=0)
all_subject_partial_matrices_fisherz = np.stack(all_subject_partial_matrices_fisherz, axis=0)

mean_full_matrix = np.nanmean(all_subject_full_matrices, axis=0)
mean_full_matrix_fisherz = np.nanmean(all_subject_full_matrices_fisherz, axis=0)
median_full_matrix = np.nanmedian(all_subject_full_matrices, axis=0)
median_full_matrix_fisherz = np.nanmedian(all_subject_full_matrices_fisherz, axis=0)
mean_partial_matrix = np.nanmean(all_subject_partial_matrices, axis=0)
mean_partial_matrix_fisherz = np.nanmean(all_subject_partial_matrices_fisherz, axis=0)
median_partial_matrix = np.nanmedian(all_subject_partial_matrices, axis=0)
median_partial_matrix_fisherz = np.nanmedian(all_subject_partial_matrices_fisherz, axis=0)

all_subject_parcel_full = np.stack(all_subject_parcel_full, axis=0)
all_subject_parcel_full_fisherz = np.stack(all_subject_parcel_full_fisherz, axis=0)
mean_parcel_full = np.nanmean(all_subject_parcel_full, axis=0)
mean_parcel_full_fisherz = np.nanmean(all_subject_parcel_full_fisherz, axis=0)
median_parcel_full = np.nanmedian(all_subject_parcel_full, axis=0)
median_parcel_full_fisherz = np.nanmedian(all_subject_parcel_full_fisherz, axis=0)

# Save results
np.savez_compressed(
    os.path.join(output_folder, "nilearn_partial_corr_task-mmp.npz"),
    mean_cluster_parcel_full=mean_full_matrix,
    mean_cluster_parcel_full_fisherz=mean_full_matrix_fisherz,
    median_cluster_parcel_full=median_full_matrix,
    median_cluster_parcel_full_fisherz=median_full_matrix_fisherz,
    mean_cluster_parcel_partial=mean_partial_matrix,
    mean_cluster_parcel_partial_fisherz=mean_partial_matrix_fisherz,
    median_cluster_parcel_partial=median_partial_matrix,
    median_cluster_parcel_partial_fisherz=median_partial_matrix_fisherz,
    mean_parcel_full=mean_parcel_full,
    mean_parcel_full_fisherz=mean_parcel_full_fisherz,
    median_parcel_full=median_parcel_full,
    median_parcel_full_fisherz=median_parcel_full_fisherz,
    subjects=np.array(subjects),
    clusters=np.array(clusters),
    parcels=np.array(parcels),
)

print(f"Cluster × parcel correlation matrices saved to {output_folder}/nilearn_partial_corr_task-mmp.npz")
