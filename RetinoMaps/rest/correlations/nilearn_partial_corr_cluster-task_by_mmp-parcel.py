#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  6 14:57:04 2026

Compute PARTIAL correlations between clusters (seeds) and parcels (targets)
using Nilearn partial correlation.

For each cluster:
  - parcels belonging to that cluster are EXCLUDED from conditioning set
  - output shape = (n_clusters, n_parcels)

Written by Marco Bedini: marco.bedini@univ-amu.fr
Refactored for clarity and debugging
"""

import os
import sys
import numpy as np
import pandas as pd
from nilearn.connectome import ConnectivityMeasure

# =========================
# PATHS
# =========================

main_data = "/scratch/mszinte/data/RetinoMaps/derivatives/pp_data"
partial_output_folder = "/scratch/mszinte/data/RetinoMaps/derivatives/pp_data/group/91k/rest/partial_corr"
os.makedirs(partial_output_folder, exist_ok=True)

USER = os.environ["USER"]
main_codes = f"/home/{USER}/projects"
utils_path = os.path.join(main_codes, "pRF_analysis/analysis_code/utils")
sys.path.append(utils_path)

from surface_utils import load_surface
from cifti_utils import from_91k_to_32k

# =========================
# SUBJECTS / ROIS
# =========================

subjects = [
    'sub-01','sub-02','sub-03','sub-04','sub-05','sub-06',
    'sub-07','sub-08','sub-09','sub-11','sub-12','sub-13',
    'sub-14','sub-17','sub-20','sub-21','sub-22','sub-23',
    'sub-24','sub-25'
]

clusters = ['mPCS','sPCS','iPCS','sIPS','iIPS','hMT+','VO','LO','V3AB','V3','V2','V1']

seed_to_parcels = {
    'mPCS': ['SCEF','p32pr','24dv'],
    'sPCS': ['FEF','i6-8','6a','6d','6mp','6ma'],
    'iPCS': ['PEF','IFJp','6v','6r','IFJa','55b'],
    'sIPS': ['VIP','LIPv','LIPd','IP2','7PC','AIP','7AL','7Am','7Pm'],
    'iIPS': ['IP0','IPS1','V7','MIP','IP1','V6A','7PL'],
    'hMT+': ['V4t','MST','MT','FST'],
    'VO': ['V8','PIT','PH','FFC','VMV1','VMV2','VMV3','VVC'],
    'LO': ['LO1','LO2','LO3'],
    'V3AB': ['V3CD','V3A','V3B'],
    'V3': ['V3','V4'],
    'V2': ['V2'],
    'V1': ['V1']
}

parcels = []
for cl in clusters:
    parcels.extend(seed_to_parcels[cl])

exclude_parcels_per_cluster = seed_to_parcels

# =========================
# STORAGE
# =========================

all_subject_partial = []
all_subject_partial_fz = []

# =========================
# SUBJECT LOOP
# =========================

for subject in subjects:

    print(f"\n=== Processing {subject} ===")

    timeseries_fn = f"{main_data}/{subject}/91k/rest/timeseries/{subject}_ses-01_task-rest_space-fsLR_den-91k_desc-denoised_bold.dtseries.nii"
    if not os.path.exists(timeseries_fn):
        print(f"  ❌ Missing timeseries: {timeseries_fn}")
        continue

    ts_img, ts_data_raw = load_surface(timeseries_fn)
    res = from_91k_to_32k(ts_img, ts_data_raw, return_concat_hemis=True)
    ts_data = res["data_concat"]
    n_time = ts_data.shape[0]

    print(f"  Loaded timeseries: {ts_data.shape}")

    # =========================
    # CLUSTER TIMESERIES
    # =========================

    cluster_ts_list = []
    cluster_names_used = []

    for roi in clusters:
        lh, rh = [
            load_surface(f"{main_data}/{subject}/91k/rest/seed/{subject}_91k_intertask_Sac_Pur_vision-pursuit-saccade_{hemi}_{roi}.shape.gii")[1]
            for hemi in ("lh", "rh")
        ]
        mask = np.hstack((lh, rh)).ravel()
        if np.any(mask):
            ts = ts_data[:, mask > 0]
            cluster_ts_list.append(ts.mean(axis=1))
            cluster_names_used.append(roi)

    cluster_ts = np.column_stack(cluster_ts_list)
    print(f"  Clusters used: {cluster_names_used}")

    # =========================
    # PARCEL TIMESERIES
    # =========================

    parcel_ts_list = []
    parcel_names_used = []

    for parcel in parcels:
        lh, rh = [
            load_surface(f"{main_codes}/pRF_analysis/RetinoMaps/rest/mmp1_clusters/parcels/{hemi}_{parcel}_ROI.shape.gii")[1]
            for hemi in ("L", "R")
        ]
        mask = np.hstack((lh, rh)).ravel()
        if np.any(mask):
            ts = ts_data[:, mask > 0]
            parcel_ts_list.append(ts.mean(axis=1))
            parcel_names_used.append(parcel)

    parcel_ts = np.column_stack(parcel_ts_list)
    print(f"  Parcels used: {len(parcel_names_used)}/{len(parcels)}")

    if cluster_ts.size == 0 or parcel_ts.size == 0:
        print("  ⚠️ Empty cluster or parcel matrix — skipping subject")
        continue

    n_clusters = cluster_ts.shape[1]

    partial_matrix = np.full((n_clusters, parcel_ts.shape[1]), np.nan)
    partial_matrix_fz = np.full_like(partial_matrix, np.nan)

    # =========================
    # PARTIAL CORRELATION
    # =========================

    for i_cl, cl_name in enumerate(cluster_names_used):

        exclude = set(exclude_parcels_per_cluster.get(cl_name, []))
        included_idx = [j for j,p in enumerate(parcel_names_used) if p not in exclude]

        if not included_idx:
            print(f"  ⚠️ No parcels left after exclusion for {cl_name}")
            continue

        combined = np.hstack([cluster_ts, parcel_ts[:, included_idx]])

        try:
            tmp = ConnectivityMeasure(kind="partial correlation").fit_transform([combined])[0]
            vals = tmp[i_cl, n_clusters:]
            partial_matrix[i_cl, included_idx] = vals
            partial_matrix_fz[i_cl, included_idx] = np.arctanh(vals)

        except Exception as e:
            print(f"  ❌ Partial corr failed for {cl_name}: {e}")

    # =========================
    # GLOBAL MATRIX FILL
    # =========================

    filled = np.full((len(clusters), len(parcels)), np.nan)
    filled_fz = filled.copy()

    for i_cl, cl in enumerate(cluster_names_used):
        gr = clusters.index(cl)
        for j_pa, pa in enumerate(parcel_names_used):
            gc = parcels.index(pa)
            filled[gr, gc] = partial_matrix[i_cl, j_pa]
            filled_fz[gr, gc] = partial_matrix_fz[i_cl, j_pa]

    sub_out = f"{main_data}/{subject}/91k/rest/corr/partial_corr"
    os.makedirs(sub_out, exist_ok=True)

    np.save(os.path.join(sub_out, "cluster_by_mmp-parcel_partial.npy"), filled)
    np.save(os.path.join(sub_out, "cluster_by_mmp-parcel_partial_fisherz.npy"), filled_fz)

    df = pd.DataFrame(filled, index=clusters, columns=parcels)
    df_fz = pd.DataFrame(filled_fz, index=clusters, columns=parcels)

    df.to_csv(os.path.join(sub_out, "cluster_by_mmp-parcel_partial.csv"))
    df_fz.to_csv(os.path.join(sub_out, "cluster_by_mmp-parcel_partial_fisherz.csv"))

    all_subject_partial.append(filled)
    all_subject_partial_fz.append(filled_fz)

# =========================
# GROUP STATS
# =========================

print("\n=== GROUP STATS ===")

all_subject_partial = np.stack(all_subject_partial, axis=0)
all_subject_partial_fz = np.stack(all_subject_partial_fz, axis=0)

mean_partial = np.nanmean(all_subject_partial, axis=0)
median_partial = np.nanmedian(all_subject_partial, axis=0)

mean_partial_fz = np.nanmean(all_subject_partial_fz, axis=0)
median_partial_fz = np.nanmedian(all_subject_partial_fz, axis=0)

df_mean = pd.DataFrame(mean_partial, index=clusters, columns=parcels)
df_median = pd.DataFrame(median_partial, index=clusters, columns=parcels)

df_mean.to_csv(os.path.join(partial_output_folder, "group_mean_cluster_by_mmp-parcel_partial.csv"))
df_median.to_csv(os.path.join(partial_output_folder, "group_median_cluster_by_mmp-parcel_partial.csv"))

np.savez_compressed(
    os.path.join(partial_output_folder, "group_partial_corr.npz"),
    mean_cluster_parcel_partial=mean_partial,
    mean_cluster_parcel_partial_fisherz=mean_partial_fz,
    median_cluster_parcel_partial=median_partial,
    median_cluster_parcel_partial_fisherz=median_partial_fz,
    subjects=np.array(subjects),
    clusters=np.array(clusters),
    parcels=np.array(parcels),
)

print("\n✅ Done.")
