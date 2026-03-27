#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mar 27, 2026

Compute INTRA-HEMISPHERIC* PARTIAL correlations between clusters (seeds)
and parcels (targets) using Nilearn partial correlation

* i.e., Left-hemisphere seeds are correlated only with left-hemisphere parcels,
and right-hemisphere seeds with right-hemisphere parcels

For each cluster:
  - parcels belonging to that cluster are EXCLUDED from the conditioning set
  - output shape per hemisphere = (n_clusters, n_parcels)

---------------------------------------------------
Written by Marco Bedini (marco.bedini@univ-amu.fr)
---------------------------------------------------
"""

import os
import sys
import numpy as np
import pandas as pd
from pathlib import Path
from nilearn.connectome import ConnectivityMeasure

# ============================================================
# Paths
# ============================================================
USER = os.environ["USER"]

main_data        = "/scratch/mszinte/data/RetinoMaps/derivatives/pp_data"
seed_folder      = main_data
atlas_folder     = "/scratch/mszinte/data/RetinoMaps/derivatives/pp_data/atlas"
partial_output_folder = (
    "/scratch/mszinte/data/RetinoMaps/derivatives/pp_data"
    "/group/91k/rest/partial_corr_by_hemi"
)
os.makedirs(partial_output_folder, exist_ok=True)

# Personal imports
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../"))
sys.path.append(os.path.abspath(os.path.join(base_dir, "analysis_code/utils")))
from settings_utils import load_settings
from surface_utils import load_surface
from cifti_utils import from_91k_to_32k

# Load settings
project_dir = "RetinoMaps"
settings_path     = os.path.join(base_dir, project_dir, "settings.yml")
prf_settings_path = os.path.join(base_dir, project_dir, "prf-analysis.yml")
settings          = load_settings([settings_path, prf_settings_path])
analysis_info     = settings[0]
subjects          = analysis_info["subjects"]

# ============================================================
# ROIs
# ============================================================

clusters = analysis_info["rois-drawn"]
seed_to_parcels = analysis_info["rois-group-mmp"]

# mPCS first
clusters.reverse()

parcels = []
for cl in clusters:
    parcels.extend(seed_to_parcels[cl])

# ============================================================
# Hemisphere configuration
#
# Each entry drives the full pipeline for one hemisphere:
#   seed_key  : suffix used in seed .shape.gii filenames  (lh / rh)
#   atlas_key : prefix used in atlas parcel filenames     (L  / R )
#   ts_key    : key returned by from_91k_to_32k for the
#               hemisphere-specific timeseries array
# ============================================================

HEMIS = [
    {"label": "LH", "seed_key": "lh", "atlas_key": "L", "ts_key": "data_L"},
    {"label": "RH", "seed_key": "rh", "atlas_key": "R", "ts_key": "data_R"},
]

# =========================
# STORAGE  (keyed by hemi label)
# =========================

all_subject_partial    = {h["label"]: [] for h in HEMIS}
all_subject_partial_fz = {h["label"]: [] for h in HEMIS}

# =========================
# SUBJECT LOOP
# =========================

for subject in subjects:

    print(f"\n=== Processing {subject} ===")

    timeseries_fn = (
        f"{main_data}/{subject}/91k/rest/timeseries/"
        f"{subject}_ses-01_task-rest_space-fsLR_den-91k_desc-denoised_bold.dtseries.nii"
    )

    ts_img, ts_data_raw = load_surface(timeseries_fn)

    # Request separate hemisphere arrays from the utils
    res = from_91k_to_32k(ts_img, ts_data_raw, return_concat_hemis=False, return_32k_mask=True) # sanity check for medial wall vertices

    for h in HEMIS:
        label     = h["label"]          # "LH" or "RH"
        seed_key  = h["seed_key"]       # "lh" or "rh"
        atlas_key = h["atlas_key"]      # "L"  or "R"
        ts_key    = h["ts_key"]         # "data_lh" or "data_rh"

        if ts_key not in res:
            raise KeyError(
                f"from_91k_to_32k result is missing key '{ts_key}'. "
                "Ensure the function returns per-hemisphere arrays."
            )

        ts_data = res[ts_key]           # shape: (n_time, n_vertices_hemi)

        print(f"\n  [{label}] Timeseries shape: {ts_data.shape}")

        # -------------------------
        # Cluster timeseries (single hemisphere)
        # -------------------------

        cluster_ts_list   = []
        cluster_names_used = []

        for roi in clusters:
            mask_img, mask_data = load_surface(
                f"{main_data}/{subject}/91k/rest/seed/"
                f"{subject}_91k_intertask_Sac_Pur_vision-pursuit-saccade"
                f"_{seed_key}_{roi}.shape.gii"
            )
            mask = mask_data.ravel()

            if not np.any(mask):
                print(f"  [{label}] ⚠️  Empty seed mask for {roi} — skipping")
                continue

            cluster_ts_list.append(ts_data[:, mask > 0].mean(axis=1))
            cluster_names_used.append(roi)

        if not cluster_ts_list:
            print(f"  [{label}] ⚠️  No valid cluster timeseries — skipping hemisphere")
            continue

        cluster_ts = np.column_stack(cluster_ts_list)
        print(f"  [{label}] Clusters used: {cluster_names_used}")

        # -------------------------
        # Parcel timeseries (single hemisphere)
        # -------------------------

        parcel_ts_list   = []
        parcel_names_used = []

        for parcel in parcels:
            mask_img, mask_data = load_surface(
                f"{atlas_folder}/mmp1_clusters/parcels/"
                f"{atlas_key}_{parcel}_ROI.shape.gii"
            )
            mask = mask_data.ravel()

            if not np.any(mask):
                continue

            parcel_ts_list.append(ts_data[:, mask > 0].mean(axis=1))
            parcel_names_used.append(parcel)

        if not parcel_ts_list:
            print(f"  [{label}] ⚠️  No valid parcel timeseries — skipping hemisphere")
            continue

        parcel_ts = np.column_stack(parcel_ts_list)
        print(f"  [{label}] Parcels used: {len(parcel_names_used)}/{len(parcels)}")

        # -------------------------
        # Partial correlation
        # -------------------------

        n_clusters_used = cluster_ts.shape[1]
        n_parcels_used  = parcel_ts.shape[1]

        partial_matrix    = np.full((n_clusters_used, n_parcels_used), np.nan)
        partial_matrix_fz = np.full_like(partial_matrix, np.nan)

        for i_cl, cl_name in enumerate(cluster_names_used):

            exclude = set(seed_to_parcels.get(cl_name, []))

            included_idx = [
                j for j, p in enumerate(parcel_names_used)
                if p not in exclude
            ]

            if not included_idx:
                print(f"  [{label}] ⚠️  No parcels left after exclusion for {cl_name}")
                continue

            X = np.column_stack([
                cluster_ts[:, i_cl],
                parcel_ts[:, included_idx],
            ])

            conn = ConnectivityMeasure(kind="partial correlation")
            C    = conn.fit_transform([X])[0]

            vals = C[0, 1:]   # seed row, skip self-correlation

            partial_matrix[i_cl, included_idx]    = vals
            partial_matrix_fz[i_cl, included_idx] = np.arctanh(vals)

        # -------------------------
        # Map back to full (clusters × parcels) grid
        # -------------------------

        filled    = np.full((len(clusters), len(parcels)), np.nan)
        filled_fz = filled.copy()

        for i_cl, cl in enumerate(cluster_names_used):
            gr = clusters.index(cl)
            for j_pa, pa in enumerate(parcel_names_used):
                gc = parcels.index(pa)
                filled[gr, gc]    = partial_matrix[i_cl, j_pa]
                filled_fz[gr, gc] = partial_matrix_fz[i_cl, j_pa]

        # -------------------------
        # Save subject-level results
        # -------------------------

        sub_out = f"{main_data}/{subject}/91k/rest/corr/partial_corr_by_hemi"
        os.makedirs(sub_out, exist_ok=True)

        tag = label.lower()   # "lh" or "rh"

        np.save(
            os.path.join(sub_out, f"cluster_by_mmp-parcel_partial_{tag}.npy"),
            filled,
        )
        np.save(
            os.path.join(sub_out, f"cluster_by_mmp-parcel_partial_fisherz_{tag}.npy"),
            filled_fz,
        )

        pd.DataFrame(filled,    index=clusters, columns=parcels).to_csv(
            os.path.join(sub_out, f"cluster_by_mmp-parcel_partial_{tag}.csv")
        )
        pd.DataFrame(filled_fz, index=clusters, columns=parcels).to_csv(
            os.path.join(sub_out, f"cluster_by_mmp-parcel_partial_fisherz_{tag}.csv")
        )

        all_subject_partial[label].append(filled)
        all_subject_partial_fz[label].append(filled_fz)

# =========================
# GROUP STATS (per hemisphere)
# =========================

print("\n=== GROUP STATS ===")

for h in HEMIS:
    label = h["label"]

    sub_list    = all_subject_partial[label]
    sub_list_fz = all_subject_partial_fz[label]

    if not sub_list:
        print(f"  [{label}] ⚠️  No subject data collected — skipping group stats")
        continue

    stacked    = np.stack(sub_list,    axis=0)   # (n_subjects, n_clusters, n_parcels)
    stacked_fz = np.stack(sub_list_fz, axis=0)

    mean_partial    = np.nanmean(stacked,    axis=0)
    median_partial  = np.nanmedian(stacked,  axis=0)
    mean_partial_fz = np.nanmean(stacked_fz, axis=0)
    median_partial_fz = np.nanmedian(stacked_fz, axis=0)

    tag = label.lower()

    pd.DataFrame(mean_partial,   index=clusters, columns=parcels).to_csv(
        os.path.join(partial_output_folder,
                     f"group_mean_cluster_by_mmp-parcel_partial_{tag}.csv")
    )
    pd.DataFrame(median_partial, index=clusters, columns=parcels).to_csv(
        os.path.join(partial_output_folder,
                     f"group_median_cluster_by_mmp-parcel_partial_{tag}.csv")
    )

    np.savez_compressed(
        os.path.join(partial_output_folder, f"group_partial_corr_{tag}.npz"),
        mean_cluster_parcel_partial         = mean_partial,
        mean_cluster_parcel_partial_fisherz = mean_partial_fz,
        median_cluster_parcel_partial       = median_partial,
        median_cluster_parcel_partial_fisherz = median_partial_fz,
        subjects  = np.array(subjects),
        clusters  = np.array(clusters),
        parcels   = np.array(parcels),
        hemi      = np.array(label),
    )

    print(f"  [{label}] Group stats saved.")

print("\nDone.")
# ============================================================