#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compute INTRA-HEMISPHERIC PARTIAL correlations between clusters (seeds)
and parcels (targets) using Nilearn partial correlation.

Left-hemisphere seeds are correlated only with left-hemisphere parcels,
and right-hemisphere seeds with right-hemisphere parcels.

Important note: the conditioning set now includes all other parcels
Whether intra or inter-hemispheric 
Matching Dawson et al. (2016) and Genç et al. (2016) 

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
    "/group/91k/rest/partial_corr/by_hemi"
)
os.makedirs(partial_output_folder, exist_ok=True)

# Personal imports
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../"))
sys.path.append(os.path.abspath(os.path.join(base_dir, "analysis_code/utils")))
from settings_utils import load_settings
from surface_utils import load_surface
from cifti_utils import from_91k_to_32k

# ============================================================
# Helpers
# ============================================================

def impute_nan_columns(X: np.ndarray, label: str = "") -> np.ndarray:
    """Replace NaN values in a (n_time × n_signals) matrix before passing to
    Nilearn's ConnectivityMeasure, which cannot handle NaNs.

    Strategy (column-wise, conservative):
      - Columns that are entirely NaN → replaced with zeros.
        These carry no information; zeroing them is equivalent to
        excluding them from the partial correlation without altering
        the matrix shape or index mapping.
      - Columns with *some* NaN timepoints → replaced with the column mean.
        This minimises distortion of the covariance structure relative to
        alternatives such as linear interpolation or global-mean imputation.

    A warning is printed for every affected column so the operator can
    identify problematic parcels/seeds in the log.

    Parameters
    ----------
    X     : array, shape (n_time, n_signals)
    label : optional string prefix for warning messages (e.g. subject + hemi)

    Returns
    -------
    X_clean : array, same shape, no NaNs
    """
    X_clean   = X.copy()
    n_signals = X_clean.shape[1]

    for j in range(n_signals):
        col      = X_clean[:, j]
        nan_mask = np.isnan(col)

        if not nan_mask.any():
            continue

        n_nan = int(nan_mask.sum())

        if nan_mask.all():
            print(
                f"  {'[' + label + '] ' if label else ''}⚠️  Column {j}: "
                f"ALL {n_nan} timepoints are NaN — replacing with zeros"
            )
            X_clean[:, j] = 0.0
        else:
            col_mean = np.nanmean(col)
            print(
                f"  {'[' + label + '] ' if label else ''}⚠️  Column {j}: "
                f"{n_nan}/{len(col)} NaN timepoints — imputing with column mean "
                f"({col_mean:.4f})"
            )
            X_clean[nan_mask, j] = col_mean

    return X_clean

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
#   ts_key    : key in from_91k_to_32k result dict when
#               return_concat_hemis=False ('data_L' / 'data_R')
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

    # return_concat_hemis=False → separate 'data_L' and 'data_R' arrays
    # return_32k_mask=True      → 'mask_32k' boolean array (True=cortex, False=medial wall)
    res = from_91k_to_32k(
        ts_img, ts_data_raw,
        return_concat_hemis=False,
        return_32k_mask=True,
    )

    mask_32k = res["mask_32k"]
    print(f"  mask_32k: {int(np.sum(~mask_32k))} medial-wall vertices masked per hemi")

    for h in HEMIS:
        label     = h["label"]     # "LH" or "RH"
        seed_key  = h["seed_key"]  # "lh" or "rh"
        atlas_key = h["atlas_key"] # "L"  or "R"
        ts_key    = h["ts_key"]    # "data_L" or "data_R"

        ts_data = res[ts_key]      # shape: (n_time, n_vertices_hemi)

        print(f"\n  [{label}] Timeseries shape: {ts_data.shape}")

        # -------------------------
        # Cluster timeseries (single hemisphere)
        # -------------------------

        cluster_ts_list   = []
        cluster_names_used = []

        for roi in clusters:
            mask_img, mask_data = load_surface(
                f"{main_data}/{subject}/91k/rest/seed/"
                f"{subject}_91k_intertask_Sac-Pur-pRF"
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
        # Parcel timeseries — BILATERAL conditioning set
        #
        # All 106 parcels (both hemispheres) are loaded into the conditioning
        # matrix X so that inter-hemispheric confounds are fully controlled.
        # We track two index sets per parcel name:
        #   bilateral_col_idx : column position in the full 106-parcel matrix
        #   ipsi_col_idx      : subset of the above that belongs to the
        #                       ipsilateral hemisphere (atlas_key), used to
        #                       read off the result values after partial corr.
        # -------------------------

        # Both hemisphere keys, ipsilateral first so the loop is explicit
        ATLAS_KEYS_BILATERAL = {"L": ["L", "R"], "R": ["R", "L"]}
        atlas_keys_ordered   = ATLAS_KEYS_BILATERAL[atlas_key]   # ipsi, contra

        parcel_ts_list    = []   # one entry per successfully loaded parcel×hemi
        parcel_names_used = []   # matching parcel name (same parcel can appear twice)
        parcel_hemi_used  = []   # "L" or "R" for each column
        ipsi_col_idx      = []   # column indices that are ipsilateral

        for ak in atlas_keys_ordered:
            ts_source = res["data_L" if ak == "L" else "data_R"]

            for parcel in parcels:
                mask_img, mask_data = load_surface(
                    f"{atlas_folder}/parcels/"
                    f"{ak}_{parcel}_ROI.shape.gii"
                )
                mask = mask_data.ravel()

                if not np.any(mask):
                    continue

                col_idx = len(parcel_ts_list)
                parcel_ts_list.append(ts_source[:, mask > 0].mean(axis=1))
                parcel_names_used.append(parcel)
                parcel_hemi_used.append(ak)

                if ak == atlas_key:
                    ipsi_col_idx.append(col_idx)

        if not parcel_ts_list:
            print(f"  [{label}] ⚠️  No valid parcel timeseries — skipping hemisphere")
            continue

        parcel_ts = np.column_stack(parcel_ts_list)

        n_ipsi   = sum(1 for h in parcel_hemi_used if h == atlas_key)
        n_contra = sum(1 for h in parcel_hemi_used if h != atlas_key)
        print(
            f"  [{label}] Loaded bilateral parcel matrix: "
            f"{parcel_ts.shape[1]} columns before exclusions "
            f"({n_ipsi} ipsi [{atlas_key}] + {n_contra} contra)"
        )

        # Impute NaNs in both timeseries matrices before entering Nilearn.
        # Cluster seeds are unlikely to have NaNs but checked defensively.
        cluster_ts = impute_nan_columns(cluster_ts, label=f"{subject} {label} seed")
        parcel_ts  = impute_nan_columns(parcel_ts,  label=f"{subject} {label} parcel")

        # -------------------------
        # Partial correlation
        #
        # Conditioning set per seed = all bilateral parcel columns MINUS:
        #   (a) the seed's own cluster parcels (ipsilateral only, by convention)
        #   (b) the target parcel being correlated (ipsilateral only)
        # Result values are read from the ipsilateral columns only.
        # -------------------------

        n_clusters_used = cluster_ts.shape[1]
        # Output matrix rows = clusters, cols = ipsilateral parcels only
        partial_matrix    = np.full((n_clusters_used, len(ipsi_col_idx)), np.nan)
        partial_matrix_fz = np.full_like(partial_matrix, np.nan)

        # Build a name→column lookup restricted to ipsilateral parcels
        ipsi_parcel_names = [parcel_names_used[j] for j in ipsi_col_idx]

        for i_cl, cl_name in enumerate(cluster_names_used):

            # (a) seed's own cluster parcels excluded from conditioning (ipsi only)
            seed_own_parcels = set(seed_to_parcels.get(cl_name, []))

            for i_target, target_parcel in enumerate(ipsi_parcel_names):

                # (b) exclude the target parcel itself from the conditioning set
                exclude_cols = set()
                for j, (pname, phemi) in enumerate(
                    zip(parcel_names_used, parcel_hemi_used)
                ):
                    if phemi == atlas_key and (
                        pname in seed_own_parcels or pname == target_parcel
                    ):
                        exclude_cols.add(j)

                conditioning_idx = [
                    j for j in range(len(parcel_names_used))
                    if j not in exclude_cols
                    and j != ipsi_col_idx[i_target]   # target itself not in X
                ]

                if not conditioning_idx:
                    print(
                        f"  [{label}] ⚠️  No conditioning parcels left "
                        f"for seed {cl_name} → target {target_parcel}"
                    )
                    continue

                # X = [seed | target_parcel | conditioning parcels]
                # C[0,1] gives the seed ↔ target partial correlation
                X = np.column_stack([
                    cluster_ts[:, i_cl],
                    parcel_ts[:, ipsi_col_idx[i_target]],
                    parcel_ts[:, conditioning_idx],
                ])

                conn = ConnectivityMeasure(kind="partial correlation")
                C    = conn.fit_transform([X])[0]

                val = C[0, 1]   # seed ↔ target, conditioned on everything else

                partial_matrix[i_cl, i_target]    = val
                partial_matrix_fz[i_cl, i_target] = np.arctanh(val)

        print(
            f"  [{label}] Partial corr matrix: "
            f"{partial_matrix.shape[0]} clusters × {partial_matrix.shape[1]} ipsi parcels"
        )

        # -------------------------
        # Map back to full (clusters × parcels) grid
        # columns correspond to the ipsilateral parcels only
        # -------------------------

        filled    = np.full((len(clusters), len(parcels)), np.nan)
        filled_fz = filled.copy()

        for i_cl, cl in enumerate(cluster_names_used):
            gr = clusters.index(cl)
            for i_target, pa in enumerate(ipsi_parcel_names):
                gc = parcels.index(pa)
                filled[gr, gc]    = partial_matrix[i_cl, i_target]
                filled_fz[gr, gc] = partial_matrix_fz[i_cl, i_target]

        # -------------------------
        # Save subject-level results
        # -------------------------

        sub_out = f"{main_data}/{subject}/91k/rest/corr/partial_corr/by_hemi"
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
# GROUP STATS  (per hemisphere)
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