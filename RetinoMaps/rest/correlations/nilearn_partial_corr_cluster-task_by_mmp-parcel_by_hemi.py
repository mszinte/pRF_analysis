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

Outputs are also computed by run so we can exclude some based on the quality check info
(Either eye tracking data or denoising or preproc bugs)

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

main_data = "/scratch/mszinte/data/RetinoMaps/derivatives/pp_data"
seed_folder = main_data
atlas_folder = "/scratch/mszinte/data/RetinoMaps/derivatives/pp_data/atlas"
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
    """
    X_clean = X.copy()
    n_signals = X_clean.shape[1]

    for j in range(n_signals):
        col = X_clean[:, j]
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
settings_path = os.path.join(base_dir, project_dir, "settings.yml")
prf_settings_path = os.path.join(base_dir, project_dir, "prf-analysis.yml")
settings = load_settings([settings_path, prf_settings_path])
analysis_info = settings[0]
subjects = analysis_info["subjects"]

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
# ============================================================

HEMIS = [
    {"label": "LH", "seed_key": "lh", "atlas_key": "L", "ts_key": "data_L"},
    {"label": "RH", "seed_key": "rh", "atlas_key": "R", "ts_key": "data_R"},
]

# ============================================================
# Runs
# ============================================================

RUNS = ["run-01", "run-02", ""]

# =========================
# STORAGE (run × hemisphere)
# =========================

all_subject_partial = {
    run: {h["label"]: [] for h in HEMIS}
    for run in RUNS
}

all_subject_partial_fz = {
    run: {h["label"]: [] for h in HEMIS}
    for run in RUNS
}

# =========================
# SUBJECT LOOP
# =========================

for subject in subjects:

    for run in RUNS:

        run_tag = f"_{run}" if run else ""

        print(f"\n=== Processing {subject}{run_tag} ===")

        timeseries_fn = (
            f"{main_data}/{subject}/91k/rest/timeseries/"
            f"{subject}_ses-01_task-rest{run_tag}_space-fsLR_den-91k_desc-denoised_bold.dtseries.nii"
        )

        ts_img, ts_data_raw = load_surface(timeseries_fn)

        res = from_91k_to_32k(
            ts_img,
            ts_data_raw,
            return_concat_hemis=False,
            return_32k_mask=True,
        )

        mask_32k = res["mask_32k"]
        print(
            f"  mask_32k: {int(np.sum(~mask_32k))} medial-wall vertices masked per hemi"
        )

        for h in HEMIS:
            label = h["label"]
            seed_key = h["seed_key"]
            atlas_key = h["atlas_key"]
            ts_key = h["ts_key"]

            ts_data = res[ts_key]

            print(f"\n  [{label}] Timeseries shape: {ts_data.shape}")

            # -------------------------
            # Cluster timeseries
            # -------------------------

            cluster_ts_list = []
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
                print(
                    f"  [{label}] ⚠️  No valid cluster timeseries — skipping hemisphere"
                )
                continue

            cluster_ts = np.column_stack(cluster_ts_list)
            print(f"  [{label}] Clusters used: {cluster_names_used}")

            # -------------------------
            # Parcel timeseries (bilateral conditioning)
            # -------------------------

            ATLAS_KEYS_BILATERAL = {"L": ["L", "R"], "R": ["R", "L"]}
            atlas_keys_ordered = ATLAS_KEYS_BILATERAL[atlas_key]

            parcel_ts_list = []
            parcel_names_used = []
            parcel_hemi_used = []
            ipsi_col_idx = []

            for ak in atlas_keys_ordered:
                ts_source = res["data_L" if ak == "L" else "data_R"]

                for parcel in parcels:
                    mask_img, mask_data = load_surface(
                        f"{atlas_folder}/parcels/{ak}_{parcel}_ROI.shape.gii"
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
                print(
                    f"  [{label}] ⚠️  No valid parcel timeseries — skipping hemisphere"
                )
                continue

            parcel_ts = np.column_stack(parcel_ts_list)

            n_ipsi = sum(1 for h_ in parcel_hemi_used if h_ == atlas_key)
            n_contra = sum(1 for h_ in parcel_hemi_used if h_ != atlas_key)

            print(
                f"  [{label}] Loaded bilateral parcel matrix: "
                f"{parcel_ts.shape[1]} columns before exclusions "
                f"({n_ipsi} ipsi [{atlas_key}] + {n_contra} contra)"
            )

            cluster_ts = impute_nan_columns(
                cluster_ts, label=f"{subject}{run_tag} {label} seed"
            )
            parcel_ts = impute_nan_columns(
                parcel_ts, label=f"{subject}{run_tag} {label} parcel"
            )

            # -------------------------
            # Partial correlation
            # -------------------------

            n_clusters_used = cluster_ts.shape[1]

            partial_matrix = np.full(
                (n_clusters_used, len(ipsi_col_idx)),
                np.nan
            )
            partial_matrix_fz = np.full_like(partial_matrix, np.nan)

            ipsi_parcel_names = [parcel_names_used[j] for j in ipsi_col_idx]

            for i_cl, cl_name in enumerate(cluster_names_used):

                seed_own_parcels = set(seed_to_parcels.get(cl_name, []))

                for i_target, target_parcel in enumerate(ipsi_parcel_names):

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
                        and j != ipsi_col_idx[i_target]
                    ]

                    if not conditioning_idx:
                        print(
                            f"  [{label}] ⚠️  No conditioning parcels left "
                            f"for seed {cl_name} → target {target_parcel}"
                        )
                        continue

                    X = np.column_stack([
                        cluster_ts[:, i_cl],
                        parcel_ts[:, ipsi_col_idx[i_target]],
                        parcel_ts[:, conditioning_idx],
                    ])

                    conn = ConnectivityMeasure(kind="partial correlation")
                    C = conn.fit_transform([X])[0]

                    val = C[0, 1]

                    partial_matrix[i_cl, i_target] = val
                    partial_matrix_fz[i_cl, i_target] = np.arctanh(val)

            print(
                f"  [{label}] Partial corr matrix: "
                f"{partial_matrix.shape[0]} clusters × "
                f"{partial_matrix.shape[1]} ipsi parcels"
            )

            # -------------------------
            # Map back to full grid
            # -------------------------

            filled = np.full((len(clusters), len(parcels)), np.nan)
            filled_fz = filled.copy()

            for i_cl, cl in enumerate(cluster_names_used):
                gr = clusters.index(cl)

                for i_target, pa in enumerate(ipsi_parcel_names):
                    gc = parcels.index(pa)

                    filled[gr, gc] = partial_matrix[i_cl, i_target]
                    filled_fz[gr, gc] = partial_matrix_fz[i_cl, i_target]

            # -------------------------
            # Save subject-level results
            # -------------------------

            sub_out = f"{main_data}/{subject}/91k/rest/corr/partial_corr/by_hemi"
            os.makedirs(sub_out, exist_ok=True)

            tag = label.lower()

            np.save(
                os.path.join(
                    sub_out,
                    f"cluster_by_mmp-parcel_partial{run_tag}_{tag}.npy"
                ),
                filled,
            )

            np.save(
                os.path.join(
                    sub_out,
                    f"cluster_by_mmp-parcel_partial_fisherz{run_tag}_{tag}.npy"
                ),
                filled_fz,
            )

            pd.DataFrame(
                filled,
                index=clusters,
                columns=parcels
            ).to_csv(
                os.path.join(
                    sub_out,
                    f"cluster_by_mmp-parcel_partial{run_tag}_{tag}.csv"
                )
            )

            pd.DataFrame(
                filled_fz,
                index=clusters,
                columns=parcels
            ).to_csv(
                os.path.join(
                    sub_out,
                    f"cluster_by_mmp-parcel_partial_fisherz{run_tag}_{tag}.csv"
                )
            )

            all_subject_partial[run][label].append(filled)
            all_subject_partial_fz[run][label].append(filled_fz)

print("\nDone.")
# ============================================================