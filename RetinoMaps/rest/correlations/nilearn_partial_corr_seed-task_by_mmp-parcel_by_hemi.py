#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compute INTRA-HEMISPHERIC PARTIAL correlations between ROI macro-regions (seeds)
and MMP parcels (targets) using Nilearn ConnectivityMeasure.

Design decisions
----------------
Seed timeseries  : ipsilateral hemisphere only (macro-region mask → mean vertex signal)
Target parcels   : ipsilateral hemisphere only (results are reported per hemisphere)
                   targets are not constrained by task results but span all MMP parcels
Conditioning set : ALL 106 parcels (both hemispheres), following Dawson et al. (2016)
                   and Genç et al. (2016).  Controlling for contralateral parcels
                   removes inter-hemispheric confounds that would otherwise inflate
                   apparent ipsilateral connectivity.

Exclusions from the conditioning set (ipsilateral only):
  (a) Parcels belonging to the seed's own macro-region — these share signal with
      the seed by definition and would otherwise dominate the partial correlation.
  (b) The target parcel itself — it cannot condition on itself.

Self-correlation masking
  A seed macro-region and a target parcel that belongs to the SAME macro-region
  will share vertices, so their partial correlation will be spuriously high
  (approaching 1).  These entries are left as NaN in the output rather than
  filled with a misleading value.  A separate macro-region-by-macro-region summary
  matrix is saved as a sanity check; off-diagonal entries should be interpretable
  connectivity estimates, diagonal entries are NaN by construction.

Standardize flag
  ConnectivityMeasure(standardize=False) because the concatenated XCP-D timeseries
  are already z-scored.  Using standardize=True on pre-standardised data would
  re-scale and distort the covariance structure.

Outputs (per subject, per run, per hemisphere)
  seed-task_by_mmp-parcel_partial_{run}_{hemi}.npy / .csv          — Pearson r
  seed-task_by_mmp-parcel_partial_fisherz_{run}_{hemi}.npy / .csv  — Fisher z
  seed-task_by_macro-region_partial_fisherz_{run}_{hemi}.npy / .csv — sanity check

Group aggregation is handled by group_stats_partial_corr.py (always in Fisher-z
space; back-transformed to r only at the reporting stage).

---------------------------------------------------
Written by Marco Bedini (marco.bedini@univ-amu.fr)
---------------------------------------------------
"""

import os
import sys
import numpy as np
import pandas as pd
from nilearn.connectome import ConnectivityMeasure

# ============================================================
# Paths
# ============================================================

main_data    = "/scratch/mszinte/data/RetinoMaps/derivatives/pp_data"
atlas_folder = "/scratch/mszinte/data/RetinoMaps/derivatives/pp_data/atlas"

base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../"))

sys.path.append(os.path.abspath(os.path.join(base_dir, "analysis_code/utils")))
from settings_utils import load_settings
from surface_utils import load_surface
from cifti_utils import from_91k_to_32k

sys.path.append(os.path.abspath(os.path.join(base_dir, "RetinoMaps/rest/utils")))
from rest_utils import impute_nan_columns

# ============================================================
# Settings
# ============================================================

project_dir       = "RetinoMaps"
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

# Reverse so mPCS is first (matches downstream visualisation scripts)
clusters.reverse()

# Flat ordered parcel list — defines the column order of all output matrices
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

RUNS = ["run-01", "run-02", ""]   # "" = full concatenated session

# ============================================================
# Subject loop
# ============================================================

for subject in subjects:

    for run in RUNS:

        run_tag = f"_{run}" if run else ""

        print(f"\n=== Processing {subject}{run_tag} ===")

        timeseries_fn = (
            f"{main_data}/{subject}/91k/rest/timeseries/"
            f"{subject}_ses-01_task-rest{run_tag}_space-fsLR_den-91k_desc-denoised_bold.dtseries.nii"
        )

        ts_img, ts_data_raw = load_surface(timeseries_fn)

        # Separate LH and RH 32k surface arrays.
        # return_32k_mask=True returns a boolean mask (True=cortex, False=medial wall)
        # for downstream QC — medial wall vertices are never included in parcel masks.
        res = from_91k_to_32k(
            ts_img, ts_data_raw,
            return_concat_hemis=False,
            return_32k_mask=True,
        )

        mask_32k = res["mask_32k"]
        print(f"  mask_32k: {int(np.sum(~mask_32k))} medial-wall vertices per hemi")

        # ----------------------------------------------------------
        # Hemisphere loop
        # ----------------------------------------------------------

        for h in HEMIS:
            label     = h["label"]     # "LH" or "RH"
            seed_key  = h["seed_key"]  # "lh" or "rh"  (seed filename suffix)
            atlas_key = h["atlas_key"] # "L"  or "R"   (atlas filename prefix)
            ts_key    = h["ts_key"]    # "data_L" or "data_R"

            ts_data = res[ts_key]   # (n_time, n_vertices_hemi)
            print(f"\n  [{label}] Timeseries shape: {ts_data.shape}")

            # ------------------------------------------------------
            # Step 1 — Seed (macro-region) timeseries, ipsilateral only
            #
            # Each macro-region has a pre-computed binary mask (.shape.gii).
            # Mean signal across all mask vertices → one timeseries per macro-region.
            # ------------------------------------------------------

            cluster_ts_list    = []
            cluster_names_used = []

            for roi in clusters:
                _, mask_data = load_surface(
                    f"{main_data}/{subject}/91k/rest/seed/"
                    f"{subject}_91k_intertask_Sac-Pur-pRF_{seed_key}_{roi}.shape.gii"
                )
                mask = mask_data.ravel()

                if not np.any(mask):
                    print(f"  [{label}] ⚠️  Empty seed mask for {roi} — skipping")
                    continue

                cluster_ts_list.append(ts_data[:, mask > 0].mean(axis=1))
                cluster_names_used.append(roi)

            if not cluster_ts_list:
                print(f"  [{label}] ⚠️  No valid seed timeseries — skipping hemisphere")
                continue

            cluster_ts = np.column_stack(cluster_ts_list)
            print(f"  [{label}] Seeds loaded: {cluster_names_used}")

            # ------------------------------------------------------
            # Step 2 — Parcel timeseries, BOTH hemispheres
            #
            # All 106 parcels (53 ipsi + 53 contra) are loaded to form the
            # bilateral conditioning set.  Three tracking lists are built:
            #   parcel_names_used : parcel name for each matrix column
            #   parcel_hemi_used  : "L" or "R" for each column
            #   ipsi_col_idx      : column indices belonging to the ipsilateral
            #                       hemisphere — these are the targets and define
            #                       the result columns in the output matrix
            #
            # Ipsilateral columns are listed first so ipsi_col_idx entries
            # are always the lowest indices (easier to audit in logs).
            # ------------------------------------------------------

            atlas_keys_ordered = [atlas_key, "R" if atlas_key == "L" else "L"]

            parcel_ts_list    = []
            parcel_names_used = []
            parcel_hemi_used  = []
            ipsi_col_idx      = []

            for ak in atlas_keys_ordered:
                ts_source = res["data_L" if ak == "L" else "data_R"]

                for parcel in parcels:
                    _, mask_data = load_surface(
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
                print(f"  [{label}] ⚠️  No valid parcel timeseries — skipping hemisphere")
                continue

            parcel_ts = np.column_stack(parcel_ts_list)

            n_ipsi   = sum(1 for hm in parcel_hemi_used if hm == atlas_key)
            n_contra = sum(1 for hm in parcel_hemi_used if hm != atlas_key)
            print(
                f"  [{label}] Bilateral parcel matrix: {parcel_ts.shape[1]} columns "
                f"({n_ipsi} ipsi [{atlas_key}] + {n_contra} contra)"
            )

            # ------------------------------------------------------
            # Step 3 — NaN imputation (must happen before Nilearn)
            #
            # Expected cases: sub-22 run-02 (bbregister failure),
            #                 sub-25 parcel 6mp (FOV truncation).
            # ------------------------------------------------------

            cluster_ts = impute_nan_columns(cluster_ts, label=f"{subject}{run_tag} {label} seed")
            parcel_ts  = impute_nan_columns(parcel_ts,  label=f"{subject}{run_tag} {label} parcel")

            # ------------------------------------------------------
            # Step 4 — Partial correlations
            #
            # For every (seed_macro-region, target_parcel) pair:
            #
            #   X = [seed | target | conditioning_parcels]
            #   C = partial_corr(X)        ← via Nilearn ConnectivityMeasure
            #   result = C[0, 1]           ← seed ↔ target, all others partialled out
            #
            # Exclusions from the conditioning set (ipsilateral only):
            #   (a) Seed's own macro-region parcels — pre-computed once per seed
            #       as seed_own_cols (constant across all targets for this seed).
            #   (b) The target parcel column itself — added per iteration.
            #
            # SELF-CORRELATION MASKING:
            #   When the target parcel belongs to the seed's own macro-region, the
            #   seed and target timeseries are derived from overlapping vertices and
            #   their partial correlation is trivially high.  These entries are left
            #   as NaN and are not computed.  They are summarised in the macro-region-
            #   by-macro-region sanity check matrix (Step 5) instead.
            # ------------------------------------------------------

            ipsi_parcel_names = [parcel_names_used[j] for j in ipsi_col_idx]
            n_ipsi_parcels    = len(ipsi_col_idx)
            n_clusters_used   = cluster_ts.shape[1]

            # Result matrices: rows = macro-regions, cols = ipsilateral target parcels
            partial_r  = np.full((n_clusters_used, n_ipsi_parcels), np.nan)
            partial_fz = np.full_like(partial_r, np.nan)

            # Contralateral column indices and names (complement of ipsi_col_idx)
            ipsi_col_set        = set(ipsi_col_idx)
            contra_col_idx      = [j for j in range(len(parcel_names_used)) if j not in ipsi_col_set]
            contra_parcel_names = [parcel_names_used[j] for j in contra_col_idx]
            n_contra_parcels    = len(contra_col_idx)

            # Contralateral result matrices
            partial_r_contra  = np.full((n_clusters_used, n_contra_parcels), np.nan)
            partial_fz_contra = np.full_like(partial_r_contra, np.nan)

            # Single ConnectivityMeasure instance reused for every (seed, target) pair
            # Standardize set to false because now XCP-D z-scores the concat runs
            conn = ConnectivityMeasure(kind="partial correlation", standardize=False)

            for i_cl, cl_name in enumerate(cluster_names_used):

                # (a) ipsilateral columns belonging to this seed's own macro-region
                seed_own_parcels = set(seed_to_parcels.get(cl_name, []))
                seed_own_cols = {
                    j for j, (pname, phemi) in enumerate(
                        zip(parcel_names_used, parcel_hemi_used)
                    )
                    if phemi == atlas_key and pname in seed_own_parcels
                }

                # --- Ipsilateral targets ---
                for i_target, target_parcel in enumerate(ipsi_parcel_names):

                    # Self-correlation: skip, leave as NaN
                    if target_parcel in seed_own_parcels:
                        continue

                    # (b) exclude the target column itself from conditioning
                    target_col   = ipsi_col_idx[i_target]
                    exclude_cols = seed_own_cols | {target_col}

                    conditioning_idx = [
                        j for j in range(len(parcel_names_used))
                        if j not in exclude_cols
                    ]

                    if not conditioning_idx:
                        print(f"  [{label}] ⚠️  No conditioning parcels left for {cl_name} → {target_parcel}")
                        continue

                    # Column 0 = seed, column 1 = target, columns 2: = conditioning parcels
                    X = np.column_stack([
                        cluster_ts[:, i_cl],
                        parcel_ts[:, target_col],
                        parcel_ts[:, conditioning_idx],
                    ])

                    C = conn.fit_transform([X])[0]

                    r = C[0, 1]
                    partial_r[i_cl,  i_target] = r
                    partial_fz[i_cl, i_target] = np.arctanh(r)

                # --- Contralateral targets ---
                # No self-correlation masking needed (seed is ipsilateral, target is
                # contralateral — they cannot share vertices).
                # Conditioning set excludes seed's own ipsilateral parcels and the
                # target column, same logic as above.
                for i_target, target_parcel in enumerate(contra_parcel_names):

                    target_col   = contra_col_idx[i_target]
                    exclude_cols = seed_own_cols | {target_col}

                    conditioning_idx = [
                        j for j in range(len(parcel_names_used))
                        if j not in exclude_cols
                    ]

                    if not conditioning_idx:
                        print(f"  [{label}] ⚠️  No conditioning parcels left for {cl_name} → contra {target_parcel}")
                        continue

                    X = np.column_stack([
                        cluster_ts[:, i_cl],
                        parcel_ts[:, target_col],
                        parcel_ts[:, conditioning_idx],
                    ])

                    C = conn.fit_transform([X])[0]

                    r = C[0, 1]
                    partial_r_contra[i_cl,  i_target] = r
                    partial_fz_contra[i_cl, i_target] = np.arctanh(r)

            print(
                f"  [{label}] Partial corr: "
                f"{n_clusters_used} seeds × {n_ipsi_parcels} ipsi + {n_contra_parcels} contra parcels"
            )

            # ------------------------------------------------------
            # Step 5 — Macro-region-by-macro-region sanity check
            #
            # Averages partial_fz over the target parcels belonging to each
            # macro-region, giving an (n_clusters × n_clusters) summary matrix.
            #   - Diagonal entries are NaN by construction (self-correlation mask).
            #   - Off-diagonal entries are mean partial connectivity between
            #     macro-region pairs — used for QC and cluster-level reporting.
            # ------------------------------------------------------

            n_cl = len(clusters)
            cluster_by_cluster_fz = np.full((n_cl, n_cl), np.nan)

            for i_cl, cl_seed in enumerate(cluster_names_used):
                gr = clusters.index(cl_seed)

                for j_cl, cl_target in enumerate(clusters):
                    target_parcel_names = seed_to_parcels.get(cl_target, [])

                    target_col_idx = [
                        i for i, pname in enumerate(ipsi_parcel_names)
                        if pname in target_parcel_names
                    ]

                    if not target_col_idx:
                        continue

                    vals = partial_fz[i_cl, target_col_idx]
                    if not np.all(np.isnan(vals)):
                        cluster_by_cluster_fz[gr, j_cl] = np.nanmean(vals)

            # ------------------------------------------------------
            # Step 6 — Map to full output grids
            #
            # Two grids are produced:
            #
            #   filled_r / filled_fz  (n_clusters × n_parcels)
            #     Ipsilateral only. Columns = parcels in canonical YAML order.
            #     This is the primary output and is what all downstream scripts
            #     (group stats, visualisation) currently consume.
            #
            #   filled_r_bilateral / filled_fz_bilateral  (n_clusters × 2*n_parcels)
            #     Columns = [ipsi_parcels | contra_parcels], both in canonical
            #     YAML order.  The ipsilateral half is always [:, :n_parcels],
            #     so any downstream script can recover it without change.
            #     Saved with a _bilateral suffix to avoid breaking existing loaders.
            #
            # cluster_names_used may be a subset of clusters if any seed masks
            # were empty — both grids are initialised to NaN so missing rows
            # are explicitly absent rather than silently zero.
            # ------------------------------------------------------

            n_parcels_total = len(parcels)

            filled_r  = np.full((len(clusters), n_parcels_total), np.nan)
            filled_fz = np.full_like(filled_r, np.nan)

            # Bilateral: [ipsi | contra], both in parcels order
            filled_r_bilateral  = np.full((len(clusters), 2 * n_parcels_total), np.nan)
            filled_fz_bilateral = np.full_like(filled_r_bilateral, np.nan)

            # Column labels for the bilateral DataFrame
            contra_key   = "R" if atlas_key == "L" else "L"
            parcels_bilateral = (
                [f"{p}_{atlas_key}" for p in parcels] +   # ipsi half
                [f"{p}_{contra_key}" for p in parcels]    # contra half
            )

            for i_cl, cl in enumerate(cluster_names_used):
                gr = clusters.index(cl)

                # Ipsilateral half
                for i_target, pa in enumerate(ipsi_parcel_names):
                    gc = parcels.index(pa)
                    filled_r[gr, gc]  = partial_r[i_cl,  i_target]
                    filled_fz[gr, gc] = partial_fz[i_cl, i_target]
                    # Bilateral ipsi half (columns 0 : n_parcels)
                    filled_r_bilateral[gr, gc]  = partial_r[i_cl,  i_target]
                    filled_fz_bilateral[gr, gc] = partial_fz[i_cl, i_target]

                # Contralateral half (columns n_parcels : 2*n_parcels)
                for i_target, pa in enumerate(contra_parcel_names):
                    gc = parcels.index(pa)
                    filled_r_bilateral[gr, n_parcels_total + gc]  = partial_r_contra[i_cl,  i_target]
                    filled_fz_bilateral[gr, n_parcels_total + gc] = partial_fz_contra[i_cl, i_target]

            # ------------------------------------------------------
            # Step 7 — Save subject-level outputs
            # ------------------------------------------------------

            sub_out = f"{main_data}/{subject}/91k/rest/corr/partial_corr/by_hemi/task-free"
            os.makedirs(sub_out, exist_ok=True)

            tag = label.lower()   # "lh" or "rh"

            # --- Ipsilateral outputs (primary; consumed by all downstream scripts) ---
            np.save(
                os.path.join(sub_out, f"seed-task_by_mmp-parcel_partial{run_tag}_{tag}.npy"),
                filled_r,
            )
            np.save(
                os.path.join(sub_out, f"seed-task_by_mmp-parcel_partial_fisherz{run_tag}_{tag}.npy"),
                filled_fz,
            )
            pd.DataFrame(filled_r,  index=clusters, columns=parcels).to_csv(
                os.path.join(sub_out, f"seed-task_by_mmp-parcel_partial{run_tag}_{tag}.csv")
            )
            pd.DataFrame(filled_fz, index=clusters, columns=parcels).to_csv(
                os.path.join(sub_out, f"seed-task_by_mmp-parcel_partial_fisherz{run_tag}_{tag}.csv")
            )

            # --- Bilateral outputs ([ipsi | contra]; ipsi half = [:, :n_parcels]) ---
            np.save(
                os.path.join(sub_out, f"seed-task_by_mmp-parcel_partial{run_tag}_{tag}_bilateral.npy"),
                filled_r_bilateral,
            )
            np.save(
                os.path.join(sub_out, f"seed-task_by_mmp-parcel_partial_fisherz{run_tag}_{tag}_bilateral.npy"),
                filled_fz_bilateral,
            )
            pd.DataFrame(filled_r_bilateral,  index=clusters, columns=parcels_bilateral).to_csv(
                os.path.join(sub_out, f"seed-task_by_mmp-parcel_partial{run_tag}_{tag}_bilateral.csv")
            )
            pd.DataFrame(filled_fz_bilateral, index=clusters, columns=parcels_bilateral).to_csv(
                os.path.join(sub_out, f"seed-task_by_mmp-parcel_partial_fisherz{run_tag}_{tag}_bilateral.csv")
            )

            # --- Macro-region-by-macro-region sanity check (Fisher-z, ipsi only) ---
            np.save(
                os.path.join(sub_out, f"seed-task_by_macro-region_partial_fisherz{run_tag}_{tag}.npy"),
                cluster_by_cluster_fz,
            )
            pd.DataFrame(cluster_by_cluster_fz, index=clusters, columns=clusters).to_csv(
                os.path.join(sub_out, f"seed-task_by_macro-region_partial_fisherz{run_tag}_{tag}.csv")
            )

            print(f"  [{label}] Saved to {sub_out}")

print("\nDone. Run group_stats_partial_corr.py to aggregate across subjects.")
# ============================================================