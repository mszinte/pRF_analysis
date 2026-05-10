#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Compute INTRA-HEMISPHERIC PARTIAL correlations between ROI clusters (seeds)
and MMP parcels (targets) using Nilearn ConnectivityMeasure

Design decisions
----------------
Seed timeseries  : ipsilateral hemisphere only (cluster mask → mean vertex signal)
Target parcels   : ipsilateral hemisphere only (results are reported per hemisphere)
                   note that in this script targets are not constrained by task results but defined by the whole MMP parcels
Conditioning set : ALL 106 parcels (both hemispheres), following Dawson et al. (2016)
                   and Genç et al. (2016).  Controlling for contralateral parcels
                   removes inter-hemispheric confounds that would otherwise inflate
                   apparent ipsilateral connectivity.

Exclusions from the conditioning set (ipsilateral only):
  (a) Parcels belonging to the seed's own cluster — these share signal with the seed
      by definition and would otherwise dominate the partial correlation.
  (b) The target parcel itself — it cannot condition on itself.

Self-correlation masking
  A seed cluster and a target parcel that belongs to the SAME cluster will share
  vertices, so their partial correlation will be spuriously high (approaching 1).
  These entries are left as NaN in the output rather than filled with a misleading
  value.  A separate cluster-by-cluster summary matrix is saved as a sanity check;
  those values should indeed be close to 1 for on-diagonal entries.

Standardize flag
  ConnectivityMeasure(standardize=False) because the concatenated XCP-D timeseries
  are already z-scored.  Using standardize=True on pre-standardised data would
  re-scale and distort the covariance structure.

Outputs (per subject, per run, per hemisphere)
  cluster_by_mmp-parcel_partial_{run}_{hemi}.npy / .csv          — Pearson r
  cluster_by_mmp-parcel_partial_fisherz_{run}_{hemi}.npy / .csv  — Fisher z
  cluster_by_cluster_partial_{run}_{hemi}.npy / .csv             — sanity check

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

# Reverse lookup: parcel name → its parent cluster
parcel_to_cluster = {
    p: cl
    for cl, plist in seed_to_parcels.items()
    for p in plist
}

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
# Helper: NaN imputation
# ============================================================

def impute_nan_columns(X: np.ndarray, label: str = "") -> np.ndarray:
    """Replace NaN values in a (n_time × n_signals) matrix.

    Nilearn's ConnectivityMeasure cannot handle NaN input.

    Strategy (column-wise):
      - All-NaN column  → replaced with zeros.  The column carries no information;
                          zeroing is equivalent to removing it from the partial
                          correlation without altering the matrix shape or index mapping.
      - Partial-NaN col → replaced with the column mean.  Minimal distortion of the
                          covariance structure relative to interpolation or global-mean
                          imputation.

    A warning is printed for every affected column so affected parcels/subjects
    can be audited.  Expected to occur only for sub-22 run-02 where bbregister failed
    And sub-25 / parcel 6mp due to FOV
    """
    X_clean = X.copy()

    for j in range(X_clean.shape[1]):
        col      = X_clean[:, j]
        nan_mask = np.isnan(col)

        if not nan_mask.any():
            continue

        n_nan = int(nan_mask.sum())
        prefix = f"  {'[' + label + '] ' if label else ''}⚠️  Column {j}: "

        if nan_mask.all():
            print(f"{prefix}ALL {n_nan} timepoints NaN — replacing with zeros")
            X_clean[:, j] = 0.0
        else:
            col_mean = float(np.nanmean(col))
            print(f"{prefix}{n_nan}/{len(col)} NaN timepoints — imputing with column mean ({col_mean:.4f})")
            X_clean[nan_mask, j] = col_mean

    return X_clean

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

        # Separate LH and RH 32k surface arrays
        # return_32k_mask=True returns a boolean mask (True=cortex, False=medial wall)
        # for downstream QC — medial wall vertices are never included in parcel masks
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
            # Step 1 — Seed (cluster) timeseries, ipsilateral only
            #
            # Each cluster has a pre-computed binary mask (.shape.gii).
            # Mean signal across all vertices in the mask → one timeseries per cluster.
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
            # bilateral conditioning set.  For each column we track:
            #   parcel_names_used : parcel name
            #   parcel_hemi_used  : "L" or "R"
            #   ipsi_col_idx      : column indices belonging to the ipsilateral hemisphere
            #                       (these are the targets and the result columns)
            # ------------------------------------------------------

            # Ipsilateral hemisphere is listed first so ipsi_col_idx entries
            # are always the lowest indices — easier to audit
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
            # ------------------------------------------------------

            cluster_ts = impute_nan_columns(cluster_ts, label=f"{subject}{run_tag} {label} seed")
            parcel_ts  = impute_nan_columns(parcel_ts,  label=f"{subject}{run_tag} {label} parcel")

            # ------------------------------------------------------
            # Step 4 — Partial correlations
            #
            # For every (seed_cluster, target_parcel) pair:
            #
            #   X = [seed | target | conditioning_parcels]
            #   C = partial_corr(X)
            #   result = C[0, 1]   ← seed ↔ target, all others partialled out
            #
            # Exclusions from the conditioning set (ipsilateral only):
            #   (a) seed's own cluster parcels — shared signal with seed
            #   (b) the target parcel itself   — cannot condition on the target
            #
            # SELF-CORRELATION MASKING:
            #   When target_parcel belongs to seed's own cluster, the seed timeseries
            #   and target timeseries are derived from overlapping vertices, so the
            #   partial correlation approaches 1 trivially.  These entries are left
            #   as NaN.  They are summarised separately in the cluster-by-cluster
            #   matrix (Step 5) as a sanity check.
            #
            # Note: seed-own parcels are pre-computed outside the target loop
            # (they are constant for a given seed) to avoid rebuilding the same
            # set on every iteration.
            # ------------------------------------------------------

            ipsi_parcel_names = [parcel_names_used[j] for j in ipsi_col_idx]
            n_ipsi_parcels    = len(ipsi_col_idx)
            n_clusters_used   = cluster_ts.shape[1]

            # Result matrices: rows = seeds, cols = ipsilateral targets
            partial_r  = np.full((n_clusters_used, n_ipsi_parcels), np.nan)
            partial_fz = np.full_like(partial_r, np.nan)

            for i_cl, cl_name in enumerate(cluster_names_used):

                # (a) parcels belonging to this seed's cluster (ipsilateral only)
                seed_own_parcels = set(seed_to_parcels.get(cl_name, []))

                # Conditioning columns that are excluded regardless of target:
                # ipsilateral columns belonging to the seed's own cluster.
                # (b) target-specific exclusion is added inside the target loop.
                seed_own_cols = {
                    j for j, (pname, phemi) in enumerate(
                        zip(parcel_names_used, parcel_hemi_used)
                    )
                    if phemi == atlas_key and pname in seed_own_parcels
                }

                for i_target, target_parcel in enumerate(ipsi_parcel_names):

                    # Self-correlation: target belongs to the seed's own cluster.
                    # Leave as NaN — reported separately in the cluster-by-cluster matrix.
                    if target_parcel in seed_own_parcels:
                        continue

                    # (b) also exclude the target parcel's own ipsilateral column
                    target_col = ipsi_col_idx[i_target]
                    exclude_cols = seed_own_cols | {target_col}

                    conditioning_idx = [
                        j for j in range(len(parcel_names_used))
                        if j not in exclude_cols
                    ]

                    if not conditioning_idx:
                        print(
                            f"  [{label}] ⚠️  No conditioning parcels left "
                            f"for {cl_name} → {target_parcel}"
                        )
                        continue

                    # X layout: column 0 = seed, column 1 = target, columns 2: = conditioning
                    X = np.column_stack([
                        cluster_ts[:, i_cl],
                        parcel_ts[:, target_col],
                        parcel_ts[:, conditioning_idx],
                    ])

                    # standardize=False: XCP-D timeseries are already z-scored
                    conn = ConnectivityMeasure(kind="partial correlation", standardize=False)
                    C    = conn.fit_transform([X])[0]

                    r  = C[0, 1]   # seed ↔ target partial correlation
                    partial_r[i_cl,  i_target] = r
                    partial_fz[i_cl, i_target] = np.arctanh(r)

            print(
                f"  [{label}] Partial corr matrix: "
                f"{n_clusters_used} clusters × {n_ipsi_parcels} ipsi parcels"
            )

            # ------------------------------------------------------
            # Step 5 — Cluster-by-cluster sanity check matrix
            #
            # Average partial_fz values over parcels belonging to the TARGET macro-region (e.g., sPCS)
            # On-diagonal entries (seed cluster == target cluster) will be NaN
            # because those pairs were masked in Step 4.
            # Off-diagonal entries give the mean partial correlation between
            # seed cluster i and target cluster j — useful for QC and for
            # reporting cluster-level connectivity alongside parcel-level detail.
            # ------------------------------------------------------

            n_cl = len(clusters)
            cluster_by_cluster_fz = np.full((n_cl, n_cl), np.nan)

            for i_cl, cl_seed in enumerate(cluster_names_used):
                gr = clusters.index(cl_seed)

                for j_cl, cl_target in enumerate(clusters):
                    target_parcel_names = seed_to_parcels.get(cl_target, [])

                    # Column indices in ipsi_parcel_names that belong to this target cluster
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
            # Step 6 — Map results to full (n_clusters × n_parcels) output grid
            #
            # The computation used only cluster_names_used (may be a subset of
            # clusters if some seed masks were empty). Map back to the full
            # ordered grid so output shape is always consistent.
            # ------------------------------------------------------

            filled_r  = np.full((len(clusters), len(parcels)), np.nan)
            filled_fz = np.full_like(filled_r, np.nan)

            for i_cl, cl in enumerate(cluster_names_used):
                gr = clusters.index(cl)
                for i_target, pa in enumerate(ipsi_parcel_names):
                    gc = parcels.index(pa)
                    filled_r[gr, gc]  = partial_r[i_cl,  i_target]
                    filled_fz[gr, gc] = partial_fz[i_cl, i_target]

            # ------------------------------------------------------
            # Step 7 — Save subject-level outputs
            # ------------------------------------------------------

            sub_out = f"{main_data}/{subject}/91k/rest/corr/partial_corr/by_hemi"
            os.makedirs(sub_out, exist_ok=True)

            tag = label.lower()   # "lh" or "rh"

            # Parcel-level partial correlation (r and Fisher-z)
            np.save(os.path.join(sub_out, f"cluster_by_mmp-parcel_partial{run_tag}_{tag}.npy"), filled_r)
            np.save(os.path.join(sub_out, f"cluster_by_mmp-parcel_partial_fisherz{run_tag}_{tag}.npy"), filled_fz)

            pd.DataFrame(filled_r,  index=clusters, columns=parcels).to_csv(
                os.path.join(sub_out, f"cluster_by_mmp-parcel_partial{run_tag}_{tag}.csv")
            )
            pd.DataFrame(filled_fz, index=clusters, columns=parcels).to_csv(
                os.path.join(sub_out, f"cluster_by_mmp-parcel_partial_fisherz{run_tag}_{tag}.csv")
            )

            # Cluster-by-cluster sanity check (Fisher-z only)
            np.save(os.path.join(sub_out, f"cluster_by_cluster_partial_fisherz{run_tag}_{tag}.npy"), cluster_by_cluster_fz)

            pd.DataFrame(cluster_by_cluster_fz, index=clusters, columns=clusters).to_csv(
                os.path.join(sub_out, f"cluster_by_cluster_partial_fisherz{run_tag}_{tag}.csv")
            )

            print(f"  [{label}] Saved to {sub_out}")

print("\nDone. Run group_stats_partial_corr.py to aggregate across subjects.")
# ============================================================