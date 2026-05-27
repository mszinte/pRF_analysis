#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Compute TASK-CONSTRAINED INTRA-HEMISPHERIC PARTIAL correlations between
ROI macro-regions (seeds) and TASK-DEFINED macro-region targets using
Nilearn ConnectivityMeasure.

Differs from the parcel-level (task-free) script in two key respects:

1. Targets are task-defined macro-regions (same binary masks as seeds,
   loaded from the /target/ folder), not atlas MMP parcels.
   Output matrix shape: (n_clusters × n_clusters).

2. Conditioning set = the OTHER macro-region timeseries within the same
   task-defined network (excluding seed and target).
   This mirrors the logic of the task-free script exactly:
     - Task-free:        condition on all MMP parcels except seed's own + target
     - Task-constrained: condition on all task macro-regions except seed + target
   Full vs partial correlations are therefore directly comparable within
   each analysis variant (task-free / task-constrained).

   Note: conditioning on atlas parcels outside the task network would change
   what is being measured and break the symmetry with the task-free analysis.

Design decisions
----------------
Seed timeseries   : ipsilateral hemisphere only (macro-region mask → mean vertex signal)
Target timeseries : ipsilateral (primary output) + contralateral (bilateral output)
Conditioning set  : remaining task macro-regions, ipsilateral only
                    (same hemisphere restriction as in the task-free script)

Self-correlation masking
  Seed == target entries are left as NaN (diagonal of the output matrix).

Standardize flag
  'zscore_sample' throughout: safe on already-z-scored XCP-D data (near-no-op
  since mean≈0, std≈1) and correct on non-z-scored individual runs.
  Forward-compatible with Nilearn ≥ 0.15.

Outputs (per subject, per run, per hemisphere)
  seed-task_by_macror-task_partial_{run}_{hemi}.npy / .csv           — Pearson r
  seed-task_by_macror-task_partial_fisherz_{run}_{hemi}.npy / .csv   — Fisher z
  seed-task_by_macror-task_partial_{run}_{hemi}_bilateral.npy / .csv — bilateral

Group aggregation: run group_stats_partial_corr.py (Fisher-z space throughout).

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

        run_tag      = f"_{run}" if run else "_concat"   # "" → "_concat" in filenames
        run_tag_ts   = f"_{run}" if run else ""          # "" for timeseries filename (no tag)

        print(f"\n=== Processing {subject}{run_tag} ===")
        
        timeseries_fn = (
            f"{main_data}/{subject}/91k/rest/timeseries/"
            f"{subject}_ses-01_task-rest{run_tag_ts}_space-fsLR_den-91k_desc-denoised_bold.dtseries.nii"
        )

        ts_img, ts_data_raw = load_surface(timeseries_fn)

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
            label      = h["label"]     # "LH" or "RH"
            seed_key   = h["seed_key"]  # "lh" or "rh"
            atlas_key  = h["atlas_key"] # "L"  or "R"
            ts_key     = h["ts_key"]    # "data_L" or "data_R"
            contra_key = "rh" if seed_key == "lh" else "lh"
            ts_key_contra = "data_R" if ts_key == "data_L" else "data_L"

            ts_data       = res[ts_key]        # (n_time, n_vertices_ipsi)
            ts_data_contra = res[ts_key_contra] # (n_time, n_vertices_contra)
            print(f"\n  [{label}] Timeseries shape: {ts_data.shape}")

            # ------------------------------------------------------
            # Step 1 — Load ALL task macro-region timeseries
            #
            # Seeds, ipsilateral targets, and contralateral targets are all
            # loaded up front into a single dict so the conditioning set
            # can be assembled cleanly in Step 4.
            #
            # Both /seed/ and /target/ use the same file naming convention;
            # they are loaded from separate folders to keep the analysis
            # intent explicit.
            # ------------------------------------------------------

            # macro_ts[roi] = ipsilateral timeseries for that macro-region
            macro_ts        = {}
            macro_ts_contra = {}   # contralateral timeseries
            loaded_rois        = []
            loaded_rois_contra = []

            for roi in clusters:

                # Ipsilateral (used as seed and ipsilateral conditioning)
                _, mask_data = load_surface(
                    f"{main_data}/{subject}/91k/rest/seed/"
                    f"{subject}_91k_intertask_Sac-Pur-pRF_{seed_key}_{roi}.shape.gii"
                )
                mask = mask_data.ravel()
                if np.any(mask):
                    macro_ts[roi] = ts_data[:, mask > 0].mean(axis=1)
                    loaded_rois.append(roi)
                else:
                    print(f"  [{label}] ⚠️  Empty ipsi mask for {roi} — skipping")

                # Contralateral (used as contra target and contralateral conditioning)
                _, mask_data = load_surface(
                    f"{main_data}/{subject}/91k/rest/target/"
                    f"{subject}_91k_intertask_Sac-Pur-pRF_{contra_key}_{roi}.shape.gii"
                )
                mask = mask_data.ravel()
                if np.any(mask):
                    macro_ts_contra[roi] = ts_data_contra[:, mask > 0].mean(axis=1)
                    loaded_rois_contra.append(roi)
                else:
                    print(f"  [{label}] ⚠️  Empty contra mask for {roi} — skipping")

            if not loaded_rois:
                print(f"  [{label}] ⚠️  No valid macro-region timeseries — skipping hemisphere")
                continue

            print(f"  [{label}] Ipsi macro-regions:  {loaded_rois}")
            print(f"  [{label}] Contra macro-regions: {loaded_rois_contra}")

            # ------------------------------------------------------
            # Step 2 — NaN imputation
            # ------------------------------------------------------

            for roi in list(macro_ts.keys()):
                clean = impute_nan_columns(
                    macro_ts[roi][:, np.newaxis],
                    label=f"{subject}{run_tag} {label} {roi}"
                )
                macro_ts[roi] = clean[:, 0]

            for roi in list(macro_ts_contra.keys()):
                clean = impute_nan_columns(
                    macro_ts_contra[roi][:, np.newaxis],
                    label=f"{subject}{run_tag} {label} contra {roi}"
                )
                macro_ts_contra[roi] = clean[:, 0]

            # ------------------------------------------------------
            # Step 3 — Partial correlations
            #
            # For every (seed, target) macro-region pair:
            #
            #   X = [seed | target | conditioning_macro-regions]
            #   C = partial_corr(X)
            #   result = C[0, 1]   ← seed ↔ target, all others partialled out
            #
            # Conditioning set = ALL OTHER macro-region timeseries, both
            # ipsilateral and contralateral, excluding seed and target.
            # This mirrors the task-free script where the bilateral MMP parcel
            # set is used as conditioning — inter-hemispheric confounds within
            # the task network are removed in the same way.
            #
            # Ipsilateral pairs:
            #   conditioning = remaining ipsi (excl. seed + target) + all contra
            #
            # Contralateral pairs:
            #   seed is ipsi, target is contra
            #   conditioning = remaining ipsi (excl. seed) + remaining contra (excl. target)
            #
            # Diagonal (seed == target) is left as NaN.
            # standardize='zscore_sample': safe on already-z-scored data
            # (near-no-op) and correct on non-z-scored individual runs.
            # ------------------------------------------------------

            n_rois            = len(clusters)
            partial_r_ipsi    = np.full((n_rois, n_rois), np.nan)
            partial_fz_ipsi   = np.full_like(partial_r_ipsi, np.nan)
            partial_r_contra  = np.full((n_rois, n_rois), np.nan)
            partial_fz_contra = np.full_like(partial_r_contra, np.nan)

            conn = ConnectivityMeasure(kind="partial correlation", standardize="zscore_sample")

            for seed_name in loaded_rois:
                i_seed = clusters.index(seed_name)

                # --- Ipsilateral targets ---
                for tgt_name in loaded_rois:

                    if seed_name == tgt_name:
                        continue   # diagonal: leave as NaN

                    i_tgt = clusters.index(tgt_name)

                    # Conditioning: remaining ipsi (excl. seed + target) + all contra
                    ipsi_cond  = [macro_ts[r]        for r in loaded_rois        if r not in (seed_name, tgt_name)]
                    contra_cond = [macro_ts_contra[r] for r in loaded_rois_contra]
                    conditioning = ipsi_cond + contra_cond

                    if not conditioning:
                        print(f"  [{label}] ⚠️  No conditioning regions for {seed_name} → {tgt_name}; computing full correlation")

                    X = np.column_stack([macro_ts[seed_name], macro_ts[tgt_name]] + conditioning)

                    C = conn.fit_transform([X])[0]
                    r = C[0, 1]
                    partial_r_ipsi[i_seed,  i_tgt] = r
                    partial_fz_ipsi[i_seed, i_tgt] = np.arctanh(r)

                # --- Contralateral targets ---
                for tgt_name in loaded_rois_contra:

                    i_tgt = clusters.index(tgt_name)

                    # Conditioning: remaining ipsi (excl. seed) + remaining contra (excl. target)
                    ipsi_cond   = [macro_ts[r]        for r in loaded_rois        if r != seed_name]
                    contra_cond = [macro_ts_contra[r] for r in loaded_rois_contra if r != tgt_name]
                    conditioning = ipsi_cond + contra_cond

                    if not conditioning:
                        print(f"  [{label}] ⚠️  No conditioning regions for {seed_name} → contra {tgt_name}; computing full correlation")

                    X = np.column_stack([macro_ts[seed_name], macro_ts_contra[tgt_name]] + conditioning)

                    C = conn.fit_transform([X])[0]
                    r = C[0, 1]
                    partial_r_contra[i_seed,  i_tgt] = r
                    partial_fz_contra[i_seed, i_tgt] = np.arctanh(r)

            print(
                f"  [{label}] Partial corr: "
                f"{len(loaded_rois)} × {len(loaded_rois)} ipsi  |  "
                f"{len(loaded_rois)} × {len(loaded_rois_contra)} contra"
            )

            # ------------------------------------------------------
            # Step 4 — Map to full output grids
            #
            # Primary:   (n_clusters × n_clusters), ipsilateral only
            # Bilateral: (n_clusters × 2*n_clusters), [ipsi | contra]
            #            ipsi half always recoverable as [:, :n_clusters]
            #            column labels: {cluster}_ipsi / {cluster}_contra
            # ------------------------------------------------------

            n_cl = len(clusters)

            filled_r  = partial_r_ipsi.copy()    # already (n_cl × n_cl)
            filled_fz = partial_fz_ipsi.copy()

            filled_r_bilateral  = np.full((n_cl, 2 * n_cl), np.nan)
            filled_fz_bilateral = np.full_like(filled_r_bilateral, np.nan)

            # Ipsi half
            filled_r_bilateral[:,  :n_cl] = partial_r_ipsi
            filled_fz_bilateral[:, :n_cl] = partial_fz_ipsi

            # Contra half
            filled_r_bilateral[:,  n_cl:] = partial_r_contra
            filled_fz_bilateral[:, n_cl:] = partial_fz_contra

            clusters_bilateral = (
                [f"{c}_ipsi"  for c in clusters] +
                [f"{c}_contra" for c in clusters]
            )

            # ------------------------------------------------------
            # Step 5 — Save subject-level outputs
            # ------------------------------------------------------

            sub_out = (
                f"{main_data}/{subject}/91k/rest/corr/partial_corr"
                f"/by_hemi/task-constrained"
            )
            os.makedirs(sub_out, exist_ok=True)

            tag = label.lower()   # "lh" or "rh"

            # Ipsilateral outputs (primary; consumed by downstream scripts)
            np.save(
                os.path.join(sub_out, f"seed-task_by_macror-task_partial{run_tag}_{tag}.npy"),
                filled_r,
            )
            np.save(
                os.path.join(sub_out, f"seed-task_by_macror-task_partial_fisherz{run_tag}_{tag}.npy"),
                filled_fz,
            )
            pd.DataFrame(filled_r,  index=clusters, columns=clusters).to_csv(
                os.path.join(sub_out, f"seed-task_by_macror-task_partial{run_tag}_{tag}.csv")
            )
            pd.DataFrame(filled_fz, index=clusters, columns=clusters).to_csv(
                os.path.join(sub_out, f"seed-task_by_macror-task_partial_fisherz{run_tag}_{tag}.csv")
            )

            # Bilateral outputs ([ipsi | contra]; ipsi half = [:, :n_clusters])
            np.save(
                os.path.join(sub_out, f"seed-task_by_macror-task_partial{run_tag}_{tag}_bilateral.npy"),
                filled_r_bilateral,
            )
            np.save(
                os.path.join(sub_out, f"seed-task_by_macror-task_partial_fisherz{run_tag}_{tag}_bilateral.npy"),
                filled_fz_bilateral,
            )
            pd.DataFrame(filled_r_bilateral,  index=clusters, columns=clusters_bilateral).to_csv(
                os.path.join(sub_out, f"seed-task_by_macror-task_partial{run_tag}_{tag}_bilateral.csv")
            )
            pd.DataFrame(filled_fz_bilateral, index=clusters, columns=clusters_bilateral).to_csv(
                os.path.join(sub_out, f"seed-task_by_macror-task_partial_fisherz{run_tag}_{tag}_bilateral.csv")
            )

            print(f"  [{label}] Saved to {sub_out}")

print("\nDone. Run group_stats_partial_corr.py to aggregate across subjects.")
# ============================================================