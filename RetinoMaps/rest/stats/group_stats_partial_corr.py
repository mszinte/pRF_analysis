#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Group-level statistics for intra-hemispheric partial correlations.

Operates on per-subject Fisher-z matrices produced by the computation script.
Averaging is always performed in Fisher-z space; Pearson r is recovered only
at the final reporting stage via tanh().

This separation is intentional:
  - Fisher-z has approximately constant variance ~1/(n-3), making it the
    correct space for averaging and any subsequent parametric tests
  - Raw Pearson r values have r-dependent variance and should never be
    averaged directly

Outputs per hemisphere × run:
  - group mean and median in Fisher-z space      (.npy, .csv)
  - group mean and median back-converted to r    (.npy, .csv)
  - compressed archive with all four arrays      (.npz)

---------------------------------------------------
Written by Marco Bedini (marco.bedini@univ-amu.fr)
---------------------------------------------------
"""

import os
import sys
import numpy as np
import pandas as pd

# ============================================================
# Paths
# ============================================================

main_data             = "/scratch/mszinte/data/RetinoMaps/derivatives/pp_data"
partial_output_folder = (
    "/scratch/mszinte/data/RetinoMaps/derivatives/pp_data"
    "/group/91k/rest/partial_corr/by_hemi"
)
os.makedirs(partial_output_folder, exist_ok=True)

# Personal imports
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../"))
sys.path.append(os.path.abspath(os.path.join(base_dir, "analysis_code/utils")))
from settings_utils import load_settings

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
# ROIs — needed for DataFrame index/columns and .npz metadata
# ============================================================

clusters = analysis_info["rois-drawn"]
seed_to_parcels = analysis_info["rois-group-mmp"]

clusters.reverse()   # mPCS first

parcels = []
for cl in clusters:
    parcels.extend(seed_to_parcels[cl])

# ============================================================
# Configuration
# ============================================================

HEMIS = [
    {"label": "LH", "tag": "lh"},
    {"label": "RH", "tag": "rh"},
]

RUNS = ["run-01", "run-02", ""]   # "" = full concatenated session

# ============================================================
# GROUP STATS
# ============================================================

print("=== GROUP STATS ===\n")

for run in RUNS:

    run_tag = f"_{run}" if run else ""

    for h in HEMIS:
        label = h["label"]
        tag   = h["tag"]

        # -------------------------
        # Load per-subject Fisher-z matrices
        # -------------------------

        fz_stack = []
        missing  = []

        for subject in subjects:
            sub_path = os.path.join(
                main_data, subject,
                "91k/rest/corr/partial_corr/by_hemi",
                f"cluster_by_mmp-parcel_partial_fisherz{run_tag}_{tag}.npy"
            )

            if not os.path.isfile(sub_path):
                print(f"  [{label}{run_tag}] ⚠️  Missing: {subject} — skipped")
                missing.append(subject)
                continue

            fz_stack.append(np.load(sub_path))

        if not fz_stack:
            print(
                f"  [{label}{run_tag}] ⚠️  No subject files found — "
                "skipping group stats"
            )
            continue

        n_loaded = len(fz_stack)
        n_total  = len(subjects)
        print(
            f"  [{label}{run_tag}] Loaded {n_loaded}/{n_total} subjects"
            + (f" (missing: {missing})" if missing else "")
        )

        stacked_fz = np.stack(fz_stack, axis=0)   # (n_subjects, n_clusters, n_parcels)

        # -------------------------
        # Aggregate in Fisher-z space
        # -------------------------

        mean_fz   = np.nanmean(stacked_fz,   axis=0)
        median_fz = np.nanmedian(stacked_fz, axis=0)

        # -------------------------
        # Convert to Pearson r only at reporting stage
        # -------------------------

        mean_r   = np.tanh(mean_fz)
        median_r = np.tanh(median_fz)

        # -------------------------
        # Save
        # -------------------------

        # Fisher-z arrays
        np.save(
            os.path.join(
                partial_output_folder,
                f"group_mean_cluster_by_mmp-parcel_partial_fisherz{run_tag}_{tag}.npy"
            ),
            mean_fz,
        )
        np.save(
            os.path.join(
                partial_output_folder,
                f"group_median_cluster_by_mmp-parcel_partial_fisherz{run_tag}_{tag}.npy"
            ),
            median_fz,
        )

        pd.DataFrame(mean_fz,   index=clusters, columns=parcels).to_csv(
            os.path.join(
                partial_output_folder,
                f"group_mean_cluster_by_mmp-parcel_partial_fisherz{run_tag}_{tag}.csv"
            )
        )
        pd.DataFrame(median_fz, index=clusters, columns=parcels).to_csv(
            os.path.join(
                partial_output_folder,
                f"group_median_cluster_by_mmp-parcel_partial_fisherz{run_tag}_{tag}.csv"
            )
        )

        # Pearson r arrays
        np.save(
            os.path.join(
                partial_output_folder,
                f"group_mean_cluster_by_mmp-parcel_partial_r{run_tag}_{tag}.npy"
            ),
            mean_r,
        )
        np.save(
            os.path.join(
                partial_output_folder,
                f"group_median_cluster_by_mmp-parcel_partial_r{run_tag}_{tag}.npy"
            ),
            median_r,
        )

        pd.DataFrame(mean_r,   index=clusters, columns=parcels).to_csv(
            os.path.join(
                partial_output_folder,
                f"group_mean_cluster_by_mmp-parcel_partial_r{run_tag}_{tag}.csv"
            )
        )
        pd.DataFrame(median_r, index=clusters, columns=parcels).to_csv(
            os.path.join(
                partial_output_folder,
                f"group_median_cluster_by_mmp-parcel_partial_r{run_tag}_{tag}.csv"
            )
        )

        # Compressed archive with all four arrays + metadata
        np.savez_compressed(
            os.path.join(
                partial_output_folder,
                f"group_partial_corr{run_tag}_{tag}.npz"
            ),
            mean_fz            = mean_fz,
            median_fz          = median_fz,
            mean_r             = mean_r,
            median_r           = median_r,
            n_subjects_loaded  = np.array(n_loaded),
            subjects_loaded    = np.array([s for s in subjects if s not in missing]),
            subjects_missing   = np.array(missing),
            clusters           = np.array(clusters),
            parcels            = np.array(parcels),
            hemi               = np.array(label),
            run                = np.array(run if run else "concatenated"),
        )

        print(f"  [{label}{run_tag}] Saved. Shape: {mean_fz.shape}\n")

print("Done.")
# ============================================================