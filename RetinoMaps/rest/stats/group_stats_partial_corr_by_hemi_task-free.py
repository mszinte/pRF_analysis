#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
group_partial_corr_by_hemi_task-free.py
------------------------------------------------------------------------------------------
Goal:
    Compute group-level Fisher-z statistics from per-subject partial-correlation
    matrices produced by Nilearn.

    Operates on per-subject .npy files (already in seeds × parcels format,
    YAML canonical parcel order). No TSV loading or atlas-key remapping needed.

    For each hemisphere × variant the script:
        1. Resolves which .npy file to load for each subject (run variant logic,
           including concat_clean fallback for RUN02_EXCLUDED subjects).
        2. Stacks the (n_clusters × n_parcels) matrices across subjects.
        3. Computes group mean and median in Fisher-z space, back-converts to r.
        4. Saves .npy, .csv and a compressed .npz archive per hemisphere × variant.

    Averaging is always in Fisher-z space; Pearson r is recovered only at the
    final reporting stage via tanh().  This is intentional:
      - Fisher-z has approximately constant variance ~1/(n-3), making it the
        correct space for averaging and any subsequent parametric tests.
      - Raw Pearson r values have r-dependent variance and must never be
        averaged directly.

    Covariance estimator (added for consistency with the task-constrained
    group stats script): ESTIMATOR_TAG (raw / ledoit-wolf / graphical-lasso)
    must match whatever was passed to nilearn_partial_corr_task_free.py at
    the subject level — subject-level output filenames are tagged with it,
    so this script needs the same tag to find the right files.
------------------------------------------------------------------------------------------
Output filename convention (harmonized with full-corr group stats):

    seed-task_by_mmp-parcel_partial-corr_{space}_{stat}_{run_label}_{hemi}_{estimator}.npy / .csv
    seed-task_by_mmp-parcel_partial-corr_{run_label}_{hemi}_{estimator}.npz

    e.g. seed-task_by_mmp-parcel_partial-corr_fisherz_median_concat_clean_rh_ledoit-wolf.npy

    Full-corr equivalent for reference:
    seed-task_by_macro_full-corr_fisherz_median_{run_label}_{hemi}.npy

    (partial corr has no mode token — there is only one parcellation variant)
------------------------------------------------------------------------------------------
Run variants:
    concat       — concatenated-run .npy, all subjects
    concat_clean — best available run per subject:
                     · RUN02_EXCLUDED subjects → run-01 .npy
                     · all other subjects      → concatenated-run .npy
    run-01       — run-01 .npy, all subjects
    run-02       — run-02 .npy, all subjects (bad subjects kept intentionally
                   to expose the registration artifact in group plots)
------------------------------------------------------------------------------------------
Inputs (sys.argv):
    1: main project directory   (e.g. /scratch/mszinte/data)
    2: project name/directory   (e.g. RetinoMaps)
    3: server group             (e.g. 327)
    4: server project           (e.g. b327)
    5: covariance estimator — one of "raw", "ledoit-wolf", "graphical-lasso"
       (optional; default: ledoit-wolf). Must match the ESTIMATOR_TAG used
       when running nilearn_partial_corr_task_free.py at the subject level.

Outputs (per hemisphere × variant):
    seed-task_by_mmp-parcel_partial-corr_fisherz_mean_{run_label}_{hemi}_{estimator}.npy / .csv
    seed-task_by_mmp-parcel_partial-corr_fisherz_median_{run_label}_{hemi}_{estimator}.npy / .csv
    seed-task_by_mmp-parcel_partial-corr_r_mean_{run_label}_{hemi}_{estimator}.npy / .csv
    seed-task_by_mmp-parcel_partial-corr_r_median_{run_label}_{hemi}_{estimator}.npy / .csv
    seed-task_by_mmp-parcel_partial-corr_{run_label}_{hemi}_{estimator}.npz

    Rows    : seed/cluster names  (n_clusters)
    Columns : parcel names in YAML-derived canonical order  (n_parcels)

To run:
    $ cd projects/pRF_analysis/RetinoMaps/rest/stats
    $ python group_partial_corr_by_hemi_task-free.py /scratch/mszinte/data RetinoMaps 327 b327 ledoit-wolf
------------------------------------------------------------------------------------------
Written by Marco Bedini (marco.bedini@univ-amu.fr)
------------------------------------------------------------------------------------------
"""

import os
import sys
import numpy as np
import pandas as pd
from pathlib import Path
from typing import List, Optional

# ============================================================
# Parse and validate arguments
# ============================================================
VALID_ESTIMATORS = ("raw", "ledoit-wolf", "graphical-lasso")

USAGE = (
    "Usage: python group_partial_corr_by_hemi_task-free.py "
    "<main_dir> <project_dir> <group> <server> [estimator]\n"
    f"  estimator: one of {', '.join(VALID_ESTIMATORS)} (default: ledoit-wolf)\n"
    "  Must match the ESTIMATOR_TAG used when running the subject-level\n"
    "  nilearn_partial_corr_task_free.py script."
)

if len(sys.argv) not in (5, 6):
    print(f"ERROR: expected 4 or 5 arguments, got {len(sys.argv) - 1}.\n{USAGE}")
    sys.exit(1)

main_dir    = sys.argv[1]
project_dir = sys.argv[2]
group       = sys.argv[3]
server      = sys.argv[4]

ESTIMATOR_TAG = sys.argv[5] if len(sys.argv) == 6 else "ledoit-wolf"
if ESTIMATOR_TAG not in VALID_ESTIMATORS:
    print(f"ERROR: unrecognised estimator '{ESTIMATOR_TAG}'.\n{USAGE}")
    sys.exit(1)

print("=" * 80)
print("GROUP PARTIAL CORRELATION (task-free) — Fisher-z statistics")
print(f"Estimator: {ESTIMATOR_TAG}")
print("=" * 80)
print(f"  main_dir    : {main_dir}")
print(f"  project_dir : {project_dir}")
print(f"  group       : {group}")
print(f"  server      : {server}")

# ============================================================
# Load settings
# ============================================================
base_dir            = os.path.abspath(os.path.join(os.getcwd(), "../../../"))
sys.path.append(os.path.abspath(os.path.join(base_dir, "analysis_code/utils")))
from settings_utils import load_settings
sys.path.append(os.path.abspath(os.path.join(base_dir, "RetinoMaps/rest/utils")))
from rest_utils import VARIANTS

settings_path       = os.path.join(base_dir, project_dir, "settings.yml")
prf_settings_path   = os.path.join(base_dir, project_dir, "prf-analysis.yml")
settings            = load_settings([settings_path, prf_settings_path])
rest_settings_path  = os.path.join(base_dir, project_dir, "rest-settings.yml")
rest_settings       = load_settings([rest_settings_path])[0]
analysis_info       = settings[0]
subjects            = analysis_info["subjects"]
RUNS                = rest_settings["runs"][0]
RUN02_EXCLUDED      = frozenset(rest_settings["run02_excluded"][0])

# ============================================================
# ROIs — canonical order from YAML config
# ============================================================
clusters        = list(analysis_info["rois-drawn"])
seed_to_parcels = analysis_info["rois-group-mmp"]
clusters.reverse()   # mPCS first

parcels: List[str] = []
for cl in clusters:
    parcels.extend(seed_to_parcels[cl])

n_clusters = len(clusters)
n_parcels  = len(parcels)

print(f"\n  Seeds   (n={n_clusters}): {clusters}")
print(f"  Parcels (n={n_parcels}): {parcels[:5]} ... {parcels[-5:]}")

# ============================================================
# Paths
# ============================================================
main_data     = Path(main_dir) / project_dir / "derivatives/pp_data"
output_folder = main_data / "group/91k/rest/partial_corr/by_hemi/task-free"
output_folder.mkdir(parents=True, exist_ok=True)

# ============================================================
# Helper: resolve per-subject .npy path
#
# Subject-level files live under:
#   {subject}/91k/rest/corr/partial_corr/by_hemi/task-free/
#       seed-task_by_mmp-parcel_partial_fisherz{run_entity}_{hemi}_{estimator}.npy
#
# Note: subject-level files use "partial" (no hyphen), not "partial-corr".
#       The hyphenated form is only used in group output filenames.
#
# CORRECTED run_entity convention: the subject-level script
# (nilearn_partial_corr_task_free.py) uses run_tag = f"_{run}" if run else ""
# — i.e. NO suffix at all for the concat variant, not "_concat". This
# matches the convention used everywhere else in the pipeline (e.g.
# rest_utils.tsv_path). An earlier version of this function used "_concat"
# as the concat fallback, which does not match any file the subject-level
# script actually produces — fixed here.
# ============================================================
def npy_path(subject: str, hemi: str, run_tag: Optional[str]) -> Path:
    run_entity = f"_{run_tag}" if run_tag is not None else ""
    fname      = (
        f"seed-task_by_mmp-parcel_partial_fisherz{run_entity}"
        f"_{hemi}_{ESTIMATOR_TAG}.npy"
    )
    return (
        main_data / subject
        / "91k/rest/corr/partial_corr/by_hemi/task-free"
        / fname
    )

# ============================================================
# Output filename stem builder
#
# Pattern: seed-task_by_mmp-parcel_partial-corr_{space}_{stat}_{run_label}_{hemi}_{estimator}
# Matches full-corr convention: seed-task_by_macro_full-corr_fisherz_{stat}_{run_label}_{hemi}
# ============================================================
def _stem(stat: str, space: str, run_label: str, hemi: str) -> str:
    return (
        f"seed-task_by_mmp-parcel_partial-corr"
        f"_{space}_{stat}_{run_label}_{hemi}_{ESTIMATOR_TAG}"
    )

# ============================================================
# Main loop — hemisphere × variant
# ============================================================

for hemi in ("lh", "rh"):
    print(f"\n{'='*80}")
    print(f"Processing hemisphere: {hemi.upper()}")
    print("=" * 80)

    for variant, (normal_tag, excluded_tag, _skip) in VARIANTS.items():
        # run-02 intentionally keeps all subjects (group QC strategy),
        # regardless of the skip_excluded flag used in the WTA script.
        print(f"\n  --- Variant: {variant} ---")

        fz_stack:    List[np.ndarray] = []
        subject_ids: List[str]        = []
        missing:     List[str]        = []

        for subject in subjects:
            is_excluded = subject in RUN02_EXCLUDED
            run_tag     = excluded_tag if is_excluded else normal_tag

            fpath = npy_path(subject, hemi, run_tag)

            if not fpath.exists():
                print(f"    WARNING [{subject}]: missing {fpath.name} — skipped")
                missing.append(subject)
                continue

            mat = np.load(fpath)

            if mat.shape != (n_clusters, n_parcels):
                raise ValueError(
                    f"[{subject} {hemi} {variant}] Unexpected shape {mat.shape} "
                    f"in {fpath.name} (expected ({n_clusters}, {n_parcels}))."
                )

            if variant == "concat_clean" and is_excluded:
                print(f"    {subject}: OK (fallback → run-01)")
            else:
                print(f"    {subject}: OK")

            fz_stack.append(mat)
            subject_ids.append(subject)

        if not fz_stack:
            print(f"    ERROR: no valid subjects for {hemi} / {variant} — skipping.")
            continue

        n_valid = len(fz_stack)
        print(f"\n    Valid subjects: {n_valid}/{len(subjects)}")
        if missing:
            print(f"    Missing       : {missing}")

        # Stack → (n_subjects × n_clusters × n_parcels)
        stacked_fz = np.stack(fz_stack, axis=0)
        assert stacked_fz.shape == (n_valid, n_clusters, n_parcels), (
            f"Unexpected stack shape {stacked_fz.shape}."
        )

        # Group statistics — always in Fisher-z space
        mean_fz   = np.nanmean(  stacked_fz, axis=0)   # (n_clusters × n_parcels)
        median_fz = np.nanmedian(stacked_fz, axis=0)

        # Back-convert to Pearson r only at reporting stage
        mean_r   = np.tanh(mean_fz)
        median_r = np.tanh(median_fz)

        print(f"    Fisher-z mean   range : [{np.nanmin(mean_fz):.4f},   {np.nanmax(mean_fz):.4f}]")
        print(f"    Fisher-z median range : [{np.nanmin(median_fz):.4f}, {np.nanmax(median_fz):.4f}]")

        # run_label: None normal_tag → variant name (e.g. "concat", "concat_clean")
        run_label = normal_tag if normal_tag is not None else variant

        # Save Fisher-z arrays
        for stat, arr in (("mean", mean_fz), ("median", median_fz)):
            stem = _stem(stat, "fisherz", run_label, hemi)
            np.save(output_folder / f"{stem}.npy", arr)
            pd.DataFrame(arr, index=clusters, columns=parcels).to_csv(
                output_folder / f"{stem}.csv"
            )
            print(f"    Saved: {stem}.npy / .csv")

        # Save Pearson r arrays
        for stat, arr in (("mean", mean_r), ("median", median_r)):
            stem = _stem(stat, "r", run_label, hemi)
            np.save(output_folder / f"{stem}.npy", arr)
            pd.DataFrame(arr, index=clusters, columns=parcels).to_csv(
                output_folder / f"{stem}.csv"
            )
            print(f"    Saved: {stem}.npy / .csv")

        # Compressed archive with all four arrays + metadata
        npz_stem = (
            f"seed-task_by_mmp-parcel_partial-corr"
            f"_{run_label}_{hemi}_{ESTIMATOR_TAG}"
        )
        np.savez_compressed(
            output_folder / f"{npz_stem}.npz",
            mean_fz           = mean_fz,
            median_fz         = median_fz,
            mean_r            = mean_r,
            median_r          = median_r,
            n_subjects_loaded = np.array(n_valid),
            subjects_loaded   = np.array(subject_ids),
            subjects_missing  = np.array(missing),
            clusters          = np.array(clusters),
            parcels           = np.array(parcels),
            hemi              = np.array(hemi),
            variant           = np.array(variant),
            cov_estimator     = np.array(ESTIMATOR_TAG),
        )
        print(f"    Saved: {npz_stem}.npz")

print("\n" + "=" * 80)
print("ALL HEMISPHERES × VARIANTS COMPLETE")
print("=" * 80)
print(f"\nOutputs written to: {output_folder}")
print(
    "\nNote: Pearson r group values are recovered back from Fisher-z outputs with np.tanh()"
    "Use these values at the final reporting or plotting stage as is standard practice."
)