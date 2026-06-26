#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
group_partial_corr_by_hemi_task-constrained.py
------------------------------------------------------------------------------------------
Goal:
    Compute group-level Fisher-z statistics from per-subject partial-correlation
    matrices produced by partial_corr_hemi_task_constrained.py

    Operates on per-subject .npy files in (n_clusters × n_clusters) format
    (macro-region × macro-region, ipsilateral only).

    For each hemisphere × variant the script:
        1. Resolves which .npy file to load per subject (run variant logic,
           including concat_clean fallback for RUN02_EXCLUDED subjects).
        2. Stacks the (n_clusters × n_clusters) matrices across subjects.
        3. Computes group mean and median in Fisher-z space, back-converts to r.
        4. Saves .npy, .csv and a compressed .npz archive per hemisphere × variant.

    Averaging is always in Fisher-z space; Pearson r is recovered only at the
    final reporting stage via tanh().

------------------------------------------------------------------------------------------
Run variants:
    concat       — concatenated-run .npy, all subjects
    concat_clean — best available run per subject:
                     · RUN02_EXCLUDED subjects → run-01 .npy
                     · all other subjects      → concatenated-run .npy
    run-01       — run-01 .npy, all subjects
    run-02       — run-02 .npy, all subjects

------------------------------------------------------------------------------------------
Inputs (sys.argv):
    1: main project directory   (e.g. /scratch/mszinte/data)
    2: project name/directory   (e.g. RetinoMaps)
    3: server group             (e.g. 327)
    4: server project           (e.g. b327)

Outputs (per hemisphere × variant):
    seed-task_by_macror-task_partial-corr_fisherz_mean_{run_label}_{hemi}.npy / .csv
    seed-task_by_macror-task_partial-corr_fisherz_median_{run_label}_{hemi}.npy / .csv
    seed-task_by_macror-task_partial-corr_r_mean_{run_label}_{hemi}.npy / .csv
    seed-task_by_macror-task_partial-corr_r_median_{run_label}_{hemi}.npy / .csv
    seed-task_by_macror-task_partial-corr_{run_label}_{hemi}.npz

    e.g. seed-task_by_macror-task_partial-corr_fisherz_median_concat_clean_rh.npy

    Full-corr equivalent for reference:
    seed-task_by_macro_full-corr_fisherz_median_{run_label}_{hemi}.npy

To run:
    $ cd projects/pRF_analysis/RetinoMaps/rest/stats
    $ python group_partial_corr_by_hemi_task-constrained.py /scratch/mszinte/data RetinoMaps 327 b327
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
# Parse arguments
# ============================================================
USAGE = (
    "Usage: python group_partial_corr_by_hemi_task-constrained.py "
    "<main_dir> <project_dir> <group> <server>"
)

if len(sys.argv) != 5:
    print(f"ERROR: expected 4 arguments, got {len(sys.argv) - 1}.\n{USAGE}")
    sys.exit(1)

main_dir    = sys.argv[1]
project_dir = sys.argv[2]
group       = sys.argv[3]
server      = sys.argv[4]

print("=" * 80)
print("GROUP PARTIAL CORRELATION (task-constrained) — Fisher-z statistics")
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
clusters = list(analysis_info["rois-drawn"])
clusters.reverse()   # mPCS first

n_clusters = len(clusters)

print(f"\n  Macro-regions (n={n_clusters}): {clusters}")

# ============================================================
# Paths
# ============================================================
main_data     = Path(main_dir) / project_dir / "derivatives/pp_data"
output_folder = main_data / "group/91k/rest/partial_corr/by_hemi/task-constrained"
output_folder.mkdir(parents=True, exist_ok=True)

# ============================================================
# Helper: resolve per-subject .npy path
#
# Subject-level files live under:
#   {subject}/91k/rest/corr/partial_corr/by_hemi/task-constrained/
#       seed-task_by_macror-task_partial_fisherz{run_tag}_{hemi}.npy
#
# run_tag mapping:
#   run-01     → "_run-01"
#   run-02     → "_run-02"
#   concat     → "_concat"      (full concatenated session)
#   concat_clean uses "_run-01" for RUN02_EXCLUDED, "_concat" for others
# ============================================================

def npy_path(subject: str, hemi: str, run_tag: str) -> Path:
    fname = f"seed-task_by_macror-task_partial_fisherz{run_tag}_{hemi}.npy"
    return (
        main_data / subject
        / "91k/rest/corr/partial_corr/by_hemi/task-constrained"
        / fname
    )

# Map VARIANTS normal_tag / excluded_tag (None = concat) to filename run_tag
def to_run_tag(tag: Optional[str]) -> str:
    """Convert VARIANTS run_tag (None = concat, str = run-XX) to filename tag."""
    return f"_{tag}" if tag is not None else "_concat"

# ============================================================
# Main loop — hemisphere × variant
# ============================================================

for hemi in ("lh", "rh"):
    print(f"\n{'='*80}")
    print(f"Processing hemisphere: {hemi.upper()}")
    print("=" * 80)

    for variant, (normal_tag, excluded_tag, _skip) in VARIANTS.items():

        print(f"\n  --- Variant: {variant} ---")

        fz_stack:    List[np.ndarray] = []
        subject_ids: List[str]        = []
        missing:     List[str]        = []

        for subject in subjects:
            is_excluded = subject in RUN02_EXCLUDED
            tag         = excluded_tag if is_excluded else normal_tag
            run_tag     = to_run_tag(tag)

            fpath = npy_path(subject, hemi, run_tag)

            if not fpath.exists():
                print(f"    WARNING [{subject}]: missing {fpath.name} — skipped")
                missing.append(subject)
                continue

            mat = np.load(fpath)

            # Expected shape: (n_clusters × n_clusters)
            if mat.shape != (n_clusters, n_clusters):
                raise ValueError(
                    f"[{subject} {hemi} {variant}] Unexpected shape {mat.shape} "
                    f"in {fpath.name} (expected ({n_clusters}, {n_clusters}))."
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

        # Stack → (n_subjects × n_clusters × n_clusters)
        stacked_fz = np.stack(fz_stack, axis=0)
        assert stacked_fz.shape == (n_valid, n_clusters, n_clusters), (
            f"Unexpected stack shape {stacked_fz.shape}."
        )

        # ── Group statistics in Fisher-z space ───────────────────────────────
        mean_fz   = np.nanmean(  stacked_fz, axis=0)   # (n_clusters × n_clusters)
        median_fz = np.nanmedian(stacked_fz, axis=0)

        # ── Back-convert to Pearson r at reporting stage only ─────────────────
        mean_r   = np.tanh(mean_fz)
        median_r = np.tanh(median_fz)

        print(f"    Fisher-z mean   range : [{np.nanmin(mean_fz):.4f},  {np.nanmax(mean_fz):.4f}]")
        print(f"    Fisher-z median range : [{np.nanmin(median_fz):.4f}, {np.nanmax(median_fz):.4f}]")

        # run_label: None normal_tag → use variant name (e.g. "concat", "concat_clean")
        run_label = normal_tag if normal_tag is not None else variant

        # Output filename pattern (harmonized with full-corr):
        #   seed-task_by_macror-task_partial-corr_{space}_{stat}_{run_label}_{hemi}
        # e.g. seed-task_by_macror-task_partial-corr_fisherz_median_concat_clean_rh.npy
        def _stem(stat: str, space: str) -> str:
            return (
                f"seed-task_by_macror-task_partial-corr"
                f"_{space}_{stat}_{run_label}_{hemi}"
            )

        # Save Fisher-z arrays
        for stat, arr in (("mean", mean_fz), ("median", median_fz)):
            stem = _stem(stat, "fisherz")
            np.save(output_folder / f"{stem}.npy", arr)
            pd.DataFrame(arr, index=clusters, columns=clusters).to_csv(
                output_folder / f"{stem}.csv"
            )
            print(f"    Saved: {stem}.npy / .csv")

        # Save Pearson r arrays
        for stat, arr in (("mean", mean_r), ("median", median_r)):
            stem = _stem(stat, "r")
            np.save(output_folder / f"{stem}.npy", arr)
            pd.DataFrame(arr, index=clusters, columns=clusters).to_csv(
                output_folder / f"{stem}.csv"
            )
            print(f"    Saved: {stem}.npy / .csv")

        # Compressed archive — keys use the same space_stat pattern for consistency
        npz_stem = f"seed-task_by_macror-task_partial-corr_{run_label}_{hemi}"
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
            hemi              = np.array(hemi),
            variant           = np.array(variant),
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