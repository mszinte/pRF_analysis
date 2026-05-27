#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
group_full_corr_by_hemi.py
------------------------------------------------------------------------------------------
Goal:
    Compute group-level Fisher-z statistics from per-subject, per-ROI, per-hemi
    full-correlation TSV files produced by the workbench parcellation step

    Each input TSV has:
        - 106 rows : parcels in atlas-key order
                     rows  0–52  → right hemisphere  (RH atlas keys 1–163)
                     rows 53–105 → left  hemisphere  (LH atlas keys 181–343)
        - 1 column : Fisher-z correlation of that seed/ROI with each parcel

    For each hemisphere × variant the script:
        1. Resolves which TSV to load for each subject (run variant logic).
        2. Loads all seed TSVs for a subject → slices the hemi-specific rows.
        3. Re-maps the parcel order from atlas-key order to the canonical order
           defined by the YAML config (rois-group-mmp), so columns are
           consistent with the partial-correlation pipeline and the WTA function.
        4. Stacks the (n_seeds × n_parcels) matrices across subjects.
        5. Computes the group median and standard deviation in Fisher-z space.
        6. Saves one .npy + one .csv per statistic per hemisphere × variant.

    Averaging is always in Fisher-z space; Pearson r is recovered only at the
    final reporting stage via tanh().  This is intentional: Fisher-z has
    approximately constant variance ~1/(n-3), making it the correct space for
    averaging and any subsequent parametric tests.  Raw Pearson r values have
    r-dependent variance and must never be averaged directly.
------------------------------------------------------------------------------------------
Run variants:
    concat       — concatenated-run TSV, all subjects
    concat_clean — best available run per subject:
                     · RUN02_EXCLUDED subjects → run-01 TSV
                     · all other subjects      → concatenated-run TSV
    run-01       — run-01 TSV, all subjects
    run-02       — run-02 TSV, all subjects (bad subjects kept intentionally
                   to expose the registration artifact in group plots)
------------------------------------------------------------------------------------------
Inputs (sys.argv):
    1: main project directory   (e.g. /scratch/mszinte/data)
    2: project name/directory   (e.g. RetinoMaps)
    3: server group             (e.g. 327)
    4: server project           (e.g. b327)
    5: parcellation mode:
           "default"     → _parcellated.tsv
           "legacy"      → _parcellated_legacy-mode.tsv
           "no_outliers" → _parcellated_no_outliers.tsv

Outputs (per hemisphere × variant):
    cluster_by_mmp-parcel_full-corr_fisherz_median_{run_label}_{hemi}_{mode}.npy / .csv
    cluster_by_mmp-parcel_full-corr_fisherz_std_{run_label}_{hemi}_{mode}.npy  / .csv

    Rows    : seed/cluster names  (n_clusters)
    Columns : parcel names in YAML-derived canonical order  (n_parcels)

To run:
    $ cd projects/pRF_analysis/RetinoMaps/rest/stats
    $ python group_stats_full_corr_by_hemi_task-free.py /scratch/mszinte/data RetinoMaps 327 b327 default
------------------------------------------------------------------------------------------
Written by Marco Bedini (marco.bedini@univ-amu.fr)
------------------------------------------------------------------------------------------
"""

import os
import sys
import warnings
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, List
warnings.filterwarnings("ignore")

# ============================================================
# Personal imports
# ============================================================
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../"))

sys.path.append(os.path.abspath(os.path.join(base_dir, "analysis_code/utils")))
from settings_utils import load_settings

sys.path.append(os.path.abspath(os.path.join(base_dir, "RetinoMaps/rest/utils")))
from rest_utils import (
    RUN02_EXCLUDED,
    VARIANTS,
    MODE_SUFFIX,
    ATLAS_ORDER,
    HEMI_ROW_SLICE,
    N_PARCELS_PER_HEMI,
    N_PARCELS_TOTAL,
    build_remap,
    load_full_corr_matrix,
)

# ============================================================
# Parse and validate arguments
# ============================================================
USAGE = (
    "Usage: python group_full_corr_by_hemi.py "
    "<main_dir> <project_dir> <group> <server> <mode>\n"
    f"  <mode> must be one of: {', '.join(MODE_SUFFIX)}"
)

if len(sys.argv) != 6:
    print(f"ERROR: expected 5 arguments, got {len(sys.argv) - 1}.\n{USAGE}")
    sys.exit(1)

main_dir    = sys.argv[1]
project_dir = sys.argv[2]
group       = sys.argv[3]
server      = sys.argv[4]
mode        = sys.argv[5]

if mode not in MODE_SUFFIX:
    print(f"ERROR: unrecognised mode '{mode}'.\n  Accepted: {', '.join(MODE_SUFFIX)}\n{USAGE}")
    sys.exit(1)

tsv_suffix = MODE_SUFFIX[mode]

print("=" * 80)
print("GROUP FULL CORRELATION — Fisher-z statistics (workbench parcellated TSVs)")
print("=" * 80)
print(f"  main_dir    : {main_dir}")
print(f"  project_dir : {project_dir}")
print(f"  group       : {group}")
print(f"  server      : {server}")
print(f"  mode        : {mode!r}  →  suffix: '_parcellated{tsv_suffix}.tsv'")

# ============================================================
# Load settings
# ============================================================
settings_path     = os.path.join(base_dir, project_dir, "settings.yml")
prf_settings_path = os.path.join(base_dir, project_dir, "prf-analysis.yml")
settings          = load_settings([settings_path, prf_settings_path])
analysis_info     = settings[0]
subjects          = analysis_info["subjects"]

# ============================================================
# ROIs — canonical order from YAML config
# ============================================================
clusters        = list(analysis_info["rois-drawn"])
seed_to_parcels = analysis_info["rois-group-mmp"]   # {seed_name: [parcel_names]}
clusters.reverse()                                   # mPCS first

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
output_folder = main_data / "group/91k/rest/full_corr/by_hemi"
output_folder.mkdir(parents=True, exist_ok=True)

# ============================================================
# Pre-compute remap indices: atlas-key order → canonical YAML order
# Done once per hemisphere; passed into load_full_corr_matrix() each call.
# ============================================================
_REMAP: Dict[str, tuple] = {}
for _hemi in ("rh", "lh"):
    _remap_idx, _present = build_remap(ATLAS_ORDER[_hemi], parcels, _hemi)
    _REMAP[_hemi] = (_remap_idx, _present)
    print(
        f"  Remap [{_hemi.upper()}]: {len(_present)}/{n_parcels} canonical parcels "
        "found in atlas key table."
    )

# ============================================================
# Main loop — hemisphere × variant
# ============================================================

for hemi in ("lh", "rh"):
    print(f"\n{'='*80}")
    print(f"Processing hemisphere: {hemi.upper()}")
    print("=" * 80)

    remap_idx, present = _REMAP[hemi]

    for variant, (normal_tag, excluded_tag, _skip) in VARIANTS.items():
        # run-02 intentionally keeps all subjects here (group QC strategy),
        # regardless of the skip_excluded flag used in the WTA script.
        print(f"\n  --- Variant: {variant} ---")

        subject_matrices: List[np.ndarray] = []
        subject_ids:      List[str]        = []

        for subject in subjects:
            is_excluded = subject in RUN02_EXCLUDED
            run_tag     = excluded_tag if is_excluded else normal_tag

            df = load_full_corr_matrix(
                subject, hemi, run_tag,
                clusters, parcels,
                remap_idx, present,
                main_data, tsv_suffix,
            )

            if df is None:
                print(f"    {subject}: SKIPPED (missing files)")
                continue

            # load_full_corr_matrix returns a DataFrame; .values gives the
            # (n_clusters × n_parcels) ndarray needed for np.stack below
            mat = df.values

            assert mat.shape == (n_clusters, n_parcels), (
                f"[{subject} {hemi} {variant}] Unexpected matrix shape {mat.shape}, "
                f"expected ({n_clusters}, {n_parcels})."
            )

            if variant == "concat_clean" and is_excluded:
                print(f"    {subject}: OK (fallback → run-01)")
            else:
                print(f"    {subject}: OK")

            subject_matrices.append(mat)
            subject_ids.append(subject)

        if not subject_matrices:
            print(f"    ERROR: no valid subjects for {hemi} / {variant} — skipping.")
            continue

        n_valid = len(subject_matrices)
        print(f"\n    Valid subjects: {n_valid}/{len(subjects)}")

        # Stack → (n_subjects × n_clusters × n_parcels)
        stack = np.stack(subject_matrices, axis=0)
        assert stack.shape == (n_valid, n_clusters, n_parcels), (
            f"Unexpected stack shape {stack.shape}."
        )

        # Group statistics in Fisher-z space.
        # nanmedian / nanstd correctly ignore NaN cells (absent parcels or |r|=1).
        group_median = np.nanmedian(stack, axis=0)   # (n_clusters × n_parcels)
        group_std    = np.nanstd(   stack, axis=0)   # (n_clusters × n_parcels)

        n_nan_med = int(np.isnan(group_median).sum())
        n_nan_std = int(np.isnan(group_std).sum())
        if n_nan_med:
            print(f"    WARNING: {n_nan_med}/{group_median.size} NaN cells in group median.")
        if n_nan_std:
            print(f"    WARNING: {n_nan_std}/{group_std.size} NaN cells in group std.")

        print(f"    Group median range : [{np.nanmin(group_median):.4f}, {np.nanmax(group_median):.4f}]")
        print(f"    Group std  range   : [{np.nanmin(group_std):.4f},  {np.nanmax(group_std):.4f}]")

        # Filename run-label: None normal_tag → use variant name (e.g. "concat")
        run_label = normal_tag if normal_tag is not None else variant

        # Save group median
        stem_median = (
            f"seed-task_by_mmp-parcel_full-corr_fisherz_median"
            f"_{run_label}_{hemi}_{mode}"
        )
        np.save(output_folder / f"{stem_median}.npy", group_median)
        pd.DataFrame(group_median, index=clusters, columns=parcels).to_csv(
            output_folder / f"{stem_median}.csv"
        )

        # Save group std
        stem_std = (
            f"seed-task_by_mmp-parcel_full-corr_fisherz_std"
            f"_{run_label}_{hemi}_{mode}"
        )
        np.save(output_folder / f"{stem_std}.npy", group_std)
        pd.DataFrame(group_std, index=clusters, columns=parcels).to_csv(
            output_folder / f"{stem_std}.csv"
        )

        print(f"    Saved: {stem_median}.npy / .csv")
        print(f"    Saved: {stem_std}.npy / .csv")

print("\n" + "=" * 80)
print("ALL HEMISPHERES × VARIANTS COMPLETE")
print("=" * 80)
print(f"\nOutputs written to: {output_folder}")
print(
    "\nNote: outputs are in Fisher-z space. Apply np.tanh() to recover Pearson r "
    "only at the final reporting or plotting stage."
)