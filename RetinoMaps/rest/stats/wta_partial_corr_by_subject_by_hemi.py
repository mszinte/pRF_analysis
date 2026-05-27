#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
wta_partial_corr_by_subject_by_hemi.py
-----------------------------------------------------------------------------------------
Goal:
    Compute winner-take-all (WTA) from per-subject partial correlation CSVs.
    Outputs one combined table per hemisphere × variant containing:
        - one row per subject   (winner seed label per parcel)
        - one GROUP row         (modal winner across subjects per parcel)
        - one CONSISTENCY row   (% of subjects matching the group winner)

    Input CSVs are already in seeds × parcels format (no slicing needed).
    Working on fisher-z outputs that stabilize variance across subjects/parcels.
    All shared logic (constants, WTA, tie-breaking) is imported from rest_utils.py.

    Tie-breaking cascade in append_group_and_consistency():
        1. Vote count   — unique modal winner across subjects → done
        2. Fisher-z     — concat_clean group-median partial-corr Fisher-z used
                          to resolve remaining ties (same reference for all variants)
        3. Seed number  — deterministic fallback: lowest seed number

    Group median file (from group_partial_corr_by_hemi.py):
        seed-task_by_mmp-parcel_partial-corr_fisherz_median_concat_clean_{hemi}.npy
-----------------------------------------------------------------------------------------
Run variants (from rest_utils.VARIANTS):
    concat       — concatenated-run CSV, all subjects
    concat_clean — best available run per subject (concat or run-01 fallback)
    run-01       — run-01 CSV, all subjects
    run-02       — run-02 CSV, RUN02_EXCLUDED subjects skipped with WARNING

Filename convention (input CSVs):
    concat / concat_clean : cluster_by_mmp-parcel_partial_fisherz_{hemi}.csv
    run-01 / run-02       : cluster_by_mmp-parcel_partial_fisherz_{run_tag}_{hemi}.csv
-----------------------------------------------------------------------------------------
Inputs (sys.argv):
    1: main project directory   (e.g. /scratch/mszinte/data)
    2: project name/directory   (e.g. RetinoMaps)
    3: server group             (e.g. 327)
    4: server project           (e.g. b327)

Output:
    Per hemisphere × variant, one CSV:
        winning_seeds_by_subject_partial_corr_{hemi}_{variant}.csv
    Rows: subjects + GROUP + CONSISTENCY_%; columns: parcel names.

To run:
    $ cd projects/pRF_analysis/RetinoMaps/rest/stats
    $ python wta_partial_corr_by_subject_by_hemi.py /scratch/mszinte/data RetinoMaps 327 b327
-----------------------------------------------------------------------------------------
Written by Marco Bedini (marco.bedini@univ-amu.fr)
-----------------------------------------------------------------------------------------
"""

import warnings
warnings.filterwarnings("ignore")

import os
import sys
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, List, Optional

# ============================================================
# Personal imports — settings utils
# ============================================================
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../"))
sys.path.append(os.path.abspath(os.path.join(base_dir, "analysis_code/utils")))
from settings_utils import load_settings

# ============================================================
# rest_utils — shared pipeline constants and functions
# ============================================================
sys.path.append(os.path.abspath(os.path.join(base_dir, "RetinoMaps/rest/utils")))
from rest_utils import (
    RUN02_EXCLUDED,
    VARIANTS,
    compute_winners,
    append_group_and_consistency,
)

# ============================================================
# Parse and validate arguments
# ============================================================
USAGE = (
    "Usage: python wta_partial_corr_by_subject_by_hemi.py "
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
print("WTA — partial correlation (Nilearn parcellated fisher-z CSVs)")
print("=" * 80)
print(f"  main_dir    : {main_dir}")
print(f"  project_dir : {project_dir}")
print(f"  group       : {group}")
print(f"  server      : {server}")

# ============================================================
# Load settings
# ============================================================
settings_path     = os.path.join(base_dir, project_dir, "settings.yml")
prf_settings_path = os.path.join(base_dir, project_dir, "prf-analysis.yml")
settings          = load_settings([settings_path, prf_settings_path])
analysis_info     = settings[0]
subjects          = analysis_info["subjects"]

# ============================================================
# Paths
# ============================================================
main_data          = Path(main_dir) / project_dir / "derivatives/pp_data"
output_folder      = main_data / "group/91k/rest/wta/nilearn"
group_corr_folder  = main_data / "group/91k/rest/partial_corr/by_hemi"
output_folder.mkdir(parents=True, exist_ok=True)

# ============================================================
# ROI / cluster settings
# ============================================================
clusters: List[str]                   = list(analysis_info["rois-drawn"])
seed_to_parcels: Dict[str, List[str]] = analysis_info["rois-group-mmp"]
clusters.reverse()                    # mPCS first

seed_to_number: Dict[str, int] = {s: i + 1 for i, s in enumerate(clusters)}

n_clusters = len(clusters)

# ============================================================
# Load concat_clean group-median Fisher-z matrices
#
# Used as the level-2 Fisher-z reference in break_ties_wta() for ALL variants.
# One matrix per hemisphere: shape (n_clusters × n_parcels), rows in clusters
# order, columns in canonical parcels order.
#
# Filename produced by group_partial_corr_by_hemi.py:
#   seed-task_by_mmp-parcel_partial-corr_fisherz_median_concat_clean_{hemi}.npy
# ============================================================
def load_group_median(hemi: str, n_parcels: int) -> Optional[np.ndarray]:
    """
    Load the concat_clean group-median partial-corr Fisher-z .npy for one hemi.
    Returns None with a WARNING if the file is not found.
    """
    fname = (
        f"seed-task_by_mmp-parcel_partial-corr_fisherz"
        f"_median_concat_clean_{hemi}.npy"
    )
    fpath = group_corr_folder / fname
    if not fpath.exists():
        print(f"  WARNING: group median not found — {fpath.name}")
        print("    Tie-breaking will fall back to level 3 (lowest seed number).")
        return None
    mat = np.load(fpath)
    if mat.shape != (n_clusters, n_parcels):
        raise ValueError(
            f"Group median shape {mat.shape} does not match expected "
            f"({n_clusters}, {n_parcels}) for {fpath.name}.\n"
            "  Check that the group stats script used the same YAML config."
        )
    print(f"  Group median [{hemi.upper()}]: loaded {fpath.name}  shape={mat.shape}")
    return mat

# ============================================================
# Partial-corr-specific I/O helpers
# (not in rest_utils: path logic and file format are Nilearn-specific)
# ============================================================

def csv_path(subject: str, hemi: str, run_tag: Optional[str]) -> Path:
    """
    Return expected path for one subject / hemi partial corr CSV.

    run_tag = None  → concatenated file:
        cluster_by_mmp-parcel_partial_fisherz_{hemi}.csv
    run_tag = str   → per-run file:
        cluster_by_mmp-parcel_partial_fisherz_{run_tag}_{hemi}.csv
    """
    subj_dir   = main_data / subject / "91k/rest/corr/partial_corr/by_hemi"
    run_entity = f"_{run_tag}" if run_tag is not None else ""
    fname      = f"cluster_by_mmp-parcel_partial_fisherz{run_entity}_{hemi}.csv"
    return subj_dir / fname


def load_partial_corr_matrix(
    subject: str,
    hemi: str,
    run_tag: Optional[str],
) -> Optional[pd.DataFrame]:
    """
    Load the seeds × parcels partial correlation matrix for one subject / hemi.

    Returns None if the file is missing.
    Raises ValueError if any expected seed is absent from the CSV index.
    """
    fpath = csv_path(subject, hemi, run_tag)

    if not fpath.exists():
        print(f"  WARNING [{subject} {hemi}]: missing {fpath.name}")
        return None

    df = pd.read_csv(fpath, index_col=0)

    missing_seeds = set(clusters) - set(df.index)
    if missing_seeds:
        raise ValueError(
            f"Seeds missing from {fpath.name}: {sorted(missing_seeds)}"
        )

    return df

# ============================================================
# Main loop — hemisphere × variant
# ============================================================

for hemi in ("lh", "rh"):
    print(f"\n{'='*80}")
    print(f"Processing hemisphere: {hemi.upper()}")
    print("=" * 80)

    # parcel_columns captured from first valid subject — consistent across variants
    # within a hemisphere (same CSV structure guaranteed by the Nilearn pipeline)
    parcel_columns: Optional[List[str]] = None

    for variant, (normal_tag, excluded_tag, skip_excluded) in VARIANTS.items():
        print(f"\n  --- Variant: {variant} ---")

        all_winners: List[np.ndarray] = []
        subject_ids: List[str]        = []

        for subject in subjects:
            is_excluded = subject in RUN02_EXCLUDED

            if is_excluded and skip_excluded:
                print(f"    WARNING [{variant}]: {subject} is in RUN02_EXCLUDED — SKIPPED")
                continue

            run_tag = excluded_tag if is_excluded else normal_tag

            df_corr = load_partial_corr_matrix(subject, hemi, run_tag)
            if df_corr is None:
                print(f"    {subject}: SKIPPED (missing file)")
                continue

            # Capture and enforce consistent parcel column order across all
            # subjects and variants within this hemisphere
            if parcel_columns is None:
                parcel_columns = list(df_corr.columns)
            elif list(df_corr.columns) != parcel_columns:
                raise ValueError(
                    f"Column order mismatch for {subject} {hemi} {variant} — "
                    "parcels are not consistent across subjects."
                )

            if variant == "concat_clean" and is_excluded:
                print(f"    {subject}: OK (fallback → run-01)")
            else:
                print(f"    {subject}: OK")

            all_winners.append(
                compute_winners(
                    df_corr         = df_corr,
                    seed_to_parcels = seed_to_parcels,
                    seed_to_number  = seed_to_number,
                )
            )
            subject_ids.append(subject)

        if not all_winners:
            print(f"    ERROR: no valid subjects for {hemi} / {variant} — skipping.")
            continue

        # Load group median now that we know n_parcels from the data
        # (done once per hemi, lazily on first variant that produces results)
        n_parcels = len(parcel_columns)
        group_median = load_group_median(hemi, n_parcels)

        subject_df = pd.DataFrame(
            all_winners, index=subject_ids, columns=parcel_columns
        )

        # Tie-breaking cascade:
        #   Level 1 — unique vote winner
        #   Level 2 — concat_clean group-median partial-corr Fisher-z
        #   Level 3 — lowest seed number (deterministic fallback)
        combined_df = append_group_and_consistency(
            subject_df     = subject_df,
            clusters       = clusters,
            seed_to_number = seed_to_number,
            group_median   = group_median,
            parcels        = parcel_columns,
        )

        out_csv = (
            output_folder
            / f"winning_seeds_by_subject_partial_corr_{hemi}_{variant}.csv"
        )
        combined_df.to_csv(out_csv)
        print(f"    Saved: {out_csv.name}")
        print(f"    Rows : {list(combined_df.index)}")

print("\n" + "=" * 80)
print("ALL HEMISPHERES × VARIANTS COMPLETE")
print("=" * 80)
print(f"\nOutputs written to: {output_folder}")