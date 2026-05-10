#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Feb 27, 2025

Compute winner take all on partial correlation results
From Nilearn (outputs are parcellated by definition)
Working on fisherz outputs that stabilize variance across subjects and parcels

wta_partial_corr_by_subject_by_hemi_per_run.py
-----------------------------------------------------------------------------------------
Goal:
    Compute winner-take-all (WTA) from per-subject partial correlation CSVs
    Outputs one combined table per hemisphere × variant containing:
        - one row per subject   (winner seed label per parcel)
        - one GROUP row         (modal winner across subjects per parcel)
        - one CONSISTENCY row   (% of subjects matching the group winner)

    Input CSVs are already in seeds × parcels format (no slicing needed)

Pipeline per hemisphere × variant:
    1. For each subject, resolve which CSV to load (varies by variant)
    2. Load the seeds × parcels correlation matrix
    3. Apply compute_winners() → winning seed label per parcel
    4. Collect across subjects → compute GROUP (mode) and CONSISTENCY rows
    5. Save combined table
-----------------------------------------------------------------------------------------
Run variants:
    concat       — concatenated-run CSV, all subjects
    concat_clean — best available run per subject:
                     · RUN02_EXCLUDED subjects → run-01 CSV
                     · all other subjects      → concatenated-run CSV
    run-01       — run-01 CSV, all subjects
    run-02       — run-02 CSV, RUN02_EXCLUDED subjects skipped with WARNING

Filename convention:
    concat / concat_clean : cluster_by_mmp-parcel_partial_{hemi}.csv
    run-01 / run-02       : cluster_by_mmp-parcel_partial_{run_tag}_{hemi}.csv
-----------------------------------------------------------------------------------------
Inputs (sys.argv):
    1: main project directory   (e.g. /scratch/mszinte/data)
    2: project name/directory   (e.g. RetinoMaps)
    3: server group             (e.g. 327)
    4: server project           (e.g. b327)

Output:
    Per hemisphere × variant, one CSV:
        winning_seeds_by_subject_partial_corr_{hemi}_{variant}.csv
    Rows: subjects + GROUP + CONSISTENCY; columns: parcel names

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
from scipy import stats as scipy_stats

# ============================================================
# Personal imports
# ============================================================
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../"))
sys.path.append(os.path.abspath(os.path.join(base_dir, "analysis_code/utils")))
from settings_utils import load_settings

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
print("WTA — partial correlation (Nilearn parcellated CSVs)")
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
main_data     = Path(main_dir) / project_dir / "derivatives/pp_data"
output_folder = main_data / "group/91k/rest/wta/nilearn"
output_folder.mkdir(parents=True, exist_ok=True)

# ============================================================
# Run variants
#
# Subjects excluded from run-02 (bad data / registration error (sub-22)).
# For concat_clean these subjects fall back to their run-01 file instead.
# ============================================================
RUN02_EXCLUDED: set = {"sub-03", "sub-04", "sub-14", "sub-21", "sub-22", "sub-23"}

# Variant definitions ─────────────────────────────────────────────────────────
# Each entry:  variant_tag → (run_tag_for_normal_subjects,
#                             run_tag_for_RUN02_EXCLUDED_subjects,
#                             skip_RUN02_EXCLUDED)
#
#   run_tag = None  → concatenated file (no run-XX in filename)
#   run_tag = str   → per-run file (run-XX inserted in filename)
#   skip_RUN02_EXCLUDED = True → those subjects are skipped with a WARNING
#
# variant        normal_tag  excluded_tag  skip_excluded
VARIANTS: Dict[str, tuple] = {
    "concat":       (None,     None,     False),  # concat file, all subjects
    "concat_clean": (None,     "run-01", False),  # concat for good, run-01 for bad
    "run-01":       ("run-01", "run-01", False),  # run-01 file, all subjects
    "run-02":       ("run-02", None,     True),   # run-02 file, skip bad subjects
}

# ============================================================
# ROI / cluster settings
# ============================================================
clusters        = list(analysis_info["rois-drawn"])
seed_to_parcels = analysis_info["rois-group-mmp"]   # {seed_name: [parcel_names]}
clusters.reverse()                                   # mPCS first, same as original

seed_to_number: Dict[str, int] = {s: i + 1 for i, s in enumerate(clusters)}

# ============================================================
# Helpers
# ============================================================

def csv_path(subject: str, hemi: str, run_tag: Optional[str]) -> Path:
    """
    Return expected path for one subject / hemi partial corr CSV.

    run_tag = None  → concatenated file:
        cluster_by_mmp-parcel_partial_{hemi}.csv
    run_tag = str   → per-run file:
        cluster_by_mmp-parcel_partial_{run_tag}_{hemi}.csv
    """
    subj_dir   = main_data / subject / "91k/rest/corr/partial_corr/by_hemi"
    run_entity = f"_{run_tag}" if run_tag is not None else ""
    fname      = f"cluster_by_mmp-parcel_partial_fisherz{run_entity}_{hemi}.csv"
    return subj_dir / fname


def load_corr_matrix(
    subject: str,
    hemi: str,
    run_tag: Optional[str],
) -> Optional[pd.DataFrame]:
    """
    Load the seeds × parcels partial correlation matrix for one subject / hemi.

    Returns None if the file is missing.
    """
    fpath = csv_path(subject, hemi, run_tag)

    if not fpath.exists():
        print(f"  WARNING [{subject} {hemi}]: missing {fpath.name}")
        return None

    df = pd.read_csv(fpath, index_col=0)

    # Defensive: verify expected seeds are present
    missing_seeds = set(clusters) - set(df.index)
    if missing_seeds:
        raise ValueError(
            f"Seeds missing from {fpath.name}: {sorted(missing_seeds)}"
        )

    return df


def compute_winners(df_corr: pd.DataFrame) -> np.ndarray:
    """
    Return a 1-D array of winning seed labels (1-based int) per parcel.
    Self-seed parcels are masked to NaN before the argmax.
    Parcels where all seeds are NaN receive NaN.
    """
    df = df_corr.copy()

    for seed, plist in seed_to_parcels.items():
        for p in plist:
            if seed in df.index and p in df.columns:
                df.loc[seed, p] = np.nan

    winners = []
    for parcel in df.columns:
        col = df[parcel]
        winners.append(np.nan if col.isna().all() else seed_to_number[col.idxmax()])

    return np.array(winners)


def append_group_and_consistency(subject_df: pd.DataFrame) -> pd.DataFrame:
    """
    Append GROUP and CONSISTENCY rows to the subject winner table.

    GROUP row:
        Modal winner per parcel across subjects (NaN ignored).
        Ties broken by lowest seed number (scipy default).

    CONSISTENCY row:
        % of subjects whose winner matches GROUP (NaN where GROUP is NaN).
    """
    group_winners = []
    for parcel in subject_df.columns:
        col = subject_df[parcel].dropna()
        if col.empty:
            group_winners.append(np.nan)
        else:
            group_winners.append(float(scipy_stats.mode(col, keepdims=True).mode[0]))

    group_series = pd.Series(group_winners, index=subject_df.columns, name="GROUP")

    consistency = []
    for parcel in subject_df.columns:
        gw = group_series[parcel]
        if pd.isna(gw):
            consistency.append(np.nan)
        else:
            consistency.append(
                100.0 * (subject_df[parcel] == gw).sum() / len(subject_df)
            )

    consistency_series = pd.Series(
        consistency, index=subject_df.columns, name="CONSISTENCY_%"
    )

    return pd.concat(
        [subject_df, group_series.to_frame().T, consistency_series.to_frame().T]
    )

# ============================================================
# Main loop — hemisphere × variant
# ============================================================

for hemi in ("lh", "rh"):
    print(f"\n{'='*80}")
    print(f"Processing hemisphere: {hemi.upper()}")
    print("=" * 80)

    for variant, (normal_tag, excluded_tag, skip_excluded) in VARIANTS.items():
        print(f"\n  --- Variant: {variant} ---")

        all_winners: List[np.ndarray] = []
        subject_ids: List[str]        = []
        parcel_columns: Optional[List[str]] = None

        for subject in subjects:
            is_excluded = subject in RUN02_EXCLUDED

            if is_excluded and skip_excluded:
                print(f"    WARNING [{variant}]: {subject} is in RUN02_EXCLUDED — SKIPPED")
                continue

            run_tag = excluded_tag if is_excluded else normal_tag

            df_corr = load_corr_matrix(subject, hemi, run_tag)
            if df_corr is None:
                print(f"    {subject}: SKIPPED (missing file)")
                continue

            # Capture parcel column order from the first valid subject
            if parcel_columns is None:
                parcel_columns = list(df_corr.columns)

            if list(df_corr.columns) != parcel_columns:
                raise ValueError(
                    f"Column order mismatch for {subject} {hemi} {variant} — "
                    "parcels are not consistent across subjects."
                )

            if variant == "concat_clean" and is_excluded:
                print(f"    {subject}: OK (fallback → run-01)")
            else:
                print(f"    {subject}: OK")

            all_winners.append(compute_winners(df_corr))
            subject_ids.append(subject)

        if not all_winners:
            print(f"    ERROR: no valid subjects for {hemi} / {variant} — skipping.")
            continue

        subject_df  = pd.DataFrame(all_winners, index=subject_ids, columns=parcel_columns)
        combined_df = append_group_and_consistency(subject_df)

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