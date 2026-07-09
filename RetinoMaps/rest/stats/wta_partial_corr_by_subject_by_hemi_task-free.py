#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
wta_partial_corr_by_subject_by_hemi_task-free.py
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
    $ python wta_partial_corr_by_subject_by_hemi_task-free.py /scratch/mszinte/data RetinoMaps 327 b327
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
    "Usage: python wta_partial_corr_by_subject_by_hemi_task-free.py "
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
group_corr_folder  = main_data / "group/91k/rest/partial_corr/by_hemi/task-free"
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
    subj_dir   = main_data / subject / "91k/rest/corr/partial_corr/by_hemi/task-free"
    run_entity = f"_{run_tag}" if run_tag is not None else ""
    fname      = f"seed-task_by_mmp-parcel_partial_fisherz{run_entity}_{hemi}.csv"
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

# ============================================================
# Consistency range report — concat_clean only
#
# For each target macro-region × hemisphere, report the min and max
# CONSISTENCY_% across parcels whose GROUP winner is one of the 5 target
# macro-regions (in-scope only; self-seed and out-of-scope parcels excluded).
#
# Output: one tidy CSV + ready-to-paste results text printed to stdout.
# ============================================================

TARGET_MACROS  = ["mPCS", "sPCS", "iPCS", "sIPS", "iIPS"]
number_to_seed = {v: k for k, v in seed_to_number.items()}


def consistency_ranges() -> None:
    """
    Load the concat_clean WTA CSVs for both hemispheres and compute
    per-macro-region consistency ranges (min/max parcel + value).
    Results are saved to a tidy CSV and printed as ready-to-paste text.
    """
    print(f"\n{'='*80}")
    print("CONSISTENCY RANGE REPORT — partial_corr / concat_clean")
    print("=" * 80)

    rows = []   # for the tidy CSV

    for hemi in ("lh", "rh"):
        wta_fname = (
            f"winning_seeds_by_subject_partial_corr_{hemi}_concat_clean.csv"
        )
        wta_path = output_folder / wta_fname
        if not wta_path.exists():
            print(f"  WARNING: {wta_fname} not found — skipping {hemi.upper()}.")
            continue

        df = pd.read_csv(wta_path, index_col=0)

        if "GROUP" not in df.index or "CONSISTENCY_%" not in df.index:
            print(f"  WARNING: GROUP or CONSISTENCY_% row missing in {wta_fname}.")
            continue

        group_row       = df.loc["GROUP"].astype(float)
        consistency_row = df.loc["CONSISTENCY_%"].astype(float)

        for target_macro in TARGET_MACROS:
            target_parcels = seed_to_parcels[target_macro]

            # Collect eligible parcels: in target macro-region, GROUP winner is
            # an in-scope seed (not the target itself, not out-of-scope)
            eligible = []
            for parcel in target_parcels:
                if parcel not in group_row.index:
                    continue
                winner_label = group_row[parcel]
                if pd.isna(winner_label):
                    continue
                winner_macro = number_to_seed[int(winner_label)]
                if winner_macro == target_macro:
                    continue              # self-seed excluded
                if winner_macro not in TARGET_MACROS:
                    continue              # out-of-scope excluded
                pct = consistency_row.get(parcel, np.nan)
                if pd.isna(pct):
                    continue
                eligible.append((parcel, winner_macro, pct))

            if not eligible:
                print(
                    f"  {hemi.upper()} {target_macro}: no eligible parcels "
                    "(all self-seed or out-of-scope)."
                )
                rows.append({
                    "hemi": hemi, "macro_region": target_macro,
                    "n_eligible_parcels": 0,
                    "min_pct": np.nan, "min_parcel": "", "min_winner": "",
                    "max_pct": np.nan, "max_parcel": "", "max_winner": "",
                })
                continue

            eligible_df = pd.DataFrame(
                eligible, columns=["parcel", "winner_macro", "consistency_pct"]
            ).sort_values("consistency_pct")

            min_row = eligible_df.iloc[0]
            max_row = eligible_df.iloc[-1]

            rows.append({
                "hemi":               hemi,
                "macro_region":       target_macro,
                "n_eligible_parcels": len(eligible_df),
                "min_pct":            round(min_row["consistency_pct"], 1),
                "min_parcel":         min_row["parcel"],
                "min_winner":         min_row["winner_macro"],
                "max_pct":            round(max_row["consistency_pct"], 1),
                "max_parcel":         max_row["parcel"],
                "max_winner":         max_row["winner_macro"],
            })

    if not rows:
        print("  No data to report.")
        return

    report_df = pd.DataFrame(rows)

    # Save tidy CSV
    report_path = output_folder / "consistency_ranges_partial_corr_concat_clean.csv"
    report_df.to_csv(report_path, index=False)
    print(f"\n  Saved tidy table: {report_path.name}")

    # Print ready-to-paste results text
    print("\n  --- Ready-to-paste results text ---\n")
    for hemi in ("lh", "rh"):
        hemi_label = "left hemisphere" if hemi == "lh" else "right hemisphere"
        hemi_rows  = report_df[report_df["hemi"] == hemi]
        parts = []
        for _, r in hemi_rows.iterrows():
            if pd.isna(r["min_pct"]):
                parts.append(f"{r['macro_region']} no eligible parcels")
            else:
                parts.append(
                    f"{r['macro_region']} {r['min_pct']:.1f}%"
                    f" ({r['min_parcel']})"
                    f"–{r['max_pct']:.1f}%"
                    f" ({r['max_parcel']})"
                )
        print(f"  {hemi_label}: {'; '.join(parts)}")

    print()


consistency_ranges()