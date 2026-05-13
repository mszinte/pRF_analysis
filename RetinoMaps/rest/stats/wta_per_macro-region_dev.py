"""
Created on May 3, 2025

wta_per_macro_region.py
-----------------------------------------------------------------------------------------
Goal:
    Compute macro-region-level winner-take-all (WTA) by aggregating the
    per-parcel WTA outputs produced by:
        - wta_full_corr_by_subject_by_hemi.py
        - wta_partial_corr_by_subject_by_hemi.py

    For each subject and each TARGET macro-region, find which SOURCE
    macro-region (seed) wins the most parcels inside it — excluding the
    target macro-region's own seed (self-exclusion).

    Two parallel scoring methods are computed:
        raw      — each parcel contributes 1 vote to its winning seed's
                   macro-region tally
        weighted — each parcel contributes its CONSISTENCY_% value (from the
                   GROUP row of the input table) as a fractional vote

    Output tables mirror the per-parcel WTA structure:
        rows    = subjects + GROUP + CONSISTENCY
        columns = target macro-regions (one winner label per target)

    Two output CSVs are saved per hemisphere × variant × corr_type:
        ..._raw.csv
        ..._weighted.csv
-----------------------------------------------------------------------------------------
Inputs (sys.argv):
    1: main project directory   (e.g. /scratch/mszinte/data)
    2: project name/directory   (e.g. RetinoMaps)
    3: server group             (e.g. 327)
    4: server project           (e.g. b327)
    5: parcellation mode for full_corr inputs:
           "default"     → files ending _default.csv
           "legacy"      → files ending _legacy.csv
           "no_outliers" → files ending _no_outliers.csv

Output:
    Per hemisphere × variant × corr_type × scoring:
        wta_macro_region_{corr_type}_{hemi}_{variant}_{mode}_raw.csv
        wta_macro_region_{corr_type}_{hemi}_{variant}_{mode}_weighted.csv
    (partial corr outputs omit the mode token)

To run:
    $ cd projects/pRF_analysis/RetinoMaps/rest/stats
    $ python wta_per_macro-region.py /scratch/mszinte/data RetinoMaps 327 b327 default
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
from typing import Dict, List, Optional, Tuple
from scipy import stats as scipy_stats

# ============================================================
# Personal imports
# ============================================================
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../"))
sys.path.append(os.path.abspath(os.path.join(base_dir, "analysis_code/utils")))
from settings_utils import load_settings

# ============================================================
# Parcellation mode → filename token (full_corr only)
# ============================================================
MODE_SUFFIX: Dict[str, str] = {
    "default":     "",
    "legacy":      "_legacy-mode",
    "no_outliers": "_no_outliers",
}

USAGE = (
    "Usage: python wta_per_macro_region.py "
    "<main_dir> <project_dir> <group> <server> <mode>\n"
    f"  <mode> must be one of: {', '.join(MODE_SUFFIX)}"
)

# ============================================================
# Parse and validate arguments
# ============================================================
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

print("=" * 80)
print("WTA — macro-region level")
print("=" * 80)
print(f"  main_dir    : {main_dir}")
print(f"  project_dir : {project_dir}")
print(f"  group       : {group}")
print(f"  server      : {server}")
print(f"  mode        : {mode!r}")

# ============================================================
# Load settings
# ============================================================
settings_path     = os.path.join(base_dir, project_dir, "settings.yml")
prf_settings_path = os.path.join(base_dir, project_dir, "prf-analysis.yml")
settings          = load_settings([settings_path, prf_settings_path])
analysis_info     = settings[0]

# ============================================================
# Paths
# ============================================================
main_data      = Path(main_dir) / project_dir / "derivatives/pp_data"
wta_wb_folder  = main_data / "group/91k/rest/wta/workbench"   # full_corr inputs
wta_nl_folder  = main_data / "group/91k/rest/wta/nilearn"     # partial_corr inputs
output_folder  = main_data / "group/91k/rest/wta/macro_region"
output_folder.mkdir(parents=True, exist_ok=True)

# ============================================================
# ROI / macro-region definitions  (from settings)
#
# rois-group-mmp: {macro_region: [mmp_parcel, ...]}
# Invert to get:  parcel_to_macro: {mmp_parcel: macro_region}
# Seeds are the macro-region names themselves.
# ============================================================
seed_to_parcels: Dict[str, List[str]] = analysis_info["rois-group-mmp"]

# Ordered macro-region list — mPCS first (same convention as parcel WTA)
macro_regions: List[str] = list(analysis_info["rois-drawn"])
macro_regions.reverse()   # mPCS first

# 1-based label for each macro-region (parallel to seed_to_number in parcel WTA)
macro_to_number: Dict[str, int] = {m: i + 1 for i, m in enumerate(macro_regions)}
number_to_macro: Dict[int, str] = {v: k for k, v in macro_to_number.items()}

# Parcel → macro-region lookup (inverted from rois-group-mmp)
parcel_to_macro: Dict[str, str] = {}
for macro, parcels in seed_to_parcels.items():
    for p in parcels:
        if p in parcel_to_macro:
            raise ValueError(
                f"Parcel '{p}' is assigned to more than one macro-region "
                f"('{parcel_to_macro[p]}' and '{macro}'). Check rois-group-mmp."
            )
        parcel_to_macro[p] = macro

# ============================================================
# Run variants (identical definition to parcel WTA scripts)
# ============================================================
VARIANTS: Dict[str, tuple] = {
    "concat":       (None,     None,     False),
    "concat_clean": (None,     "run-01", False),
    "run-01":       ("run-01", "run-01", False),
    "run-02":       ("run-02", None,     True),
}

# ============================================================
# Input file locators
# ============================================================

def full_corr_path(hemi: str, variant: str) -> Path:
    """Path to per-parcel WTA CSV for full_corr."""
    return (
        wta_wb_folder
        / f"winning_seeds_by_subject_full_corr_{hemi}_{variant}_{mode}.csv"
    )


def partial_corr_path(hemi: str, variant: str) -> Path:
    """Path to per-parcel WTA CSV for partial_corr."""
    return (
        wta_nl_folder
        / f"winning_seeds_by_subject_partial_corr_{hemi}_{variant}.csv"
    )


# ============================================================
# Core: macro-region WTA for one subject row
#
# Given:
#   parcel_winners  — Series, index=parcel names, values=seed macro-region
#                     labels (1-based int or NaN); one row from the input table
#   consistency_row — Series, index=parcel names, values=CONSISTENCY_%
#                     (used as weights for the weighted scoring)
#
# Returns:
#   raw_winner      — 1-based int label of the winning source macro-region
#                     for each TARGET macro-region (raw count argmax)
#   weighted_winner — same but argmax of consistency-weighted tally
# ============================================================

def macro_winners_for_subject(
    parcel_winners: pd.Series,
    consistency_row: pd.Series,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    For each target macro-region, find the source macro-region whose parcels
    win the most parcels inside it (raw and consistency-weighted).

    Self-exclusion: the target macro-region's own seed is excluded from
    competition (its parcels are masked before the argmax, mirroring the
    per-parcel WTA self-exclusion logic).

    Parameters
    ----------
    parcel_winners  : Series (index=parcel, value=1-based seed label or NaN)
    consistency_row : Series (index=parcel, value=CONSISTENCY_% float or NaN)

    Returns
    -------
    raw_winners      : ndarray, shape (n_macro_regions,), 1-based int or NaN
    weighted_winners : ndarray, shape (n_macro_regions,), 1-based int or NaN
    """
    raw_winners      = []
    weighted_winners = []

    for target_macro in macro_regions:
        target_parcels = seed_to_parcels[target_macro]  # parcels in this target

        # Accumulate votes per source macro-region
        raw_tally      = {m: 0.0 for m in macro_regions}
        weighted_tally = {m: 0.0 for m in macro_regions}

        for parcel in target_parcels:
            if parcel not in parcel_winners.index:
                continue

            winner_label = parcel_winners[parcel]
            if pd.isna(winner_label):
                continue

            source_macro = number_to_macro[int(winner_label)]

            # Self-exclusion: skip parcels whose winner is the target itself
            if source_macro == target_macro:
                continue

            weight = consistency_row.get(parcel, np.nan)
            raw_tally[source_macro]      += 1.0
            weighted_tally[source_macro] += 0.0 if pd.isna(weight) else weight

        # Argmax — NaN if no valid votes
        def tally_argmax(tally: Dict[str, float]) -> float:
            total = sum(tally.values())
            if total == 0:
                return np.nan
            best = max(tally, key=tally.__getitem__)
            return float(macro_to_number[best])

        raw_winners.append(tally_argmax(raw_tally))
        weighted_winners.append(tally_argmax(weighted_tally))

    return np.array(raw_winners), np.array(weighted_winners)


# ============================================================
# Group + consistency rows  (same logic as parcel WTA)
# ============================================================

def append_group_and_consistency(subject_df: pd.DataFrame) -> pd.DataFrame:
    """
    Append GROUP (modal winner) and CONSISTENCY_% rows to subject table.
    Ties broken by lowest macro-region number (scipy default).
    """
    group_winners = []
    for col in subject_df.columns:
        vals = subject_df[col].dropna()
        if vals.empty:
            group_winners.append(np.nan)
        else:
            group_winners.append(float(scipy_stats.mode(vals, keepdims=True).mode[0]))

    group_series = pd.Series(group_winners, index=subject_df.columns, name="GROUP")

    consistency = []
    for col in subject_df.columns:
        gw = group_series[col]
        if pd.isna(gw):
            consistency.append(np.nan)
        else:
            consistency.append(
                100.0 * (subject_df[col] == gw).sum() / len(subject_df)
            )

    consistency_series = pd.Series(
        consistency, index=subject_df.columns, name="CONSISTENCY_%"
    )

    return pd.concat(
        [subject_df, group_series.to_frame().T, consistency_series.to_frame().T]
    )


# ============================================================
# Process one input WTA table → macro-region WTA tables
# ============================================================

def process_wta_table(
    wta_csv: Path,
    hemi: str,
    variant: str,
    corr_type: str,
    mode_token: str,
) -> None:
    """
    Load a per-parcel WTA table, aggregate to macro-region level,
    and save raw + weighted output CSVs.

    Parameters
    ----------
    wta_csv    : path to the per-parcel WTA CSV (subjects + GROUP + CONSISTENCY rows)
    hemi       : 'lh' or 'rh'
    variant    : one of VARIANTS keys
    corr_type  : 'full_corr' or 'partial_corr'
    mode_token : mode string to embed in output filenames (empty str for partial)
    """
    if not wta_csv.exists():
        print(f"    WARNING: input not found — {wta_csv.name}")
        return

    df_all = pd.read_csv(wta_csv, index_col=0)

    # Separate subject rows from the two summary rows
    if "CONSISTENCY_%" not in df_all.index or "GROUP" not in df_all.index:
        raise ValueError(
            f"Expected GROUP and CONSISTENCY_% rows in {wta_csv.name}. "
            "Was this file produced by the current WTA pipeline?"
        )

    consistency_row = df_all.loc["CONSISTENCY_%"]
    subject_df_in   = df_all.drop(index=["GROUP", "CONSISTENCY_%"])

    raw_rows      = []
    weighted_rows = []
    subject_ids   = []

    for subject in subject_df_in.index:
        parcel_winners = subject_df_in.loc[subject]
        raw_w, weighted_w = macro_winners_for_subject(parcel_winners, consistency_row)

        raw_rows.append(raw_w)
        weighted_rows.append(weighted_w)
        subject_ids.append(subject)

    if not raw_rows:
        print(f"    ERROR: no subjects processed for {corr_type} {hemi} {variant}.")
        return

    raw_df      = pd.DataFrame(raw_rows,      index=subject_ids, columns=macro_regions)
    weighted_df = pd.DataFrame(weighted_rows, index=subject_ids, columns=macro_regions)

    raw_combined      = append_group_and_consistency(raw_df)
    weighted_combined = append_group_and_consistency(weighted_df)

    # Build output filename tokens
    mode_part = f"_{mode_token}" if mode_token else ""
    stem = f"wta_macro_region_{corr_type}_{hemi}_{variant}{mode_part}"

    raw_out      = output_folder / f"{stem}_raw.csv"
    weighted_out = output_folder / f"{stem}_weighted.csv"

    raw_combined.to_csv(raw_out)
    weighted_combined.to_csv(weighted_out)

    print(f"    Saved raw      : {raw_out.name}")
    print(f"    Saved weighted : {weighted_out.name}")


# ============================================================
# Main loop — corr_type × hemisphere × variant
# ============================================================

CORR_CONFIGS = [
    # (label,          path_fn,            mode_token_in_output)
    ("full_corr",    full_corr_path,     mode),
    ("partial_corr", partial_corr_path,  ""),
]

for corr_type, path_fn, mode_token in CORR_CONFIGS:
    print(f"\n{'='*80}")
    print(f"Correlation type: {corr_type}")
    print("=" * 80)

    for hemi in ("lh", "rh"):
        print(f"\n  Hemisphere: {hemi.upper()}")

        for variant in VARIANTS:
            print(f"\n  --- Variant: {variant} ---")
            wta_csv = path_fn(hemi, variant)
            process_wta_table(wta_csv, hemi, variant, corr_type, mode_token)

print("\n" + "=" * 80)
print("ALL DONE")
print("=" * 80)
print(f"\nOutputs written to: {output_folder}")