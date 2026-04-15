"""
Created on April 12, 2025

wta_full_corr_by_subject_by_hemi.py
-----------------------------------------------------------------------------------------
Goal of the script:
Compute winner take all on full correlation results from workbench outputs
Written to run the same pipeline as with Nilearn outputs

-----------------------------------------------------------------------------------------
Goal:
    Compute winner-take-all (WTA) from per-subject per-ROI TSV files
    Outputs one winner table and one consistency table per hemisphere

    Each TSV has:
        - 106 rows : parcels in atlas-key order
                     rows  0–52  → right hemisphere (RH_PARCEL_ORDER, keys 1–163)
                     rows 53–105 → left  hemisphere (LH_PARCEL_ORDER, keys 181–343)
        - 1 column : fisher-z correlation of that seed/ROI with each parcel

Pipeline per hemisphere:
    1. For each subject, load one TSV per seed/ROI and slice the hemi-specific rows.
    2. Stack into a (n_seeds × n_parcels_hemi) DataFrame.
    3. Apply compute_winners() → winning seed label per parcel.
    4. Collect across subjects + group → save winner table and consistency table.
-----------------------------------------------------------------------------------------

Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: server group (e.g. 327)
sys.argv[4]: server project (eg b327)
sys.argv[5]: grab legacy outputs or default
-----------------------------------------------------------------------------------------
Output(s):
TSV to import into the generate workbench dlabel file scripts
-----------------------------------------------------------------------------------------
To run:
1. cd to function
$ cd projects/pRF_analysis/RetinoMaps/rest/stats
2. run python command

-----------------------------------------------------------------------------------------
Examples:

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
# Personal imports — adjust base_dir as needed
# ============================================================
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../"))
sys.path.append(os.path.abspath(os.path.join(base_dir, "analysis_code/utils")))
from settings_utils import load_settings

# ============================================================
# Inputs
# ============================================================
main_dir = sys.argv[1]
project_dir = sys.argv[2]
group = sys.argv[3]
server = sys.argv[4]
mode = sys.argv[5] # legacy, default

# ============================================================
# Load project settings
# ============================================================
project_dir       = "RetinoMaps"
settings_path     = os.path.join(base_dir, project_dir, "settings.yml")
prf_settings_path = os.path.join(base_dir, project_dir, "prf-analysis.yml")
settings          = load_settings([settings_path, prf_settings_path])
analysis_info     = settings[0]
subjects          = analysis_info["subjects"]

# ============================================================
# Paths
# ============================================================
main_data     = Path("/scratch/mszinte/data/RetinoMaps/derivatives/pp_data")
output_folder = main_data / "group/91k/rest/wta/workbench"
output_folder.mkdir(parents=True, exist_ok=True)

# ============================================================
# Atlas parcel keys
# Both hemispheres share the same 53 parcel *names* in key-sorted order.
# RH keys: 1–163   → TSV rows  0–52
# LH keys: 181–343 → TSV rows 53–105
# ============================================================

R_KEYS: Dict[str, int] = {
    "V1":1,    "MST":2,   "V2":4,    "V3":5,    "V4":6,
    "V8":7,    "FEF":10,  "PEF":11,  "55b":12,  "V3A":13,
    "V7":16,   "IPS1":17, "FFC":18,  "V3B":19,  "LO1":20,
    "LO2":21,  "PIT":22,  "MT":23,   "7Pm":29,  "24dv":41,
    "7AL":42,  "SCEF":43, "6ma":44,  "7Am":45,  "7PL":46,
    "7PC":47,  "LIPv":48, "VIP":49,  "MIP":50,  "6d":54,
    "6mp":55,  "6v":56,   "p32pr":60,"6r":78,   "IFJa":79,
    "IFJp":80, "LIPd":95, "6a":96,   "i6-8":97, "AIP":117,
    "PH":138,  "IP2":144, "IP1":145, "IP0":146, "V6A":152,
    "VMV1":153,"VMV3":154,"V4t":156, "FST":157, "V3CD":158,
    "LO3":159, "VMV2":160,"VVC":163,
}

L_KEYS: Dict[str, int] = {
    "V1":181,  "MST":182, "V2":184,  "V3":185,  "V4":186,
    "V8":187,  "FEF":190, "PEF":191, "55b":192, "V3A":193,
    "V7":196,  "IPS1":197,"FFC":198, "V3B":199, "LO1":200,
    "LO2":201, "PIT":202, "MT":203,  "7Pm":209, "24dv":221,
    "7AL":222, "SCEF":223,"6ma":224, "7Am":225, "7PL":226,
    "7PC":227, "LIPv":228,"VIP":229, "MIP":230, "6d":234,
    "6mp":235, "6v":236,  "p32pr":240,"6r":258,  "IFJa":259,
    "IFJp":260,"LIPd":275,"6a":276,  "i6-8":277,"AIP":297,
    "PH":318,  "IP2":324, "IP1":325, "IP0":326, "V6A":332,
    "VMV1":333,"VMV3":334,"V4t":336, "FST":337, "V3CD":338,
    "LO3":339, "VMV2":340,"VVC":343,
}

# Parcel names in key-sorted order — same sequence for both hemispheres
RH_PARCEL_ORDER: List[str] = [name for name, _ in sorted(R_KEYS.items(), key=lambda x: x[1])]
LH_PARCEL_ORDER: List[str] = [name for name, _ in sorted(L_KEYS.items(), key=lambda x: x[1])]

N_PARCELS_PER_HEMI = 53
N_PARCELS_TOTAL    = 106  # 53 RH + 53 LH per TSV

assert len(RH_PARCEL_ORDER) == N_PARCELS_PER_HEMI
assert len(LH_PARCEL_ORDER) == N_PARCELS_PER_HEMI
assert RH_PARCEL_ORDER == LH_PARCEL_ORDER, (
    "Parcel name order differs between hemispheres — check key tables."
)

# Row slices within each 106-row TSV
HEMI_ROW_SLICE: Dict[str, slice] = {
    "rh": slice(0,  N_PARCELS_PER_HEMI),
    "lh": slice(N_PARCELS_PER_HEMI, N_PARCELS_TOTAL),
}
HEMI_PARCEL_ORDER: Dict[str, List[str]] = {
    "rh": RH_PARCEL_ORDER,
    "lh": LH_PARCEL_ORDER,
}

# ============================================================
# ROI / cluster settings (from original script)
# ============================================================
clusters        = list(analysis_info["rois-drawn"])  # copy before reversing
seed_to_parcels = analysis_info["rois-group-mmp"]    # {seed_name: [parcel_names]}

# Keep the same seed ordering as the original script (mPCS first)
clusters.reverse()

# 1-based integer label for each seed
seed_to_number: Dict[str, int] = {s: i + 1 for i, s in enumerate(clusters)}
number_to_seed: Dict[int, str] = {v: k for k, v in seed_to_number.items()}

# ============================================================
# Helper: build TSV file path
# ============================================================

def tsv_path(subject: str, hemi: str, roi: str) -> Path:
    """Return expected path for one subject / hemi / ROI TSV file."""
    subj_dir = (
        main_data / subject
        / "91k/rest/corr/full_corr/workbench_full_corr/by_hemi"
    )
    fname = (
        f"{subject}_task-rest_space-fsLR_den-91k"
        f"_desc-fisher-z_{hemi}_{roi}_parcellated_legacy-mode.tsv"
    )
    return subj_dir / fname

# ============================================================
# Helper: load seeds × parcels correlation matrix for one subject / hemi
# ============================================================

def load_corr_matrix(subject: str, hemi: str) -> Optional[pd.DataFrame]:
    """
    Load one TSV per seed/ROI, slice the hemi-specific rows (53 of 106),
    and stack into a DataFrame of shape (n_seeds × n_parcels_hemi).

    Rows   = seed names  (from `clusters`)
    Columns = parcel names (from HEMI_PARCEL_ORDER[hemi], key-sorted)

    Returns None if any seed file is missing for this subject/hemi.
    All missing files are reported before returning.
    """
    row_slice    = HEMI_ROW_SLICE[hemi]
    parcel_order = HEMI_PARCEL_ORDER[hemi]

    seed_series: Dict[str, pd.Series] = {}
    missing_files: List[Path] = []

    for seed in clusters:
        fpath = tsv_path(subject, hemi, seed)

        if not fpath.exists():
            missing_files.append(fpath)
            continue

        raw = pd.read_csv(fpath, header=None, sep="\t")

        if raw.shape != (N_PARCELS_TOTAL, 1):
            raise ValueError(
                f"Unexpected shape {raw.shape} in {fpath.name} "
                f"(expected ({N_PARCELS_TOTAL}, 1))."
            )

        hemi_values = raw.iloc[row_slice, 0].values  # 53 correlation values

        seed_series[seed] = pd.Series(
            hemi_values,
            index=parcel_order,
            name=seed,
            dtype=float,
        )

    if missing_files:
        for f in missing_files:
            print(f"  WARNING [{subject} {hemi}]: missing {f.name}")
        return None

    # (n_seeds × n_parcels) — rows are seeds, columns are parcels
    df_corr = pd.DataFrame(seed_series).T
    assert df_corr.shape == (len(clusters), N_PARCELS_PER_HEMI), (
        f"Unexpected matrix shape {df_corr.shape} for {subject} {hemi}."
    )
    return df_corr

# ============================================================
# WTA: compute winning seed per parcel
# ============================================================

def compute_winners(df_corr, seed_to_parcels, seed_to_number):
    """
    df_corr : DataFrame (seeds x parcels)
    """
    df = df_corr.copy()

    # exclude self-seed parcels
    for seed, plist in seed_to_parcels.items():
        for p in plist:
            if seed in df.index and p in df.columns:
                df.loc[seed, p] = np.nan

    winners = []
    for parcel in df.columns:
        col = df[parcel]

        # if entire column is NaN → no winner
        if col.isna().all():
            winners.append(np.nan)
        else:
            winners.append(seed_to_number[col.idxmax()])

    return np.array(winners)

# ============================================================
# Main loop
# ============================================================

for hemi in ("lh", "rh"):
    print(f"\n{'='*80}")
    print(f"Processing hemisphere: {hemi.upper()}")
    print("="*80)

    parcel_order = HEMI_PARCEL_ORDER[hemi]
    all_winners: List[pd.Series] = []
    subject_ids: List[str]       = []

    # ------------------------------------------------------------------
    # Subject-level WTA
    # ------------------------------------------------------------------
    print("Processing subjects...")
    for subject in subjects:
        df_corr = load_corr_matrix(subject, hemi)
        if df_corr is None:
            print(f"  {subject}: SKIPPED")
            continue

        winners = compute_winners(df_corr, seed_to_parcels, seed_to_number)
        all_winners.append(winners)
        subject_ids.append(subject)
        print(f"  {subject}: OK")

    if not all_winners:
        print(f"  ERROR: No valid subjects for {hemi} — skipping hemisphere.")
        continue

    # ------------------------------------------------------------------
    # Group-level WTA
    # ------------------------------------------------------------------
    group_result_path = main_data / "group/91k/rest/full_corr/by_hemi"
    group_csv         = (
        group_result_path
        / f"group_median_cluster_by_mmp-parcel_full_corr_{hemi}_by_hemi.csv"
    )

    if not group_csv.exists():
        print(f"  WARNING: Group file not found: {group_csv}")
        group_winners = pd.Series(np.nan, index=parcel_order)
    else:
        df_group      = pd.read_csv(group_csv, index_col=0)
        group_winners = compute_winners(df_group, seed_to_parcels, seed_to_number)

    all_winners.append(group_winners)
    subject_ids.append("GROUP")

    # ------------------------------------------------------------------
    # Assemble winner table: rows = subjects + GROUP, columns = parcels
    # ------------------------------------------------------------------
    winners_df = pd.DataFrame(all_winners, index=subject_ids)

    # Defensive: flag any parcels that dropped out
    missing_parcels = set(parcel_order) - set(winners_df.columns)
    if missing_parcels:
        print(f"  WARNING: parcels absent from winner table: {sorted(missing_parcels)}")

    out_csv = output_folder / f"winning_seeds_by_subject_full_corr_{hemi}.csv"
    winners_df.to_csv(out_csv)
    print(f"\n  Saved winner table : {out_csv}")

    # ------------------------------------------------------------------
    # Consistency: % of subjects matching the group winner, per parcel
    # ------------------------------------------------------------------
    subject_df = winners_df.drop(index="GROUP", errors="ignore")

    consistency_rows = []
    for parcel in winners_df.columns:
        gw = (
            winners_df.loc["GROUP", parcel]
            if "GROUP" in winners_df.index
            else np.nan
        )
        if pd.isna(gw):
            pct = np.nan
        else:
            pct = 100.0 * (subject_df[parcel] == gw).sum() / len(subject_df)

        consistency_rows.append({
            "Parcel":        parcel,
            "Group_Winner":  gw,
            "Consistency_%": pct,
        })

    consistency_df = (
        pd.DataFrame(consistency_rows)
        .sort_values("Consistency_%", ascending=False)
    )

    consistency_file = output_folder / f"winner_consistency_by_parcel_full_corr_{hemi}.csv"
    consistency_df.to_csv(consistency_file, index=False)
    print(f"  Saved consistency  : {consistency_file}")

print("\n" + "="*80)
print("ALL HEMISPHERES COMPLETE")
print("="*80)
print(f"\nOutputs written to: {output_folder}")