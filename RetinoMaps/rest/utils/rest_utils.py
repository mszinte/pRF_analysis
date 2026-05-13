#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
rest_utils.py
-----------------------------------------------------------------------------
Shared utilities for the RetinoMaps resting-state pipeline

Example import:

    import os, sys
    base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../"))
    sys.path.append(os.path.abspath(os.path.join(base_dir, "RetinoMaps/rest/utils")))
    from rest_utils import (
        RUN02_EXCLUDED,
        VARIANTS,
        MODE_SUFFIX,
        ATLAS_KEY_TABLES,
        build_remap,
        tsv_path,
        load_full_corr_matrix,
        compute_winners,
        append_group_and_consistency,
    )

Design principles
-----------------
- Functions that require pipeline-level state (paths, mode suffix, parcel
  lists) accept those as explicit arguments rather than relying on module-
  level globals — this keeps them testable in isolation and safe to call
  from scripts that set up their own argument parsing
- Constants that are project-wide and never change (excluded subjects, atlas
  key tables) are module-level and exported directly
- Averaging (Fisher-z group mean/median/std) lives in group_full_corr_by_hemi.py
  and group_partial_corr_by_hemi.py — those scripts own the decision about
  which statistic to compute
-----------------------------------------------------------------------------
Written by Marco Bedini (marco.bedini@univ-amu.fr)
-----------------------------------------------------------------------------
"""

import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, List, Optional, Tuple
from scipy import stats as scipy_stats

# ============================================================
# Run variant configuration
# ============================================================

# Subjects excluded from run-02 (bad data / registration error).
# Used across all WTA and group-stats scripts.
RUN02_EXCLUDED: frozenset = frozenset(
    {"sub-03", "sub-04", "sub-14", "sub-21", "sub-22", "sub-23"}
)

# Canonical four-variant table.
#   Keys   : variant name (used in output filenames)
#   Values : (normal_tag, excluded_tag, skip_excluded)
#     normal_tag   : run_tag for non-excluded subjects (None = concat file)
#     excluded_tag : run_tag for RUN02_EXCLUDED subjects (None = concat file)
#     skip_excluded: if True, RUN02_EXCLUDED subjects are skipped with WARNING
#                    (run-02 only — their TSV data is genuinely unusable)
VARIANTS: Dict[str, Tuple[Optional[str], Optional[str], bool]] = {
    "concat":       (None,     None,     False),
    "concat_clean": (None,     "run-01", False),
    "run-01":       ("run-01", "run-01", False),
    "run-02":       ("run-02", None,     True),
}

# Parcellation mode → TSV filename suffix
MODE_SUFFIX: Dict[str, str] = {
    "default":     "",
    "legacy":      "_legacy-mode",
    "no_outliers": "_no_outliers",
}

# ============================================================
# Atlas parcel key tables
#
# Define the row ordering inside workbench-parcellated TSV files.
# Each TSV has 106 rows total:
#   rows  0–52  → RH parcels sorted by atlas key (1–163)
#   rows 53–105 → LH parcels sorted by atlas key (181–343)
# ============================================================
ATLAS_KEY_TABLES: Dict[str, Dict[str, int]] = {
    "rh": {
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
    },
    "lh": {
        "V1":181,  "MST":182, "V2":184,  "V3":185,  "V4":186,
        "V8":187,  "FEF":190, "PEF":191, "55b":192, "V3A":193,
        "V7":196,  "IPS1":197,"FFC":198, "V3B":199, "LO1":200,
        "LO2":201, "PIT":202, "MT":203,  "7Pm":209, "24dv":221,
        "7AL":222, "SCEF":223,"6ma":224, "7Am":225, "7PL":226,
        "7PC":227, "LIPv":228,"VIP":229, "MIP":230, "6d":234,
        "6mp":235, "6v":236,  "p32pr":240,"6r":258, "IFJa":259,
        "IFJp":260,"LIPd":275,"6a":276,  "i6-8":277,"AIP":297,
        "PH":318,  "IP2":324, "IP1":325, "IP0":326, "V6A":332,
        "VMV1":333,"VMV3":334,"V4t":336, "FST":337, "V3CD":338,
        "LO3":339, "VMV2":340,"VVC":343,
    },
}

N_PARCELS_PER_HEMI: int = 53
N_PARCELS_TOTAL:    int = 106

# Atlas-key-sorted parcel name list (actual TSV row order) — derived once
ATLAS_ORDER: Dict[str, List[str]] = {
    hemi: [n for n, _ in sorted(keys.items(), key=lambda x: x[1])]
    for hemi, keys in ATLAS_KEY_TABLES.items()
}

assert len(ATLAS_ORDER["rh"]) == N_PARCELS_PER_HEMI
assert len(ATLAS_ORDER["lh"]) == N_PARCELS_PER_HEMI
assert ATLAS_ORDER["rh"] == ATLAS_ORDER["lh"], (
    "Parcel name order differs between RH and LH atlas key tables — "
    "check ATLAS_KEY_TABLES."
)

HEMI_ROW_SLICE: Dict[str, slice] = {
    "rh": slice(0, N_PARCELS_PER_HEMI),
    "lh": slice(N_PARCELS_PER_HEMI, N_PARCELS_TOTAL),
}

# ============================================================
# Remap: atlas-key order → canonical YAML parcel order
# ============================================================

def build_remap(
    atlas_order: List[str],
    canonical_order: List[str],
    hemi: str,
) -> Tuple[List[int], List[str]]:
    """
    Compute the index array that converts from atlas-key-sorted TSV rows to
    canonical YAML parcel order (rois-group-mmp flattened in cluster order).

    Parameters
    ----------
    atlas_order    : parcel names in atlas-key order (from ATLAS_ORDER[hemi])
    canonical_order: parcel names in YAML canonical order (parcels list)
    hemi           : "rh" or "lh" — used only for error/warning messages

    Returns
    -------
    remap_idx : indices into atlas_order such that
                ``atlas_values[remap_idx]`` gives values in canonical_order
    present   : subset of canonical_order present in the TSV;
                absent parcels receive NaN in the output

    Raises
    ------
    ValueError if the atlas key table contains parcels absent from the YAML.
    """
    atlas_set     = set(atlas_order)
    canonical_set = set(canonical_order)

    extra_in_tsv = atlas_set - canonical_set
    if extra_in_tsv:
        raise ValueError(
            f"[{hemi.upper()}] ATLAS_KEY_TABLES contains parcels absent from "
            f"YAML rois-group-mmp: {sorted(extra_in_tsv)}\n"
            "  → Update ATLAS_KEY_TABLES in rest_utils.py or the YAML."
        )

    missing_in_tsv = canonical_set - atlas_set
    if missing_in_tsv:
        print(
            f"  WARNING [{hemi.upper()}]: {len(missing_in_tsv)} YAML parcel(s) absent "
            f"from atlas key table (will be NaN in output): {sorted(missing_in_tsv)}"
        )

    atlas_index = {name: i for i, name in enumerate(atlas_order)}
    remap_idx: List[int] = []
    present:   List[str] = []
    for name in canonical_order:
        if name in atlas_index:
            remap_idx.append(atlas_index[name])
            present.append(name)

    return remap_idx, present


# ============================================================
# TSV I/O  (full correlation — workbench parcellated files)
# ============================================================

def tsv_path(
    subject: str,
    hemi: str,
    roi: str,
    run_tag: Optional[str],
    main_data: Path,
    tsv_suffix: str,
) -> Path:
    """
    Return the expected path for one subject / hemi / seed TSV file.

    Parameters
    ----------
    subject    : e.g. "sub-01"
    hemi       : "lh" or "rh"
    roi        : seed/cluster name (e.g. "mPCS")
    run_tag    : None → concat file; str → per-run file (e.g. "run-01")
    main_data  : Path to {main_dir}/{project_dir}/derivatives/pp_data
    tsv_suffix : from MODE_SUFFIX[mode], e.g. "" or "_legacy-mode"
    """
    subj_dir   = main_data / subject / "91k/rest/corr/full_corr/by_hemi"
    run_entity = f"_{run_tag}" if run_tag is not None else ""
    fname = (
        f"{subject}_task-rest{run_entity}_space-fsLR_den-91k"
        f"_desc-fisher-z_{hemi}_{roi}"
        f"_parcellated{tsv_suffix}.tsv"
    )
    return subj_dir / fname


def load_full_corr_matrix(
    subject: str,
    hemi: str,
    run_tag: Optional[str],
    clusters: List[str],
    parcels: List[str],
    remap_idx: List[int],
    present: List[str],
    main_data: Path,
    tsv_suffix: str,
) -> Optional[pd.DataFrame]:
    """
    Load all seed TSVs for one subject/hemi/run_tag, remap to canonical YAML
    parcel order, and return a (n_seeds × n_parcels) DataFrame.

    Returns None if any seed file is missing; all missing paths are printed.
    Raises ValueError on unexpected TSV shape.
    """
    row_slice       = HEMI_ROW_SLICE[hemi]
    n_parcels       = len(parcels)
    n_clusters      = len(clusters)
    present_col_idx = [parcels.index(p) for p in present]

    seed_rows: Dict[str, np.ndarray] = {}
    missing:   List[Path]            = []

    for seed in clusters:
        fpath = tsv_path(subject, hemi, seed, run_tag, main_data, tsv_suffix)

        if not fpath.exists():
            missing.append(fpath)
            continue

        raw = pd.read_csv(fpath, header=None, sep="\t")

        if raw.shape != (N_PARCELS_TOTAL, 1):
            raise ValueError(
                f"[{subject} {hemi}] Unexpected shape {raw.shape} in {fpath.name} "
                f"(expected ({N_PARCELS_TOTAL}, 1))."
            )

        hemi_values = raw.iloc[row_slice, 0].values.astype(float)

        row = np.full(n_parcels, np.nan, dtype=float)
        row[present_col_idx] = hemi_values[remap_idx]
        seed_rows[seed] = row

    if missing:
        for f in missing:
            print(f"  WARNING [{subject} {hemi}]: missing {f.name}")
        return None

    df = pd.DataFrame(seed_rows, index=parcels).T   # (n_seeds × n_parcels)
    df.index.name   = None
    df.columns.name = None

    assert df.shape == (n_clusters, n_parcels), (
        f"[{subject} {hemi}] Unexpected matrix shape {df.shape}."
    )
    return df


# ============================================================
# Winner-take-all
# ============================================================

def compute_winners(
    df_corr: pd.DataFrame,
    seed_to_parcels: Dict[str, List[str]],
    seed_to_number: Dict[str, int],
) -> np.ndarray:
    """
    Return a 1-D array of winning seed labels (1-based int) per parcel.

    Self-seed parcels are masked to NaN before argmax to avoid circularity.
    Parcels where all seeds are NaN receive NaN.

    Parameters
    ----------
    df_corr         : (n_seeds × n_parcels) DataFrame; index = seed names,
                      columns = parcel names in canonical order.
    seed_to_parcels : {seed_name: [parcel_names]} from YAML rois-group-mmp.
    seed_to_number  : {seed_name: 1-based integer label}.

    Notes
    -----
    col.idxmax() raises on an all-NaN column — the col.isna().all() guard
    is essential and must not be removed.
    """
    df = df_corr.copy()

    for seed, plist in seed_to_parcels.items():
        for p in plist:
            if seed in df.index and p in df.columns:
                df.loc[seed, p] = np.nan

    winners = []
    for parcel in df.columns:
        col = df[parcel]
        winners.append(
            np.nan if col.isna().all() else seed_to_number[col.idxmax()]
        )

    return np.array(winners)


def append_group_and_consistency(subject_df: pd.DataFrame) -> pd.DataFrame:
    """
    Append GROUP and CONSISTENCY_% rows to a per-subject winner table.

    Parameters
    ----------
    subject_df : DataFrame with one row per subject, one column per parcel,
                 values are 1-based winner seed integers (or NaN).

    Returns
    -------
    DataFrame: rows = subjects + GROUP + CONSISTENCY_%

    GROUP
        Modal winner per parcel across subjects (NaN ignored).
        Ties broken by lowest seed number (scipy.stats.mode default).

    CONSISTENCY_%
        Percentage of subjects whose winner matches GROUP.
        NaN where GROUP is NaN.
    """
    group_winners = []
    for parcel in subject_df.columns:
        col = subject_df[parcel].dropna()
        if col.empty:
            group_winners.append(np.nan)
        else:
            group_winners.append(
                float(scipy_stats.mode(col, keepdims=True).mode[0])
            )

    group_series = pd.Series(
        group_winners, index=subject_df.columns, name="GROUP"
    )

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