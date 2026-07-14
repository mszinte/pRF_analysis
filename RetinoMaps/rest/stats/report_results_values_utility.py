#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
report_results_values_utility.py
------------------------------------------------------------------------------------------
Reads group-level WTA and correlation CSV files and writes a ready-to-paste
manuscript results paragraph covering:

    1. WTA in-scope parcel counts  (X out of Y parcels per hemisphere per corr type)
    2. Consistency ranges per macro-region (min% (parcel) – max% (parcel))
    3. Task-constrained correlation values (median r [p25; p75] per seed × target)

Output paragraph order matches the manuscript:
    - IPEF first (sIPS, iIPS), then PCEF (mPCS, sPCS, iPCS)
    - Full corr block, then partial corr block
    - Within each block: LH then RH

Reads:
    WTA consistency CSVs (from wta_*_by_subject_by_hemi.py):
        {wta_folder}/consistency_ranges_{corr_type}_concat_clean{_mode}.csv

    Task-constrained correlation CSVs (from group stats scripts):
        {group_folder}/seed-task_by_macror-task_{corr_type}_r_{stat}_{variant}_{hemi}{_mode}.csv

Usage:
    python report_results_values_utility.py
    (edit the CONFIG block at the bottom)
------------------------------------------------------------------------------------------
Written by Marco Bedini (marco.bedini@univ-amu.fr)
------------------------------------------------------------------------------------------
"""

import os
import sys
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Optional, List, Tuple, Dict

# ============================================================
# Region order constants
# ============================================================

# IPEF first, then PCEF — matches manuscript paragraph order
IPEF = ["sIPS", "iIPS"]
PCEF = ["mPCS", "sPCS", "iPCS"]
TARGET_MACROS = IPEF + PCEF   # display order for consistency block

EYE_FIELDS = ["mPCS", "sPCS", "iPCS", "sIPS", "iIPS"]   # for corr values block

# Total parcels in the 5 macro-regions (denominator for the counts sentence)
# Computed from rois-group-mmp in settings; hardcoded here for standalone use.
# mPCS(3) + sPCS(6) + iPCS(6) + sIPS(9) + iIPS(7) = 31
# NOTE: update if YAML changes
TOTAL_PARCELS_IN_SCOPE = 31

# ============================================================
# Section 1: WTA in-scope counts from consistency CSV
# ============================================================

def load_consistency_csv(
    wta_folder:  Path,
    corr_type:   str,          # "full_corr" or "partial_corr"
    mode_suffix: str = "",     # e.g. "_default" or ""
) -> Optional[pd.DataFrame]:
    """
    Load the consistency ranges CSV saved by the WTA scripts.

    Expected columns: hemi, macro_region, n_eligible_parcels,
                      min_pct, min_parcel, min_winner,
                      max_pct, max_parcel, max_winner
    """
    fname = f"consistency_ranges_{corr_type}_concat_clean{mode_suffix}.csv"
    fpath = wta_folder / fname
    if not fpath.exists():
        raise FileNotFoundError(
            f"Consistency CSV not found: {fpath}\n"
            "Run the WTA script first to generate it."
        )
    return pd.read_csv(fpath)


def count_inscope_parcels(df: pd.DataFrame, hemi: str) -> int:
    """
    Sum n_eligible_parcels across all macro-regions for one hemisphere.
    This gives the count of parcels won by an in-scope seed.
    """
    sub = df[df["hemi"] == hemi]
    return int(sub["n_eligible_parcels"].sum())


# ============================================================
# Section 2: Consistency range text block
# ============================================================

def _format_range(row: pd.Series) -> str:
    """Format one macro-region range as: XX% (parcel) – XX% (parcel)"""
    if pd.isna(row["min_pct"]):
        return "no eligible parcels"
    return (
        f"{row['min_pct']:.0f}% ({row['min_parcel']})"
        f" \u2013 "
        f"{row['max_pct']:.0f}% ({row['max_parcel']})"
    )


def build_consistency_block(
    df:         pd.DataFrame,
    hemi:       str,
    hemi_label: str,
) -> str:
    """
    Build the consistency sub-clause for one hemisphere:
        left hemisphere: sIPS, XX% (P) – XX% (P); iIPS, XX% (P) – XX% (P);
                         mPCS, ...; sPCS, ...; iPCS, ...
    """
    sub   = df[df["hemi"] == hemi].set_index("macro_region")
    parts = []
    for macro in TARGET_MACROS:
        if macro not in sub.index:
            parts.append(f"{macro}, no data")
            continue
        parts.append(f"{macro}, {_format_range(sub.loc[macro])}")
    return f"{hemi_label}: " + "; ".join(parts)


# ============================================================
# Section 3: Task-constrained correlation values
# ============================================================

def _load_corr_csv(
    folder:     Path,
    corr_type:  str,
    stat:       str,
    variant:    str,
    hemi:       str,
    mode_suffix: str = "",
) -> pd.DataFrame:
    fname = (
        f"seed-task_by_macror-task_{corr_type}"
        f"_r_{stat}_{variant}_{hemi}{mode_suffix}.csv"
    )
    fpath = folder / fname
    if not fpath.exists():
        raise FileNotFoundError(
            f"Correlation CSV not found: {fpath}"
        )
    return pd.read_csv(fpath, index_col=0)


def build_corr_values_block(
    folder:      Path,
    corr_type:   str,
    variant:     str,
    hemi:        str,
    hemi_label:  str,
    pct_lo_tag:  str = "p25",
    pct_hi_tag:  str = "p75",
    mode_suffix: str = "",
) -> str:
    """
    Build the correlation values sub-clause for one hemisphere:
        left hemisphere: mPCS seed: r = X.XX [X.XX; X.XX], sPCS seed: ...
    For each target in TARGET_MACROS, list all non-self seeds.
    """
    df_med = _load_corr_csv(folder, corr_type, "median", variant, hemi, mode_suffix)
    df_plo = _load_corr_csv(folder, corr_type, pct_lo_tag, variant, hemi, mode_suffix)
    df_phi = _load_corr_csv(folder, corr_type, pct_hi_tag, variant, hemi, mode_suffix)

    target_parts = []
    for target in TARGET_MACROS:
        col = f"{target}_ipsi"
        seed_parts = []
        for seed in EYE_FIELDS:
            if seed == target:
                continue
            med = df_med.loc[seed, col]
            plo = df_plo.loc[seed, col]
            phi = df_phi.loc[seed, col]
            seed_parts.append(
                f"{seed}: r\u00a0=\u00a0{med:.2f} [{plo:.2f};\u00a0{phi:.2f}]"
            )
        target_parts.append(
            f"{target} target ({', '.join(seed_parts)})"
        )

    return f"{hemi_label}: " + "; ".join(target_parts)


# ============================================================
# Master paragraph builder
# ============================================================

def generate_results_paragraph(
    # WTA folders (contain consistency_ranges_*.csv)
    wta_full_corr_folder:    Path,
    wta_partial_corr_folder: Path,
    # Group correlation folders (contain seed-task_by_macror-task_*.csv)
    group_full_corr_folder:    Path,
    group_partial_corr_folder: Path,
    # Options
    variant:        str = "concat_clean",
    pct_lo_tag:     str = "p25",
    pct_hi_tag:     str = "p75",
    full_mode_suffix:    str = "_default",   # e.g. "_default" or "_legacy"
    partial_mode_suffix: str = "",           # partial corr has no mode
    total_parcels:  int = TOTAL_PARCELS_IN_SCOPE,
) -> str:
    """
    Generate the full manuscript results paragraph.

    Returns a single string ready to paste into the manuscript.
    """

    # ── Load consistency CSVs ──────────────────────────────────────────────
    df_cons_full    = load_consistency_csv(
        wta_full_corr_folder,    "full_corr",    full_mode_suffix
    )
    df_cons_partial = load_consistency_csv(
        wta_partial_corr_folder, "partial_corr", partial_mode_suffix
    )

    # ── In-scope parcel counts ─────────────────────────────────────────────
    counts = {}
    for label, df in [("full", df_cons_full), ("partial", df_cons_partial)]:
        for hemi in ("lh", "rh"):
            counts[(label, hemi)] = count_inscope_parcels(df, hemi)

    counts_sentence = (
        "This winner-take-all (WTA) analysis showed that PCEF and IPEF were "
        "preferentially connected to parcels contained within the five frontal "
        "and parietal macro-regions, in both the left hemisphere "
        f"(full correlation: {counts[('full','lh')]} out of the "
        f"{total_parcels} parcels in the macro-regions; "
        f"partial correlation: {counts[('partial','lh')]}/{total_parcels}) "
        f"and right hemisphere "
        f"(full correlation: {counts[('full','rh')]}/{total_parcels}; "
        f"partial correlation: {counts[('partial','rh')]}/{total_parcels})."
    )

    # ── Consistency range paragraphs ───────────────────────────────────────
    def _cons_paragraph(df: pd.DataFrame, corr_label: str, supp_ref: str) -> str:
        lh_block = build_consistency_block(df, "lh", "left hemisphere")
        rh_block = build_consistency_block(df, "rh", "right hemisphere")

        # Split into IPEF and PCEF sub-clauses for each hemisphere
        def _split_hemi(block: str) -> Tuple[str, str]:
            # block is e.g. "left hemisphere: sIPS, ...; iIPS, ...; mPCS, ...; ..."
            hemi_prefix, content = block.split(": ", 1)
            items = [x.strip() for x in content.split(";")]
            ipef_items = [x for x in items if any(x.startswith(m) for m in IPEF)]
            pcef_items = [x for x in items if any(x.startswith(m) for m in PCEF)]
            return (
                f"{hemi_prefix}: " + "; ".join(ipef_items),
                "; ".join(pcef_items),
            )

        lh_ipef, lh_pcef = _split_hemi(lh_block)
        rh_ipef, rh_pcef = _split_hemi(rh_block)

        return (
            f"Winner assignments for {corr_label} results were quite consistent "
            f"in the IPEF ({lh_ipef}; {rh_ipef}) "
            f"and somewhat more variable in the PCEF "
            f"({lh_pcef}; right hemisphere: {rh_pcef}; "
            f"see {supp_ref} for individual maps)."
        )

    cons_full    = _cons_paragraph(df_cons_full,    "full correlation",    "Supplementary Figure 3a")
    cons_partial = _cons_paragraph(df_cons_partial, "partial correlation", "Supplementary Figure 3a")

    # ── Assemble final paragraph ───────────────────────────────────────────
    paragraph = "\n\n".join([
        counts_sentence,
        cons_full,
        cons_partial,
    ])

    return paragraph


# ============================================================
# CONFIG — edit before running
# ============================================================

if __name__ == "__main__":

    MAIN_DATA = Path("/scratch/mszinte/data/RetinoMaps/derivatives/pp_data")

    # WTA output folders (where consistency_ranges_*.csv files live)
    WTA_FULL_CORR_FOLDER    = MAIN_DATA / "group/91k/rest/wta/workbench"
    WTA_PARTIAL_CORR_FOLDER = MAIN_DATA / "group/91k/rest/wta/nilearn"

    # Group correlation folders (task-constrained, for the r-value block)
    GROUP_FULL_CORR_FOLDER    = MAIN_DATA / "group/91k/rest/full_corr/by_hemi/task-constrained"
    GROUP_PARTIAL_CORR_FOLDER = MAIN_DATA / "group/91k/rest/partial_corr/by_hemi/task-constrained"

    VARIANT          = "concat_clean"
    PCT_LO_TAG       = "p25"
    PCT_HI_TAG       = "p75"
    FULL_MODE_SUFFIX = "_default"   # appended to full-corr filenames
    PARTIAL_MODE_SUFFIX = ""        # partial corr has no mode token

    paragraph = generate_results_paragraph(
        wta_full_corr_folder    = WTA_FULL_CORR_FOLDER,
        wta_partial_corr_folder = WTA_PARTIAL_CORR_FOLDER,
        group_full_corr_folder    = GROUP_FULL_CORR_FOLDER,
        group_partial_corr_folder = GROUP_PARTIAL_CORR_FOLDER,
        variant             = VARIANT,
        pct_lo_tag          = PCT_LO_TAG,
        pct_hi_tag          = PCT_HI_TAG,
        full_mode_suffix    = FULL_MODE_SUFFIX,
        partial_mode_suffix = PARTIAL_MODE_SUFFIX,
    )

    print(paragraph)
    print()

    # Save to results_reports folders
    for label, folder in [
        ("full-corr",    WTA_FULL_CORR_FOLDER    / "results_reports"),
        ("partial-corr", WTA_PARTIAL_CORR_FOLDER / "results_reports"),
    ]:
        folder.mkdir(parents=True, exist_ok=True)
        out = folder / f"results_paragraph_wta_{label}_{VARIANT}.txt"
        out.write_text(paragraph, encoding="utf-8")
        print(f"Saved: {out}")