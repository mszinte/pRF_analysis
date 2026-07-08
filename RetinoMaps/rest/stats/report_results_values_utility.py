#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
report_results_values_utility.py
------------------------------------------------------------------------------------------
Reads group-level correlation CSV files and writes a ready-to-paste manuscript
results paragraph for the task-constrained full and partial correlation analyses.

Output format (matching document 25):
    {target} target (left hemisphere: {seed}: r = X.XX [X.XX; X.XX], ...;
                      right hemisphere: {seed}: r = X.XX [X.XX; X.XX], ...) ...

Reads three CSV files per correlation type × variant:
    seed-task_by_macror-task_{corr_type}_r_median_{variant}_{hemi}.csv
    seed-task_by_macror-task_{corr_type}_r_p25_{variant}_{hemi}.csv
    seed-task_by_macror-task_{corr_type}_r_p75_{variant}_{hemi}.csv

Rows    = EYE_FIELDS seeds  (mPCS, sPCS, iPCS, sIPS, iIPS)
Columns = TARGET_COLUMNS    (mPCS_ipsi, sPCS_ipsi, ..., iIPS_contra)

Usage:
    python generate_results_text.py
    
    (edit the CONFIG block at the bottom to point to your group folders)
------------------------------------------------------------------------------------------
Written by Marco Bedini (marco.bedini@univ-amu.fr)
------------------------------------------------------------------------------------------
"""

import os
import pandas as pd
from pathlib import Path
from typing import Optional

# ============================================================
# Core region order — must match EYE_FIELDS in the group stats scripts
# ============================================================
EYE_FIELDS = ["mPCS", "sPCS", "iPCS", "sIPS", "iIPS"]

# ============================================================
# Helper: load the three r-space CSVs for one hemi
# ============================================================

def _load_hemi_csvs(
    folder:     Path,
    corr_type:  str,   # "full-corr" or "partial-corr"
    variant:    str,   # e.g. "concat_clean"
    hemi:       str,   # "lh" or "rh"
    pct_lo_tag: str,   # e.g. "p25"
    pct_hi_tag: str,   # e.g. "p75"
    mode_label: Optional[str] = None,  # "legacy" for full-corr, None for partial
) -> tuple:
    """Load median, p_lo, p_hi DataFrames for one hemisphere.

    Returns (df_median, df_plo, df_phi) with index=seeds, columns=TARGET_COLUMNS.
    Raises FileNotFoundError with a clear message if any file is missing.
    """
    suffix = f"_{mode_label}" if mode_label else ""

    def _path(stat: str) -> Path:
        fname = (
            f"seed-task_by_macror-task_{corr_type}"
            f"_r_{stat}_{variant}_{hemi}{suffix}.csv"
        )
        return folder / fname

    dfs = []
    for stat in ("median", pct_lo_tag, pct_hi_tag):
        p = _path(stat)
        if not p.exists():
            raise FileNotFoundError(
                f"Expected CSV not found: {p}\n"
                f"Check folder, corr_type, variant, hemi, and pct tags."
            )
        dfs.append(pd.read_csv(p, index_col=0))

    return tuple(dfs)   # (df_median, df_plo, df_phi)


# ============================================================
# Main function: generate results text
# ============================================================

def generate_results_text(
    full_corr_folder:  Path,
    partial_corr_folder: Path,
    variant:    str  = "concat_clean",
    pct_lo_tag: str  = "p25",
    pct_hi_tag: str  = "p75",
    mode_label: str  = "legacy",      # suffix for full-corr files
) -> str:
    """
    Generate a manuscript-ready results paragraph for full and partial
    correlations, matching the format in the results section.

    Parameters
    ----------
    full_corr_folder    : Path to the group full-corr output directory
    partial_corr_folder : Path to the group partial-corr output directory
    variant             : run variant label (default "concat_clean")
    pct_lo_tag          : lower percentile file label (default "p25")
    pct_hi_tag          : upper percentile file label (default "p75")
    mode_label          : mode suffix for full-corr files (default "legacy")

    Returns
    -------
    str : formatted results text ready to paste into the manuscript
    """
    lines = []

    for corr_type, folder, ml, label in [
        ("full-corr",    full_corr_folder,    mode_label, "full correlations"),
        ("partial-corr", partial_corr_folder, None,       "partial correlations"),
    ]:
        # Load LH and RH DataFrames
        df_med_lh, df_plo_lh, df_phi_lh = _load_hemi_csvs(
            folder, corr_type, variant, "lh", pct_lo_tag, pct_hi_tag, ml
        )
        df_med_rh, df_plo_rh, df_phi_rh = _load_hemi_csvs(
            folder, corr_type, variant, "rh", pct_lo_tag, pct_hi_tag, ml
        )

        # Build one sentence per target, looping over EYE_FIELDS as targets
        target_phrases = []

        for target in EYE_FIELDS:

            hemi_phrases = []

            for hemi_label, df_med, df_plo, df_phi, col_suffix in [
                ("left hemisphere",  df_med_lh, df_plo_lh, df_phi_lh, "_ipsi"),
                ("right hemisphere", df_med_rh, df_plo_rh, df_phi_rh, "_ipsi"),
            ]:
                # Column name in the CSV for this target in the ipsilateral half
                col = f"{target}{col_suffix}"

                seed_phrases = []

                for seed in EYE_FIELDS:
                    if seed == target:
                        continue   # diagonal: no self-correlation reported

                    med = df_med.loc[seed, col]
                    plo = df_plo.loc[seed, col]
                    phi = df_phi.loc[seed, col]

                    seed_phrases.append(
                        f"{seed} seed: r\u00a0=\u00a0{med:.2f},"
                        f" [{plo:.2f};\u00a0{phi:.2f}]"
                    )

                hemi_phrases.append(
                    f"{hemi_label}: " + ", ".join(seed_phrases)
                )

            target_phrases.append(
                f"{target} target ({'; '.join(hemi_phrases)})"
            )

        # Join all targets into the paragraph sentence
        if len(target_phrases) > 1:
            body = (
                ", ".join(target_phrases[:-1])
                + ", and "
                + target_phrases[-1]
            )
        else:
            body = target_phrases[0]

        pct_note = (
            "reported as a function of the target, "
            "with median r and 25th/75th percentiles across participants"
        )

        lines.append(
            f"{label.capitalize()} ({pct_note}) were as follows: {body}."
        )
        lines.append("")   # blank line between full and partial

    return "\n".join(lines).strip()


# ============================================================
# CONFIG — edit these paths before running
# ============================================================

if __name__ == "__main__":

    MAIN_DATA = Path(
        "/scratch/mszinte/data/RetinoMaps/derivatives/pp_data"
    )

    FULL_CORR_FOLDER    = MAIN_DATA / "group/91k/rest/full_corr/by_hemi/task-constrained"
    PARTIAL_CORR_FOLDER = MAIN_DATA / "group/91k/rest/partial_corr/by_hemi/task-constrained"

    VARIANT    = "concat_clean"
    PCT_LO_TAG = "p25"
    PCT_HI_TAG = "p75"
    MODE_LABEL = "legacy"

    text = generate_results_text(
        full_corr_folder    = FULL_CORR_FOLDER,
        partial_corr_folder = PARTIAL_CORR_FOLDER,
        variant             = VARIANT,
        pct_lo_tag          = PCT_LO_TAG,
        pct_hi_tag          = PCT_HI_TAG,
        mode_label          = MODE_LABEL,
    )

    print(text)
    print()

    # Save full-corr and partial-corr reports to separate results_reports folders
    for report_label, base_folder in [
        ("full-corr_task-constrained",    MAIN_DATA / "group/91k/rest/full_corr/results_reports"),
        ("partial-corr_task-constrained", MAIN_DATA / "group/91k/rest/partial_corr/results_reports"),
    ]:
        base_folder.mkdir(parents=True, exist_ok=True)
        out_path = base_folder / f"results_text_{report_label}_{VARIANT}.txt"
        out_path.write_text(text, encoding="utf-8")
        print(f"Saved: {out_path}")