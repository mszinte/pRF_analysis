#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
group_stats_partial_corr_by_hemi_task-constrained.py
------------------------------------------------------------------------------------------
Goal:
    Compute group-level Fisher-z statistics from per-subject partial-correlation
    produced by nilearn_partial_corr_seed-task_by_macror-task_by_hemi.py

    Operates on per-subject bilateral .npy files with shape
    (n_clusters × 2*n_clusters), where:
        columns  0 : n_clusters   → ipsilateral targets  (same hemi as seed)
        columns n_clusters : end  → contralateral targets

    REPORTING MATRIX FORMAT  (5 seeds × 10 targets):
    Both seed rows and target columns are restricted to the 5 core eye-field
    ROIs (mPCS, sPCS, iPCS, sIPS, iIPS), matching the full-corr output format
    exactly.  Target columns are split into ipsilateral and contralateral halves
    with self-documenting labels that require no downstream slicing:
        mPCS_ipsi, sPCS_ipsi, iPCS_ipsi, sIPS_ipsi, iIPS_ipsi,
        mPCS_contra, sPCS_contra, iPCS_contra, sIPS_contra, iIPS_contra

    For each hemisphere × variant the script:
        1. Resolves which .npy file to load per subject (run variant logic,
           including concat_clean fallback for RUN02_EXCLUDED subjects).
        2. Slices ipsi and contra eye-field columns from the bilateral matrix,
           assembles a (5 × 10) Fisher-z matrix per subject.
        3. Stacks matrices across subjects → (n_subjects × 5 × 10).
        4. Computes group mean and median in Fisher-z space.
        5. Back-converts to Pearson r via tanh() at the reporting stage only.
        6. Saves .npy + .csv for both spaces (fisherz and r) per hemisphere ×
           variant, plus a compressed .npz archive with full metadata.
        7. For concat_clean only: saves two long-format reporting TSVs
           (one for ipsi targets, one for contra targets) with one row per
           subject × seed, values in Pearson r, plus a GROUP row at the
           bottom whose values come from tanh(nanmedian(Fisher-z)) — the
           same group median already computed for the .npy outputs.

    Averaging is always in Fisher-z space; Pearson r is recovered only at the
    final reporting stage via tanh().  This is intentional:
      - Fisher-z has approximately constant variance ~1/(n-3), making it the
        correct space for averaging and any subsequent parametric tests.
      - Raw Pearson r values have r-dependent variance and must never be
        averaged directly.
------------------------------------------------------------------------------------------
Output filename convention (harmonized with full-corr group stats):

    seed-task_by_macror-task_partial-corr_{space}_{stat}_{run_label}_{hemi}.npy / .csv
    seed-task_by_macror-task_partial-corr_{run_label}_{hemi}.npz

    Reporting TSVs (concat_clean only):
    seed-task_by_macror-task_partial-corr_r_report_ipsi_{hemi}.tsv
    seed-task_by_macror-task_partial-corr_r_report_contra_{hemi}.tsv

    Reporting TSV structure:
        subject   : subject ID (e.g. "sub-01") or "GROUP"
        seed      : eye-field seed name (mPCS / sPCS / iPCS / sIPS / iIPS)
        mPCS / sPCS / iPCS / sIPS / iIPS : Pearson r for each target region
            (ipsi half in the ipsi table, contra half in the contra table)
        GROUP row values = tanh(nanmedian(Fisher-z)) across subjects — the
            same quantity saved in the _r_median_ .npy / .csv outputs.

    Full-corr equivalent for reference:
    seed-task_by_macror-task_full-corr_fisherz_median_{run_label}_{hemi}_legacy.npy
------------------------------------------------------------------------------------------
Run variants:
    concat       — concatenated-run .npy, all subjects
    concat_clean — best available run per subject:
                     · RUN02_EXCLUDED subjects → run-01 .npy
                     · all other subjects      → concatenated-run .npy
    run-01       — run-01 .npy, all subjects
    run-02       — run-02 .npy, all subjects (bad subjects retained intentionally
                   to expose the registration artifact in group plots)
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
    seed-task_by_macror-task_partial-corr_r_p25_{run_label}_{hemi}.npy / .csv
    seed-task_by_macror-task_partial-corr_r_p75_{run_label}_{hemi}.npy / .csv
    seed-task_by_macror-task_partial-corr_{run_label}_{hemi}.npz

    Reporting TSVs (concat_clean only, per hemisphere, written to
    {pp_data}/group/91k/rest/partial_corr/tables/):
    seed-task_by_macror-task_partial-corr_r_report_ipsi_{hemi}.tsv
    seed-task_by_macror-task_partial-corr_r_report_contra_{hemi}.tsv

    Rows    : 5 eye-field seeds in canonical order (mPCS first)
    Columns : 5 eye-field regions × 2 hemispheres (ipsi then contra)  (n = 10)
              mPCS_ipsi, sPCS_ipsi, iPCS_ipsi, sIPS_ipsi, iIPS_ipsi,
              mPCS_contra, sPCS_contra, iPCS_contra, sIPS_contra, iIPS_contra

To run:
    $ cd projects/pRF_analysis/RetinoMaps/rest/stats
    $ python group_stats_partial_corr_by_hemi_task-constrained.py /scratch/mszinte/data RetinoMaps 327 b327
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
from typing import Dict, List, Optional

# Suppress only RuntimeWarnings (e.g. nanmean of all-NaN slice on diagonal);
# keep all other warning categories visible.
warnings.filterwarnings("ignore", category=RuntimeWarning)

# ============================================================
# Personal imports
# ============================================================
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../"))

sys.path.append(os.path.abspath(os.path.join(base_dir, "analysis_code/utils")))
from settings_utils import load_settings

sys.path.append(os.path.abspath(os.path.join(base_dir, "RetinoMaps/rest/utils")))
from rest_utils import RUN02_EXCLUDED, VARIANTS

# ============================================================
# Parse and validate arguments
# ============================================================
USAGE = (
    "Usage: python group_stats_partial_corr_by_hemi_task-constrained.py "
    "<main_dir> <project_dir> <group> <server>"
)

if len(sys.argv) != 5:
    print(f"ERROR: expected 4 arguments, got {len(sys.argv) - 1}.\n{USAGE}")
    sys.exit(1)

main_dir    = sys.argv[1]
project_dir = sys.argv[2]
group       = sys.argv[3]
server      = sys.argv[4]

# Percentile bounds for the subject distribution saved alongside group stats.
PCT_LO, PCT_HI = 25.0, 75.0

print("=" * 80)
print("GROUP PARTIAL CORRELATION (TASK-CONSTRAINED) — Fisher-z statistics")
print("Eye-field seeds × eye-field targets (5 × 10)")
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
subjects          = analysis_info["subjects"]  # type: List[str]

# ============================================================
# Macro-regions — full canonical list (mPCS first)
# ============================================================
macro_regions = list(analysis_info["rois-drawn"])
macro_regions.reverse()   # mPCS first

N_MACRO = len(macro_regions)   # expected: 12

# ============================================================
# Eye-field regions — the 5 core ROIs used for all outputs
# ============================================================
N_EYE_FIELDS = 5
EYE_FIELDS   = macro_regions[:N_EYE_FIELDS]  # type: List[str]

assert EYE_FIELDS == ["mPCS", "sPCS", "iPCS", "sIPS", "iIPS"], (
    f"EYE_FIELDS resolved to {EYE_FIELDS}, expected the 5 eye-field ROIs "
    "in mPCS-first order. Check rois-drawn ordering in settings.yml."
)

# Row/column indices of EYE_FIELDS within the macro_regions list
EYE_FIELDS_IDX = [macro_regions.index(r) for r in EYE_FIELDS]

# Target column labels: [5 ipsi | 5 contra], self-documenting for co-authors
TARGET_COLUMNS = (
    [f"{r}_ipsi"   for r in EYE_FIELDS] +
    [f"{r}_contra" for r in EYE_FIELDS]
)  # type: List[str]

print(f"\n  All macro-regions (n={N_MACRO}): {macro_regions}")
print(f"  Eye-field regions (n={N_EYE_FIELDS}): {EYE_FIELDS}")
print(f"  Output shape : ({N_EYE_FIELDS} seeds × {2*N_EYE_FIELDS} targets)")
print(f"  Target cols  : {TARGET_COLUMNS}")

# ============================================================
# Paths
# ============================================================
main_data      = Path(main_dir) / project_dir / "derivatives/pp_data"
output_folder  = main_data / "group/91k/rest/partial_corr/by_hemi/task-constrained"
tables_folder  = main_data / "group/91k/rest/partial_corr/tables"

output_folder.mkdir(parents=True, exist_ok=True)
tables_folder.mkdir(parents=True, exist_ok=True)

# ============================================================
# Helper: resolve per-subject bilateral .npy path
# ============================================================
def npy_path(subject: str, hemi: str, run_tag: Optional[str]) -> Path:
    # Actual BIDS-style filename produced by the subject-level script:
    #   sub-05_task-rest_run-01_space-fsLR_den-91k_desc-fisher-z_lh_task-constrained_bilateral.npy  (per-run)
    #   sub-05_task-rest_space-fsLR_den-91k_desc-fisher-z_lh_task-constrained_bilateral.npy          (concat)
    run_entity = f"_{run_tag}" if run_tag is not None else ""
    fname = (
        f"{subject}_task-rest{run_entity}_space-fsLR_den-91k"
        f"_desc-fisher-z_{hemi}_task-constrained_bilateral.npy"
    )
    return (
        main_data / subject
        / "91k/rest/corr/partial_corr/by_hemi/task-constrained"
        / fname
    )

# ============================================================
# Output filename stem builder
# ============================================================
def _stem(stat: str, space: str, run_label: str, hemi: str) -> str:
    return (
        f"seed-task_by_macror-task_partial-corr"
        f"_{space}_{stat}_{run_label}_{hemi}"
    )

# ============================================================
# Reporting TSV builder (concat_clean only)
#
# Produces two long-format tables — one for ipsi targets, one for contra —
# with the structure:
#   subject | seed | mPCS | sPCS | iPCS | sIPS | iIPS
#
# Subject rows: raw Pearson r = tanh applied per-subject to their Fisher-z
#   slice, never averaged across subjects.
# GROUP row: tanh(nanmedian(Fisher-z)) = median_r passed in — identical to
#   the value saved in the _r_median_ .npy / .csv outputs.
# ============================================================
def _save_reporting_tsvs(
    stacked_fz: np.ndarray,
    median_r:   np.ndarray,
    subject_ids: List[str],
    hemi: str,
) -> None:
    ipsi_cols   = list(range(N_EYE_FIELDS))
    contra_cols = list(range(N_EYE_FIELDS, 2 * N_EYE_FIELDS))

    for side, col_idx in (("ipsi", ipsi_cols), ("contra", contra_cols)):
        rows = []  # type: List[Dict]

        # ── Per-subject rows ──────────────────────────────────────────────
        for s_idx, subj in enumerate(subject_ids):
            subj_r = np.tanh(stacked_fz[s_idx])   # (N_EYE_FIELDS × 10)
            for seed_idx, seed in enumerate(EYE_FIELDS):
                row = {"subject": subj, "seed": seed}
                for t_idx, t_col in enumerate(col_idx):
                    row[EYE_FIELDS[t_idx]] = subj_r[seed_idx, t_col]
                rows.append(row)

        # ── GROUP row ─────────────────────────────────────────────────────
        for seed_idx, seed in enumerate(EYE_FIELDS):
            row = {"subject": "GROUP", "seed": seed}
            for t_idx, t_col in enumerate(col_idx):
                row[EYE_FIELDS[t_idx]] = median_r[seed_idx, t_col]
            rows.append(row)

        col_order = ["subject", "seed"] + EYE_FIELDS
        df = pd.DataFrame(rows, columns=col_order)

        fname = (
            f"seed-task_by_macror-task_partial-corr"
            f"_r_report_{side}_{hemi}.tsv"
        )
        df.to_csv(tables_folder / fname, sep="\t", index=False,
                  float_format="%.4f")
        print(f"    Saved reporting TSV: {fname}")

# ============================================================
# Main loop — hemisphere × variant
# ============================================================

for hemi in ("lh", "rh"):
    print(f"\n{'='*80}")
    print(f"Processing hemisphere: {hemi.upper()}")
    print("=" * 80)

    for variant, (normal_tag, excluded_tag, _skip) in VARIANTS.items():
        print(f"\n  --- Variant: {variant} ---")

        subject_matrices = []  # type: List[np.ndarray]
        subject_ids      = []  # type: List[str]
        missing_subjects = []  # type: List[str]

        for subject in subjects:
            is_excluded = subject in RUN02_EXCLUDED
            run_tag     = excluded_tag if is_excluded else normal_tag

            fpath = npy_path(subject, hemi, run_tag)

            if not fpath.exists():
                print(f"    WARNING [{subject}]: missing {fpath.name} — skipped")
                missing_subjects.append(subject)
                continue

            bilateral = np.load(fpath)

            if bilateral.shape != (N_MACRO, 2 * N_MACRO):
                raise ValueError(
                    f"[{subject} {hemi} {variant}] Unexpected shape "
                    f"{bilateral.shape} in {fpath.name} "
                    f"(expected ({N_MACRO}, {2 * N_MACRO}))."
                )

            ipsi_half   = bilateral[:, :N_MACRO]
            contra_half = bilateral[:, N_MACRO:]

            ipsi_ef   = ipsi_half[np.ix_(EYE_FIELDS_IDX, EYE_FIELDS_IDX)]    # (5×5)
            contra_ef = contra_half[np.ix_(EYE_FIELDS_IDX, EYE_FIELDS_IDX)]  # (5×5)

            mat = np.concatenate([ipsi_ef, contra_ef], axis=1)  # (5×10)

            if mat.shape != (N_EYE_FIELDS, 2 * N_EYE_FIELDS):
                raise ValueError(
                    f"[{subject} {hemi} {variant}] Unexpected assembled shape "
                    f"{mat.shape}, expected ({N_EYE_FIELDS}, {2 * N_EYE_FIELDS})."
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
        if missing_subjects:
            print(f"    Missing       : {missing_subjects}")

        stacked_fz = np.stack(subject_matrices, axis=0)

        if stacked_fz.shape != (n_valid, N_EYE_FIELDS, 2 * N_EYE_FIELDS):
            raise ValueError(
                f"Unexpected stack shape {stacked_fz.shape} for "
                f"{hemi} / {variant}."
            )

        # ── Group statistics in Fisher-z space ───────────────────────────
        mean_fz   = np.nanmean(  stacked_fz, axis=0)
        median_fz = np.nanmedian(stacked_fz, axis=0)
        pct_lo_fz = np.nanpercentile(stacked_fz, PCT_LO, axis=0)
        pct_hi_fz = np.nanpercentile(stacked_fz, PCT_HI, axis=0)

        # ── Back-convert to Pearson r at reporting stage only ─────────────
        mean_r   = np.tanh(mean_fz)
        median_r = np.tanh(median_fz)
        pct_lo_r = np.tanh(pct_lo_fz)
        pct_hi_r = np.tanh(pct_hi_fz)

        print(f"    Fisher-z mean   range : [{np.nanmin(mean_fz):.4f},   {np.nanmax(mean_fz):.4f}]")
        print(f"    Fisher-z median range : [{np.nanmin(median_fz):.4f}, {np.nanmax(median_fz):.4f}]")

        run_label  = normal_tag if normal_tag is not None else variant
        pct_lo_tag = f"p{int(PCT_LO):02d}"
        pct_hi_tag = f"p{int(PCT_HI):02d}"

        # ── Save Fisher-z and r arrays (.npy + .csv) ─────────────────────
        for space, arrays in (
            ("fisherz", (("mean",     mean_fz),
                         ("median",   median_fz),
                         (pct_lo_tag, pct_lo_fz),
                         (pct_hi_tag, pct_hi_fz))),
            ("r",       (("mean",     mean_r),
                         ("median",   median_r),
                         (pct_lo_tag, pct_lo_r),
                         (pct_hi_tag, pct_hi_r))),
        ):
            for stat, arr in arrays:
                stem = _stem(stat, space, run_label, hemi)
                np.save(output_folder / f"{stem}.npy", arr)
                pd.DataFrame(
                    arr, index=EYE_FIELDS, columns=TARGET_COLUMNS
                ).to_csv(
                    output_folder / f"{stem}.csv", float_format="%.4f"
                )
                print(f"    Saved: {stem}.npy / .csv")

        # ── Compressed archive with all arrays + metadata ─────────────────
        npz_stem = f"seed-task_by_macror-task_partial-corr_{run_label}_{hemi}"
        np.savez_compressed(
            output_folder / f"{npz_stem}.npz",
            mean_fz           = mean_fz,
            median_fz         = median_fz,
            pct_lo_fz         = pct_lo_fz,
            pct_hi_fz         = pct_hi_fz,
            mean_r            = mean_r,
            median_r          = median_r,
            pct_lo_r          = pct_lo_r,
            pct_hi_r          = pct_hi_r,
            n_subjects_loaded = np.array(n_valid),
            subjects_loaded   = np.array(subject_ids),
            subjects_missing  = np.array(missing_subjects),
            eye_fields        = np.array(EYE_FIELDS),
            target_columns    = np.array(TARGET_COLUMNS),
            hemi              = np.array(hemi),
            variant           = np.array(variant),
        )
        print(f"    Saved: {npz_stem}.npz")

        # ── Reporting TSVs — concat_clean only ───────────────────────────
        if variant == "concat_clean":
            _save_reporting_tsvs(stacked_fz, median_r, subject_ids, hemi)

print("\n" + "=" * 80)
print("ALL HEMISPHERES × VARIANTS COMPLETE")
print("=" * 80)
print(f"\nStats outputs : {output_folder}")
print(f"Reporting TSVs: {tables_folder}")
print(
    "\nNote: Fisher-z outputs are in z-space. Apply np.tanh() to recover Pearson r "
    "only at the final reporting or plotting stage."
)