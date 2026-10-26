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

    TWO REPORTING SCOPES are produced (EXTENDED — mirrors the equivalent
    extension in group_stats_full_corr_by_hemi_task-constrained.py):

    1. EYE-FIELD SCOPE  (5 seeds × 10 targets) — original/primary scope,
       filenames UNCHANGED from before this extension so nothing downstream
       (violin plots, co-author's work) breaks.
       Both seed rows and target columns are restricted to the 5 core
       eye-field ROIs (mPCS, sPCS, iPCS, sIPS, iIPS), matching the full-corr
       output format exactly. Target columns split ipsi/contra:
           mPCS_ipsi, sPCS_ipsi, iPCS_ipsi, sIPS_ipsi, iIPS_ipsi,
           mPCS_contra, sPCS_contra, iPCS_contra, sIPS_contra, iIPS_contra

    2. ALL-MACRO SCOPE  (12 seeds × 24 targets) — NEW, added for a
       supplementary figure. Seed rows and target columns span ALL 12
       macro-regions. Filenames carry an extra "_all-macro" token so they
       never collide with the eye-field-scope files.

    Both scopes are slices of the SAME per-subject bilateral matrix: the
    subject-level .npy already contains the full (12 × 24) Fisher-z matrix
    (all macro-region seeds × all macro-region ipsi+contra targets) — the
    eye-field scope was always just a (5 × 10) slice of it, with the rest
    of the matrix loaded and then discarded. This extension simply keeps
    the full matrix instead of discarding it, so eye-field-scope numbers
    are unchanged (same slice, same source array) and no extra file I/O is
    needed for the all-macro scope.

    For each hemisphere × variant the script:
        1. Resolves which .npy file to load per subject (run variant logic,
           including concat_clean fallback for RUN02_EXCLUDED subjects).
        2. Loads the full (12 × 24) bilateral Fisher-z matrix per subject.
        3. Slices out the eye-field (5 × 10) sub-matrix for the original scope.
        4. Stacks matrices across subjects for BOTH scopes independently
           → (n_subjects × 5 × 10) and (n_subjects × 12 × 24).
        5. Computes group mean/median/p25/p75 in Fisher-z space, per scope.
        6. Back-converts to Pearson r via tanh() at the reporting stage only.
        7. Saves .npy + .csv (both spaces) + a compressed .npz archive per
           hemisphere × variant, EYE-FIELD SCOPE ONLY (as before — the
           all-macro scope only produces the reporting TSVs described next;
           ask if the full stat-array set is also wanted for that scope).
        8. For concat_clean only: saves TWO PAIRS of long-format reporting
           TSVs (ipsi/contra) — one pair for eye-field scope (unchanged
           filenames), one pair for all-macro scope (new). Each row is one
           subject × seed, values in Pearson r, plus THREE summary rows per
           seed: GROUP (median), GROUP_p25, GROUP_p75 — all three computed
           as tanh(percentile(Fisher-z)) across subjects, never
           percentile(r). Added so each reporting TSV is self-contained for
           a median + IQR heatmap without needing the separate stat-array
           files.

    Averaging is always in Fisher-z space; Pearson r is recovered only at the
    final reporting stage via tanh().  This is intentional:
      - Fisher-z has approximately constant variance ~1/(n-3), making it the
        correct space for averaging and any subsequent parametric tests.
      - Raw Pearson r values have r-dependent variance and must never be
        averaged directly.
------------------------------------------------------------------------------------------
Output filename convention (harmonized with full-corr group stats):

    Eye-field scope (UNCHANGED filenames):
    seed-task_by_macror-task_partial-corr_{space}_{stat}_{run_label}_{hemi}_{estimator}.npy / .csv
    seed-task_by_macror-task_partial-corr_{run_label}_{hemi}_{estimator}.npz
    seed-task_by_macror-task_partial-corr_r_report_{side}_{hemi}_{estimator}.tsv

    All-macro scope (NEW — reporting TSVs only):
    seed-task_by_macror-task_partial-corr_r_report_{side}_all-macro_{hemi}_{estimator}.tsv

    Reporting TSV structure (both scopes):
        subject   : subject ID (e.g. "sub-01"), or one of three summary
                    labels: "GROUP", "GROUP_p25", "GROUP_p75"
        seed      : seed macro-region name (5 eye-fields, or all 12 for
                    all-macro scope)
        <region columns> : Pearson r for each target region (ipsi half in
            the ipsi table, contra half in the contra table) — 5 columns
            for eye-field scope, 12 for all-macro scope
        GROUP / GROUP_p25 / GROUP_p75 row values = tanh(nanmedian /
            nanpercentile(25) / nanpercentile(75) of Fisher-z) across
            subjects — for eye-field scope these are the same quantities
            saved in the _r_median_ / _r_p25_ / _r_p75_ .npy / .csv outputs.

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
    5: covariance estimator — one of "raw", "ledoit-wolf", "graphical-lasso"
       (optional; default: ledoit-wolf). Must match the ESTIMATOR_TAG used
       at the subject level.

Outputs (per hemisphere × variant):
    Eye-field scope stat arrays (unchanged):
    seed-task_by_macror-task_partial-corr_fisherz_mean_{run_label}_{hemi}_{estimator}.npy / .csv
    seed-task_by_macror-task_partial-corr_fisherz_median_{run_label}_{hemi}_{estimator}.npy / .csv
    seed-task_by_macror-task_partial-corr_r_mean_{run_label}_{hemi}_{estimator}.npy / .csv
    seed-task_by_macror-task_partial-corr_r_median_{run_label}_{hemi}_{estimator}.npy / .csv
    seed-task_by_macror-task_partial-corr_r_p25_{run_label}_{hemi}_{estimator}.npy / .csv
    seed-task_by_macror-task_partial-corr_r_p75_{run_label}_{hemi}_{estimator}.npy / .csv
    seed-task_by_macror-task_partial-corr_{run_label}_{hemi}_{estimator}.npz

    Reporting TSVs (concat_clean only, per hemisphere, written to
    {pp_data}/group/91k/rest/partial_corr/tables/):
    Eye-field scope   : seed-task_by_macror-task_partial-corr_r_report_ipsi_{hemi}_{estimator}.tsv
                        seed-task_by_macror-task_partial-corr_r_report_contra_{hemi}_{estimator}.tsv
    All-macro scope   : seed-task_by_macror-task_partial-corr_r_report_ipsi_all-macro_{hemi}_{estimator}.tsv
                        seed-task_by_macror-task_partial-corr_r_report_contra_all-macro_{hemi}_{estimator}.tsv

    Rows    : 5 eye-field seeds (or 12 macro-regions for all-macro scope),
              canonical order (mPCS first)
    Columns : 5 (or 12) regions × 2 hemispheres (ipsi then contra)

To run:
    $ cd projects/pRF_analysis/RetinoMaps/rest/stats
    $ python group_stats_partial_corr_by_hemi_task-constrained.py /scratch/mszinte/data RetinoMaps 327 b327 ledoit-wolf
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
VALID_ESTIMATORS = ("raw", "ledoit-wolf", "graphical-lasso")

USAGE = (
    "Usage: python group_stats_partial_corr_by_hemi_task-constrained.py "
    "<main_dir> <project_dir> <group> <server> [estimator]\n"
    f"  estimator: one of {', '.join(VALID_ESTIMATORS)} (default: ledoit-wolf)\n"
    "  Must match the ESTIMATOR_TAG used when running the subject-level\n"
    "  nilearn_partial_corr_task_constrained.py script."
)

if len(sys.argv) not in (5, 6):
    print(f"ERROR: expected 4 or 5 arguments, got {len(sys.argv) - 1}.\n{USAGE}")
    sys.exit(1)

main_dir    = sys.argv[1]
project_dir = sys.argv[2]
group       = sys.argv[3]
server      = sys.argv[4]

ESTIMATOR_TAG = sys.argv[5] if len(sys.argv) == 6 else "ledoit-wolf"
if ESTIMATOR_TAG not in VALID_ESTIMATORS:
    print(f"ERROR: unrecognised estimator '{ESTIMATOR_TAG}'.\n{USAGE}")
    sys.exit(1)

# Percentile bounds for the subject distribution saved alongside group stats.
PCT_LO, PCT_HI = 25.0, 75.0

print("=" * 80)
print("GROUP PARTIAL CORRELATION (TASK-CONSTRAINED) — Fisher-z statistics")
print(f"Eye-field scope (5x10) + all-macro scope (12x24) — estimator: {ESTIMATOR_TAG}")
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
# Eye-field regions — the 5 core ROIs used for the original scope
# ============================================================
N_EYE_FIELDS = 5
EYE_FIELDS   = macro_regions[:N_EYE_FIELDS]  # type: List[str]

assert EYE_FIELDS == ["mPCS", "sPCS", "iPCS", "sIPS", "iIPS"], (
    f"EYE_FIELDS resolved to {EYE_FIELDS}, expected the 5 eye-field ROIs "
    "in mPCS-first order. Check rois-drawn ordering in settings.yml."
)

# Row/column indices of EYE_FIELDS within the macro_regions list (== [0..4])
EYE_FIELDS_IDX = [macro_regions.index(r) for r in EYE_FIELDS]

# Target column labels — eye-field scope (5 seeds x 10 targets)
TARGET_COLUMNS = (
    [f"{r}_ipsi"   for r in EYE_FIELDS] +
    [f"{r}_contra" for r in EYE_FIELDS]
)  # type: List[str]

print(f"\n  All macro-regions (n={N_MACRO}): {macro_regions}")
print(f"  Eye-field regions (n={N_EYE_FIELDS}): {EYE_FIELDS}")
print(f"  Eye-field output shape : ({N_EYE_FIELDS} seeds × {2*N_EYE_FIELDS} targets)")
print(f"  All-macro  output shape : ({N_MACRO} seeds × {2*N_MACRO} targets)")

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
    #   sub-05_task-rest_run-01_space-fsLR_den-91k_desc-fisher-z_lh_task-constrained_ledoit-wolf_bilateral.npy  (per-run)
    #   sub-05_task-rest_space-fsLR_den-91k_desc-fisher-z_lh_task-constrained_ledoit-wolf_bilateral.npy          (concat)
    # ESTIMATOR_TAG must match what was used at the subject level (see USAGE).
    run_entity = f"_{run_tag}" if run_tag is not None else ""
    fname = (
        f"{subject}_task-rest{run_entity}_space-fsLR_den-91k"
        f"_desc-fisher-z_{hemi}_task-constrained_{ESTIMATOR_TAG}_bilateral.npy"
    )
    return (
        main_data / subject
        / "91k/rest/corr/partial_corr/by_hemi/task-constrained"
        / fname
    )

# ============================================================
# Filename stem builders — pure functions, no Path objects, no I/O.
# Directory joining is done inline at each call site instead of being
# wrapped in a function.
# ============================================================

def _stem(stat: str, space: str, run_label: str, hemi: str) -> str:
    """Eye-field scope stat-array filename stem (unchanged from before)."""
    return (
        f"seed-task_by_macror-task_partial-corr"
        f"_{space}_{stat}_{run_label}_{hemi}_{ESTIMATOR_TAG}"
    )


def _report_stem(side: str, hemi: str, scope_token: str = "") -> str:
    """
    Reporting-TSV filename stem.

    scope_token: "" for the eye-field scope (matches the original filename
    exactly); "_all-macro" for the new all-macro scope.
    """
    return (
        f"seed-task_by_macror-task_partial-corr"
        f"_r_report_{side}{scope_token}_{hemi}_{ESTIMATOR_TAG}"
    )


# ============================================================
# Reporting TSV builder — generalized over scope (eye-field or all-macro)
#
# Produces two long-format tables — one for ipsi targets, one for contra —
# with the structure:
#   subject | seed | <region columns...>
#
# Subject rows: raw Pearson r = tanh applied per-subject to their Fisher-z
#   slice, never averaged across subjects. Three summary rows per seed:
#   GROUP / GROUP_p25 / GROUP_p75 = tanh(nanmedian / nanpercentile(25/75)
#   of Fisher-z) across subjects — never percentile(r).
#
# Parameters
# ----------
# stacked_fz    : (n_subjects x n_seeds x n_targets) Fisher-z array
# median_r      : (n_seeds x n_targets) group median in Pearson r space
# pct_lo_r      : (n_seeds x n_targets) group 25th percentile in Pearson r
# pct_hi_r      : (n_seeds x n_targets) group 75th percentile in Pearson r
# subject_ids   : list of subject ID strings, length n_subjects
# hemi          : "lh" or "rh"
# seeds         : seed names in row order (EYE_FIELDS or macro_regions)
# region_cols   : region names used as BOTH the per-side column set and the
#                 seed set — n_targets == 2 * len(region_cols)
# scope_token   : "" or "_all-macro", passed through to _report_stem()
# ============================================================
def _save_reporting_tsvs(
    stacked_fz: np.ndarray,
    median_r:   np.ndarray,
    pct_lo_r:   np.ndarray,
    pct_hi_r:   np.ndarray,
    subject_ids: List[str],
    hemi: str,
    seeds: List[str],
    region_cols: List[str],
    scope_token: str = "",
) -> None:
    n_regions = len(region_cols)
    ipsi_cols   = list(range(n_regions))
    contra_cols = list(range(n_regions, 2 * n_regions))

    for side, col_idx in (("ipsi", ipsi_cols), ("contra", contra_cols)):
        rows = []  # type: List[Dict]

        # ── Per-subject rows ──────────────────────────────────────────────
        for s_idx, subj in enumerate(subject_ids):
            subj_r = np.tanh(stacked_fz[s_idx])   # (n_seeds x n_targets)
            for seed_idx, seed in enumerate(seeds):
                row = {"subject": subj, "seed": seed}
                for t_idx, t_col in enumerate(col_idx):
                    row[region_cols[t_idx]] = subj_r[seed_idx, t_col]
                rows.append(row)

        # ── GROUP summary rows ────────────────────────────────────────────
        for label, arr in (("GROUP", median_r),
                            ("GROUP_p25", pct_lo_r),
                            ("GROUP_p75", pct_hi_r)):
            for seed_idx, seed in enumerate(seeds):
                row = {"subject": label, "seed": seed}
                for t_idx, t_col in enumerate(col_idx):
                    row[region_cols[t_idx]] = arr[seed_idx, t_col]
                rows.append(row)

        col_order = ["subject", "seed"] + region_cols
        df = pd.DataFrame(rows, columns=col_order)

        fname = _report_stem(side, hemi, scope_token) + ".tsv"
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

        # Per-subject full (12 x 24) bilateral matrices — the single
        # source both scopes are derived from.
        full_matrices    = []  # type: List[np.ndarray]
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

            if variant == "concat_clean" and is_excluded:
                print(f"    {subject}: OK (fallback → run-01)")
            else:
                print(f"    {subject}: OK")

            # Keep the FULL (12 x 24) matrix — no eye-field restriction here.
            # This is the same array that was previously sliced down
            # immediately and the rest discarded; nothing about how it's
            # loaded has changed.
            full_matrices.append(bilateral)
            subject_ids.append(subject)

        if not full_matrices:
            print(f"    ERROR: no valid subjects for {hemi} / {variant} — skipping.")
            continue

        n_valid = len(full_matrices)
        print(f"\n    Valid subjects: {n_valid}/{len(subjects)}")
        if missing_subjects:
            print(f"    Missing       : {missing_subjects}")

        # Stack -> (n_subjects x 12 x 24) in Fisher-z space — all-macro scope
        stacked_fz_allmacro = np.stack(full_matrices, axis=0)
        if stacked_fz_allmacro.shape != (n_valid, N_MACRO, 2 * N_MACRO):
            raise ValueError(
                f"Unexpected all-macro stack shape {stacked_fz_allmacro.shape} "
                f"for {hemi} / {variant}."
            )

        # Slice out eye-field scope: rows = EYE_FIELDS_IDX (ipsi+contra
        # targets selected via np.ix_ semantics, applied per-half exactly
        # as before) -> (n_subjects x 5 x 10), bit-for-bit identical to the
        # pre-extension slicing (same source array, same indices).
        ipsi_idx_full   = EYE_FIELDS_IDX
        contra_idx_full = [N_MACRO + i for i in EYE_FIELDS_IDX]
        eye_field_col_idx = ipsi_idx_full + contra_idx_full

        stacked_fz = (
            stacked_fz_allmacro[:, EYE_FIELDS_IDX, :][:, :, eye_field_col_idx]
        )
        if stacked_fz.shape != (n_valid, N_EYE_FIELDS, 2 * N_EYE_FIELDS):
            raise ValueError(
                f"Unexpected eye-field stack shape {stacked_fz.shape} for "
                f"{hemi} / {variant}."
            )

        # ── Group statistics in Fisher-z space — EYE-FIELD SCOPE ─────────
        # (unchanged from before the extension: stat arrays are only saved
        # for this scope, matching the original script's outputs exactly)
        mean_fz   = np.nanmean(  stacked_fz, axis=0)
        median_fz = np.nanmedian(stacked_fz, axis=0)
        pct_lo_fz = np.nanpercentile(stacked_fz, PCT_LO, axis=0)
        pct_hi_fz = np.nanpercentile(stacked_fz, PCT_HI, axis=0)

        mean_r   = np.tanh(mean_fz)
        median_r = np.tanh(median_fz)
        pct_lo_r = np.tanh(pct_lo_fz)
        pct_hi_r = np.tanh(pct_hi_fz)

        print(f"    [eye-field] Fisher-z mean   range : [{np.nanmin(mean_fz):.4f},   {np.nanmax(mean_fz):.4f}]")
        print(f"    [eye-field] Fisher-z median range : [{np.nanmin(median_fz):.4f}, {np.nanmax(median_fz):.4f}]")

        # ── Group statistics in Fisher-z space — ALL-MACRO SCOPE ─────────
        # median + p25/p75 needed for the reporting TSV's GROUP,
        # GROUP_p25, GROUP_p75 rows — not saved as standalone .npy/.csv
        # stat-array files; ask if those are also wanted for this scope.
        median_fz_allmacro = np.nanmedian(stacked_fz_allmacro, axis=0)
        pct_lo_fz_allmacro  = np.nanpercentile(stacked_fz_allmacro, PCT_LO, axis=0)
        pct_hi_fz_allmacro  = np.nanpercentile(stacked_fz_allmacro, PCT_HI, axis=0)

        median_r_allmacro  = np.tanh(median_fz_allmacro)
        pct_lo_r_allmacro  = np.tanh(pct_lo_fz_allmacro)
        pct_hi_r_allmacro  = np.tanh(pct_hi_fz_allmacro)

        print(f"    [all-macro] Fisher-z median range : [{np.nanmin(median_fz_allmacro):.4f}, {np.nanmax(median_fz_allmacro):.4f}]")

        run_label  = normal_tag if normal_tag is not None else variant
        pct_lo_tag = f"p{int(PCT_LO):02d}"
        pct_hi_tag = f"p{int(PCT_HI):02d}"

        # ── Save Fisher-z and r arrays (.npy + .csv) — eye-field scope ───
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

        # ── Compressed archive with all arrays + metadata — eye-field scope
        npz_stem = (
            f"seed-task_by_macror-task_partial-corr"
            f"_{run_label}_{hemi}_{ESTIMATOR_TAG}"
        )
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
            cov_estimator     = np.array(ESTIMATOR_TAG),
        )
        print(f"    Saved: {npz_stem}.npz")

        # ── Reporting TSVs — concat_clean only, BOTH scopes ──────────────
        if variant == "concat_clean":
            _save_reporting_tsvs(
                stacked_fz, median_r, pct_lo_r, pct_hi_r, subject_ids, hemi,
                seeds=EYE_FIELDS, region_cols=list(EYE_FIELDS),
                scope_token="",
            )
            _save_reporting_tsvs(
                stacked_fz_allmacro, median_r_allmacro,
                pct_lo_r_allmacro, pct_hi_r_allmacro, subject_ids, hemi,
                seeds=macro_regions, region_cols=list(macro_regions),
                scope_token="_all-macro",
            )

print("\n" + "=" * 80)
print("ALL HEMISPHERES × VARIANTS COMPLETE")
print("=" * 80)
print(f"\nStats outputs (eye-field scope) : {output_folder}")
print(f"Reporting TSVs (both scopes)     : {tables_folder}")
print(
    "\nNote: Fisher-z outputs are in z-space. Apply np.tanh() to recover Pearson r "
    "only at the final reporting or plotting stage."
)