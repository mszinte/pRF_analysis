#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
group_stats_full_corr_by_hemi_task-constrained.py
------------------------------------------------------------------------------------------
Goal:

    Compute group-level Fisher-z statistics from per-subject, per-seed,
    task-constrained full-correlation TSV files parcellated by macro-region.

    "Task-constrained" means both the seed and the target regions come from
    task (pRF/saccade) results rather than from the full Glasser MMP atlas.
    The counterpart script (group_full_corr_by_hemi.py) handles the
    "task-free" case (53 MMP atlas parcels per hemisphere).

    Each input TSV has:
        - 24 rows : macro-regions in rois-drawn canonical order, both
                    hemispheres concatenated in a single file:
                      rows  0–11  → LH macro-regions (mPCS = row 0, V1 = row 11)
                      rows 12–23  → RH macro-regions (mPCS = row 12, V1 = row 23)
        - 1 column : Fisher-z correlation of that seed with each macro-region

    Files are always in legacy mode (hardcoded); no mode argument is needed
    since this was the only error-free mode in -cifti-parcellate when
    parcellating by macro-regions (some files have missing vertices).

    TWO REPORTING SCOPES are produced (EXTENDED — see below):

    1. EYE-FIELD SCOPE  (5 seeds × 10 targets) — original/primary scope,
       filenames UNCHANGED from before this extension so nothing downstream
       (violin plots, co-author's work) breaks.
       Both seed rows and target columns are restricted to the 5 core
       eye-field ROIs (mPCS, sPCS, iPCS, sIPS, iIPS). Target columns split
       into ipsi/contra halves:
           mPCS_ipsi, sPCS_ipsi, iPCS_ipsi, sIPS_ipsi, iIPS_ipsi,
           mPCS_contra, sPCS_contra, iPCS_contra, sIPS_contra, iIPS_contra

    2. ALL-MACRO SCOPE  (12 seeds × 24 targets) — NEW, added for a
       supplementary figure. Seed rows and target columns span ALL 12
       macro-regions (not just the 5 eye-fields). Target columns split
       into ipsi/contra halves, same principle as above but full breadth.
       Filenames carry an extra "_all-macro" token so they never collide
       with the eye-field-scope files.

    Both scopes are derived from the SAME per-subject data: each subject's
    12 seed TSVs (all macro-regions, not just eye-fields) are loaded ONCE
    into a (12 × 24) Fisher-z matrix; the eye-field (5 × 10) matrix is a
    slice of that full matrix, not a separately-loaded quantity — this
    guarantees the eye-field-scope numbers are bit-for-bit identical to
    before the extension (loading fewer TSVs than the full set would have
    given identical values; loading the full set and slicing is equivalent
    and avoids reading each seed TSV twice).

    For each hemisphere × variant the script:
        1. Resolves which TSV to load per subject (run variant logic, including
           concat_clean fallback for RUN02_EXCLUDED subjects).
        2. Loads ALL 12 macro-region seed TSVs per subject, assembles a
           (12 × 24) Fisher-z matrix (ipsi + contra targets, all macro-regions).
        3. Slices out the eye-field (5 × 10) sub-matrix for the original scope.
        4. Stacks matrices across subjects for BOTH scopes independently
           → (n_subjects × 5 × 10) and (n_subjects × 12 × 24).
        5. Computes group mean/median/p25/p75 in Fisher-z space, per scope.
        6. Back-converts to Pearson r via tanh() at the reporting stage only.
        7. Saves .npy + .csv (both spaces) + a compressed .npz archive per
           hemisphere × variant, EYE-FIELD SCOPE ONLY (as before — the
           all-macro scope only produces the reporting TSVs described next,
           not the full stat-array set; see note below if that's also needed).
        8. For concat_clean only: saves TWO PAIRS of long-format reporting
           TSVs (ipsi/contra) — one pair for eye-field scope (unchanged
           filenames), one pair for all-macro scope (new). Each has one row
           per subject × seed, values in Pearson r, plus a GROUP row from
           tanh(nanmedian(Fisher-z)) — the same group median used for the
           eye-field scope's saved .npy/.csv outputs.

    Averaging is always in Fisher-z space; Pearson r is recovered only at the
    final reporting stage via tanh().
------------------------------------------------------------------------------------------
Output filename convention (harmonized with partial-corr group stats):

    Eye-field scope (UNCHANGED filenames):
    seed-task_by_macror-task_full-corr_{space}_{stat}_{run_label}_{hemi}_legacy.npy / .csv
    seed-task_by_macror-task_full-corr_{run_label}_{hemi}_legacy.npz
    seed-task_by_macror-task_full-corr_r_report_{side}_{hemi}_legacy.tsv

    All-macro scope (NEW — reporting TSVs only):
    seed-task_by_macror-task_full-corr_r_report_{side}_all-macro_{hemi}_legacy.tsv

    Partial-corr equivalent for reference:
    seed-task_by_macror-task_partial-corr_fisherz_median_{run_label}_{hemi}.npy
------------------------------------------------------------------------------------------
Run variants:
    concat       — concatenated-run TSV, all subjects
    concat_clean — best available run per subject:
                     · RUN02_EXCLUDED subjects → run-01 TSV
                     · all other subjects      → concatenated-run TSV
    run-01       — run-01 TSV, all subjects
    run-02       — run-02 TSV, all subjects (bad subjects retained intentionally
                   to expose the registration artifact in group plots)
------------------------------------------------------------------------------------------
Inputs (sys.argv):
    1: main project directory   (e.g. /scratch/mszinte/data)
    2: project name/directory   (e.g. RetinoMaps)
    3: server group             (e.g. 327)
    4: server project           (e.g. b327)

Outputs (per hemisphere × variant):
    Eye-field scope stat arrays (unchanged):
    seed-task_by_macror-task_full-corr_fisherz_mean_{run_label}_{hemi}_legacy.npy / .csv
    seed-task_by_macror-task_full-corr_fisherz_median_{run_label}_{hemi}_legacy.npy / .csv
    seed-task_by_macror-task_full-corr_r_mean_{run_label}_{hemi}_legacy.npy / .csv
    seed-task_by_macror-task_full-corr_r_median_{run_label}_{hemi}_legacy.npy / .csv
    seed-task_by_macror-task_full-corr_r_p25_{run_label}_{hemi}_legacy.npy / .csv
    seed-task_by_macror-task_full-corr_r_p75_{run_label}_{hemi}_legacy.npy / .csv
    seed-task_by_macror-task_full-corr_{run_label}_{hemi}_legacy.npz

    Reporting TSVs (concat_clean only, per hemisphere):
    Eye-field scope   : seed-task_by_macror-task_full-corr_r_report_ipsi_{hemi}_legacy.tsv
                        seed-task_by_macror-task_full-corr_r_report_contra_{hemi}_legacy.tsv
    All-macro scope   : seed-task_by_macror-task_full-corr_r_report_ipsi_all-macro_{hemi}_legacy.tsv
                        seed-task_by_macror-task_full-corr_r_report_contra_all-macro_{hemi}_legacy.tsv

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
            saved in the _r_median_ / _r_p25_ / _r_p75_ .npy / .csv
            outputs. Added so each reporting TSV is self-contained for a
            median + IQR heatmap without needing the separate stat-array
            files.

Filename example (input TSV, run-01, lh seed hMT+):
    sub-05_task-rest_run-01_space-fsLR_den-91k_desc-fisher-z_lh_hMT+
        _task-constrained_parcellated_by_macro_legacy-mode.tsv

To run:
    $ cd projects/pRF_analysis/RetinoMaps/rest/stats
    $ python group_stats_full_corr_by_hemi_task-constrained.py /scratch/mszinte/data RetinoMaps 327 b327
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
from typing import Dict, List

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
    "Usage: python group_stats_full_corr_by_hemi_task-constrained.py "
    "<main_dir> <project_dir> <group> <server>"
)

if len(sys.argv) != 5:
    print("ERROR: expected 4 arguments, got {0}.\n{1}".format(
        len(sys.argv) - 1, USAGE))
    sys.exit(1)

main_dir    = sys.argv[1]
project_dir = sys.argv[2]
group       = sys.argv[3]
server      = sys.argv[4]

# Hardcoded: input TSVs are always produced in legacy mode.
TSV_SUFFIX = "_legacy-mode"
MODE_LABEL = "legacy"

# Percentile bounds for the subject distribution saved alongside group stats.
PCT_LO, PCT_HI = 25.0, 75.0

print("=" * 80)
print("GROUP FULL CORRELATION (TASK-CONSTRAINED) — Fisher-z statistics")
print("Eye-field scope (5x10) + all-macro scope (12x24), legacy mode")
print("=" * 80)
print("  main_dir    : {0}".format(main_dir))
print("  project_dir : {0}".format(project_dir))
print("  group       : {0}".format(group))
print("  server      : {0}".format(server))

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
#
# The full list defines the TSV row order (12 per hemisphere block) and is
# used both as the ALL-MACRO scope's seed/target set, and to locate
# EYE_FIELDS within it via EYE_FIELDS_IDX for the eye-field scope.
# ============================================================
macro_regions = list(analysis_info["rois-drawn"])
macro_regions.reverse()   # mPCS first

N_MACRO      = len(macro_regions)   # expected: 12
N_ROWS_TOTAL = N_MACRO * 2          # 24 rows per TSV

HEMI_ROW_SLICE = {
    "lh": slice(0,       N_MACRO),
    "rh": slice(N_MACRO, N_ROWS_TOTAL),
}  # type: Dict[str, slice]

# ============================================================
# Eye-field regions — the 5 core ROIs used for the original scope
#
# Derived positionally (first 5 of macro_regions after reversing rois-drawn)
# so it stays correct as long as rois-drawn ordering is maintained.
# The assertion below guards against any silent mismatch, and is also what
# makes EYE_FIELDS_IDX == [0, 1, 2, 3, 4] exactly — relied on below when
# slicing the eye-field sub-matrix out of the full (12 x 24) matrix.
# ============================================================
N_EYE_FIELDS = 5
EYE_FIELDS   = macro_regions[:N_EYE_FIELDS]  # type: List[str]

assert EYE_FIELDS == ["mPCS", "sPCS", "iPCS", "sIPS", "iIPS"], (
    "EYE_FIELDS resolved to {0}, expected the 5 eye-field ROIs "
    "in mPCS-first order. Check rois-drawn ordering in settings.yml.".format(
        EYE_FIELDS)
)

# Row indices of EYE_FIELDS within a 12-row hemisphere block (== [0..4])
EYE_FIELDS_IDX = [macro_regions.index(r) for r in EYE_FIELDS]

# Column labels — eye-field scope (5 seeds x 10 targets)
TARGET_COLUMNS = (
    ["{0}_ipsi".format(r)  for r in EYE_FIELDS] +
    ["{0}_contra".format(r) for r in EYE_FIELDS]
)  # type: List[str]

# Column indices of the eye-field targets within the full (12 x 24) matrix's
# 24 columns, which are laid out as [all-ipsi(12) | all-contra(12)] in
# macro_regions order. Since EYE_FIELDS_IDX == [0..4] (asserted above),
# eye-field ipsi columns are 0-4 and eye-field contra columns are 12-16.
EYE_FIELD_COL_IDX = EYE_FIELDS_IDX + [N_MACRO + i for i in EYE_FIELDS_IDX]

# Column labels — all-macro scope (12 seeds x 24 targets)
ALL_TARGET_COLUMNS = (
    ["{0}_ipsi".format(r)  for r in macro_regions] +
    ["{0}_contra".format(r) for r in macro_regions]
)  # type: List[str]

print("\n  All macro-regions (n={0}): {1}".format(N_MACRO, macro_regions))
print("  Eye-field regions (n={0}): {1}".format(N_EYE_FIELDS, EYE_FIELDS))
print("  TSV layout   : {0} rows — LH rows 0-{1}, RH rows {2}-{3}".format(
    N_ROWS_TOTAL, N_MACRO - 1, N_MACRO, N_ROWS_TOTAL - 1))
print("  Eye-field output shape : ({0} seeds x {1} targets)".format(
    N_EYE_FIELDS, 2 * N_EYE_FIELDS))
print("  All-macro  output shape : ({0} seeds x {1} targets)".format(
    N_MACRO, 2 * N_MACRO))

# ============================================================
# Paths
# ============================================================
main_data     = Path(main_dir) / project_dir / "derivatives/pp_data"
output_folder = main_data / "group/91k/rest/full_corr/by_hemi/task-constrained"
tables_folder = main_data / "group/91k/rest/full_corr/tables"

output_folder.mkdir(parents=True, exist_ok=True)
tables_folder.mkdir(parents=True, exist_ok=True)

# ============================================================
# Filename stem builders — pure functions, no Path objects, no I/O.
# Directory joining is done inline at each call site instead of being
# wrapped in a function, since the directory structure itself is static
# and doesn't need testing the way filename construction does.
# ============================================================

def _stem(stat, space, run_label, hemi):
    # type: (str, str, str, str) -> str
    """Eye-field scope stat-array filename stem (unchanged from before)."""
    return (
        "seed-task_by_macror-task_full-corr"
        "_{space}_{stat}_{run_label}_{hemi}_{mode}".format(
            space=space, stat=stat, run_label=run_label,
            hemi=hemi, mode=MODE_LABEL)
    )


def _report_stem(side, hemi, scope_token=""):
    # type: (str, str, str) -> str
    """
    Reporting-TSV filename stem.

    scope_token: "" for the eye-field scope (matches the original filename
    exactly, so nothing downstream breaks); "_all-macro" for the new
    all-macro scope.
    """
    return (
        "seed-task_by_macror-task_full-corr"
        "_r_report_{side}{scope}_{hemi}_{mode}".format(
            side=side, scope=scope_token, hemi=hemi, mode=MODE_LABEL)
    )


# ============================================================
# Reporting TSV builder — generalized over scope (eye-field or all-macro)
#
# Produces two long-format tables — one for ipsi targets, one for contra —
# with the structure:
#   subject | seed | <region columns...>
#
# Subject rows contain raw Pearson r values (tanh of individual Fisher-z,
# never averaged). Three summary rows are appended per seed:
#   GROUP      = tanh(nanmedian(Fisher-z across subjects))
#   GROUP_p25  = tanh(nanpercentile(Fisher-z across subjects, 25))
#   GROUP_p75  = tanh(nanpercentile(Fisher-z across subjects, 75))
#
# These three rows make the table self-contained for a median + IQR
# heatmap without needing to cross-reference the separate .npy/.csv stat
# files — added specifically because the reporting TSVs previously only
# carried the median (as the plain "GROUP" row), which wasn't enough for
# that purpose.
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
#                 seed set (EYE_FIELDS for eye-field scope, macro_regions
#                 for all-macro scope) — n_targets == 2 * len(region_cols)
# scope_token   : "" or "_all-macro", passed through to _report_stem()
# ============================================================
def _save_reporting_tsvs(stacked_fz, median_r, pct_lo_r, pct_hi_r,
                          subject_ids, hemi, seeds, region_cols,
                          scope_token=""):
    # type: (np.ndarray, np.ndarray, np.ndarray, np.ndarray, List[str], str, List[str], List[str], str) -> None

    n_regions = len(region_cols)
    ipsi_cols   = list(range(n_regions))
    contra_cols = list(range(n_regions, 2 * n_regions))

    for side, col_idx in (("ipsi", ipsi_cols), ("contra", contra_cols)):
        rows = []  # type: List[Dict]

        # ── Per-subject rows ──────────────────────────────────────────────
        # Convert each subject's Fisher-z slice to Pearson r individually
        # (tanh applied per subject, not to an average)
        for s_idx, subj in enumerate(subject_ids):
            subj_fz = stacked_fz[s_idx]          # (n_seeds x n_targets)
            subj_r  = np.tanh(subj_fz)           # Pearson r, same shape

            for seed_idx, seed in enumerate(seeds):
                row = {"subject": subj, "seed": seed}
                for t_idx, t_col in enumerate(col_idx):
                    row[region_cols[t_idx]] = subj_r[seed_idx, t_col]
                rows.append(row)

        # ── GROUP summary rows ────────────────────────────────────────────
        # median / p25 / p75 all use their respective pre-computed r arrays
        # (each already tanh(percentile(Fisher-z)), never percentile(r))
        for label, arr in (("GROUP", median_r),
                            ("GROUP_p25", pct_lo_r),
                            ("GROUP_p75", pct_hi_r)):
            for seed_idx, seed in enumerate(seeds):
                row = {"subject": label, "seed": seed}
                for t_idx, t_col in enumerate(col_idx):
                    row[region_cols[t_idx]] = arr[seed_idx, t_col]
                rows.append(row)

        # ── Save ─────────────────────────────────────────────────────────
        col_order = ["subject", "seed"] + region_cols
        df = pd.DataFrame(rows, columns=col_order)

        fname = _report_stem(side, hemi, scope_token) + ".tsv"
        df.to_csv(tables_folder / fname, sep="\t", index=False,
                  float_format="%.4f")
        print("    Saved reporting TSV: {0}".format(fname))


# ============================================================
# Main loop — hemisphere x variant
# ============================================================

for hemi in ("lh", "rh"):
    print("\n" + "=" * 80)
    print("Processing hemisphere: {0}".format(hemi.upper()))
    print("=" * 80)

    contra_hemi      = "rh" if hemi == "lh" else "lh"
    row_slice_ipsi   = HEMI_ROW_SLICE[hemi]
    row_slice_contra = HEMI_ROW_SLICE[contra_hemi]

    for variant, (normal_tag, excluded_tag, _skip) in VARIANTS.items():
        # run-02: intentionally retains all subjects for group QC (artifact
        # visibility), regardless of the skip_excluded flag used in WTA scripts.
        print("\n  --- Variant: {0} ---".format(variant))

        # Per-subject full (12 x 24) matrices — the single source both
        # scopes are derived from.
        full_matrices    = []  # type: List[np.ndarray]
        subject_ids      = []  # type: List[str]
        missing_subjects = []  # type: List[str]

        for subject in subjects:
            is_excluded = subject in RUN02_EXCLUDED
            run_tag     = excluded_tag if is_excluded else normal_tag
            run_entity  = "_{0}".format(run_tag) if run_tag is not None else ""

            subj_dir = (
                main_data / subject
                / "91k/rest/corr/full_corr/by_hemi/task-constrained"
            )

            seed_rows     = {}   # type: Dict[str, np.ndarray]
            missing_files = []   # type: List[str]

            # Load ALL 12 macro-region seed TSVs (not just eye-fields) —
            # needed for the all-macro scope, and the eye-field scope is
            # sliced out of the same data rather than loaded separately.
            for seed in macro_regions:
                fname = (
                    "{subject}_task-rest{run}_space-fsLR_den-91k"
                    "_desc-fisher-z_{hemi}_{seed}"
                    "_task-constrained_parcellated_by_macro{suffix}.tsv".format(
                        subject=subject, run=run_entity, hemi=hemi,
                        seed=seed, suffix=TSV_SUFFIX)
                )
                fpath = subj_dir / fname

                if not fpath.exists():
                    missing_files.append(fname)
                    continue

                raw = pd.read_csv(fpath, header=None, sep="\t")

                if raw.shape != (N_ROWS_TOTAL, 1):
                    raise ValueError(
                        "[{subject} {hemi}] Unexpected shape {shape} in "
                        "{fname} (expected ({rows}, 1)).".format(
                            subject=subject, hemi=hemi,
                            shape=raw.shape, fname=fname,
                            rows=N_ROWS_TOTAL)
                    )

                ipsi_block   = raw.iloc[row_slice_ipsi,   0].values.astype(float)
                contra_block = raw.iloc[row_slice_contra, 0].values.astype(float)

                # Full-breadth row: all 12 macro-regions ipsi, then all 12
                # contra, in macro_regions order — no eye-field restriction
                # applied here (that happens later, positionally, when
                # slicing out the eye-field sub-matrix).
                seed_rows[seed] = np.concatenate([ipsi_block, contra_block])  # (24,)

            if missing_files:
                for f in missing_files:
                    print("    WARNING [{0} {1}]: missing {2}".format(
                        subject, hemi, f))
                print("    {0}: SKIPPED".format(subject))
                missing_subjects.append(subject)
                continue

            mat_full = np.stack([seed_rows[s] for s in macro_regions], axis=0)

            if mat_full.shape != (N_MACRO, 2 * N_MACRO):
                raise ValueError(
                    "[{0} {1} {2}] Unexpected full matrix shape {3}, "
                    "expected ({4}, {5}).".format(
                        subject, hemi, variant, mat_full.shape,
                        N_MACRO, 2 * N_MACRO)
                )

            if variant == "concat_clean" and is_excluded:
                print("    {0}: OK (fallback -> run-01)".format(subject))
            else:
                print("    {0}: OK".format(subject))

            full_matrices.append(mat_full)
            subject_ids.append(subject)

        if not full_matrices:
            print("    ERROR: no valid subjects for {0} / {1} — skipping.".format(
                hemi, variant))
            continue

        n_valid = len(full_matrices)
        print("\n    Valid subjects: {0}/{1}".format(n_valid, len(subjects)))
        if missing_subjects:
            print("    Missing       : {0}".format(missing_subjects))

        # Stack -> (n_subjects x 12 x 24) in Fisher-z space — all-macro scope
        stacked_fz_allmacro = np.stack(full_matrices, axis=0)
        if stacked_fz_allmacro.shape != (n_valid, N_MACRO, 2 * N_MACRO):
            raise ValueError(
                "Unexpected all-macro stack shape {0} for {1} / {2}.".format(
                    stacked_fz_allmacro.shape, hemi, variant)
            )

        # Slice out eye-field scope: rows = EYE_FIELDS_IDX, cols = EYE_FIELD_COL_IDX
        # -> (n_subjects x 5 x 10), bit-for-bit identical to loading only the
        # 5 eye-field seed TSVs directly (same underlying files, same indices).
        stacked_fz = stacked_fz_allmacro[:, EYE_FIELDS_IDX, :][:, :, EYE_FIELD_COL_IDX]
        if stacked_fz.shape != (n_valid, N_EYE_FIELDS, 2 * N_EYE_FIELDS):
            raise ValueError(
                "Unexpected eye-field stack shape {0} for {1} / {2}.".format(
                    stacked_fz.shape, hemi, variant)
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

        print("    [eye-field] Fisher-z mean   range : [{0:.4f}, {1:.4f}]".format(
            np.nanmin(mean_fz), np.nanmax(mean_fz)))
        print("    [eye-field] Fisher-z median range : [{0:.4f}, {1:.4f}]".format(
            np.nanmin(median_fz), np.nanmax(median_fz)))

        # ── Group statistics in Fisher-z space — ALL-MACRO SCOPE ─────────
        # median AND p25/p75 are needed for the reporting TSV's GROUP,
        # GROUP_p25, GROUP_p75 rows (see _save_reporting_tsvs docstring) —
        # not saved as standalone .npy/.csv stat-array files; ask if those
        # are also wanted.
        median_fz_allmacro = np.nanmedian(stacked_fz_allmacro, axis=0)
        pct_lo_fz_allmacro  = np.nanpercentile(stacked_fz_allmacro, PCT_LO, axis=0)
        pct_hi_fz_allmacro  = np.nanpercentile(stacked_fz_allmacro, PCT_HI, axis=0)

        median_r_allmacro  = np.tanh(median_fz_allmacro)
        pct_lo_r_allmacro  = np.tanh(pct_lo_fz_allmacro)
        pct_hi_r_allmacro  = np.tanh(pct_hi_fz_allmacro)

        print("    [all-macro] Fisher-z median range : [{0:.4f}, {1:.4f}]".format(
            np.nanmin(median_fz_allmacro), np.nanmax(median_fz_allmacro)))

        run_label = normal_tag if normal_tag is not None else variant

        pct_lo_tag = "p{0:02d}".format(int(PCT_LO))
        pct_hi_tag = "p{0:02d}".format(int(PCT_HI))

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
                np.save(output_folder / "{0}.npy".format(stem), arr)
                pd.DataFrame(
                    arr, index=EYE_FIELDS, columns=TARGET_COLUMNS
                ).to_csv(
                    output_folder / "{0}.csv".format(stem),
                    float_format="%.4f"
                )
                print("    Saved: {0}.npy / .csv".format(stem))

        # ── Compressed archive with all arrays + metadata — eye-field scope
        npz_stem = (
            "seed-task_by_macror-task_full-corr"
            "_{run_label}_{hemi}_{mode}".format(
                run_label=run_label, hemi=hemi, mode=MODE_LABEL)
        )
        np.savez_compressed(
            output_folder / "{0}.npz".format(npz_stem),
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
            mode              = np.array(MODE_LABEL),
        )
        print("    Saved: {0}.npz".format(npz_stem))

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
print("ALL HEMISPHERES x VARIANTS COMPLETE")
print("=" * 80)
print("\nStats outputs (eye-field scope) : {0}".format(output_folder))
print("Reporting TSVs (both scopes)     : {0}".format(tables_folder))
print(
    "\nNote: Fisher-z outputs are in z-space. Apply np.tanh() to recover "
    "Pearson r only at the final reporting or plotting stage."
)