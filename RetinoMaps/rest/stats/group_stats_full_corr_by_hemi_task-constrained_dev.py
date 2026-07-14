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

    REPORTING MATRIX FORMAT  (5 seeds × 10 targets):
    Both seed rows and target columns are restricted to the 5 core eye-field
    ROIs (mPCS, sPCS, iPCS, sIPS, iIPS).  Target columns are further split
    into ipsilateral and contralateral halves, producing explicitly labelled
    columns that require no downstream slicing:
        mPCS_ipsi, sPCS_ipsi, iPCS_ipsi, sIPS_ipsi, iIPS_ipsi,
        mPCS_contra, sPCS_contra, iPCS_contra, sIPS_contra, iIPS_contra

    For each hemisphere × variant the script:
        1. Resolves which TSV to load per subject (run variant logic, including
           concat_clean fallback for RUN02_EXCLUDED subjects).
        2. Loads all seed TSVs for each subject in EYE_FIELDS order, extracts
           ipsi + contra eye-field values, assembles a (5 × 10) Fisher-z matrix.
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
    final reporting stage via tanh().
------------------------------------------------------------------------------------------
Output filename convention (harmonized with partial-corr group stats):

    seed-task_by_macror-task_full-corr_{space}_{stat}_{run_label}_{hemi}_legacy.npy / .csv
    seed-task_by_macror-task_full-corr_{run_label}_{hemi}_legacy.npz

    Reporting TSVs (concat_clean only):
    seed-task_by_macror-task_full-corr_r_report_ipsi_{hemi}_legacy.tsv
    seed-task_by_macror-task_full-corr_r_report_contra_{hemi}_legacy.tsv

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
    seed-task_by_macror-task_full-corr_fisherz_mean_{run_label}_{hemi}_legacy.npy / .csv
    seed-task_by_macror-task_full-corr_fisherz_median_{run_label}_{hemi}_legacy.npy / .csv
    seed-task_by_macror-task_full-corr_r_mean_{run_label}_{hemi}_legacy.npy / .csv
    seed-task_by_macror-task_full-corr_r_median_{run_label}_{hemi}_legacy.npy / .csv
    seed-task_by_macror-task_full-corr_r_p25_{run_label}_{hemi}_legacy.npy / .csv
    seed-task_by_macror-task_full-corr_r_p75_{run_label}_{hemi}_legacy.npy / .csv
    seed-task_by_macror-task_full-corr_{run_label}_{hemi}_legacy.npz

    Reporting TSVs (concat_clean only, per hemisphere):
    seed-task_by_macror-task_full-corr_r_report_ipsi_{hemi}_legacy.tsv
    seed-task_by_macror-task_full-corr_r_report_contra_{hemi}_legacy.tsv

    Reporting TSV structure:
        subject   : subject ID (e.g. "sub-01") or "GROUP"
        seed      : eye-field seed name (mPCS / sPCS / iPCS / sIPS / iIPS)
        mPCS / sPCS / iPCS / sIPS / iIPS : Pearson r for each target region
            (ipsi half in the ipsi table, contra half in the contra table)
        GROUP row values = tanh(nanmedian(Fisher-z)) across subjects — the
            same quantity saved in the _r_median_ .npy / .csv outputs.

    Rows    : 5 eye-field seeds in canonical order (mPCS first)
    Columns : 5 eye-field regions × 2 hemispheres (ipsi then contra)  (n = 10)
              mPCS_ipsi, sPCS_ipsi, iPCS_ipsi, sIPS_ipsi, iIPS_ipsi,
              mPCS_contra, sPCS_contra, iPCS_contra, sIPS_contra, iIPS_contra

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
print("Eye-field seeds x eye-field targets (5 x 10), legacy mode")
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
# used to locate EYE_FIELDS within each block via EYE_FIELDS_IDX.
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
# Eye-field regions — the 5 core ROIs used for all outputs
#
# Derived positionally (first 5 of macro_regions after reversing rois-drawn)
# so it stays correct as long as rois-drawn ordering is maintained.
# The assertion below guards against any silent mismatch.
# ============================================================
N_EYE_FIELDS = 5
EYE_FIELDS   = macro_regions[:N_EYE_FIELDS]  # type: List[str]

assert EYE_FIELDS == ["mPCS", "sPCS", "iPCS", "sIPS", "iIPS"], (
    "EYE_FIELDS resolved to {0}, expected the 5 eye-field ROIs "
    "in mPCS-first order. Check rois-drawn ordering in settings.yml.".format(
        EYE_FIELDS)
)

# Row indices of EYE_FIELDS within a 12-row hemisphere block
EYE_FIELDS_IDX = [macro_regions.index(r) for r in EYE_FIELDS]

# Column labels for the full (5 x 10) output matrices
TARGET_COLUMNS = (
    ["{0}_ipsi".format(r)  for r in EYE_FIELDS] +
    ["{0}_contra".format(r) for r in EYE_FIELDS]
)  # type: List[str]

print("\n  All macro-regions (n={0}): {1}".format(N_MACRO, macro_regions))
print("  Eye-field regions (n={0}): {1}".format(N_EYE_FIELDS, EYE_FIELDS))
print("  TSV layout   : {0} rows — LH rows 0-{1}, RH rows {2}-{3}".format(
    N_ROWS_TOTAL, N_MACRO - 1, N_MACRO, N_ROWS_TOTAL - 1))
print("  Output shape : ({0} seeds x {1} targets)".format(
    N_EYE_FIELDS, 2 * N_EYE_FIELDS))
print("  Target cols  : {0}".format(TARGET_COLUMNS))

# ============================================================
# Paths
# ============================================================
main_data     = Path(main_dir) / project_dir / "derivatives/pp_data"
output_folder = (
    main_data / "group/91k/rest/full_corr/by_hemi/task-constrained"
)
output_folder.mkdir(parents=True, exist_ok=True)

# ============================================================
# Output filename stem builder
# ============================================================
def _stem(stat, space, run_label, hemi):
    # type: (str, str, str, str) -> str
    return (
        "seed-task_by_macror-task_full-corr"
        "_{space}_{stat}_{run_label}_{hemi}_{mode}".format(
            space=space, stat=stat, run_label=run_label,
            hemi=hemi, mode=MODE_LABEL)
    )


# ============================================================
# Reporting TSV builder (concat_clean only)
#
# Produces two long-format tables — one for ipsi targets, one for contra —
# with the structure:
#   subject | seed | mPCS | sPCS | iPCS | sIPS | iIPS
#
# Subject rows contain raw Pearson r values (tanh of individual Fisher-z,
# never averaged).  The GROUP row at the bottom contains the group median
# in Pearson r space: tanh(nanmedian(Fisher-z across subjects)).
#
# Parameters
# ----------
# stacked_fz  : (n_subjects x N_EYE_FIELDS x 10) Fisher-z array
# median_r    : (N_EYE_FIELDS x 10) group median in Pearson r space
# subject_ids : list of subject ID strings, length n_subjects
# hemi        : "lh" or "rh"
# ============================================================
def _save_reporting_tsvs(stacked_fz, median_r, subject_ids, hemi):
    # type: (np.ndarray, np.ndarray, List[str], str) -> None

    # Column indices for ipsi (first 5) and contra (last 5) within the 10-col
    # matrices — defined positionally to match TARGET_COLUMNS construction.
    ipsi_cols  = list(range(N_EYE_FIELDS))
    contra_cols = list(range(N_EYE_FIELDS, 2 * N_EYE_FIELDS))

    for side, col_idx in (("ipsi", ipsi_cols), ("contra", contra_cols)):
        rows = []  # type: List[Dict]

        # ── Per-subject rows ──────────────────────────────────────────────
        # Convert each subject's Fisher-z slice to Pearson r individually
        # (tanh applied per subject, not to an average).
        for s_idx, subj in enumerate(subject_ids):
            subj_fz = stacked_fz[s_idx]          # (N_EYE_FIELDS x 10)
            subj_r  = np.tanh(subj_fz)           # Pearson r, same shape

            for seed_idx, seed in enumerate(EYE_FIELDS):
                row = {"subject": subj, "seed": seed}
                for t_idx, t_col in enumerate(col_idx):
                    row[EYE_FIELDS[t_idx]] = subj_r[seed_idx, t_col]
                rows.append(row)

        # ── GROUP row ─────────────────────────────────────────────────────
        # Uses the pre-computed median_r which is tanh(nanmedian(Fisher-z))
        # — never averaged Pearson r values.
        for seed_idx, seed in enumerate(EYE_FIELDS):
            row = {"subject": "GROUP", "seed": seed}
            for t_idx, t_col in enumerate(col_idx):
                row[EYE_FIELDS[t_idx]] = median_r[seed_idx, t_col]
            rows.append(row)

        # ── Save ─────────────────────────────────────────────────────────
        col_order = ["subject", "seed"] + EYE_FIELDS
        df = pd.DataFrame(rows, columns=col_order)

        fname = (
            "seed-task_by_macror-task_full-corr"
            "_r_report_{side}_{hemi}_{mode}.tsv".format(
                side=side, hemi=hemi, mode=MODE_LABEL)
        )
        df.to_csv(output_folder / fname, sep="\t", index=False,
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

        subject_matrices = []  # type: List[np.ndarray]
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

            for seed in EYE_FIELDS:
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

                ipsi_ef   = ipsi_block[EYE_FIELDS_IDX]    # (5,)
                contra_ef = contra_block[EYE_FIELDS_IDX]  # (5,)

                seed_rows[seed] = np.concatenate([ipsi_ef, contra_ef])  # (10,)

            if missing_files:
                for f in missing_files:
                    print("    WARNING [{0} {1}]: missing {2}".format(
                        subject, hemi, f))
                print("    {0}: SKIPPED".format(subject))
                missing_subjects.append(subject)
                continue

            mat = np.stack([seed_rows[s] for s in EYE_FIELDS], axis=0)

            if mat.shape != (N_EYE_FIELDS, 2 * N_EYE_FIELDS):
                raise ValueError(
                    "[{0} {1} {2}] Unexpected matrix shape {3}, "
                    "expected ({4}, {5}).".format(
                        subject, hemi, variant, mat.shape,
                        N_EYE_FIELDS, 2 * N_EYE_FIELDS)
                )

            if variant == "concat_clean" and is_excluded:
                print("    {0}: OK (fallback -> run-01)".format(subject))
            else:
                print("    {0}: OK".format(subject))

            subject_matrices.append(mat)
            subject_ids.append(subject)

        if not subject_matrices:
            print("    ERROR: no valid subjects for {0} / {1} — skipping.".format(
                hemi, variant))
            continue

        n_valid = len(subject_matrices)
        print("\n    Valid subjects: {0}/{1}".format(n_valid, len(subjects)))
        if missing_subjects:
            print("    Missing       : {0}".format(missing_subjects))

        # Stack -> (n_subjects x N_EYE_FIELDS x 10) in Fisher-z space
        stacked_fz = np.stack(subject_matrices, axis=0)

        if stacked_fz.shape != (n_valid, N_EYE_FIELDS, 2 * N_EYE_FIELDS):
            raise ValueError(
                "Unexpected stack shape {0} for {1} / {2}.".format(
                    stacked_fz.shape, hemi, variant)
            )

        # ── Group statistics in Fisher-z space ───────────────────────────
        mean_fz   = np.nanmean(  stacked_fz, axis=0)
        median_fz = np.nanmedian(stacked_fz, axis=0)
        pct_lo_fz = np.nanpercentile(stacked_fz, PCT_LO, axis=0)
        pct_hi_fz = np.nanpercentile(stacked_fz, PCT_HI, axis=0)

        # ── Back-convert to Pearson r at reporting stage only ────────────
        mean_r   = np.tanh(mean_fz)
        median_r = np.tanh(median_fz)
        pct_lo_r = np.tanh(pct_lo_fz)
        pct_hi_r = np.tanh(pct_hi_fz)

        print("    Fisher-z mean   range : [{0:.4f}, {1:.4f}]".format(
            np.nanmin(mean_fz), np.nanmax(mean_fz)))
        print("    Fisher-z median range : [{0:.4f}, {1:.4f}]".format(
            np.nanmin(median_fz), np.nanmax(median_fz)))

        run_label = normal_tag if normal_tag is not None else variant

        pct_lo_tag = "p{0:02d}".format(int(PCT_LO))
        pct_hi_tag = "p{0:02d}".format(int(PCT_HI))

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
                np.save(output_folder / "{0}.npy".format(stem), arr)
                pd.DataFrame(
                    arr, index=EYE_FIELDS, columns=TARGET_COLUMNS
                ).to_csv(
                    output_folder / "{0}.csv".format(stem),
                    float_format="%.4f"
                )
                print("    Saved: {0}.npy / .csv".format(stem))

        # ── Compressed archive with all arrays + metadata ─────────────────
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

        # ── Reporting TSVs — concat_clean only ───────────────────────────
        if variant == "concat_clean":
            _save_reporting_tsvs(stacked_fz, median_r, subject_ids, hemi)

print("\n" + "=" * 80)
print("ALL HEMISPHERES x VARIANTS COMPLETE")
print("=" * 80)
print("\nOutputs written to: {0}".format(output_folder))
print(
    "\nNote: Fisher-z outputs are in z-space. Apply np.tanh() to recover "
    "Pearson r only at the final reporting or plotting stage."
)