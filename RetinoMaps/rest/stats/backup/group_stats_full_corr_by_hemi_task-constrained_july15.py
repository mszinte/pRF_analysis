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

    Files are always in legacy mode (hardcoded): no mode argument is needed
    since this was the only error-free mode in -cifti-parcellate when
    parcellating by macro-regions (some files have missing vertices).

    REPORTING MATRIX FORMAT  (5 seeds × 10 targets):
    Both seed rows and target columns are restricted to the 5 core eye-field
    ROIs (mPCS, sPCS, iPCS, sIPS, iIPS).  Target columns are further split
    into ipsilateral and contralateral halves, producing explicitly labelled
    columns that require no downstream slicing:
        mPCS_ipsi, sPCS_ipsi, iPCS_ipsi, sIPS_ipsi, iIPS_ipsi,
        mPCS_contra, sPCS_contra, iPCS_contra, sIPS_contra, iIPS_contra

    The TSV already contains both hemispheres' target values (24 rows).
    Earlier versions of this script discarded the contralateral half; this
    version extracts both halves and restricts each to EYE_FIELDS.

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

    Averaging is always in Fisher-z space; Pearson r is recovered only at the
    final reporting stage via tanh().  This is intentional:
      - Fisher-z has approximately constant variance ~1/(n-3), making it the
        correct space for averaging and any subsequent parametric tests.
      - Raw Pearson r values have r-dependent variance and must never be
        averaged directly.
------------------------------------------------------------------------------------------
Output filename convention (harmonized with partial-corr group stats):

    seed-task_by_macror-task_full-corr_{space}_{stat}_{run_label}_{hemi}_legacy.npy / .csv
    seed-task_by_macror-task_full-corr_{run_label}_{hemi}_legacy.npz

    e.g. seed-task_by_macror-task_full-corr_fisherz_median_concat_clean_rh_legacy.npy
         seed-task_by_macror-task_full-corr_r_median_concat_clean_rh_legacy.npy

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
    seed-task_by_macror-task_full-corr_{run_label}_{hemi}_legacy.npz

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
    "Usage: python group_full_corr_task_constrained_by_hemi.py "
    "<main_dir> <project_dir> <group> <server> [--percentiles LO HI]\n"
    "  --percentiles  lower and upper percentile bounds for the subject\n"
    "                 distribution saved alongside group mean/median\n"
    "                 (default: 25 75)"
)

if len(sys.argv) not in (5, 8):
    print(f"ERROR: unexpected number of arguments.\n{USAGE}")
    sys.exit(1)

main_dir    = sys.argv[1]
project_dir = sys.argv[2]
group       = sys.argv[3]
server      = sys.argv[4]

# Optional percentile bounds (default: lower and upper quartiles)
pct_lo, pct_hi = 25.0, 75.0
if len(sys.argv) == 8:
    if sys.argv[5] != "--percentiles":
        print(f"ERROR: unknown option '{sys.argv[5]}'.\n{USAGE}")
        sys.exit(1)
    try:
        pct_lo = float(sys.argv[6])
        pct_hi = float(sys.argv[7])
    except ValueError:
        print("ERROR: --percentiles values must be floats (e.g. 2.5 97.5)")
        sys.exit(1)
    if not (0 <= pct_lo < pct_hi <= 100):
        print("ERROR: percentile values must satisfy 0 ≤ lo < hi ≤ 100")
        sys.exit(1)

# Hardcoded: input TSVs are always produced in legacy mode.
TSV_SUFFIX = "_legacy-mode"
MODE_LABEL = "legacy"

print("=" * 80)
print("GROUP FULL CORRELATION (TASK-CONSTRAINED) — Fisher-z statistics")
print("Eye-field seeds × eye-field targets (5 × 10), legacy mode")
print("=" * 80)
print(f"  main_dir    : {main_dir}")
print(f"  project_dir : {project_dir}")
print(f"  group       : {group}")
print(f"  server      : {server}")
print(f"  percentiles : {pct_lo} / {pct_hi}")

# ============================================================
# Load settings
# ============================================================
settings_path     = os.path.join(base_dir, project_dir, "settings.yml")
prf_settings_path = os.path.join(base_dir, project_dir, "prf-analysis.yml")
settings          = load_settings([settings_path, prf_settings_path])
analysis_info     = settings[0]
subjects          = analysis_info["subjects"]

# ============================================================
# Macro-regions — full canonical list (mPCS first)
#
# The full list defines the TSV row order (12 per hemisphere block) and is
# used to locate EYE_FIELDS within each block via EYE_FIELDS_IDX.
# ============================================================
macro_regions: List[str] = list(analysis_info["rois-drawn"])
macro_regions.reverse()   # mPCS first

N_MACRO = len(macro_regions)   # expected: 12

# TSV layout: 24 rows total, LH block first then RH block.
N_ROWS_TOTAL = N_MACRO * 2
HEMI_ROW_SLICE: Dict[str, slice] = {
    "lh": slice(0,       N_MACRO),
    "rh": slice(N_MACRO, N_ROWS_TOTAL),
}

# ============================================================
# Eye-field regions — the 5 core ROIs used for all outputs
#
# EYE_FIELDS is derived positionally (first 5 of macro_regions) so it stays
# correct as long as rois-drawn ordering is maintained.  The assertion below
# guards against any silent mismatch — it must pass before processing begins.
#
# Seed rows AND target columns are both restricted to EYE_FIELDS.
# Target columns are further split into ipsi and contra halves (10 total).
# ============================================================
N_EYE_FIELDS = 5
EYE_FIELDS   = macro_regions[:N_EYE_FIELDS]

assert EYE_FIELDS == ["mPCS", "sPCS", "iPCS", "sIPS", "iIPS"], (
    f"EYE_FIELDS resolved to {EYE_FIELDS}, expected the 5 eye-field ROIs "
    "in mPCS-first order. Check rois-drawn ordering in settings.yml."
)

# Row indices of EYE_FIELDS within macro_regions (used to locate seeds in TSV)
# — identical to positional order here but kept explicit for clarity
EYE_FIELDS_IDX = [macro_regions.index(r) for r in EYE_FIELDS]

# Target column labels: [5 ipsi | 5 contra], self-documenting for co-authors
TARGET_COLUMNS = (
    [f"{r}_ipsi"  for r in EYE_FIELDS] +
    [f"{r}_contra" for r in EYE_FIELDS]
)

print(f"\n  All macro-regions (n={N_MACRO}): {macro_regions}")
print(f"  Eye-field regions (n={N_EYE_FIELDS}): {EYE_FIELDS}")
print(f"  TSV layout   : {N_ROWS_TOTAL} rows — LH rows 0–{N_MACRO-1}, "
      f"RH rows {N_MACRO}–{N_ROWS_TOTAL-1}")
print(f"  Output shape : ({N_EYE_FIELDS} seeds × {2*N_EYE_FIELDS} targets)")
print(f"  Target cols  : {TARGET_COLUMNS}")

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
#
# Pattern: seed-task_by_macror-task_full-corr_{space}_{stat}_{run_label}_{hemi}_{mode}
# Matches partial-corr: seed-task_by_macror-task_partial-corr_{space}_{stat}_{run_label}_{hemi}
#
# run_label is passed explicitly so the function remains pure and testable
# without relying on loop-scoped variables.
# ============================================================
def _stem(stat: str, space: str, run_label: str, hemi: str) -> str:
    return (
        f"seed-task_by_macror-task_full-corr"
        f"_{space}_{stat}_{run_label}_{hemi}_{MODE_LABEL}"
    )

# ============================================================
# Main loop — hemisphere × variant
# ============================================================

for hemi in ("lh", "rh"):
    print(f"\n{'='*80}")
    print(f"Processing hemisphere: {hemi.upper()}")
    print("=" * 80)

    # Row slices into the 24-row TSV for ipsilateral and contralateral blocks
    contra_hemi      = "rh" if hemi == "lh" else "lh"
    row_slice_ipsi   = HEMI_ROW_SLICE[hemi]
    row_slice_contra = HEMI_ROW_SLICE[contra_hemi]

    for variant, (normal_tag, excluded_tag, _skip) in VARIANTS.items():
        # run-02: intentionally retains all subjects for group QC (artifact
        # visibility), regardless of the skip_excluded flag used in WTA scripts.
        print(f"\n  --- Variant: {variant} ---")

        subject_matrices: List[np.ndarray] = []
        subject_ids:      List[str]        = []
        missing_subjects: List[str]        = []

        for subject in subjects:
            is_excluded = subject in RUN02_EXCLUDED
            run_tag     = excluded_tag if is_excluded else normal_tag
            run_entity  = f"_{run_tag}" if run_tag is not None else ""

            subj_dir = (
                main_data / subject
                / "91k/rest/corr/full_corr/by_hemi/task-constrained"
            )

            # Load one TSV per eye-field seed; accumulate into seed_rows dict
            seed_rows:     Dict[str, np.ndarray] = {}
            missing_files: List[str]             = []

            for seed in EYE_FIELDS:
                fname = (
                    f"{subject}_task-rest{run_entity}_space-fsLR_den-91k"
                    f"_desc-fisher-z_{hemi}_{seed}"
                    f"_task-constrained_parcellated_by_macro{TSV_SUFFIX}.tsv"
                )
                fpath = subj_dir / fname

                if not fpath.exists():
                    missing_files.append(fname)
                    continue

                raw = pd.read_csv(fpath, header=None, sep="\t")

                if raw.shape != (N_ROWS_TOTAL, 1):
                    raise ValueError(
                        f"[{subject} {hemi}] Unexpected shape {raw.shape} in "
                        f"{fname} (expected ({N_ROWS_TOTAL}, 1))."
                    )

                # Extract both hemisphere blocks from the TSV, then restrict
                # each to the EYE_FIELDS positions within that block.
                # The TSV has 24 rows: rows 0–11 are LH targets, 12–23 are RH.
                # Earlier pipeline versions only read the ipsi block and silently
                # discarded the contra half — both are now retained.
                ipsi_block   = raw.iloc[row_slice_ipsi,   0].values.astype(float)
                contra_block = raw.iloc[row_slice_contra, 0].values.astype(float)

                ipsi_ef   = ipsi_block[EYE_FIELDS_IDX]    # (5,)
                contra_ef = contra_block[EYE_FIELDS_IDX]  # (5,)

                # Concatenate [ipsi | contra] — matches TARGET_COLUMNS order
                seed_rows[seed] = np.concatenate([ipsi_ef, contra_ef])  # (10,)

            if missing_files:
                for f in missing_files:
                    print(f"    WARNING [{subject} {hemi}]: missing {f}")
                print(f"    {subject}: SKIPPED")
                missing_subjects.append(subject)
                continue

            # Assemble (N_EYE_FIELDS × 2*N_EYE_FIELDS) Fisher-z matrix.
            # Rows follow EYE_FIELDS order; columns follow TARGET_COLUMNS order.
            mat = np.stack([seed_rows[s] for s in EYE_FIELDS], axis=0)

            if mat.shape != (N_EYE_FIELDS, 2 * N_EYE_FIELDS):
                raise ValueError(
                    f"[{subject} {hemi} {variant}] Unexpected matrix shape "
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

        # Stack → (n_subjects × N_EYE_FIELDS × 2*N_EYE_FIELDS)
        # Input matrices are already in Fisher-z space (workbench output).
        stacked_fz = np.stack(subject_matrices, axis=0)

        if stacked_fz.shape != (n_valid, N_EYE_FIELDS, 2 * N_EYE_FIELDS):
            raise ValueError(
                f"Unexpected stack shape {stacked_fz.shape} for "
                f"{hemi} / {variant}."
            )

        # ── Group statistics in Fisher-z space ───────────────────────────────
        mean_fz   = np.nanmean(  stacked_fz, axis=0)   # (N_EYE_FIELDS × 10)
        median_fz = np.nanmedian(stacked_fz, axis=0)

        # Inter-subject percentile bounds in Fisher-z space.
        # These describe the spread of subject values around the group statistic,
        # not uncertainty of the estimate itself (for that, use bootstrap CIs).
        pct_lo_fz = np.nanpercentile(stacked_fz, pct_lo, axis=0)
        pct_hi_fz = np.nanpercentile(stacked_fz, pct_hi, axis=0)

        # ── Back-convert to Pearson r at reporting stage only ─────────────────
        mean_r     = np.tanh(mean_fz)
        median_r   = np.tanh(median_fz)
        pct_lo_r   = np.tanh(pct_lo_fz)
        pct_hi_r   = np.tanh(pct_hi_fz)

        print(f"    Fisher-z mean   range : [{np.nanmin(mean_fz):.4f},   {np.nanmax(mean_fz):.4f}]")
        print(f"    Fisher-z median range : [{np.nanmin(median_fz):.4f}, {np.nanmax(median_fz):.4f}]")

        # run_label: None normal_tag → variant name ("concat", "concat_clean")
        run_label = normal_tag if normal_tag is not None else variant

        # Percentile label strings used in filenames (e.g. "p25", "p75")
        pct_lo_tag = f"p{int(pct_lo):02d}"
        pct_hi_tag = f"p{int(pct_hi):02d}"

        # ── Save Fisher-z arrays (full precision .npy; 2 d.p. .csv) ──────────
        for stat, arr in (
            ("mean",     mean_fz),
            ("median",   median_fz),
            (pct_lo_tag, pct_lo_fz),
            (pct_hi_tag, pct_hi_fz),
        ):
            stem = _stem(stat, "fisherz", run_label, hemi)
            np.save(output_folder / f"{stem}.npy", arr)
            pd.DataFrame(arr, index=EYE_FIELDS, columns=TARGET_COLUMNS).to_csv(
                output_folder / f"{stem}.csv", float_format="%.2f"
            )
            print(f"    Saved: {stem}.npy / .csv")

        # ── Save Pearson r arrays (full precision .npy; 2 d.p. .csv) ─────────
        for stat, arr in (
            ("mean",     mean_r),
            ("median",   median_r),
            (pct_lo_tag, pct_lo_r),
            (pct_hi_tag, pct_hi_r),
        ):
            stem = _stem(stat, "r", run_label, hemi)
            np.save(output_folder / f"{stem}.npy", arr)
            pd.DataFrame(arr, index=EYE_FIELDS, columns=TARGET_COLUMNS).to_csv(
                output_folder / f"{stem}.csv", float_format="%.2f"
            )
            print(f"    Saved: {stem}.npy / .csv")

        # ── Compressed archive with all arrays + metadata ─────────────────────
        npz_stem = (
            f"seed-task_by_macror-task_full-corr_{run_label}_{hemi}_{MODE_LABEL}"
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
            mode              = np.array(MODE_LABEL),
            percentiles       = np.array([pct_lo, pct_hi]),
        )
        print(f"    Saved: {npz_stem}.npz")

print("\n" + "=" * 80)
print("ALL HEMISPHERES × VARIANTS COMPLETE")
print("=" * 80)
print(f"\nOutputs written to: {output_folder}")
print(
    "\nNote: Fisher-z outputs are in z-space. Apply np.tanh() to recover Pearson r "
    "only at the final reporting or plotting stage."
)