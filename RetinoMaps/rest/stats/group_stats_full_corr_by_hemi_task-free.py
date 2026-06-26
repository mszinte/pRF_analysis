#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
group_full_corr_task_constrained_by_hemi.py
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

    Files are always in legacy mode (hardcoded):
    no mode argument is needed since this was the only error-free mode in
    -cifti-parcellate when parcellating by macro-regions (some files have
    missing vertices).

    For each hemisphere × variant the script:
        1. Resolves which TSV to load per subject (same run variant logic as
           the task-free pipeline).
        2. Loads all seed TSVs for a subject, slices the hemisphere-specific
           12 rows, assembles a (n_seeds × 12) Fisher-z matrix.
        3. Stacks matrices across subjects.
        4. Computes group mean and median in Fisher-z space.
        5. Back-converts to Pearson r via tanh().
        6. Saves .npy + .csv for both spaces (fisherz and r) per hemisphere × variant,
           plus a compressed .npz archive.

    Averaging is always in Fisher-z space; Pearson r is recovered only at the
    final reporting stage via tanh().  This is intentional:
      - Fisher-z has approximately constant variance ~1/(n-3), making it the
        correct space for averaging and any subsequent parametric tests.
      - Raw Pearson r values have r-dependent variance and must never be
        averaged directly.

------------------------------------------------------------------------------------------
Run variants (identical to task-free pipeline):
    concat       — concatenated-run TSV, all subjects
    concat_clean — best available run per subject:
                     · RUN02_EXCLUDED subjects → run-01 TSV
                     · all other subjects      → concatenated-run TSV
    run-01       — run-01 TSV, all subjects
    run-02       — run-02 TSV, all subjects (bad subjects kept intentionally
                   to expose the registration artifact in group plots)
------------------------------------------------------------------------------------------
Inputs (sys.argv):
    1: main project directory   (e.g. /scratch/mszinte/data)
    2: project name/directory   (e.g. RetinoMaps)
    3: server group             (e.g. 327)
    4: server project           (e.g. b327)

Outputs (per hemisphere × variant, written to
         {pp_data}/group/91k/rest/full_corr/by_hemi/task-constrained/):
    seed-task_by_mmp-parcel_full-corr_fisherz_mean_{run_label}_{hemi}_legacy.npy / .csv
    seed-task_by_mmp-parcel_full-corr_fisherz_median_{run_label}_{hemi}_legacy.npy / .csv
    seed-task_by_mmp-parcel_full-corr_r_mean_{run_label}_{hemi}_legacy.npy / .csv
    seed-task_by_mmp-parcel_full-corr_r_median_{run_label}_{hemi}_legacy.npy / .csv
    seed-task_by_mmp-parcel_full-corr_{run_label}_{hemi}_legacy.npz

    Rows    : seed/cluster names in canonical order (mPCS first)  (n = 12)
    Columns : macro-region names in canonical order (mPCS first)  (n = 12)

Filename example (input TSV, run-01, lh seed hMT+):
    sub-05_task-rest_run-01_space-fsLR_den-91k_desc-fisher-z_lh_hMT+
        _task-free_parcellated_by_macro_legacy-mode.tsv

To run:
    $ cd projects/pRF_analysis/RetinoMaps/rest/stats
    $ python group_stats_full_corr_by_hemi_task-free.py /scratch/mszinte/data RetinoMaps 327 b327
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
warnings.filterwarnings("ignore")

# ============================================================
# Parse and validate arguments
# ============================================================
USAGE = (
    "Usage: python group_full_corr_task_constrained_by_hemi.py "
    "<main_dir> <project_dir> <group> <server>"
)

if len(sys.argv) != 5:
    print(f"ERROR: expected 4 arguments, got {len(sys.argv) - 1}.\n{USAGE}")
    sys.exit(1)

main_dir    = sys.argv[1]
project_dir = sys.argv[2]
group       = sys.argv[3]
server      = sys.argv[4]

# Hardcoded: these files are always produced in legacy mode.
TSV_SUFFIX = "_legacy-mode"
MODE_LABEL = "legacy"

print("=" * 80)
print("GROUP FULL CORRELATION (TASK-CONSTRAINED) — Fisher-z statistics")
print("Parcellated by macro-region (12 seeds × 12 targets), legacy mode")
print("=" * 80)
print(f"  main_dir    : {main_dir}")
print(f"  project_dir : {project_dir}")
print(f"  group       : {group}")
print(f"  server      : {server}")

# ============================================================
# Load settings
# ============================================================
base_dir            = os.path.abspath(os.path.join(os.getcwd(), "../../../"))
sys.path.append(os.path.abspath(os.path.join(base_dir, "analysis_code/utils")))
from settings_utils import load_settings
sys.path.append(os.path.abspath(os.path.join(base_dir, "RetinoMaps/rest/utils")))
from rest_utils import VARIANTS

settings_path       = os.path.join(base_dir, project_dir, "settings.yml")
prf_settings_path   = os.path.join(base_dir, project_dir, "prf-analysis.yml")
settings            = load_settings([settings_path, prf_settings_path])
rest_settings_path  = os.path.join(base_dir, project_dir, "rest-settings.yml")
rest_settings       = load_settings([rest_settings_path])[0]
analysis_info       = settings[0]
subjects            = analysis_info["subjects"]
RUNS                = rest_settings["runs"][0]
RUN02_EXCLUDED      = frozenset(rest_settings["run02_excluded"][0])

# ============================================================
# Macro-regions — canonical order from YAML rois-drawn (reversed)
#
# Both the seed (row) and target (column) dimensions use the same list:
# mPCS first, V1 last.  This matches the row order inside the subject TSVs
# and the display convention used in all figures.
# ============================================================
macro_regions: List[str] = list(analysis_info["rois-drawn"])
macro_regions.reverse()   # mPCS first

N_MACRO = len(macro_regions)   # expected: 12

# TSV layout: 24 rows total, LH first then RH.
N_ROWS_TOTAL = N_MACRO * 2
HEMI_ROW_SLICE: Dict[str, slice] = {
    "lh": slice(0,       N_MACRO),
    "rh": slice(N_MACRO, N_ROWS_TOTAL),
}

print(f"\n  Macro-regions (n={N_MACRO}): {macro_regions}")
print(f"  TSV layout   : {N_ROWS_TOTAL} rows — LH rows 0–{N_MACRO-1}, RH rows {N_MACRO}–{N_ROWS_TOTAL-1}")

# ============================================================
# Paths
# ============================================================
main_data     = Path(main_dir) / project_dir / "derivatives/pp_data"
output_folder = main_data / "group/91k/rest/full_corr/by_hemi/task-free"
output_folder.mkdir(parents=True, exist_ok=True)

# ============================================================
# Output filename stem builder
#
# Pattern: seed-task_by_macro_full-corr_{space}_{stat}_{run_label}_{hemi}_{mode}
# Matches partial-corr convention: seed-task_by_macror-task_partial-corr_{space}_{stat}_{run_label}_{hemi}
# ============================================================
def _stem(stat: str, space: str, run_label: str, hemi: str) -> str:
    return (
        f"seed-task_by_mmp-parcel_full-corr"
        f"_{space}_{stat}_{run_label}_{hemi}_{MODE_LABEL}"
    )

# ============================================================
# Main loop — hemisphere × variant
# ============================================================

for hemi in ("lh", "rh"):
    print(f"\n{'='*80}")
    print(f"Processing hemisphere: {hemi.upper()}")
    print("=" * 80)

    row_slice = HEMI_ROW_SLICE[hemi]

    for variant, (normal_tag, excluded_tag, _skip) in VARIANTS.items():
        # run-02: intentionally keeps all subjects (group QC strategy).
        print(f"\n  --- Variant: {variant} ---")

        subject_matrices: List[np.ndarray] = []
        subject_ids:      List[str]        = []
        missing_subjects: List[str]        = []

        for subject in subjects:
            is_excluded = subject in RUN02_EXCLUDED
            run_tag     = excluded_tag if is_excluded else normal_tag
            run_entity  = f"_{run_tag}" if run_tag is not None else ""

            subj_dir = (
                main_data / subject / "91k/rest/corr/full_corr/by_hemi/task-free"
            )

            seed_rows: Dict[str, np.ndarray] = {}
            missing_files: List[str]         = []

            for seed in macro_regions:
                fname = (
                    f"{subject}_task-rest{run_entity}_space-fsLR_den-91k"
                    f"_desc-fisher-z_{hemi}_{seed}"
                    f"_task-free_parcellated_{TSV_SUFFIX}.tsv"
                )
                fpath = subj_dir / fname

                if not fpath.exists():
                    missing_files.append(fname)
                    continue

                raw = pd.read_csv(fpath, header=None, sep="\t")

                if raw.shape != (N_ROWS_TOTAL, 1):
                    raise ValueError(
                        f"[{subject} {hemi}] Unexpected shape {raw.shape} "
                        f"in {fname} (expected ({N_ROWS_TOTAL}, 1))."
                    )

                seed_rows[seed] = raw.iloc[row_slice, 0].values.astype(float)

            if missing_files:
                for f in missing_files:
                    print(f"    WARNING [{subject} {hemi}]: missing {f}")
                print(f"    {subject}: SKIPPED")
                missing_subjects.append(subject)
                continue

            # Assemble (n_seeds × n_macro) Fisher-z matrix, seeds in canonical order
            mat = np.stack([seed_rows[s] for s in macro_regions], axis=0)

            assert mat.shape == (N_MACRO, N_MACRO), (
                f"[{subject} {hemi} {variant}] Unexpected matrix shape "
                f"{mat.shape}, expected ({N_MACRO}, {N_MACRO})."
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

        # Stack → (n_subjects × n_seeds × n_macro); inputs are already Fisher-z
        stacked_fz = np.stack(subject_matrices, axis=0)
        assert stacked_fz.shape == (n_valid, N_MACRO, N_MACRO), (
            f"Unexpected stack shape {stacked_fz.shape}."
        )

        # Group statistics in Fisher-z space
        mean_fz   = np.nanmean(  stacked_fz, axis=0)   # (n_seeds × n_macro)
        median_fz = np.nanmedian(stacked_fz, axis=0)

        # Back-convert to Pearson r at reporting stage only
        mean_r   = np.tanh(mean_fz)
        median_r = np.tanh(median_fz)

        print(f"    Fisher-z mean   range : [{np.nanmin(mean_fz):.4f},   {np.nanmax(mean_fz):.4f}]")
        print(f"    Fisher-z median range : [{np.nanmin(median_fz):.4f}, {np.nanmax(median_fz):.4f}]")

        # run_label: None normal_tag → use variant name (e.g. "concat", "concat_clean")
        run_label = normal_tag if normal_tag is not None else variant

        # Save Fisher-z arrays
        for stat, arr in (("mean", mean_fz), ("median", median_fz)):
            stem = _stem(stat, "fisherz", run_label, hemi)
            np.save(output_folder / f"{stem}.npy", arr)
            pd.DataFrame(arr, index=macro_regions, columns=macro_regions).to_csv(
                output_folder / f"{stem}.csv"
            )
            print(f"    Saved: {stem}.npy / .csv")

        # Save Pearson r arrays
        for stat, arr in (("mean", mean_r), ("median", median_r)):
            stem = _stem(stat, "r", run_label, hemi)
            np.save(output_folder / f"{stem}.npy", arr)
            pd.DataFrame(arr, index=macro_regions, columns=macro_regions).to_csv(
                output_folder / f"{stem}.csv"
            )
            print(f"    Saved: {stem}.npy / .csv")

        # Compressed archive with all four arrays + metadata
        npz_stem = f"seed-task_by_mmp-parcel_full-corr_{run_label}_{hemi}_{MODE_LABEL}"
        np.savez_compressed(
            output_folder / f"{npz_stem}.npz",
            mean_fz           = mean_fz,
            median_fz         = median_fz,
            mean_r            = mean_r,
            median_r          = median_r,
            n_subjects_loaded = np.array(n_valid),
            subjects_loaded   = np.array(subject_ids),
            subjects_missing  = np.array(missing_subjects),
            macro_regions     = np.array(macro_regions),
            hemi              = np.array(hemi),
            variant           = np.array(variant),
            mode              = np.array(MODE_LABEL),
        )
        print(f"    Saved: {npz_stem}.npz")

print("\n" + "=" * 80)
print("ALL HEMISPHERES × VARIANTS COMPLETE")
print("=" * 80)
print(f"\nOutputs written to: {output_folder}")
print(
    "\nNote: Pearson r group values are recovered back from Fisher-z outputs with np.tanh()"
    "Use these values at the final reporting or plotting stage as is standard practice."
)