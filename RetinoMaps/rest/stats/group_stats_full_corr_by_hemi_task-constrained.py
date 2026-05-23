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

    Files are always in legacy mode (hardcoded) 
    no mode argument is needed since this was the only error free mode in -cifti-parcellate

    For each hemisphere × variant the script:
        1. Resolves which TSV to load per subject (same run variant logic as
           the task-free pipeline).
        2. Loads all seed TSVs for a subject, slices the hemisphere-specific
           12 rows, assembles a (n_seeds × 12) Fisher-z matrix.
        3. Stacks matrices across subjects.
        4. Computes group median and std in Fisher-z space.
        5. Saves one .npy + one .csv per statistic per hemisphere × variant.

    Averaging is always in Fisher-z space; apply np.tanh() only at the final
    reporting or plotting stage.
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
         {pp_data}/group/91k/rest/full_corr/by_hemi/task_constrained/):
    seed-task_by_macro_full-corr_fisherz_median_{run_label}_{hemi}_legacy.npy / .csv
    seed-task_by_macro_full-corr_fisherz_std_{run_label}_{hemi}_legacy.npy  / .csv

    Rows    : seed/cluster names in canonical order (mPCS first)  (n = 12)
    Columns : macro-region names in canonical order (mPCS first)  (n = 12)

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
from typing import Dict, List, Optional
warnings.filterwarnings("ignore")

# ============================================================
# Personal imports
# ============================================================
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../"))

sys.path.append(os.path.abspath(os.path.join(base_dir, "analysis_code/utils")))
from settings_utils import load_settings

sys.path.append(os.path.abspath(os.path.join(base_dir, "RetinoMaps/rest/utils")))
from rest_utils import (
    RUN02_EXCLUDED,
    VARIANTS,
)

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
settings_path     = os.path.join(base_dir, project_dir, "settings.yml")
prf_settings_path = os.path.join(base_dir, project_dir, "prf-analysis.yml")
settings          = load_settings([settings_path, prf_settings_path])
analysis_info     = settings[0]
subjects          = analysis_info["subjects"]

# ============================================================
# Macro-regions — canonical order from YAML rois-drawn (reversed)
#
# Both the seed (row) and target (column) dimensions use the same list:
# mPCS first, V1 last.  This matches the row order inside the subject TSVs
# and the display convention used in all figures.
# ============================================================
macro_regions: List[str] = list(analysis_info["rois-drawn"])
macro_regions.reverse()   # mPCS first — intentional display ordering

N_MACRO = len(macro_regions)   # expected: 12

# TSV layout: 24 rows total, LH first then RH.
N_ROWS_TOTAL = N_MACRO * 2
HEMI_ROW_SLICE: Dict[str, slice] = {
    "lh": slice(0,      N_MACRO),
    "rh": slice(N_MACRO, N_ROWS_TOTAL),
}

print(f"\n  Macro-regions (n={N_MACRO}): {macro_regions}")
print(f"  TSV layout   : {N_ROWS_TOTAL} rows — LH rows 0–{N_MACRO-1}, RH rows {N_MACRO}–{N_ROWS_TOTAL-1}")

# ============================================================
# Paths
# ============================================================
main_data     = Path(main_dir) / project_dir / "derivatives/pp_data"
output_folder = (
    main_data / "group/91k/rest/full_corr/by_hemi/task-constrained"
)
output_folder.mkdir(parents=True, exist_ok=True)

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

        for subject in subjects:
            is_excluded = subject in RUN02_EXCLUDED
            run_tag     = excluded_tag if is_excluded else normal_tag
            run_entity  = f"_{run_tag}" if run_tag is not None else ""

            subj_dir = (
                main_data / subject / "91k/rest/corr/full_corr/by_hemi/task-constrained"
            )

            seed_rows: Dict[str, np.ndarray] = {}
            missing:   List[str]             = []

            for seed in macro_regions:
                fname = (
                    f"{subject}_task-rest{run_entity}_space-fsLR_den-91k"
                    f"_desc-fisher-z_{hemi}_{seed}"
                    f"_task-constrained_parcellated_by_macro{TSV_SUFFIX}.tsv"
                )
                fpath = subj_dir / fname

                if not fpath.exists():
                    missing.append(fname)
                    continue

                raw = pd.read_csv(fpath, header=None, sep="\t")

                if raw.shape != (N_ROWS_TOTAL, 1):
                    raise ValueError(
                        f"[{subject} {hemi}] Unexpected shape {raw.shape} "
                        f"in {fname} (expected ({N_ROWS_TOTAL}, 1))."
                    )

                seed_rows[seed] = raw.iloc[row_slice, 0].values.astype(float)

            if missing:
                for f in missing:
                    print(f"    WARNING [{subject} {hemi}]: missing {f}")
                print(f"    {subject}: SKIPPED")
                continue

            # Assemble (n_seeds × n_macro) matrix, seeds in canonical order
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

        # Stack → (n_subjects × n_seeds × n_macro)
        stack = np.stack(subject_matrices, axis=0)
        assert stack.shape == (n_valid, N_MACRO, N_MACRO), (
            f"Unexpected stack shape {stack.shape}."
        )

        group_median = np.nanmedian(stack, axis=0)   # (n_seeds × n_macro)
        group_std    = np.nanstd(   stack, axis=0)

        n_nan_med = int(np.isnan(group_median).sum())
        n_nan_std = int(np.isnan(group_std).sum())
        if n_nan_med:
            print(f"    WARNING: {n_nan_med}/{group_median.size} NaN cells in group median.")
        if n_nan_std:
            print(f"    WARNING: {n_nan_std}/{group_std.size} NaN cells in group std.")

        print(
            f"    Group median range : "
            f"[{np.nanmin(group_median):.4f}, {np.nanmax(group_median):.4f}]"
        )
        print(
            f"    Group std  range   : "
            f"[{np.nanmin(group_std):.4f},  {np.nanmax(group_std):.4f}]"
        )

        # run_label: None normal_tag → use variant name (e.g. "concat")
        run_label = normal_tag if normal_tag is not None else variant

        stem_median = (
            f"seed-task_by_macro_full-corr_fisherz_median"
            f"_{run_label}_{hemi}_{MODE_LABEL}"
        )
        np.save(output_folder / f"{stem_median}.npy", group_median)
        pd.DataFrame(
            group_median, index=macro_regions, columns=macro_regions
        ).to_csv(output_folder / f"{stem_median}.csv")

        stem_std = (
            f"seed-task_by_macro_full-corr_fisherz_std"
            f"_{run_label}_{hemi}_{MODE_LABEL}"
        )
        np.save(output_folder / f"{stem_std}.npy", group_std)
        pd.DataFrame(
            group_std, index=macro_regions, columns=macro_regions
        ).to_csv(output_folder / f"{stem_std}.csv")

        print(f"    Saved: {stem_median}.npy / .csv")
        print(f"    Saved: {stem_std}.npy / .csv")

print("\n" + "=" * 80)
print("ALL HEMISPHERES × VARIANTS COMPLETE")
print("=" * 80)
print(f"\nOutputs written to: {output_folder}")
print(
    "\nNote: outputs are in Fisher-z space. Apply np.tanh() to recover "
    "Pearson r only at the final reporting or plotting stage."
)