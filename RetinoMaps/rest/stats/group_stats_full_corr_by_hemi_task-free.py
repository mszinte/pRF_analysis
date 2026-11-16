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
# [CORRECTION #50] optional 5th argument <mode> (default|legacy|no_outliers, default legacy): the Workbench script
# [CORRECTION #50] writes all 3 modes and the full-correlation WTA reads ..._{mode}.npy -> run this script with the same mode as the WTA
# BEFORE: USAGE = (
# BEFORE:     "Usage: python group_full_corr_task_constrained_by_hemi.py "
# BEFORE:     "<main_dir> <project_dir> <group> <server>"
# BEFORE: )
# BEFORE:
# BEFORE: if len(sys.argv) != 5:
# BEFORE:     print(f"ERROR: expected 4 arguments, got {len(sys.argv) - 1}.\n{USAGE}")
# BEFORE:     sys.exit(1)
USAGE = (
    "Usage: python group_stats_full_corr_by_hemi_task-free.py "
    "<main_dir> <project_dir> <group> <server> [mode: default|legacy|no_outliers]"
)

if len(sys.argv) not in (5, 6):
    print(f"ERROR: expected 4 or 5 arguments, got {len(sys.argv) - 1}.\n{USAGE}")
    sys.exit(1)

main_dir    = sys.argv[1]
project_dir = sys.argv[2]
group       = sys.argv[3]
server      = sys.argv[4]

# [CORRECTION #50] suffix taken from rest_utils.MODE_SUFFIX (no extra '_'): "_parcellated_{TSV_SUFFIX}" gave
# [CORRECTION #50] "_parcellated__legacy-mode.tsv" (double underscore) -> no file found; MODE_LABEL = name of the WTA mode
# BEFORE: # Hardcoded: these files are always produced in legacy mode.
# BEFORE: TSV_SUFFIX = "_legacy-mode"
# BEFORE: MODE_LABEL = "legacy"
MODE_LABEL = sys.argv[5] if len(sys.argv) == 6 else "legacy"
if MODE_LABEL not in ("default", "legacy", "no_outliers"):
    print(f"ERROR: unrecognised mode '{MODE_LABEL}'.\n{USAGE}")
    sys.exit(1)

print("=" * 80)
# [CORRECTION #50] labels: task-free script (12 seeds x 53 MMP parcels per hemisphere)
# BEFORE: print("GROUP FULL CORRELATION (TASK-CONSTRAINED) — Fisher-z statistics")
# BEFORE: print("Parcellated by macro-region (12 seeds × 12 targets), legacy mode")
print("GROUP FULL CORRELATION (TASK-FREE) — Fisher-z statistics")
print(f"Parcellated by MMP parcel (12 seeds × 53 parcels), {MODE_LABEL} mode")
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
# [CORRECTION #50] same utilities as wta_full_corr_by_subject_by_hemi_task-free.py (suffix, atlas order, remap, loading)
# BEFORE: from rest_utils import VARIANTS
from rest_utils import VARIANTS, MODE_SUFFIX, ATLAS_ORDER, build_remap, load_full_corr_matrix
TSV_SUFFIX = MODE_SUFFIX[MODE_LABEL]

settings_path       = os.path.join(base_dir, project_dir, "settings.yml")
prf_settings_path   = os.path.join(base_dir, project_dir, "prf-analysis.yml")
settings            = load_settings([settings_path, prf_settings_path])
rest_settings_path  = os.path.join(base_dir, project_dir, "rest-settings.yml")
rest_settings       = load_settings([rest_settings_path])[0]
analysis_info       = settings[0]
subjects            = analysis_info["subjects"]
# [CORRECTION #26] load_settings(...)[0] already returns the values: the extra [0] took the first element of the list
# [CORRECTION #26] ('sub-03') and frozenset turned it into a set of letters -> no subject was excluded from run-02 (concat_clean)
# BEFORE: RUNS                = rest_settings["runs"][0]
# BEFORE: RUN02_EXCLUDED      = frozenset(rest_settings["run02_excluded"][0])
RUNS                = rest_settings["runs"]
RUN02_EXCLUDED      = frozenset(rest_settings["run02_excluded"])

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

# [CORRECTION #50] task-free TSVs have 106 rows (53 RH parcels then 53 LH, atlas key order), not 24:
# [CORRECTION #50] columns = 53 MMP parcels in YAML order (rois-group-mmp), atlas -> YAML remap as in the full-correlation WTA
# BEFORE: # TSV layout: 24 rows total, LH first then RH.
# BEFORE: N_ROWS_TOTAL = N_MACRO * 2
# BEFORE: HEMI_ROW_SLICE: Dict[str, slice] = {
# BEFORE:     "lh": slice(0,       N_MACRO),
# BEFORE:     "rh": slice(N_MACRO, N_ROWS_TOTAL),
# BEFORE: }
# BEFORE:
# BEFORE: print(f"\n  Macro-regions (n={N_MACRO}): {macro_regions}")
# BEFORE: print(f"  TSV layout   : {N_ROWS_TOTAL} rows — LH rows 0–{N_MACRO-1}, RH rows {N_MACRO}–{N_ROWS_TOTAL-1}")
seed_to_parcels: Dict[str, List[str]] = analysis_info["rois-group-mmp"]
parcels: List[str] = []
for cl in macro_regions:
    parcels.extend(seed_to_parcels[cl])
N_PARCELS = len(parcels)       # expected: 53

REMAP: Dict[str, tuple] = {
    h: build_remap(ATLAS_ORDER[h], parcels, h) for h in ("rh", "lh")
}

print(f"\n  Macro-regions (n={N_MACRO}): {macro_regions}")
# [CORRECTION #50] printout adapted to the targets (53 MMP parcels) instead of the macro-regions
# BEFORE: print(f"  TSV layout   : {N_ROWS_TOTAL} rows — LH rows 0–{N_MACRO-1}, RH rows {N_MACRO}–{N_ROWS_TOTAL-1}")
print(f"  MMP parcels  (n={N_PARCELS}), YAML order; TSV suffix: '_parcellated{TSV_SUFFIX}.tsv'")

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

    # [CORRECTION #50] atlas -> YAML order remap (instead of simply slicing 12 rows)
    # BEFORE:     row_slice = HEMI_ROW_SLICE[hemi]
    remap_idx, present = REMAP[hemi]

    for variant, (normal_tag, excluded_tag, _skip) in VARIANTS.items():
        # run-02: intentionally keeps all subjects (group QC strategy).
        print(f"\n  --- Variant: {variant} ---")

        subject_matrices: List[np.ndarray] = []
        subject_ids:      List[str]        = []
        missing_subjects: List[str]        = []

        for subject in subjects:
            is_excluded = subject in RUN02_EXCLUDED
            run_tag     = excluded_tag if is_excluded else normal_tag
            # [CORRECTION #50] loading through rest_utils.load_full_corr_matrix (exact names written by the Workbench script,
            # [CORRECTION #50] 106 rows RH then LH, atlas -> YAML remap) -> matrix (12 seeds x 53 parcels) in Fisher z
            # BEFORE:             run_entity  = f"_{run_tag}" if run_tag is not None else ""
            # BEFORE: 
            # BEFORE:             subj_dir = (
            # BEFORE:                 main_data / subject / "91k/rest/corr/full_corr/by_hemi/task-free"
            # BEFORE:             )
            # BEFORE: 
            # BEFORE:             seed_rows: Dict[str, np.ndarray] = {}
            # BEFORE:             missing_files: List[str]         = []
            # BEFORE: 
            # BEFORE:             for seed in macro_regions:
            # BEFORE:                 fname = (
            # BEFORE:                     f"{subject}_task-rest{run_entity}_space-fsLR_den-91k"
            # BEFORE:                     f"_desc-fisher-z_{hemi}_{seed}"
            # BEFORE:                     f"_task-free_parcellated_{TSV_SUFFIX}.tsv"
            # BEFORE:                 )
            # BEFORE:                 fpath = subj_dir / fname
            # BEFORE: 
            # BEFORE:                 if not fpath.exists():
            # BEFORE:                     missing_files.append(fname)
            # BEFORE:                     continue
            # BEFORE: 
            # BEFORE:                 raw = pd.read_csv(fpath, header=None, sep="\t")
            # BEFORE: 
            # BEFORE:                 if raw.shape != (N_ROWS_TOTAL, 1):
            # BEFORE:                     raise ValueError(
            # BEFORE:                         f"[{subject} {hemi}] Unexpected shape {raw.shape} "
            # BEFORE:                         f"in {fname} (expected ({N_ROWS_TOTAL}, 1))."
            # BEFORE:                     )
            # BEFORE: 
            # BEFORE:                 seed_rows[seed] = raw.iloc[row_slice, 0].values.astype(float)
            # BEFORE: 
            # BEFORE:             if missing_files:
            # BEFORE:                 for f in missing_files:
            # BEFORE:                     print(f"    WARNING [{subject} {hemi}]: missing {f}")
            # BEFORE:                 print(f"    {subject}: SKIPPED")
            # BEFORE:                 missing_subjects.append(subject)
            # BEFORE:                 continue
            # BEFORE: 
            # BEFORE:             # Assemble (n_seeds × n_macro) Fisher-z matrix, seeds in canonical order
            # BEFORE:             mat = np.stack([seed_rows[s] for s in macro_regions], axis=0)
            # BEFORE: 
            # BEFORE:             assert mat.shape == (N_MACRO, N_MACRO), (
            # BEFORE:                 f"[{subject} {hemi} {variant}] Unexpected matrix shape "
            # BEFORE:                 f"{mat.shape}, expected ({N_MACRO}, {N_MACRO})."
            # BEFORE:             )
            df_corr = load_full_corr_matrix(
                subject    = subject,
                hemi       = hemi,
                run_tag    = run_tag,
                clusters   = macro_regions,
                parcels    = parcels,
                remap_idx  = remap_idx,
                present    = present,
                main_data  = main_data,
                tsv_suffix = TSV_SUFFIX,
                task       = "task-free",
            )

            if df_corr is None:
                print(f"    {subject}: SKIPPED")
                missing_subjects.append(subject)
                continue

            # Assemble (n_seeds × n_parcels) Fisher-z matrix, seeds in canonical order
            mat = df_corr.loc[macro_regions, parcels].values.astype(float)

            assert mat.shape == (N_MACRO, N_PARCELS), (
                f"[{subject} {hemi} {variant}] Unexpected matrix shape "
                f"{mat.shape}, expected ({N_MACRO}, {N_PARCELS})."
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
        # [CORRECTION #50] expected shape (n_subjects x 12 seeds x 53 parcels)
        # BEFORE:         assert stacked_fz.shape == (n_valid, N_MACRO, N_MACRO), (
        assert stacked_fz.shape == (n_valid, N_MACRO, N_PARCELS), (
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
            # [CORRECTION #50] columns = MMP parcels (YAML order), rows = seeds
            # BEFORE:             pd.DataFrame(arr, index=macro_regions, columns=macro_regions).to_csv(
            pd.DataFrame(arr, index=macro_regions, columns=parcels).to_csv(
                output_folder / f"{stem}.csv"
            )
            print(f"    Saved: {stem}.npy / .csv")

        # Save Pearson r arrays
        for stat, arr in (("mean", mean_r), ("median", median_r)):
            stem = _stem(stat, "r", run_label, hemi)
            np.save(output_folder / f"{stem}.npy", arr)
            # [CORRECTION #50] columns = MMP parcels (YAML order), rows = seeds
            # BEFORE:             pd.DataFrame(arr, index=macro_regions, columns=macro_regions).to_csv(
            pd.DataFrame(arr, index=macro_regions, columns=parcels).to_csv(
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
            # [CORRECTION #50] column names (MMP parcels, YAML order) added to the metadata
            # BEFORE: (line(s) added)
            parcels           = np.array(parcels),
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