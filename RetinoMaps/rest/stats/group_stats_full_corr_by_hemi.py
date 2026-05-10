#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

group_full_corr_by_hemi.py
------------------------------------------------------------------------------------------
Goal:
    Compute group-level Fisher-z averages from per-subject, per-ROI, per-hemi
    full-correlation TSV files produced by the workbench parcellation step

    Each input TSV has:
        - 106 rows : parcels in atlas-key order
                     rows  0–52  → right hemisphere  (RH atlas keys 1–163)
                     rows 53–105 → left  hemisphere  (LH atlas keys 181–343)
        - 1 column : Fisher-z correlation of that seed/ROI with each parcel

    For each hemisphere × variant the script:
        1. Resolves which TSV to load for each subject (run variant logic).
        2. Loads all seed TSVs for a subject → slices the hemi-specific rows.
        3. Re-maps the parcel order from atlas-key order to the canonical order
           defined by the YAML config (rois-group-mmp), so columns are
           consistent with the partial-correlation pipeline and the WTA function.
        4. Stacks the (n_seeds × n_parcels) matrices across subjects.
        5. Computes the group median and standard deviation in Fisher-z space.
        6. Saves one .npy + one .csv per statistic per hemisphere × variant.

    Averaging is always in Fisher-z space; Pearson r is recovered only at the
    final reporting stage via tanh().  This is intentional: Fisher-z has
    approximately constant variance ~1/(n-3), making it the correct space for
    averaging and any subsequent parametric tests.  Raw Pearson r values have
    r-dependent variance and must never be averaged directly.

------------------------------------------------------------------------------------------
Run variants:
    concat       — concatenated-run TSV, all subjects
    concat_clean — best available run per subject:
                     · RUN02_EXCLUDED subjects → run-01 TSV
                     · all other subjects      → concatenated-run TSV
    run-01       — run-01 TSV, all subjects
    run-02       — run-02 TSV, all subjects (bad subjects kept intentionally
                   to expose the registration artifact in group plots)

Filename convention (mirrors wta_full_corr_by_subject_by_hemi.py):
    concat / concat_clean : sub-XX_task-rest_space-fsLR_den-91k_desc-fisher-z_...
    run-01 / run-02       : sub-XX_task-rest_run-XX_space-fsLR_den-91k_desc-fisher-z_...
------------------------------------------------------------------------------------------
Inputs (sys.argv):
    1: main project directory   (e.g. /scratch/mszinte/data)
    2: project name/directory   (e.g. RetinoMaps)
    3: server group             (e.g. 327)
    4: server project           (e.g. b327)
    5: parcellation mode:
           "default"     → _parcellated.tsv
           "legacy"      → _parcellated_legacy-mode.tsv
           "no_outliers" → _parcellated_no_outliers.tsv

Outputs (per hemisphere × variant):
    cluster_by_mmp-parcel_full-corr_fisherz_median_{run_tag}_{hemi}_{mode}.npy / .csv
    cluster_by_mmp-parcel_full-corr_fisherz_std_{run_tag}_{hemi}_{mode}.npy  / .csv

    Rows    : seed/cluster names  (n_clusters)
    Columns : parcel names in YAML-derived canonical order  (n_parcels)

To run:
    $ cd projects/pRF_analysis/RetinoMaps/rest/stats
    $ python group_stats_full_corr_by_hemi.py /scratch/mszinte/data RetinoMaps 327 b327 default
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
from typing import Dict, List, Optional, Tuple

warnings.filterwarnings("ignore")

# ============================================================
# Personal imports
# ============================================================
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../"))
sys.path.append(os.path.abspath(os.path.join(base_dir, "analysis_code/utils")))
from settings_utils import load_settings

# ============================================================
# Parcellation mode → TSV filename suffix
# ============================================================
MODE_SUFFIX: Dict[str, str] = {
    "default":     "",
    "legacy":      "_legacy-mode",
    "no_outliers": "_no_outliers",
}

USAGE = (
    "Usage: python group_full_corr_by_hemi.py "
    "<main_dir> <project_dir> <group> <server> <mode>\n"
    f"  <mode> must be one of: {', '.join(MODE_SUFFIX)}"
)

# ============================================================
# Parse and validate arguments
# ============================================================
if len(sys.argv) != 6:
    print(f"ERROR: expected 5 arguments, got {len(sys.argv) - 1}.\n{USAGE}")
    sys.exit(1)

main_dir    = sys.argv[1]
project_dir = sys.argv[2]
group       = sys.argv[3]
server      = sys.argv[4]
mode        = sys.argv[5]

if mode not in MODE_SUFFIX:
    print(f"ERROR: unrecognised mode '{mode}'.\n  Accepted: {', '.join(MODE_SUFFIX)}\n{USAGE}")
    sys.exit(1)

tsv_suffix = MODE_SUFFIX[mode]

print("=" * 80)
print("GROUP FULL CORRELATION — Fisher-z averaging (workbench parcellated TSVs)")
print("=" * 80)
print(f"  main_dir    : {main_dir}")
print(f"  project_dir : {project_dir}")
print(f"  group       : {group}")
print(f"  server      : {server}")
print(f"  mode        : {mode!r}  →  suffix: '_parcellated{tsv_suffix}.tsv'")

# ============================================================
# Load settings
# ============================================================
settings_path     = os.path.join(base_dir, project_dir, "settings.yml")
prf_settings_path = os.path.join(base_dir, project_dir, "prf-analysis.yml")
settings          = load_settings([settings_path, prf_settings_path])
analysis_info     = settings[0]
subjects          = analysis_info["subjects"]

# ============================================================
# ROIs — canonical order from YAML config
# ============================================================
clusters        = list(analysis_info["rois-drawn"])
seed_to_parcels = analysis_info["rois-group-mmp"]   # {seed_name: [parcel_names]}
clusters.reverse()                                   # mPCS first, consistent with WTA script

# Canonical parcel list: flatten in cluster order (same as partial-corr pipeline)
parcels: List[str] = []
for cl in clusters:
    parcels.extend(seed_to_parcels[cl])

n_clusters = len(clusters)
n_parcels  = len(parcels)

print(f"\n  Seeds (n={n_clusters}): {clusters}")
print(f"  Parcels (n={n_parcels}): {parcels[:5]} ... {parcels[-5:]}")

# ============================================================
# Paths
# ============================================================
main_data     = Path(main_dir) / project_dir / "derivatives/pp_data"
output_folder = main_data / "group/91k/rest/full_corr/by_hemi"
output_folder.mkdir(parents=True, exist_ok=True)

# ============================================================
# Run variants and subject exclusions
#
# RUN02_EXCLUDED: subjects with bad data / registration error in run-02.
#
# run-02 keeps ALL subjects intentionally — the bad subjects are left in
# to expose the registration artifact in group plots (QC strategy).
# concat_clean falls back to run-01 for excluded subjects.
# ============================================================
RUN02_EXCLUDED: frozenset = frozenset(
    {"sub-03", "sub-04", "sub-14", "sub-21", "sub-22", "sub-23"}
)

# Variant table:
#   normal_tag   : run_tag for non-excluded subjects  (None = concatenated file)
#   excluded_tag : run_tag for RUN02_EXCLUDED subjects (None = concatenated file)
#   note         : human-readable description for log output
VARIANTS: Dict[str, Dict] = {
    "concat": {
        "normal_tag":   None,
        "excluded_tag": None,
        "note":         "concatenated run, all subjects",
    },
    "concat_clean": {
        "normal_tag":   None,
        "excluded_tag": "run-01",
        "note":         "concatenated for good subjects, run-01 fallback for excluded",
    },
    "run-01": {
        "normal_tag":   "run-01",
        "excluded_tag": "run-01",
        "note":         "run-01, all subjects",
    },
    "run-02": {
        "normal_tag":   "run-02",
        "excluded_tag": "run-02",
        "note":         "run-02, all subjects (bad subjects kept for QC)",
    },
}

# ============================================================
# Atlas parcel keys — define the row order inside the TSV files.
#
# Each TSV has 106 rows total:
#   rows  0–52  → RH parcels, sorted by atlas key (1–163)
#   rows 53–105 → LH parcels, sorted by atlas key (181–343)
#
# These dicts map parcel name → atlas key for each hemisphere.
# Sorted by atlas key value they give the exact row ordering in the TSV.
# ============================================================
_R_KEYS: Dict[str, int] = {
    "V1":1,    "MST":2,   "V2":4,    "V3":5,    "V4":6,
    "V8":7,    "FEF":10,  "PEF":11,  "55b":12,  "V3A":13,
    "V7":16,   "IPS1":17, "FFC":18,  "V3B":19,  "LO1":20,
    "LO2":21,  "PIT":22,  "MT":23,   "7Pm":29,  "24dv":41,
    "7AL":42,  "SCEF":43, "6ma":44,  "7Am":45,  "7PL":46,
    "7PC":47,  "LIPv":48, "VIP":49,  "MIP":50,  "6d":54,
    "6mp":55,  "6v":56,   "p32pr":60,"6r":78,   "IFJa":79,
    "IFJp":80, "LIPd":95, "6a":96,   "i6-8":97, "AIP":117,
    "PH":138,  "IP2":144, "IP1":145, "IP0":146, "V6A":152,
    "VMV1":153,"VMV3":154,"V4t":156, "FST":157, "V3CD":158,
    "LO3":159, "VMV2":160,"VVC":163,
}

_L_KEYS: Dict[str, int] = {
    "V1":181,  "MST":182, "V2":184,  "V3":185,  "V4":186,
    "V8":187,  "FEF":190, "PEF":191, "55b":192, "V3A":193,
    "V7":196,  "IPS1":197,"FFC":198, "V3B":199, "LO1":200,
    "LO2":201, "PIT":202, "MT":203,  "7Pm":209, "24dv":221,
    "7AL":222, "SCEF":223,"6ma":224, "7Am":225, "7PL":226,
    "7PC":227, "LIPv":228,"VIP":229, "MIP":230, "6d":234,
    "6mp":235, "6v":236,  "p32pr":240,"6r":258, "IFJa":259,
    "IFJp":260,"LIPd":275,"6a":276,  "i6-8":277,"AIP":297,
    "PH":318,  "IP2":324, "IP1":325, "IP0":326, "V6A":332,
    "VMV1":333,"VMV3":334,"V4t":336, "FST":337, "V3CD":338,
    "LO3":339, "VMV2":340,"VVC":343,
}

N_PARCELS_PER_HEMI = 53
N_PARCELS_TOTAL    = 106

# Atlas-key order: the actual row ordering inside the TSV
_RH_ATLAS_ORDER: List[str] = [
    n for n, _ in sorted(_R_KEYS.items(), key=lambda x: x[1])
]
_LH_ATLAS_ORDER: List[str] = [
    n for n, _ in sorted(_L_KEYS.items(), key=lambda x: x[1])
]


# Row slices within the 106-row TSV
_HEMI_ROW_SLICE: Dict[str, slice] = {
    "rh": slice(0, N_PARCELS_PER_HEMI),
    "lh": slice(N_PARCELS_PER_HEMI, N_PARCELS_TOTAL),
}

# Map each hemisphere to its atlas-key-sorted parcel name list
_HEMI_ATLAS_ORDER: Dict[str, List[str]] = {
    "rh": _RH_ATLAS_ORDER,
    "lh": _LH_ATLAS_ORDER,
}

# ============================================================
# Build re-mapping index: atlas-key order → canonical YAML order
#
# The TSV rows arrive in atlas-key order (_HEMI_ATLAS_ORDER).
# The WTA function and the partial-corr pipeline expect columns in
# canonical YAML order (parcels list built from rois-group-mmp).
#
# For each hemisphere we compute a boolean mask and an integer index
# array so that:
#   tsv_values[atlas_order][remap_idx]  → values in canonical order
#
# Parcels present in the TSV but absent from the YAML (or vice-versa)
# are detected here and reported — they must never be silently dropped.
# ============================================================
def _build_remap(
    atlas_order: List[str],
    canonical_order: List[str],
    hemi: str,
) -> Tuple[List[int], List[str]]:
    """
    Return (remap_idx, canon_parcels_in_tsv) where:
        remap_idx            : indices into atlas_order that give canonical_order
        canon_parcels_in_tsv : subset of canonical_order that is also in the TSV

    A WARNING is printed for any parcel in canonical_order that is absent from
    the TSV (these are filled with NaN in the output matrix).
    An ERROR is raised if the TSV contains parcels absent from canonical_order,
    because that indicates a mismatch between the atlas key tables and the YAML.
    """
    atlas_set     = set(atlas_order)
    canonical_set = set(canonical_order)

    extra_in_tsv  = atlas_set - canonical_set
    missing_in_tsv = canonical_set - atlas_set

    if extra_in_tsv:
        raise ValueError(
            f"[{hemi.upper()}] TSV atlas-key table contains parcels not in YAML "
            f"rois-group-mmp: {sorted(extra_in_tsv)}\n"
            "  → Update the atlas key dicts or the YAML."
        )
    if missing_in_tsv:
        print(
            f"  WARNING [{hemi.upper()}]: {len(missing_in_tsv)} YAML parcel(s) absent "
            f"from TSV atlas-key table (will be NaN): {sorted(missing_in_tsv)}"
        )

    atlas_index = {name: i for i, name in enumerate(atlas_order)}
    remap_idx: List[int]   = []
    present:   List[str]   = []
    for name in canonical_order:
        if name in atlas_index:
            remap_idx.append(atlas_index[name])
            present.append(name)

    return remap_idx, present


# Pre-compute remap indices for both hemispheres once
_REMAP: Dict[str, Tuple[List[int], List[str]]] = {}
for _hemi in ("rh", "lh"):
    _remap_idx, _present = _build_remap(
        _HEMI_ATLAS_ORDER[_hemi],
        parcels,
        _hemi,
    )
    _REMAP[_hemi] = (_remap_idx, _present)
    print(
        f"  Remap [{_hemi.upper()}]: {len(_present)}/{n_parcels} canonical parcels "
        f"found in TSV atlas-key table."
    )

# ============================================================
# Helper functions
# ============================================================

def tsv_path(
    subject: str,
    hemi: str,
    roi: str,
    run_tag: Optional[str],
) -> Path:
    """
    Return the expected path for one subject / hemi / ROI TSV file.

    run_tag = None → concatenated file (no run entity in filename):
        sub-XX_task-rest_space-fsLR_den-91k_desc-fisher-z_{hemi}_{roi}_parcellated{suffix}.tsv
    run_tag = str  → per-run file (run entity inserted after task-rest):
        sub-XX_task-rest_{run_tag}_space-fsLR_den-91k_desc-fisher-z_{hemi}_{roi}_parcellated{suffix}.tsv
    """
    subj_dir   = main_data / subject / "91k/rest/corr/full_corr/by_hemi"
    run_entity = f"_{run_tag}" if run_tag is not None else ""
    fname = (
        f"{subject}_task-rest{run_entity}_space-fsLR_den-91k"
        f"_desc-fisher-z_{hemi}_{roi}"
        f"_parcellated{tsv_suffix}.tsv"
    )
    return subj_dir / fname


def load_subject_matrix(
    subject: str,
    hemi: str,
    run_tag: Optional[str],
) -> Optional[np.ndarray]:
    """
    Load all seed TSVs for one subject/hemi/run_tag and return a
    (n_clusters × n_parcels) array in canonical YAML parcel order.

    Steps:
        1. For each seed, load its TSV and slice the hemi-specific rows
           (53 of 106) — values are already in Fisher-z space.
        2. Re-order the 53 atlas-key-sorted values to canonical YAML order
           using the pre-computed remap index (_REMAP).
        3. Stack seeds → (n_clusters × n_parcels) array.
           Parcels absent from the TSV are filled with NaN.

    Returns None if any seed file is missing; all missing paths are reported
    before returning so the caller can skip the subject.
    """
    row_slice              = _HEMI_ROW_SLICE[hemi]
    atlas_order            = _HEMI_ATLAS_ORDER[hemi]
    remap_idx, present     = _REMAP[hemi]

    # Pre-allocate output with NaN — missing/absent parcels stay NaN
    matrix = np.full((n_clusters, n_parcels), np.nan, dtype=float)

    # Column positions in the output that correspond to TSV-present parcels
    present_col_idx = [parcels.index(p) for p in present]

    missing_files: List[Path] = []

    for i_seed, seed in enumerate(clusters):
        fpath = tsv_path(subject, hemi, seed, run_tag)

        if not fpath.exists():
            missing_files.append(fpath)
            continue

        raw = pd.read_csv(fpath, header=None, sep="\t")

        if raw.shape != (N_PARCELS_TOTAL, 1):
            raise ValueError(
                f"[{subject} {hemi}] Unexpected shape {raw.shape} in {fpath.name} "
                f"(expected ({N_PARCELS_TOTAL}, 1))."
            )

        # Slice hemi rows → 1-D array in atlas-key order
        hemi_values: np.ndarray = raw.iloc[row_slice, 0].values.astype(float)

        assert hemi_values.shape == (N_PARCELS_PER_HEMI,), (
            f"[{subject} {hemi} {seed}] Unexpected slice length {hemi_values.shape}."
        )

        # Check for infinite values (|r|=1 edge case → arctanh(±1) = ±inf)
        n_inf = int(np.isinf(hemi_values).sum())
        if n_inf:
            print(
                f"  WARNING [{subject} {hemi} {seed}]: {n_inf} infinite Fisher-z "
                f"value(s) (|r|=1) — will propagate as NaN in group median."
            )
            hemi_values[np.isinf(hemi_values)] = np.nan

        # Re-map from atlas-key order to canonical YAML order
        # hemi_values[remap_idx] gives values at the TSV-present parcel positions,
        # in the same order as `present` (and therefore `present_col_idx`).
        matrix[i_seed, present_col_idx] = hemi_values[remap_idx]

    if missing_files:
        for f in missing_files:
            print(f"  WARNING [{subject} {hemi}]: missing {f.name}")
        return None

    return matrix

# ============================================================
# Main loop — hemisphere × variant
# ============================================================

for hemi in ("lh", "rh"):
    print(f"\n{'='*80}")
    print(f"Processing hemisphere: {hemi.upper()}")
    print("=" * 80)

    for variant, vcfg in VARIANTS.items():
        normal_tag   = vcfg["normal_tag"]
        excluded_tag = vcfg["excluded_tag"]
        note         = vcfg["note"]

        print(f"\n  --- Variant: {variant}  ({note}) ---")

        subject_matrices: List[np.ndarray] = []
        subject_ids:      List[str]        = []

        for subject in subjects:
            is_excluded = subject in RUN02_EXCLUDED
            run_tag     = excluded_tag if is_excluded else normal_tag

            mat = load_subject_matrix(subject, hemi, run_tag)

            if mat is None:
                print(f"    {subject}: SKIPPED (missing files)")
                continue

            # Defensive shape check before appending
            assert mat.shape == (n_clusters, n_parcels), (
                f"[{subject} {hemi} {variant}] Unexpected matrix shape {mat.shape}, "
                f"expected ({n_clusters}, {n_parcels})."
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

        # Stack to (n_subjects × n_clusters × n_parcels)
        stack = np.stack(subject_matrices, axis=0)
        assert stack.shape == (n_valid, n_clusters, n_parcels), (
            f"Unexpected stack shape {stack.shape}."
        )

        # Group statistics in Fisher-z space (nanmedian / nanstd ignore NaN cells)
        group_median = np.nanmedian(stack, axis=0)   # (n_clusters × n_parcels)
        group_std  = np.nanstd( stack, axis=0)   # (n_clusters × n_parcels)

        n_nan_median = int(np.isnan(group_median).sum())
        n_nan_std  = int(np.isnan(group_std).sum())
        if n_nan_median:
            print(f"    WARNING: {n_nan_median}/{group_median.size} NaN cells in group median.")
        if n_nan_std:
            print(f"    WARNING: {n_nan_std}/{group_std.size} NaN cells in group std.")

        print(f"    Group median range : [{np.nanmin(group_median):.4f}, {np.nanmax(group_median):.4f}]")
        print(f"    Group std  range : [{np.nanmin(group_std):.4f},  {np.nanmax(group_std):.4f}]")

        # Build variant run-tag label for filename (None → "concat")
        run_label = normal_tag if normal_tag is not None else variant

        # Save group median
        stem_median = (
            f"cluster_by_mmp-parcel_full-corr_fisherz_median"
            f"_{run_label}_{hemi}_{mode}"
        )
        np.save(output_folder / f"{stem_median}.npy", group_median)
        pd.DataFrame(group_median, index=clusters, columns=parcels).to_csv(
            output_folder / f"{stem_median}.csv"
        )

        # Save group std
        stem_std = (
            f"cluster_by_mmp-parcel_full-corr_fisherz_std"
            f"_{run_label}_{hemi}_{mode}"
        )
        np.save(output_folder / f"{stem_std}.npy", group_std)
        pd.DataFrame(group_std, index=clusters, columns=parcels).to_csv(
            output_folder / f"{stem_std}.csv"
        )

        print(f"    Saved: {stem_median}.npy / .csv")
        print(f"    Saved: {stem_std}.npy / .csv")

print("\n" + "=" * 80)
print("ALL HEMISPHERES × VARIANTS COMPLETE")
print("=" * 80)
print(f"\nOutputs written to: {output_folder}")
print(
    "\nNote: outputs are in Fisher-z space. Apply np.tanh() to recover Pearson r "
    "only at the final reporting or plotting stage."
)