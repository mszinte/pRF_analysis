"""
Created on May 3, 2025

wta_full_corr_by_subject_by_hemi.py
-----------------------------------------------------------------------------------------
Goal:
    Compute winner-take-all (WTA) from per-subject full-correlation TSV files.
    Outputs one combined table per hemisphere × variant containing:
        - one row per subject   (winner seed label per parcel)
        - one GROUP row         (modal winner across subjects per parcel)
        - one CONSISTENCY row   (% of subjects matching the group winner)

    Input TSVs are already parcellated; each file covers one seed × one hemi.
    The atlas-key row order is remapped to canonical YAML parcel order on load.
    All shared logic (constants, I/O, WTA, tie-breaking) is imported from
    rest_utils.py — do not duplicate those definitions here.

Pipeline per hemisphere × variant:
    1. For each subject, resolve which TSV files to load (varies by variant).
    2. Load one TSV per seed via load_full_corr_matrix(), which slices hemi
       rows and remaps to canonical parcel order.
    3. Apply compute_winners() → winning seed label per parcel.
    4. Collect across subjects → append GROUP and CONSISTENCY via
       append_group_and_consistency() with the three-level tie-break cascade
       (votes → Fisher-z → lowest seed number).  Fisher-z tiebreaking requires
       a group-median matrix; pass None here to fall through to level 3.
    5. Save combined table.
-----------------------------------------------------------------------------------------
Run variants (from rest_utils.VARIANTS):
    concat       — concatenated-run TSV, all subjects
    concat_clean — best available run per subject:
                     · RUN02_EXCLUDED subjects → run-01 TSV
                     · all other subjects      → concatenated-run TSV
    run-01       — run-01 TSV, all subjects
    run-02       — run-02 TSV, RUN02_EXCLUDED subjects skipped with WARNING
-----------------------------------------------------------------------------------------
Inputs (sys.argv):
    1: main project directory   (e.g. /scratch/mszinte/data)
    2: project name/directory   (e.g. RetinoMaps)
    3: server group             (e.g. 327)
    4: server project           (e.g. b327)
    5: parcellation mode:
           "default"     → _parcellated.tsv
           "legacy"      → _parcellated_legacy-mode.tsv
           "no_outliers" → _parcellated_no_outliers.tsv

Output:
    Per hemisphere × variant, one CSV:
        winning_seeds_by_subject_full_corr_{hemi}_{variant}_{mode}.csv
    Rows: subjects + GROUP + CONSISTENCY_%; columns: parcel names (YAML order).

To run:
    $ cd projects/pRF_analysis/RetinoMaps/rest/stats
    $ python wta_full_corr_by_subject_by_hemi.py /scratch/mszinte/data RetinoMaps 327 b327 default
-----------------------------------------------------------------------------------------
Written by Marco Bedini (marco.bedini@univ-amu.fr)
-----------------------------------------------------------------------------------------
"""

import warnings
warnings.filterwarnings("ignore")

import os
import sys
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, List, Optional

# ============================================================
# Personal imports — settings utils
# ============================================================
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../"))
sys.path.append(os.path.abspath(os.path.join(base_dir, "analysis_code/utils")))
from settings_utils import load_settings

# ============================================================
# rest_utils — shared pipeline constants and functions
# ============================================================
sys.path.append(os.path.abspath(os.path.join(base_dir, "RetinoMaps/rest/utils")))
from rest_utils import (
    RUN02_EXCLUDED,
    VARIANTS,
    MODE_SUFFIX,
    ATLAS_KEY_TABLES,
    ATLAS_ORDER,
    HEMI_ROW_SLICE,
    N_PARCELS_PER_HEMI,
    N_PARCELS_TOTAL,
    build_remap,
    tsv_path,
    load_full_corr_matrix,
    compute_winners,
    append_group_and_consistency,
)

# ============================================================
# Parse and validate arguments
# ============================================================
USAGE = (
    "Usage: python wta_full_corr_by_subject_by_hemi.py "
    "<main_dir> <project_dir> <group> <server> <mode>\n"
    f"  <mode> must be one of: {', '.join(MODE_SUFFIX)}"
)

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
print("WTA — full correlation (workbench parcellated TSVs)")
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
clusters: List[str]             = list(analysis_info["rois-drawn"])
seed_to_parcels: Dict[str, List[str]] = analysis_info["rois-group-mmp"]
clusters.reverse()                      # mPCS first

parcels: List[str] = []
for cl in clusters:
    parcels.extend(seed_to_parcels[cl])

seed_to_number: Dict[str, int] = {s: i + 1 for i, s in enumerate(clusters)}

n_clusters = len(clusters)
n_parcels  = len(parcels)

# ============================================================
# Paths
# ============================================================
main_data     = Path(main_dir) / project_dir / "derivatives/pp_data"
output_folder = main_data / "group/91k/rest/wta/workbench"
output_folder.mkdir(parents=True, exist_ok=True)

# ============================================================
# Build remap index: atlas-key order → canonical YAML parcel order
#
# Computed once at startup via rest_utils.build_remap().
# Raises immediately if the atlas key tables and the YAML are inconsistent.
# ============================================================
_REMAP: Dict[str, tuple] = {}
for _hemi in ("rh", "lh"):
    _remap_idx, _present = build_remap(ATLAS_ORDER[_hemi], parcels, _hemi)
    _REMAP[_hemi] = (_remap_idx, _present)
    print(
        f"  Remap [{_hemi.upper()}]: {len(_present)}/{n_parcels} canonical parcels "
        "found in atlas key table."
    )

# ============================================================
# Main loop — hemisphere × variant
# ============================================================

for hemi in ("lh", "rh"):
    print(f"\n{'='*80}")
    print(f"Processing hemisphere: {hemi.upper()}")
    print("=" * 80)

    remap_idx, present = _REMAP[hemi]

    for variant, (normal_tag, excluded_tag, skip_excluded) in VARIANTS.items():
        print(f"\n  --- Variant: {variant} ---")

        all_winners: List[np.ndarray] = []
        subject_ids: List[str]        = []

        for subject in subjects:
            is_excluded = subject in RUN02_EXCLUDED

            if is_excluded and skip_excluded:
                print(f"    WARNING [{variant}]: {subject} is in RUN02_EXCLUDED — SKIPPED")
                continue

            run_tag = excluded_tag if is_excluded else normal_tag

            # load_full_corr_matrix() from rest_utils handles TSV loading,
            # hemi-row slicing, and remap to canonical parcel order.
            df_corr = load_full_corr_matrix(
                subject     = subject,
                hemi        = hemi,
                run_tag     = run_tag,
                clusters    = clusters,
                parcels     = parcels,
                remap_idx   = remap_idx,
                present     = present,
                main_data   = main_data,
                tsv_suffix  = tsv_suffix,
            )

            if df_corr is None:
                print(f"    {subject}: SKIPPED (missing files)")
                continue

            if variant == "concat_clean" and is_excluded:
                print(f"    {subject}: OK (fallback → run-01)")
            else:
                print(f"    {subject}: OK")

            # compute_winners() from rest_utils masks self-seed parcels then argmax.
            all_winners.append(
                compute_winners(
                    df_corr         = df_corr,
                    seed_to_parcels = seed_to_parcels,
                    seed_to_number  = seed_to_number,
                )
            )
            subject_ids.append(subject)

        if not all_winners:
            print(f"    ERROR: no valid subjects for {hemi} / {variant} — skipping.")
            continue

        subject_df = pd.DataFrame(all_winners, index=subject_ids, columns=parcels)

        # append_group_and_consistency() from rest_utils computes GROUP via the
        # three-level tie-break cascade (votes → Fisher-z → lowest seed number).
        # group_median=None here: no group-level Fisher-z matrix is available at
        # this stage, so ties fall through to level 3 (lowest seed number).
        combined_df = append_group_and_consistency(
            subject_df     = subject_df,
            clusters       = clusters,
            seed_to_number = seed_to_number,
            group_median   = None,
            parcels        = parcels,
        )

        out_csv = (
            output_folder
            / f"winning_seeds_by_subject_full_corr_{hemi}_{variant}_{mode}.csv"
        )
        combined_df.to_csv(out_csv)
        print(f"    Saved: {out_csv.name}")
        print(f"    Rows : {list(combined_df.index)}")

print("\n" + "=" * 80)
print("ALL HEMISPHERES × VARIANTS COMPLETE")
print("=" * 80)
print(f"\nOutputs written to: {output_folder}")