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

Pipeline per hemisphere × variant:
    1. For each subject, resolve which TSV files to load (varies by variant).
    2. Load one TSV per seed, slice hemi rows, remap to canonical parcel order.
    3. Stack into a (n_seeds × n_parcels) DataFrame.
    4. Apply compute_winners() → winning seed label per parcel.
    5. Collect across subjects → compute GROUP (mode) and CONSISTENCY rows.
    6. Save combined table.
-----------------------------------------------------------------------------------------
Run variants:
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
from typing import Dict, List, Optional, Tuple
from scipy import stats as scipy_stats

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
    "Usage: python wta_full_corr_by_subject_by_hemi.py "
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
clusters        = list(analysis_info["rois-drawn"])
seed_to_parcels = analysis_info["rois-group-mmp"]   # {seed_name: [parcel_names]}
clusters.reverse()                                   # mPCS first

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
# Run variants
#
# Subjects excluded from run-02 (bad data / registration error).
# For concat_clean these subjects fall back to their run-01 file instead.
# run-02 skips excluded subjects with a WARNING (TSV data is unusable for them).
# ============================================================
RUN02_EXCLUDED: set = {"sub-03", "sub-04", "sub-14", "sub-21", "sub-22", "sub-23"}

# variant → (normal_tag, excluded_tag, skip_excluded)
#   run_tag = None → concatenated file (no run-XX in filename)
#   run_tag = str  → per-run file (run-XX inserted after task-rest)
#   skip_excluded  → excluded subjects skipped with WARNING (run-02 only)
VARIANTS: Dict[str, tuple] = {
    "concat":       (None,     None,     False),
    "concat_clean": (None,     "run-01", False),
    "run-01":       ("run-01", "run-01", False),
    "run-02":       ("run-02", None,     True),
}

# ============================================================
# Atlas parcel keys — define the row order inside the TSV files.
#
# Each TSV has 106 rows total:
#   rows  0–52  → RH parcels sorted by atlas key (1–163)
#   rows 53–105 → LH parcels sorted by atlas key (181–343)
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

_RH_ATLAS_ORDER: List[str] = [
    n for n, _ in sorted(_R_KEYS.items(), key=lambda x: x[1])
]
_LH_ATLAS_ORDER: List[str] = [
    n for n, _ in sorted(_L_KEYS.items(), key=lambda x: x[1])
]

assert len(_RH_ATLAS_ORDER) == N_PARCELS_PER_HEMI
assert len(_LH_ATLAS_ORDER) == N_PARCELS_PER_HEMI
assert _RH_ATLAS_ORDER == _LH_ATLAS_ORDER, (
    "Parcel name order differs between hemispheres — check atlas key dicts."
)

_HEMI_ROW_SLICE: Dict[str, slice] = {
    "rh": slice(0, N_PARCELS_PER_HEMI),
    "lh": slice(N_PARCELS_PER_HEMI, N_PARCELS_TOTAL),
}
_HEMI_ATLAS_ORDER: Dict[str, List[str]] = {
    "rh": _RH_ATLAS_ORDER,
    "lh": _LH_ATLAS_ORDER,
}

# ============================================================
# Build remap index: atlas-key order → canonical YAML parcel order
#
# Computed once at startup; used in load_corr_matrix() for every subject.
# Raises immediately if the atlas key tables and the YAML are inconsistent.
# ============================================================
def _build_remap(
    atlas_order: List[str],
    canonical_order: List[str],
    hemi: str,
) -> Tuple[List[int], List[str]]:
    """
    Return (remap_idx, present) where:
        remap_idx : indices into atlas_order that select and reorder values
                    to match canonical_order
        present   : the subset of canonical_order that exists in the TSV
    """
    atlas_set     = set(atlas_order)
    canonical_set = set(canonical_order)

    extra_in_tsv = atlas_set - canonical_set
    if extra_in_tsv:
        raise ValueError(
            f"[{hemi.upper()}] Atlas key table contains parcels absent from YAML "
            f"rois-group-mmp: {sorted(extra_in_tsv)}\n"
            "  → Update the atlas key dicts or the YAML."
        )

    missing_in_tsv = canonical_set - atlas_set
    if missing_in_tsv:
        print(
            f"  WARNING [{hemi.upper()}]: {len(missing_in_tsv)} YAML parcel(s) absent "
            f"from atlas key table (will be NaN in output): {sorted(missing_in_tsv)}"
        )

    atlas_index = {name: i for i, name in enumerate(atlas_order)}
    remap_idx: List[int] = []
    present:   List[str] = []
    for name in canonical_order:
        if name in atlas_index:
            remap_idx.append(atlas_index[name])
            present.append(name)

    return remap_idx, present


_REMAP: Dict[str, Tuple[List[int], List[str]]] = {}
for _hemi in ("rh", "lh"):
    _remap_idx, _present = _build_remap(_HEMI_ATLAS_ORDER[_hemi], parcels, _hemi)
    _REMAP[_hemi] = (_remap_idx, _present)
    print(
        f"  Remap [{_hemi.upper()}]: {len(_present)}/{n_parcels} canonical parcels "
        "found in atlas key table."
    )

# ============================================================
# Helpers
# ============================================================

def tsv_path(subject: str, hemi: str, roi: str, run_tag: Optional[str]) -> Path:
    """
    Return expected path for one subject / hemi / seed TSV file.

    run_tag = None → concatenated file (no run entity):
        sub-XX_task-rest_space-fsLR_den-91k_desc-fisher-z_{hemi}_{roi}_parcellated{suffix}.tsv
    run_tag = str  → per-run file (run entity after task-rest):
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


def load_corr_matrix(
    subject: str,
    hemi: str,
    run_tag: Optional[str],
) -> Optional[pd.DataFrame]:
    """
    Load one TSV per seed, slice hemi rows, remap to canonical YAML parcel
    order, and return a (n_seeds × n_parcels) DataFrame.

    Returns None if any seed file is missing (all missing paths reported).
    Raises ValueError on unexpected TSV shape.
    """
    row_slice          = _HEMI_ROW_SLICE[hemi]
    remap_idx, present = _REMAP[hemi]
    present_col_idx    = [parcels.index(p) for p in present]

    seed_rows: Dict[str, np.ndarray] = {}
    missing:   List[Path]            = []

    for seed in clusters:
        fpath = tsv_path(subject, hemi, seed, run_tag)

        if not fpath.exists():
            missing.append(fpath)
            continue

        raw = pd.read_csv(fpath, header=None, sep="\t")

        if raw.shape != (N_PARCELS_TOTAL, 1):
            raise ValueError(
                f"[{subject} {hemi}] Unexpected shape {raw.shape} in {fpath.name} "
                f"(expected ({N_PARCELS_TOTAL}, 1))."
            )

        hemi_values = raw.iloc[row_slice, 0].values.astype(float)

        # Remap from atlas-key order to canonical YAML parcel order
        row = np.full(n_parcels, np.nan, dtype=float)
        row[present_col_idx] = hemi_values[remap_idx]
        seed_rows[seed] = row

    if missing:
        for f in missing:
            print(f"  WARNING [{subject} {hemi}]: missing {f.name}")
        return None

    df = pd.DataFrame(seed_rows, index=parcels).T   # (n_seeds × n_parcels)
    df.index.name   = None
    df.columns.name = None

    assert df.shape == (n_clusters, n_parcels), (
        f"[{subject} {hemi}] Unexpected matrix shape {df.shape}."
    )
    return df


def compute_winners(df_corr: pd.DataFrame) -> np.ndarray:
    """
    Return a 1-D array of winning seed labels (1-based int) per parcel.
    Self-seed parcels are masked to NaN before argmax.
    Parcels where all seeds are NaN receive NaN.
    """
    df = df_corr.copy()

    for seed, plist in seed_to_parcels.items():
        for p in plist:
            if seed in df.index and p in df.columns:
                df.loc[seed, p] = np.nan

    winners = []
    for parcel in df.columns:
        col = df[parcel]
        # col.isna().all() guard is essential — col.idxmax() raises on all-NaN
        winners.append(
            np.nan if col.isna().all() else seed_to_number[col.idxmax()]
        )

    return np.array(winners)


def append_group_and_consistency(subject_df: pd.DataFrame) -> pd.DataFrame:
    """
    Append GROUP and CONSISTENCY rows to the subject winner table.

    GROUP row:
        Modal winner per parcel across subjects (NaN ignored).
        Ties broken by lowest seed number (scipy default).

    CONSISTENCY row:
        % of subjects whose winner matches GROUP (NaN where GROUP is NaN).
    """
    group_winners = []
    for parcel in subject_df.columns:
        col = subject_df[parcel].dropna()
        if col.empty:
            group_winners.append(np.nan)
        else:
            group_winners.append(float(scipy_stats.mode(col, keepdims=True).mode[0]))

    group_series = pd.Series(group_winners, index=subject_df.columns, name="GROUP")

    consistency = []
    for parcel in subject_df.columns:
        gw = group_series[parcel]
        if pd.isna(gw):
            consistency.append(np.nan)
        else:
            consistency.append(
                100.0 * (subject_df[parcel] == gw).sum() / len(subject_df)
            )

    consistency_series = pd.Series(
        consistency, index=subject_df.columns, name="CONSISTENCY_%"
    )

    return pd.concat(
        [subject_df, group_series.to_frame().T, consistency_series.to_frame().T]
    )

# ============================================================
# Main loop — hemisphere × variant
# ============================================================

for hemi in ("lh", "rh"):
    print(f"\n{'='*80}")
    print(f"Processing hemisphere: {hemi.upper()}")
    print("=" * 80)

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

            df_corr = load_corr_matrix(subject, hemi, run_tag)
            if df_corr is None:
                print(f"    {subject}: SKIPPED (missing files)")
                continue

            if variant == "concat_clean" and is_excluded:
                print(f"    {subject}: OK (fallback → run-01)")
            else:
                print(f"    {subject}: OK")

            all_winners.append(compute_winners(df_corr))
            subject_ids.append(subject)

        if not all_winners:
            print(f"    ERROR: no valid subjects for {hemi} / {variant} — skipping.")
            continue

        subject_df  = pd.DataFrame(all_winners, index=subject_ids, columns=parcels)
        combined_df = append_group_and_consistency(subject_df)

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