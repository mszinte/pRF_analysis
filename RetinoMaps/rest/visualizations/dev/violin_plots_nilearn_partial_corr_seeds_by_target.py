#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Goal of the script:
Partial correlation violin plots split by hemisphere.
Seeds are on the y-axis, one figure per TARGET cluster.
This makes it easy to compare correlation strength across seeds
for a given target (winner-take-all inspection).

The target's own row is muted (self-correlation), preserving
the same logic as the seed-figure version.

Produces figures for four data variants:
  - concat        : both runs concatenated, all subjects
  - concat_clean  : both runs concatenated, excluding subjects with bad run-02
  - run-01        : run-01 only, all subjects
  - run-02        : run-02 only, ALL subjects (including bad ones, to show artifact)

Missing files are always skipped with a warning; the script runs to completion.

---------------------------------------------------
Written by Marco Bedini (marco.bedini@univ-amu.fr)
---------------------------------------------------
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
import pandas as pd
from pathlib import Path
from typing import Optional

# ============================================================
# Paths
# ============================================================
USER = os.environ["USER"]

main_data     = Path("/scratch/mszinte/data/RetinoMaps/derivatives/pp_data")
output_folder = main_data / "group/91k/rest/partial_corr"
output_folder.mkdir(parents=True, exist_ok=True)

base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../"))
sys.path.append(os.path.abspath(os.path.join(base_dir, "analysis_code/utils")))
from settings_utils import load_settings

project_dir       = "RetinoMaps"
settings_path     = os.path.join(base_dir, project_dir, "settings.yml")
prf_settings_path = os.path.join(base_dir, project_dir, "prf-analysis.yml")
settings          = load_settings([settings_path, prf_settings_path])
analysis_info     = settings[0]
subjects          = analysis_info["subjects"]

fig_dir = os.path.join(output_folder, "figures/violin_plots")
os.makedirs(fig_dir, exist_ok=True)

# ============================================================
# Run variants & exclusions
# ============================================================

# Subjects excluded from run-02 (bad data / registration error)
RUN02_EXCLUDED = {"sub-03", "sub-04", "sub-14", "sub-21", "sub-22", "sub-23"}

# Four variants: (run_tag_in_filename, subject_filter_fn)
# run_tag=None -> concatenated file (no run-XX in filename)
# run-02 uses _all intentionally: bad subjects kept to show the artifact
def _all(sub):      return True
def _no_run02(sub): return sub not in RUN02_EXCLUDED

VARIANTS = {
    "concat":       (None,     _all),
    "concat_clean": (None,     _no_run02),
    "run-01":       ("run-01", _all),
    "run-02":       ("run-02", _all),
}

# =========================
# Cluster <-> parcel mapping
# =========================

clusters_all = analysis_info["rois-drawn"]
clusters_all.reverse()   # mPCS, sPCS, iPCS, sIPS, iIPS, hMT+, ...

cluster_to_parcels = analysis_info["rois-group-mmp"]

seed_clusters   = clusters_all[:5]
target_clusters = seed_clusters.copy()

for cl in seed_clusters:
    if cl not in cluster_to_parcels:
        raise KeyError(
            "Seed cluster '{}' has no entry in 'rois-group-mmp' -- "
            "check settings.yml".format(cl)
        )
    if len(cluster_to_parcels[cl]) == 0:
        raise ValueError("Seed cluster '{}' maps to an empty parcel list".format(cl))

# =============================
# Colors
# =============================

macro_colors = {
    "mPCS": "#FF6F00",
    "sPCS": "#FFEA00",
    "iPCS": "#97FF00",
    "sIPS": "#2CFF96",
    "iIPS": "#0098FF",
}

SHADE_FACTORS = {"LH": 0.75, "RH": 1.25}
SEED_COLOR    = "#F5F5F5"

def derive_hemi_color(base_hex, hemi):
    # type: (str, str) -> np.ndarray
    """Return RGB array for hemi ('LH'/'RH') derived from base_hex.

    LH darkened (factor < 1), RH lightened (factor > 1), clipped to [0, 1].
    """
    factor = SHADE_FACTORS[hemi]
    rgb    = np.array(mcolors.to_rgb(base_hex))
    return np.clip(rgb * factor, 0, 1)

# =========================
# Build parcel indices
# =========================

clusters = clusters_all   # alias -- row order matches .npy matrix files

parcels = []
for cl in clusters_all:
    parcels.extend(cluster_to_parcels[cl])

cluster_to_parcel_idx = {
    cl: [parcels.index(p) for p in plist]
    for cl, plist in cluster_to_parcels.items()
}

print("\nParcels from YAML: {} | first: {}".format(len(parcels), parcels[0]))
print("Clusters from YAML: {}".format(clusters))
print("\nCluster -> parcel indices:")
for cl in seed_clusters:
    print("  {:5s} -> {}".format(cl, [parcels[i] for i in cluster_to_parcel_idx[cl]]))
    if len(cluster_to_parcel_idx[cl]) == 0:
        raise RuntimeError("No parcels found for cluster '{}'".format(cl))

# =========================
# .npy loader (robust)
# =========================

def load_npy_hemi(filepath):
    # type: (str) -> Optional[np.ndarray]
    """Load a partial-corr .npy matrix, returning None on any problem.

    Returns None (with a warning) on file-not-found or shape/read errors
    so the caller can skip gracefully without crashing.
    """
    if not os.path.isfile(filepath):
        print("  WARNING: missing -- {}".format(filepath))
        return None
    try:
        arr = np.load(filepath)
    except Exception as e:
        print("  WARNING: could not load {} -- {}".format(filepath, e))
        return None
    if arr.ndim != 2:
        print("  WARNING: expected 2D array, got shape {} in {} -- skipping".format(
            arr.shape, filepath))
        return None
    return arr

# =========================
# Load data for all variants
# =========================

hemis = ["lh", "rh"]

# results[variant][seed][tc][hemi] = list of per-subject values
results = {
    variant: {
        seed: {
            tc: {hemi: [] for hemi in hemis}
            for tc in target_clusters if tc != seed
        }
        for seed in seed_clusters
    }
    for variant in VARIANTS
}

missing_files = []   # global log -- reported at end

# Loading loop — replace the existing one in the partial corr script
for variant, (run_tag, subject_filter) in VARIANTS.items():
 
    variant_subjects = [s for s in subjects if subject_filter(s)]
    print("\n=== Loading variant: {} | {} subjects ===".format(
        variant, len(variant_subjects)))
 
    for sub in variant_subjects:
 
        sub_path = os.path.join(
            main_data,
            "{}/91k/rest/corr/partial_corr/by_hemi".format(sub)
        )
 
        for hemi in hemis:
 
            fname = get_fname(sub_path, sub, variant, run_tag, hemi)
 
            cluster_partial = load_npy_hemi(fname)
 
            if cluster_partial is None:
                missing_files.append((variant, fname))
                continue
 
            print("  {} ({}, {}): {}".format(sub, variant, hemi,
                                              cluster_partial.shape))
 
            for seed in seed_clusters:
                seed_idx = clusters.index(seed)
 
                for tc in target_clusters:
                    if tc == seed:
                        continue
 
                    parcel_idx = cluster_to_parcel_idx[tc]
                    vals = cluster_partial[seed_idx, parcel_idx]
                    results[variant][seed][tc][hemi].append(float(np.nanmean(vals)))


# Summary of all missing files
if missing_files:
    print("\n" + "=" * 60)
    print("WARNING: {} file(s) could not be loaded:".format(len(missing_files)))
    for v, f in missing_files:
        print("  [{}] {}".format(v, f))
    print("=" * 60 + "\n")

# =========================
# Plotting helper
# =========================

MEDIAN_COLORS   = {"lh": "#222222", "rh": "#888888"}
MEDIAN_HALF_LEN = 0.08


def plot_target_figure(target, variant, res):
    # type: (str, str, dict) -> None
    """Draw and save one violin figure for target x variant."""

    target_row_idx = seed_clusters.index(target)

    plot_rows = []
    for seed in seed_clusters:
        if seed == target:
            for hemi_label in ["LH", "RH"]:
                plot_rows.append({
                    "Seed":        seed,
                    "Correlation": np.nan,
                    "Hemisphere":  hemi_label,
                })
            continue
        for hemi in hemis:
            for val in res[seed][target][hemi]:
                plot_rows.append({
                    "Seed":        seed,
                    "Correlation": val,
                    "Hemisphere":  hemi.upper(),
                })

    df = pd.DataFrame(plot_rows)

    # Skip entirely if no data at all for this target x variant
    if df.loc[df["Seed"] != target, "Correlation"].dropna().empty:
        print("  SKIP (no data): target={}, variant={}".format(target, variant))
        return

    fig, ax = plt.subplots(figsize=(8, 6))
    n_before = len(ax.collections)

    sns.violinplot(
        data=df,
        y="Seed",
        x="Correlation",
        hue="Hemisphere",
        order=seed_clusters,
        split=True,
        inner=None,
        linewidth=1,
        palette={"LH": "white", "RH": "white"},
        ax=ax,
        legend=False,
    )

    new_collections = ax.collections[n_before:]

    # Validate collection count -- seaborn skips all-NaN rows (the target row)
    n_plottable_rows     = len(seed_clusters) - 1
    expected_collections = n_plottable_rows * 2

    assert len(new_collections) == expected_collections, (
        "[{}, target={}] Expected {} PolyCollections but got {}. "
        "seaborn layout may have changed -- check version.".format(
            variant, target, expected_collections, len(new_collections))
    )

    # Colour violin halves by seed
    collection_idx = 0
    for seed in seed_clusters:
        if seed == target:
            continue
        for hemi_label in ["LH", "RH"]:
            artist = new_collections[collection_idx]
            rgb    = derive_hemi_color(macro_colors[seed], hemi_label)
            artist.set_facecolor(rgb)
            artist.set_edgecolor(rgb)
            artist.set_alpha(0.85)
            collection_idx += 1

    # Muted band for target row
    ax.axhspan(
        target_row_idx - 0.45,
        target_row_idx + 0.45,
        color=SEED_COLOR,
        zorder=0,
        label="_nolegend_",
    )

    # Median lines
    for i, seed in enumerate(seed_clusters):
        if seed == target:
            continue
        for hemi, y_offset in zip(hemis, [-0.12, 0.12]):
            vals = np.array(res[seed][target][hemi])
            if len(vals) == 0:
                continue
            ax.plot(
                [np.nanmedian(vals)] * 2,
                [i + y_offset - MEDIAN_HALF_LEN,
                 i + y_offset + MEDIAN_HALF_LEN],
                color=MEDIAN_COLORS[hemi],
                linewidth=1.5,
                zorder=4,
            )

    # Subject dots
    rng = np.random.default_rng(seed=42)
    for i, seed in enumerate(seed_clusters):
        if seed == target:
            continue
        for hemi, y_offset in zip(hemis, [-0.1, 0.1]):
            x = np.array(res[seed][target][hemi])
            if len(x) == 0:
                continue
            y_jitter = (
                np.ones(len(x)) * i + y_offset
                + rng.normal(0, 0.03, size=len(x))
            )
            ax.scatter(x, y_jitter,
                       color="black", s=12, alpha=0.7,
                       edgecolor="none", zorder=3)

    # legend deliberately omitted -- hemisphere distinction is LH (dark) / RH (light)
    # within each violin body; a manual legend can be added externally if needed.

    ax.axvline(0, color="black", linestyle="-", alpha=0.2)
    ax.set_xlabel("Partial Correlation", fontsize=18, fontweight="bold")
    ax.set_ylabel("Seed Cluster",        fontsize=18, fontweight="bold")
    ax.set_title("{} target -- {}".format(target, variant),
                 fontsize=18, fontweight="bold")
    ax.tick_params(axis="both", which="major", labelsize=16)
    ax.grid(axis="x", alpha=0.5, linestyle=":", linewidth=0.5)
    ax.set_xlim(-0.25, 0.55)

    plt.tight_layout()

    outname = os.path.join(
        fig_dir,
        "violin_target-{}_partial_corr_{}.png".format(target, variant)
    )
    plt.savefig(outname, dpi=300, bbox_inches="tight")
    print("  Saved: {}".format(outname))
    plt.show()
    plt.close(fig)

# =========================
# Main plotting loop
# =========================

for variant in VARIANTS:
    print("\n=== Plotting variant: {} ===".format(variant))
    for target in target_clusters:
        plot_target_figure(target, variant, results[variant])