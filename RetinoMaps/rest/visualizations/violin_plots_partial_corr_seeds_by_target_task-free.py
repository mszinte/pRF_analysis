#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Goal of the script:
Partial correlation violin plots split by hemisphere.
Seeds are on the y-axis, one figure per TARGET cluster.

Produces figures for four data variants:
  - concat        : both runs concatenated, all subjects
  - concat_clean  : both runs concatenated, excluding subjects with bad run-02
  - run-01        : run-01 only, all subjects
  - run-02        : run-02 only, ALL subjects (including bad ones, to show artifact)

Missing files are always skipped with a warning

---------------------------------------------------
Written by Marco Bedini (marco.bedini@univ-amu.fr)
---------------------------------------------------
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from pathlib import Path
from typing import Optional

# ============================================================
# Paths & settings
# ============================================================
USER = os.environ["USER"]

main_data     = Path("/scratch/mszinte/data/RetinoMaps/derivatives/pp_data")
output_folder = main_data / "group/91k/rest/partial_corr"
output_folder.mkdir(parents=True, exist_ok=True)

base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../"))
sys.path.append(os.path.abspath(os.path.join(base_dir, "analysis_code/utils")))
from settings_utils import load_settings

sys.path.append(os.path.abspath(os.path.join(base_dir, "RetinoMaps/rest/utils")))
from rest_utils import (
    RUN02_EXCLUDED,
    VARIANTS,
    MACRO_COLORS,
    SHADE_FACTORS,
    SEED_COLOR,
    MEDIAN_COLORS,
    MEDIAN_HALF_LEN,
    derive_hemi_color,
    load_npy_hemi,
)

project_dir       = "RetinoMaps"
settings_path     = os.path.join(base_dir, project_dir, "settings.yml")
prf_settings_path = os.path.join(base_dir, project_dir, "prf-analysis.yml")
settings          = load_settings([settings_path, prf_settings_path])
analysis_info     = settings[0]
subjects          = analysis_info["subjects"]

fig_dir = os.path.join(output_folder, "figures/violin_plots")
os.makedirs(fig_dir, exist_ok=True)

# ============================================================
# Cluster ↔ parcel mapping
# ============================================================

clusters_all = analysis_info["rois-drawn"]
clusters_all.reverse()   # mPCS, sPCS, iPCS, sIPS, iIPS, hMT+, ...

cluster_to_parcels = analysis_info["rois-group-mmp"]

seed_clusters   = clusters_all[:5]
target_clusters = seed_clusters.copy()

for cl in seed_clusters:
    if cl not in cluster_to_parcels:
        raise KeyError(
            f"Seed cluster '{cl}' has no entry in 'rois-group-mmp' — "
            "check settings.yml"
        )
    if len(cluster_to_parcels[cl]) == 0:
        raise ValueError(f"Seed cluster '{cl}' maps to an empty parcel list")

# ============================================================
# Build parcel indices
# ============================================================

clusters = clusters_all   # alias — row order matches .npy matrix files

parcels = []
for cl in clusters_all:
    parcels.extend(cluster_to_parcels[cl])

cluster_to_parcel_idx = {
    cl: [parcels.index(p) for p in plist]
    for cl, plist in cluster_to_parcels.items()
}

print(f"\nParcels from YAML: {len(parcels)} | first: {parcels[0]}")
print(f"Clusters from YAML: {clusters}")
print("\nCluster → parcel indices:")
for cl in seed_clusters:
    print(f"  {cl:5s} → {[parcels[i] for i in cluster_to_parcel_idx[cl]]}")
    if len(cluster_to_parcel_idx[cl]) == 0:
        raise RuntimeError(f"No parcels found for cluster '{cl}'")

# ============================================================
# BIDS path builder
# ============================================================

def get_fname(
    sub_path: str,
    sub: str,
    variant: str,
    normal_tag: Optional[str],
    excluded_tag: Optional[str],
    hemi: str,
) -> str:
    """Return the partial-correlation .npy filepath for one subject/hemi.

    For concat_clean, subjects with bad run-02 fall back to excluded_tag
    ("run-01") instead of being dropped.  For all other variants every
    subject uses normal_tag — including run-02 bad subjects, so their
    missing file produces a WARNING and they are skipped, preserving the
    artifact-visibility intent of that variant.
    """
    if variant == "concat_clean" and sub in RUN02_EXCLUDED:
        effective_tag = excluded_tag
    else:
        effective_tag = normal_tag

    if effective_tag is None:
        return os.path.join(
            sub_path,
            f"cluster_by_mmp-parcel_partial_{hemi}.npy"
        )
    return os.path.join(
        sub_path,
        f"cluster_by_mmp-parcel_partial_{effective_tag}_{hemi}.npy"
    )

# ============================================================
# Load data for all variants
# ============================================================

hemis = ["lh", "rh"]

# results[variant][seed][tc][hemi] = list of per-subject scalar values
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

## Need to add this to get pearsonr not fisher-z
#
#     results[variant][seed][tc][hemi].append(
#         float(np.tanh(np.nanmean(values[idx])))
#     )

missing_files = []   # global log — reported at end

for variant, (normal_tag, excluded_tag, skip_excluded) in VARIANTS.items():

    print(f"\n=== Loading variant: {variant} | {len(subjects)} subjects ===")

    for sub in subjects:

        if skip_excluded and sub in RUN02_EXCLUDED:
            print(f"  SKIP (excluded): {sub} [{variant}]")
            continue

        sub_path = os.path.join(
            main_data, f"{sub}/91k/rest/corr/partial_corr/by_hemi"
        )

        for hemi in hemis:

            fname = get_fname(
                sub_path=sub_path,
                sub=sub,
                variant=variant,
                normal_tag=normal_tag,
                excluded_tag=excluded_tag,
                hemi=hemi,
            )

            cluster_partial = load_npy_hemi(fname)

            if cluster_partial is None:
                missing_files.append((variant, fname))
                continue

            print(f"  {sub} ({variant}, {hemi}): {cluster_partial.shape}")

            for seed in seed_clusters:
                seed_idx = clusters.index(seed)

                for tc in target_clusters:
                    if tc == seed:
                        continue
                    parcel_idx = cluster_to_parcel_idx[tc]
                    vals = cluster_partial[seed_idx, parcel_idx]
                    results[variant][seed][tc][hemi].append(
                        float(np.nanmean(vals))
                    )

if missing_files:
    print(f"\n{'='*60}")
    print(f"WARNING: {len(missing_files)} file(s) could not be loaded:")
    for variant, f in missing_files:
        print(f"  [{variant}] {f}")
    print(f"{'='*60}\n")

# ============================================================
# Plotting
# ============================================================

def plot_target_figure(target: str, variant: str, res: dict) -> None:
    """Draw and save one violin figure for *target* × *variant*."""

    target_row_idx = seed_clusters.index(target)

    plot_rows = []
    for seed in seed_clusters:
        if seed == target:
            for hemi_label in ["LH", "RH"]:
                plot_rows.append(
                    {"Seed": seed, "Correlation": np.nan, "Hemisphere": hemi_label}
                )
            continue
        for hemi in hemis:
            for val in res[seed][target][hemi]:
                plot_rows.append(
                    {"Seed": seed, "Correlation": val, "Hemisphere": hemi.upper()}
                )

    df = pd.DataFrame(plot_rows)

    if df.loc[df["Seed"] != target, "Correlation"].dropna().empty:
        print(f"  SKIP (no data): target={target}, variant={variant}")
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

    # seaborn skips all-NaN rows (the target row), so expect one pair per
    # non-target seed
    n_plottable_rows     = len(seed_clusters) - 1
    expected_collections = n_plottable_rows * 2

    assert len(new_collections) == expected_collections, (
        f"[{variant}, target={target}] Expected {expected_collections} "
        f"PolyCollections but got {len(new_collections)}. "
        "seaborn layout may have changed — check version."
    )

    # Colour violin halves
    collection_idx = 0
    for seed in seed_clusters:
        if seed == target:
            continue
        for hemi_label in ["LH", "RH"]:
            artist = new_collections[collection_idx]
            rgb    = derive_hemi_color(MACRO_COLORS[seed], hemi_label)
            artist.set_facecolor(rgb)
            artist.set_edgecolor(rgb)
            artist.set_alpha(0.85)
            collection_idx += 1

    # Muted band for target (self-seed) row
    ax.axhspan(
        target_row_idx - 0.45,
        target_row_idx + 0.45,
        color=SEED_COLOR,
        zorder=0,
        label="_nolegend_",
    )

    # Median ticks
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

    # Per-subject dots
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

    # hemisphere distinction is LH (dark) / RH (light) within each violin
    # body; a manual legend can be added externally if needed
    ax.axvline(0, color="black", linestyle="-", alpha=0.2)
    ax.set_xlabel("Partial Correlation", fontsize=18, fontweight="bold")
    ax.set_ylabel("Seed Cluster",        fontsize=18, fontweight="bold")
    ax.set_title(f"{target} target — {variant}", fontsize=18, fontweight="bold")
    ax.tick_params(axis="both", which="major", labelsize=16)
    ax.grid(axis="x", alpha=0.5, linestyle=":", linewidth=0.5)
    ax.set_xlim(-0.25, 0.55)

    plt.tight_layout()

    outname = os.path.join(
        fig_dir,
        f"violin_target-{target}_partial_corr_{variant}.png"
    )
    plt.savefig(outname, dpi=300, bbox_inches="tight")
    print(f"  Saved: {outname}")
    plt.show()
    plt.close(fig)

# ============================================================
# Main loop
# ============================================================

for variant in VARIANTS:
    print(f"\n=== Plotting variant: {variant} ===")
    for target in target_clusters:
        plot_target_figure(target, variant, results[variant])