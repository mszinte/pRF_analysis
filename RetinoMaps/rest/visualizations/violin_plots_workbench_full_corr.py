#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Full correlation violin plots split by hemisphere
From workbench outputs

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

# ============================================================
# Paths
# ============================================================
USER = os.environ["USER"]

# Main folders
main_data = Path("/scratch/mszinte/data/RetinoMaps/derivatives/pp_data")
seed_folder = main_data
atlas_folder = main_data / "atlas"
output_folder = main_data / "group/91k/rest/full_corr"
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

# =========================
# Cluster ↔ parcel mapping
# =========================

# clusters: full ordered list, reversed so mPCS is first
clusters_all = analysis_info["rois-drawn"]
clusters_all.reverse()   # mPCS, sPCS, iPCS, sIPS, iIPS, hMT+, ...

# cluster_to_parcels
cluster_to_parcels = analysis_info["rois-group-mmp"]

# seed_clusters: first 5 entries after reverse (the precentral and intraparietal eye field ROIs)
seed_clusters   = clusters_all[:5]
target_clusters = seed_clusters.copy()

# Validate: every seed cluster must have a non-empty parcel list
for cl in seed_clusters:
    if cl not in cluster_to_parcels:
        raise KeyError(
            f"Seed cluster '{cl}' has no entry in 'rois-group-mmp' — "
            "check settings.yml"
        )
    if len(cluster_to_parcels[cl]) == 0:
        raise ValueError(f"Seed cluster '{cl}' maps to an empty parcel list")

# =============================
# Colors for plotting clusters
# =============================

macro_colors = {
    "mPCS": "#FF6F00",
    "sPCS": "#FFEA00",
    "iPCS": "#97FF00",
    "sIPS": "#2CFF96",
    "iIPS": "#0098FF",
}

# Derived hemisphere shades: LH = darker, RH = lighter
SHADE_FACTORS = {"LH": 0.75, "RH": 1.25}
SEED_COLOR    = "#F5F5F5"   # near-white silver for the self-seed row


def derive_hemi_color(base_hex: str, hemi: str) -> np.ndarray:
    """Return an RGB array for *hemi* ('LH' or 'RH') derived from *base_hex*.

    LH is darkened (factor < 1), RH is lightened (factor > 1).
    Both are hard-clipped to [0, 1].
    """
    factor = SHADE_FACTORS[hemi]
    rgb = np.array(mcolors.to_rgb(base_hex))
    return np.clip(rgb * factor, 0, 1)

# =========================
# Build parcel indices
# =========================
 
# clusters_all is already ordered mPCS-first (reversed above)
# its order matches the row order of the .npy matrix files
clusters = clusters_all   # alias — used as clusters.index(seed) below
 
# Flat ordered parcel list, built from the YAML mapping in clusters_all order
# This matches the column order of the .npy matrix files (53 parcels, SCEF first)
parcels = []
for cl in clusters_all:
    parcels.extend(cluster_to_parcels[cl])
 
cluster_to_parcel_idx = {
    cl: [parcels.index(p) for p in plist]
    for cl, plist in cluster_to_parcels.items()
}
 
print("\nParcels from YAML:", len(parcels), "| first:", parcels[0])
print("Clusters from YAML:", clusters)
print("\nCluster → parcel indices:")

for cl in seed_clusters:
    print(f"  {cl:5s} → {[parcels[i] for i in cluster_to_parcel_idx[cl]]}")
    if len(cluster_to_parcel_idx[cl]) == 0:
        raise RuntimeError(f"No parcels found for cluster '{cl}'")

# =========================
# Compute BOTH hemispheres
# =========================

hemis = ["lh", "rh"]

results = {
    seed: {
        tc: {hemi: [] for hemi in hemis}
        for tc in target_clusters if tc != seed
    }
    for seed in seed_clusters
}

for sub in subjects:

    sub_path = os.path.join(main_data, f"{sub}/91k/rest/corr/full_corr")

    for seed in seed_clusters:
        for hemi in hemis:

            fname = (
                f"{sub_path}/{sub}_task-rest_space-fsLR_den-91k"
                f"_desc-fisher-z_{hemi}_{seed}_parcellated.tsv"
            )

            values_full = pd.read_csv(fname, header=None, sep='\t').squeeze().to_numpy(dtype=float)

            if len(values_full) != 106:
                raise ValueError(f"{fname} has {len(values_full)} rows, expected 106")

            values = values_full[:53] if hemi == "lh" else values_full[53:]

            print(f"{sub} ({hemi}, seed={seed}): loaded {len(values)} parcel values")

            for tc in target_clusters:
                if tc == seed:
                    continue

                idx = cluster_to_parcel_idx[tc]
                results[seed][tc][hemi].append(float(np.nanmean(values[idx])))

# =========================
# Plotting (split violins)
# =========================

for seed in seed_clusters:

    seed_row_idx = seed_clusters.index(seed)   # integer position in y-axis

    # ------------------------------------------------------------------
    # Build DataFrame — seed row gets NaN so seaborn still draws the
    # category slot (preserving y-axis ordering) but produces an empty
    # violin body that we will fully cover with a grey band.
    # Using NaN rather than zeros avoids a spurious spike at x = 0.
    # ------------------------------------------------------------------
    plot_rows = []

    for tc in seed_clusters:
        if tc == seed:
            # One NaN sentinel per hemisphere keeps the category present
            # but renders no visible violin body.
            for hemi_label in ["LH", "RH"]:
                plot_rows.append({
                    "Target":      tc,
                    "Correlation": np.nan,
                    "Hemisphere":  hemi_label,
                })
            continue

        for hemi in hemis:
            for val in results[seed][tc][hemi]:
                plot_rows.append({
                    "Target":      tc,
                    "Correlation": val,
                    "Hemisphere":  hemi.upper(),
                })

    df = pd.DataFrame(plot_rows)

    # ------------------------------------------------------------------
    # Draw violins — all white initially; we colour them right after.
    # Snapshot ax.collections BEFORE the call so we know exactly which
    # artists are newly added by violinplot.
    # ------------------------------------------------------------------

    fig, ax = plt.subplots(figsize=(8, 6))

    n_collections_before = len(ax.collections)

    sns.violinplot(
        data=df,
        y="Target",
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

    new_collections = ax.collections[n_collections_before:]

    # ------------------------------------------------------------------
    # Colour the violin halves
    #
    # seaborn (split=True) adds exactly 2 PolyCollection objects per
    # category row that has plottable data, in category order:
    # first the LH half (hue[0]), then the RH half (hue[1]).
    #
    # IMPORTANT: seaborn silently skips rows where ALL values are NaN —
    # it adds 0 collections for that row, not 2.  The seed row uses NaN
    # sentinels so it contributes 0 collections.  The expected count is
    # therefore (n_rows - 1) * 2, and we must skip the seed row when
    # mapping row index → collection index.
    # ------------------------------------------------------------------

    n_plottable_rows     = len(seed_clusters) - 1
    expected_collections = n_plottable_rows * 2

    assert len(new_collections) == expected_collections, (
        f"Expected {expected_collections} violin PolyCollections "
        f"(skipping the seed row '{seed}') but got {len(new_collections)}. "
        "seaborn's internal layout may have changed — check version."
    )

    hemi_labels    = ["LH", "RH"]
    collection_idx = 0

    for tc in seed_clusters:
        if tc == seed:
            continue   # no collections were added for this row

        for hemi_j, hemi_label in enumerate(hemi_labels):
            artist   = new_collections[collection_idx]
            base_hex = macro_colors[tc]
            rgb      = derive_hemi_color(base_hex, hemi_label)
            artist.set_facecolor(rgb)
            artist.set_edgecolor(rgb)
            artist.set_alpha(0.85)
            collection_idx += 1

    # ------------------------------------------------------------------
    # Grey band over the seed row — drawn at zorder=0 so it sits behind
    # the axes grid and all data artists
    # ------------------------------------------------------------------

    ax.axhspan(
        seed_row_idx - 0.45,
        seed_row_idx + 0.45,
        color=SEED_COLOR,
        zorder=0,
        label="_nolegend_",
    )

    # ------------------------------------------------------------------
    # Median lines — narrower span to avoid overlap; distinct grey shades
    # so LH and RH medians are visually separable without a legend.
    # LH: dark charcoal  |  RH: mid grey
    # ------------------------------------------------------------------

    MEDIAN_COLORS   = {"lh": "#222222", "rh": "#888888"}
    MEDIAN_HALF_LEN = 0.08

    for i, tc in enumerate(seed_clusters):
        if tc == seed:
            continue

        for hemi, y_offset in zip(hemis, [-0.12, 0.12]):
            vals = np.array(results[seed][tc][hemi])
            med  = np.nanmedian(vals)

            ax.plot(
                [med, med],
                [i + y_offset - MEDIAN_HALF_LEN,
                 i + y_offset + MEDIAN_HALF_LEN],
                color=MEDIAN_COLORS[hemi],
                linewidth=1.5,
                zorder=4,
            )

    # ------------------------------------------------------------------
    # Subject dots
    # ------------------------------------------------------------------

    rng = np.random.default_rng(seed=42)   # fixed seed → reproducible jitter

    for i, tc in enumerate(seed_clusters):
        if tc == seed:
            continue

        for hemi, y_offset in zip(hemis, [-0.1, 0.1]):
            x = np.array(results[seed][tc][hemi])
            y_jitter = (np.ones(len(x)) * i + y_offset
                        + rng.normal(0, 0.03, size=len(x)))

            ax.scatter(
                x,
                y_jitter,
                color="black",
                s=12,
                alpha=0.7,
                edgecolor="none",
                zorder=3,
            )

    # legend deliberately omitted — hemisphere distinction is LH (dark) / RH (light)
    # within each violin body; a manual legend can be added externally if needed.

    # ------------------------------------------------------------------
    # Styling
    # ------------------------------------------------------------------

    ax.axvline(0, color="black", linestyle="-", alpha=0.2)
    ax.set_xlabel("Partial Correlation (r)", fontsize=18, fontweight="bold")
    ax.set_ylabel("Target Cluster",          fontsize=18, fontweight="bold")
    ax.set_title(f"{seed} seed",             fontsize=18, fontweight="bold")
    ax.tick_params(axis='both', which='major', labelsize=16)

    ax.grid(axis="x", alpha=0.5, linestyle=":", linewidth=0.5)
    ax.set_xlim(-0.25, 0.55)

    plt.tight_layout()

    outname = os.path.join(
        fig_dir,
        f"violin_seed-{seed}_full_corr_by_hemi_fisher-z.png"
    )
    plt.savefig(outname, dpi=300, bbox_inches="tight")
    print(f"Saved figure: {outname}")
    plt.show()
    plt.close(fig)
