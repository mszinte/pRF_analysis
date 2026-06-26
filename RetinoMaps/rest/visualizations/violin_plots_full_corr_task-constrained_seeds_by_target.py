#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Goal of the script:
Full correlation violin plots — TASK-CONSTRAINED version
Seeds on the y-axis, one figure per TARGET macro-region

Input files are parcellated_by_macro (24 rows: 12 LH + 12 RH macro-regions,
ordered mPCS → V1 per hemisphere), so no parcel-level averaging is needed —
we index directly by cluster position

Produces figures for four data variants:
  - concat        : both runs concatenated, all subjects
  - concat_clean  : bad run-02 subjects fall back to run-01
  - run-01        : run-01 only, all subjects
  - run-02        : run-02 only, all subjects (bad subjects show artifact)

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
output_folder = main_data / "group/91k/rest/full_corr"
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
)

project_dir       = "RetinoMaps"
settings_path     = os.path.join(base_dir, project_dir, "settings.yml")
prf_settings_path = os.path.join(base_dir, project_dir, "prf-analysis.yml")
settings          = load_settings([settings_path, prf_settings_path])
analysis_info     = settings[0]
subjects          = analysis_info["subjects"]

fig_dir = os.path.join(output_folder, "figures/violin_plots/task-constrained")
os.makedirs(fig_dir, exist_ok=True)

# ============================================================
# Cluster ordering
# ============================================================

# Full ordered list, reversed so mPCS is first — this IS the row order
# inside the _parcellated_by_macro TSVs (12 LH rows then 12 RH rows).
clusters_all = analysis_info["rois-drawn"]
clusters_all.reverse()   # mPCS, sPCS, iPCS, sIPS, iIPS, hMT+, ..., V1

# The TSV has one row per macro-region per hemisphere:
#   rows 0–11  → LH clusters in clusters_all order
#   rows 12–23 → RH clusters in clusters_all order
N_CLUSTERS = len(clusters_all)   # expected 12

assert N_CLUSTERS == 12, (
    f"Expected 12 clusters after reverse, got {N_CLUSTERS}. "
    "Check 'rois-drawn' in settings.yml."
)

HEMI_ROW_SLICE = {
    "lh": slice(0,          N_CLUSTERS),
    "rh": slice(N_CLUSTERS, N_CLUSTERS * 2),
}

# Seed clusters: first 5 (mPCS → iIPS)
seed_clusters   = clusters_all[:5]
target_clusters = seed_clusters.copy()

# Validate
cluster_to_parcels = analysis_info["rois-group-mmp"]
for cl in seed_clusters:
    if cl not in cluster_to_parcels:
        raise KeyError(
            f"Seed cluster '{cl}' has no entry in 'rois-group-mmp' — "
            "check settings.yml"
        )

# ============================================================
# TSV loader — task-constrained _parcellated_by_macro files
# ============================================================

def load_macro_tsv_hemi(filepath, hemi):
    # type: (str, str) -> Optional[np.ndarray]
    """Load a 24-row parcellated_by_macro TSV and return the 12 values for hemi.

    Row layout: rows 0-11 = LH clusters, rows 12-23 = RH clusters,
    both in clusters_all order (mPCS first).

    Returns None (with a warning) on any file or format problem.
    """
    if not os.path.isfile(filepath):
        print("  WARNING: missing -- {}".format(filepath))
        return None
    try:
        values = (
            pd.read_csv(filepath, header=None, sep="\t")
            .squeeze()
            .to_numpy(dtype=float)
        )
    except Exception as e:
        print("  WARNING: could not read {} -- {}".format(filepath, e))
        return None

    expected_rows = N_CLUSTERS * 2   # 24
    if values.ndim != 1 or len(values) != expected_rows:
        print(
            "  WARNING: unexpected shape {} in {} "
            "(expected ({},)) -- skipping".format(values.shape, filepath, expected_rows)
        )
        return None

    return values[HEMI_ROW_SLICE[hemi]]

# ============================================================
# BIDS path builder
# ============================================================

def get_fname(sub_path, sub, variant, normal_tag, excluded_tag, hemi, seed):
    # type: (str, str, str, Optional[str], Optional[str], str, str) -> str
    """Return the task-constrained TSV filepath for one subject/hemi/seed.

    For concat_clean, subjects in RUN02_EXCLUDED use excluded_tag (run-01)
    so all subjects contribute one value rather than being dropped.
    """
    if variant == "concat_clean" and sub in RUN02_EXCLUDED:
        effective_tag = excluded_tag
    else:
        effective_tag = normal_tag

    task_constrained_subdir = os.path.join(sub_path, "task-constrained")

    if effective_tag is None:
        fname = (
            "{sub}_task-rest_space-fsLR_den-91k"
            "_desc-fisher-z_{hemi}_{seed}"
            "_task-constrained_parcellated_by_macro_legacy-mode.tsv"
        ).format(sub=sub, hemi=hemi, seed=seed)
    else:
        fname = (
            "{sub}_task-rest_{tag}_space-fsLR_den-91k"
            "_desc-fisher-z_{hemi}_{seed}"
            "_task-constrained_parcellated_by_macro_legacy-mode.tsv"
        ).format(sub=sub, tag=effective_tag, hemi=hemi, seed=seed)

    return os.path.join(task_constrained_subdir, fname)

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

missing_files = []

for variant, (normal_tag, excluded_tag, skip_excluded) in VARIANTS.items():

    print("\n=== Loading variant: {} | {} subjects ===".format(
        variant, len(subjects)))

    for sub in subjects:

        if skip_excluded and sub in RUN02_EXCLUDED:
            print("  SKIP (excluded): {} [{}]".format(sub, variant))
            continue

        sub_path = os.path.join(
            main_data,
            "{}/91k/rest/corr/full_corr/by_hemi".format(sub)
        )

        for seed in seed_clusters:
            for hemi in hemis:

                fname = get_fname(
                    sub_path    = sub_path,
                    sub         = sub,
                    variant     = variant,
                    normal_tag  = normal_tag,
                    excluded_tag= excluded_tag,
                    hemi        = hemi,
                    seed        = seed,
                )

                values = load_macro_tsv_hemi(fname, hemi)

                if values is None:
                    missing_files.append((variant, fname))
                    continue

                print("  {} ({}, {}, seed={}): loaded {} cluster values".format(
                    sub, variant, hemi, seed, len(values)))

                # values is now a 12-element array: one fisher-z per cluster
                # for this seed x hemi.  Index by cluster position in clusters_all.
                seed_row_in_file = clusters_all.index(seed)

                for tc in target_clusters:
                    if tc == seed:
                        continue
                    tc_row_in_file = clusters_all.index(tc)
                    results[variant][seed][tc][hemi].append(
                        float(values[tc_row_in_file])
                    )

if missing_files:
    print("\n" + "=" * 60)
    print("WARNING: {} file(s) could not be loaded:".format(len(missing_files)))
    for v, f in missing_files:
        print("  [{}] {}".format(v, f))
    print("=" * 60 + "\n")

# ============================================================
# Plotting
# ============================================================

def plot_target_figure(target, variant, res):
    # type: (str, str, dict) -> None
    """Draw and save one violin figure for target x variant."""

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

    n_plottable_rows     = len(seed_clusters) - 1
    expected_collections = n_plottable_rows * 2

    assert len(new_collections) == expected_collections, (
        "[{}, target={}] Expected {} PolyCollections but got {}. "
        "seaborn layout may have changed -- check version.".format(
            variant, target, expected_collections, len(new_collections))
    )

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

    ax.axhspan(
        target_row_idx - 0.45,
        target_row_idx + 0.45,
        color=SEED_COLOR,
        zorder=0,
        label="_nolegend_",
    )

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

    ax.axvline(0, color="black", linestyle="-", alpha=0.2)
    ax.set_xlabel("Full Correlation (task-constrained)", fontsize=18, fontweight="bold")
    ax.set_ylabel("Seed Cluster",                        fontsize=18, fontweight="bold")
    ax.set_title("{} target -- {} (task-constrained)".format(target, variant),
                 fontsize=18, fontweight="bold")
    ax.tick_params(axis="both", which="major", labelsize=16)
    ax.grid(axis="x", alpha=0.5, linestyle=":", linewidth=0.5)
    ax.set_xlim(-0.25, 0.55)

    plt.tight_layout()

    outname = os.path.join(
        fig_dir,
        "violin_target-{}_full_corr_{}_task-constrained.png".format(target, variant)
    )
    plt.savefig(outname, dpi=300, bbox_inches="tight")
    print("  Saved: {}".format(outname))
    plt.show()
    plt.close(fig)

# ============================================================
# Main loop
# ============================================================

for variant in VARIANTS:
    print("\n=== Plotting variant: {} ===".format(variant))
    for target in target_clusters:
        plot_target_figure(target, variant, results[variant])