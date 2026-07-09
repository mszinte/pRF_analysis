#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Goal of the script:
Partial correlation violin plots — TASK-CONSTRAINED version
Seeds on the y-axis, one figure per TARGET macro-region

Input files: seed-task_by_macror-task_partial_fisherz_{variant}_{hemi}.csv
These contain fisher-z values which are converted to Pearson r for plotting

Currently only concat_clean is produced (the only reliable variant for
task-constrained partial correlations)

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
    MACRO_COLORS,
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

clusters_all = analysis_info["rois-drawn"]
clusters_all.reverse()   # mPCS, sPCS, iPCS, sIPS, iIPS, hMT+, ..., V1

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
        raise ValueError(
            "Seed cluster '{}' maps to an empty parcel list".format(cl)
        )

N_CLUSTERS = len(clusters_all)   # expected 12

# ============================================================
# Variants — task-constrained only uses concat_clean for now
# ============================================================

# variant_name : (csv_tag, subject_filter_fn)
# csv_tag is the token that appears in the filename (e.g. "concat")
# For concat_clean we load the concat file for clean subjects and the
# run-01 file for RUN02_EXCLUDED subjects (same logic as task-free).

def _all(sub):      return True
def _no_run02(sub): return sub not in RUN02_EXCLUDED

TC_VARIANTS = {
    "concat_clean": ("concat", _no_run02, "run-01"),
    # Add further variants here when the files exist, e.g.:
    # "run-01": ("run-01", _all,      None),
    # "run-02": ("run-02", _all,      None),
    # "concat": ("concat", _all,      None),
}

# ============================================================
# CSV loader with fisher-z → Pearson r conversion
# ============================================================

def load_partial_csv_hemi(filepath, hemi, n_clusters):
    # type: (str, str, int) -> Optional[np.ndarray]
    """Load a task-constrained partial corr CSV and return a (n_seeds x n_targets)
    matrix in Pearson r for the requested hemisphere.

    The CSV is expected to be a (n_seeds x n_clusters) matrix of fisher-z
    values for one hemisphere.  Rows = seeds (clusters_all order),
    Columns = target macro-regions (clusters_all order).

    Returns None (with a warning) on any file or format problem.
    """
    if not os.path.isfile(filepath):
        print("  WARNING: missing -- {}".format(filepath))
        return None
    try:
        arr = pd.read_csv(filepath, header=None).to_numpy(dtype=float)
    except Exception as e:
        print("  WARNING: could not read {} -- {}".format(filepath, e))
        return None

    if arr.ndim != 2:
        print("  WARNING: expected 2D array, got shape {} in {} -- skipping".format(
            arr.shape, filepath))
        return None

    if arr.shape != (n_clusters, n_clusters):
        print(
            "  WARNING: expected ({n}x{n}) matrix, got {s} in {f} -- skipping".format(
                n=n_clusters, s=arr.shape, f=filepath)
        )
        return None

    # Convert fisher-z to Pearson r
    return np.tanh(arr)


# ============================================================
# Path builder
# ============================================================

def get_fname(sub_path, sub, csv_tag, excluded_run_tag, hemi):
    # type: (str, str, str, Optional[str], str) -> str
    """Return the task-constrained partial corr CSV path for one subject/hemi.

    For concat_clean, subjects in RUN02_EXCLUDED use excluded_run_tag
    (run-01) instead of the concat file.
    """
    task_constrained_subdir = os.path.join(sub_path, "task-constrained")

    if sub in RUN02_EXCLUDED and excluded_run_tag is not None:
        effective_tag = excluded_run_tag
    else:
        effective_tag = csv_tag

    fname = "seed-task_by_macror-task_partial_fisherz_{}_{}.csv".format(
        effective_tag, hemi
    )
    return os.path.join(task_constrained_subdir, fname)

# ============================================================
# Load data
# ============================================================

hemis = ["lh", "rh"]

# results[variant][seed][tc][hemi] = list of per-subject Pearson r values
results = {
    variant: {
        seed: {
            tc: {hemi: [] for hemi in hemis}
            for tc in target_clusters if tc != seed
        }
        for seed in seed_clusters
    }
    for variant in TC_VARIANTS
}

missing_files = []

for variant, (csv_tag, subject_filter, excluded_run_tag) in TC_VARIANTS.items():

    variant_subjects = [s for s in subjects if subject_filter(s)]
    print("\n=== Loading variant: {} | {} subjects ===".format(
        variant, len(variant_subjects)))

    for sub in variant_subjects:

        sub_path = os.path.join(
            main_data,
            "{}/91k/rest/corr/partial_corr/by_hemi".format(sub)
        )

        for hemi in hemis:

            fname = get_fname(
                sub_path         = sub_path,
                sub              = sub,
                csv_tag          = csv_tag,
                excluded_run_tag = excluded_run_tag,
                hemi             = hemi,
            )

            matrix_r = load_partial_csv_hemi(fname, hemi, N_CLUSTERS)

            if matrix_r is None:
                missing_files.append((variant, fname))
                continue

            print("  {} ({}, {}): loaded {}x{} matrix (Pearson r)".format(
                sub, variant, hemi, matrix_r.shape[0], matrix_r.shape[1]))

            # matrix_r[seed_row, tc_col] = partial corr of seed with tc
            for seed in seed_clusters:
                seed_idx = clusters_all.index(seed)
                for tc in target_clusters:
                    if tc == seed:
                        continue
                    tc_idx = clusters_all.index(tc)
                    results[variant][seed][tc][hemi].append(
                        float(matrix_r[seed_idx, tc_idx])
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

MEDIAN_COLORS_LOCAL = {"lh": "#222222", "rh": "#888888"}
MEDIAN_HALF_LEN_LOCAL = 0.08


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
                [i + y_offset - MEDIAN_HALF_LEN_LOCAL,
                 i + y_offset + MEDIAN_HALF_LEN_LOCAL],
                color=MEDIAN_COLORS_LOCAL[hemi],
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
    ax.set_xlabel("Partial Correlation (r)", fontsize=18, fontweight="bold")
    ax.set_ylabel("Seed Cluster",            fontsize=18, fontweight="bold")
    ax.set_title(
        "{} target -- {} (task-constrained)".format(target, variant),
        fontsize=18, fontweight="bold"
    )
    ax.tick_params(axis="both", which="major", labelsize=16)
    ax.grid(axis="x", alpha=0.5, linestyle=":", linewidth=0.5)
    ax.set_xlim(-0.25, 0.55)

    plt.tight_layout()

    outname = os.path.join(
        fig_dir,
        "violin_target-{}_partial_corr_{}_task-constrained.png".format(
            target, variant)
    )
    plt.savefig(outname, dpi=300, bbox_inches="tight")
    print("  Saved: {}".format(outname))
    plt.show()
    plt.close(fig)

# ============================================================
# Main loop
# ============================================================

for variant in TC_VARIANTS:
    print("\n=== Plotting variant: {} ===".format(variant))
    for target in target_clusters:
        plot_target_figure(target, variant, results[variant])