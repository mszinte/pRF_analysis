#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
violin_partial_corr_task_constrained.py
---------------------------------------------------
Partial correlation violin plots — TASK-CONSTRAINED version
Seeds on the y-axis, one figure per TARGET macro-region

Reads from the concat_clean reporting TSVs produced by
group_stats_partial_corr_by_hemi_task-constrained.py:
    partial_corr/tables/seed-task_by_macror-task_partial-corr_r_report_ipsi_{hemi}_{estimator}.tsv

TSV structure:
    subject | seed | mPCS | sPCS | iPCS | sIPS | iIPS
    (one row per subject × seed; GROUP row at the bottom)

Only ipsilateral values are plotted (one table per hemisphere).
The GROUP row is plotted as a vertical line spanning each violin half.

Covariance estimator: must match the ESTIMATOR_TAG used when running
group_stats_partial_corr_by_hemi_task-constrained.py, which in turn must
match what was used at the subject level (nilearn_partial_corr_task_constrained.py).
Pass as sys.argv[1] (one of "raw", "ledoit-wolf", "graphical-lasso");
defaults to "ledoit-wolf" if omitted. Run once per estimator to get
directly comparable, separately-saved figure sets:

    $ python violin_partial_corr_task_constrained.py raw
    $ python violin_partial_corr_task_constrained.py ledoit-wolf
    $ python violin_partial_corr_task_constrained.py graphical-lasso

NOTE (filenames): this script's table/figure naming conventions still
differ slightly from the full-corr violin script's — to be harmonized
in a later pass.

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

# ============================================================
# Covariance estimator toggle
# ============================================================
VALID_ESTIMATORS = ("raw", "ledoit-wolf", "graphical-lasso")

if len(sys.argv) > 1:
    ESTIMATOR_TAG = sys.argv[1]
    if ESTIMATOR_TAG not in VALID_ESTIMATORS:
        print(
            "ERROR: unrecognised estimator '{}'.\n"
            "  Accepted: {}".format(ESTIMATOR_TAG, ", ".join(VALID_ESTIMATORS))
        )
        sys.exit(1)
else:
    ESTIMATOR_TAG = "ledoit-wolf"

print("Plotting estimator: {}".format(ESTIMATOR_TAG))

# ============================================================
# Paths & settings
# ============================================================
main_data     = Path("/scratch/mszinte/data/RetinoMaps/derivatives/pp_data")
tables_folder = main_data / "group/91k/rest/partial_corr/tables"
output_folder = main_data / "group/91k/rest/partial_corr"

base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../"))
sys.path.append(os.path.abspath(os.path.join(base_dir, "analysis_code/utils")))
from settings_utils import load_settings

sys.path.append(os.path.abspath(os.path.join(base_dir, "RetinoMaps/rest/utils")))
from rest_utils import (
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

# Figures tagged by estimator so raw / ledoit-wolf / graphical-lasso
# outputs never overwrite each other on disk.
fig_dir = os.path.join(
    output_folder, "figures/violin_plots/task-constrained", ESTIMATOR_TAG
)
os.makedirs(fig_dir, exist_ok=True)

# ============================================================
# Cluster ordering
# ============================================================
clusters_all = list(analysis_info["rois-drawn"])
clusters_all.reverse()   # mPCS first

seed_clusters   = clusters_all[:5]   # mPCS, sPCS, iPCS, sIPS, iIPS
target_clusters = seed_clusters.copy()

# ============================================================
# Load reporting TSVs — one per hemisphere
#
# Each TSV has columns: subject | seed | mPCS | sPCS | iPCS | sIPS | iIPS
# Rows: one per subject × seed, plus GROUP rows at the bottom.
# Values are already in Pearson r (conversion done at the stats stage).
# ============================================================
hemis = ["lh", "rh"]

tables = {}
missing_tables = []

for hemi in hemis:
    fname = (
        "seed-task_by_macror-task_partial-corr"
        "_r_report_ipsi_{hemi}_{estimator}.tsv".format(
            hemi=hemi, estimator=ESTIMATOR_TAG)
    )
    fpath = tables_folder / fname
    if not fpath.exists():
        print("WARNING: missing table -- {}".format(fpath))
        missing_tables.append(str(fpath))
        continue

    df = pd.read_csv(fpath, sep="\t")
    tables[hemi] = {
        "subjects": df[df["subject"] != "GROUP"].copy(),
        "group":    df[df["subject"] == "GROUP"].copy(),
    }
    n_subj = tables[hemi]["subjects"]["subject"].nunique()
    print("Loaded [{hemi}]: {n} subjects + GROUP".format(
        hemi=hemi.upper(), n=n_subj))

if missing_tables:
    print("\nERROR: {} reporting table(s) missing -- cannot plot.".format(
        len(missing_tables)))
    sys.exit(1)

# ============================================================
# Plotting
# ============================================================

def plot_target_figure(target):
    # type: (str) -> None
    """Draw and save one violin figure (ipsi, concat_clean) for one target."""

    target_row_idx = seed_clusters.index(target)

    plot_rows = []
    for hemi in hemis:
        subj_df = tables[hemi]["subjects"]
        for seed in seed_clusters:
            seed_rows = subj_df[subj_df["seed"] == seed]
            if seed == target:
                for _ in range(len(seed_rows)):
                    plot_rows.append({
                        "Seed": seed,
                        "Correlation": np.nan,
                        "Hemisphere": hemi.upper(),
                    })
                continue
            for val in seed_rows[target].values:
                plot_rows.append({
                    "Seed": seed,
                    "Correlation": float(val),
                    "Hemisphere": hemi.upper(),
                })

    df_plot = pd.DataFrame(plot_rows)

    if df_plot.loc[df_plot["Seed"] != target, "Correlation"].dropna().empty:
        print("  SKIP (no data): target={}".format(target))
        return

    fig, ax = plt.subplots(figsize=(8, 6))
    n_before = len(ax.collections)

    sns.violinplot(
        data=df_plot,
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
        "[target={}] Expected {} PolyCollections but got {}. "
        "seaborn layout may have changed -- check version.".format(
            target, expected_collections, len(new_collections))
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

    # Median tick lines (computed from subject rows, not GROUP)
    for i, seed in enumerate(seed_clusters):
        if seed == target:
            continue
        for hemi, y_offset in zip(hemis, [-0.12, 0.12]):
            vals = tables[hemi]["subjects"]
            vals = vals[vals["seed"] == seed][target].dropna().values
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

    # Individual subject dots with jitter
    rng = np.random.default_rng(seed=42)
    for i, seed in enumerate(seed_clusters):
        if seed == target:
            continue
        for hemi, y_offset in zip(hemis, [-0.1, 0.1]):
            vals = tables[hemi]["subjects"]
            x = vals[vals["seed"] == seed][target].dropna().values.astype(float)
            if len(x) == 0:
                continue
            y_jitter = (
                np.ones(len(x)) * i + y_offset
                + rng.normal(0, 0.03, size=len(x))
            )
            ax.scatter(x, y_jitter,
                       color="black", s=12, alpha=0.7,
                       edgecolor="none", zorder=3)

    # GROUP line markers — vertical line spanning the full violin half,
    # drawn over individual subject dots so the group value is clearly visible.
    for i, seed in enumerate(seed_clusters):
        if seed == target:
            continue
        for hemi, y_offset in zip(hemis, [-0.1, 0.1]):
            group_rows = tables[hemi]["group"]
            group_val  = group_rows[group_rows["seed"] == seed][target].values
            if len(group_val) == 0 or np.isnan(group_val[0]):
                continue
            ax.plot(
                [group_val[0]] * 2,
                [i + y_offset - 0.35, i + y_offset + 0.35],
                color=MEDIAN_COLORS[hemi],
                linewidth=2.5,
                zorder=6,
            )

    ax.axvline(0, color="black", linestyle="-", alpha=0.2)
    ax.set_xlabel("Partial Correlation (r, task-constrained)",
                  fontsize=18, fontweight="bold")
    ax.set_ylabel("Seed Cluster", fontsize=18, fontweight="bold")
    ax.set_title(
        "{} target -- concat_clean (task-constrained, {})".format(
            target, ESTIMATOR_TAG),
        fontsize=18, fontweight="bold"
    )
    ax.tick_params(axis="both", which="major", labelsize=16)
    ax.grid(axis="x", alpha=0.5, linestyle=":", linewidth=0.5)
    ax.set_xlim(-1.0, 1.0)

    plt.tight_layout()

    outname = os.path.join(
        fig_dir,
        "violin_target-{}_partial_corr_concat_clean_task-constrained_{}.png".format(
            target, ESTIMATOR_TAG)
    )
    plt.savefig(outname, dpi=300, bbox_inches="tight")
    print("  Saved: {}".format(outname))
    plt.show()
    plt.close(fig)

# ============================================================
# Main loop
# ============================================================

print("\n=== Plotting partial corr task-constrained (concat_clean, ipsi, {}) ===".format(
    ESTIMATOR_TAG))
for target in target_clusters:
    plot_target_figure(target)