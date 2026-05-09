#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

histograms_plot_wta_consistency.py
-----------------------------------------------------------------------------------------
Goal:
    Plot WTA consistency figures from per-parcel WTA CSV outputs.

    Layout: 10 subplots — 5 macro-regions (columns) × 2 hemispheres (rows).
    Each subplot shows the parcels belonging to that macro-region × hemisphere:
        - bar length  = number of subjects with that GROUP winner
        - bar color   = GROUP winner seed color, full opacity if winner is one
                        of the 5 target macro-regions (mPCS/sPCS/iPCS/sIPS/iIPS),
                        low opacity (20%) otherwise
        - y-axis      = parcel names sorted by subject count descending (top = most consistent)

    Self-seed parcels (GROUP winner == target macro-region) are excluded from display.

-----------------------------------------------------------------------------------------
Inputs (sys.argv):
    1: main project directory   (e.g. /scratch/mszinte/data)
    2: project name/directory   (e.g. RetinoMaps)
    3: server group             (e.g. 327)
    4: corr_type                (partial_corr | full_corr)
    5: variant                  (concat | concat_clean | run-01 | run-02)
    6: mode                     (default | legacy | no_outliers | none)
                                (pass "none" for partial_corr)

Output:
    PDF: {fig_dir}/{corr_type}_{variant}_{mode}_wta_consistency.pdf

To run:
    $ cd projects/pRF_analysis/RetinoMaps/rest/visualizations
    $ python histograms_plot_wta_consistency.py /scratch/mszinte/data RetinoMaps 327 partial_corr concat_clean none
-----------------------------------------------------------------------------------------
Written by Marco Bedini (marco.bedini@univ-amu.fr)
Feeding Uriel's code into Claude to get consistent style
-----------------------------------------------------------------------------------------
"""

import warnings
warnings.filterwarnings("ignore")

import os
import sys
import re
import numpy as np
import pandas as pd
from pathlib import Path

import plotly.graph_objects as go
from plotly.subplots import make_subplots

# ============================================================
# Personal imports
# ============================================================
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../"))
sys.path.append(os.path.abspath(os.path.join(base_dir, "analysis_code/utils")))
from plot_utils import plotly_template
from settings_utils import load_settings

# ============================================================
# Parse and validate arguments
# ============================================================
USAGE = (
    "Usage: python wta_consistency_fig.py "
    "<main_dir> <project_dir> <group> <corr_type> <variant> <mode>\n"
    "  corr_type : partial_corr | full_corr\n"
    "  variant   : concat | concat_clean | run-01 | run-02\n"
    "  mode      : default | legacy | no_outliers | none  (none for partial_corr)"
)

VALID_CORR    = {"partial_corr", "full_corr"}
VALID_VARIANT = {"concat", "concat_clean", "run-01", "run-02"}
VALID_MODE    = {"default", "legacy", "no_outliers", "none"}

if len(sys.argv) != 7:
    print(f"ERROR: expected 6 arguments, got {len(sys.argv) - 1}.\n{USAGE}")
    sys.exit(1)

main_dir    = sys.argv[1]
project_dir = sys.argv[2]
group       = sys.argv[3]
corr_type   = sys.argv[4]
variant     = sys.argv[5]
mode        = sys.argv[6]

for val, valid, name in [
    (corr_type, VALID_CORR,    "corr_type"),
    (variant,   VALID_VARIANT, "variant"),
    (mode,      VALID_MODE,    "mode"),
]:
    if val not in valid:
        print(f"ERROR: invalid {name} '{val}'.\n  Accepted: {', '.join(sorted(valid))}\n{USAGE}")
        sys.exit(1)

if corr_type == "partial_corr" and mode != "none":
    print(f"WARNING: mode='{mode}' is ignored for partial_corr — set to 'none'.")
    mode = "none"

# ============================================================
# Load settings
# ============================================================
settings_path     = os.path.join(base_dir, project_dir, "settings.yml")
prf_settings_path = os.path.join(base_dir, project_dir, "prf-analysis.yml")
settings          = load_settings([settings_path, prf_settings_path])
analysis_info     = settings[0]

# ============================================================
# Paths
# ============================================================
main_data     = Path(main_dir) / project_dir / "derivatives/pp_data"
wta_wb_folder = main_data / "group/91k/rest/wta/workbench"
wta_nl_folder = main_data / "group/91k/rest/wta/nilearn"
fig_dir       = main_data / "group/91k/rest/wta/figures"
fig_dir.mkdir(parents=True, exist_ok=True)

# ============================================================
# ROI / macro-region definitions
# ============================================================
seed_to_parcels = analysis_info["rois-group-mmp"]   # {macro_region: [mmp_parcel, ...]}

# Ordered macro-region list — mPCS first (same as WTA pipeline)
macro_regions = list(analysis_info["rois-drawn"])
macro_regions.reverse()

# 1-based seed label → macro-region name
seed_to_number = {s: i + 1 for i, s in enumerate(macro_regions)}
number_to_seed = {v: k for k, v in seed_to_number.items()}

# The 5 target (in-scope) macro-regions
TARGET_MACROS = ["mPCS", "sPCS", "iPCS", "sIPS", "iIPS"]

# ============================================================
# Seed colors (RGB strings, matching the rest of the project)
# Full-opacity used for in-scope winners; low-opacity for out-of-scope.
# ============================================================
SEED_COLORS_RGB = {
    "mPCS": (255, 111,   0),
    "sPCS": (255, 234,   0),
    "iPCS": (151, 255,   0),
    "sIPS": ( 44, 255, 150),
    "iIPS": (  0, 152, 255),
    "hMT+": (  0,  25, 255),
    "VO":   (  0,   0, 200),
    "LO":   (150,   0,  90),
    "V3AB": (235, 127, 134),
    "V3":   (248, 160, 126),
    "V2":   (250, 196, 132),
    "V1":   (243, 231, 155),
}

OUT_OF_SCOPE_OPACITY = 0.2   # bars for winners outside the 5 target macro-regions

def seed_bar_color(winner_macro: str) -> str:
    """Return an rgba() string for the bar: full opacity if in-scope, muted if not."""
    r, g, b = SEED_COLORS_RGB.get(winner_macro, (180, 180, 180))
    if winner_macro in TARGET_MACROS:
        return f"rgba({r},{g},{b},1.0)"
    else:
        return f"rgba({r},{g},{b},{OUT_OF_SCOPE_OPACITY})"

# Legend colors always full opacity (for the legend chips to be readable)
def seed_legend_color(winner_macro: str) -> str:
    r, g, b = SEED_COLORS_RGB.get(winner_macro, (180, 180, 180))
    return f"rgb({r},{g},{b})"

# ============================================================
# Figure template
# ============================================================
template_specs = dict(
    axes_color="rgba(0, 0, 0, 1)",
    axes_width=2,
    axes_font_size=13,
    bg_col="rgba(255, 255, 255, 1)",
    font="Arial",
    title_font_size=14,
    rois_plot_width=1.5,
)
fig_template = plotly_template(template_specs)

# ============================================================
# Input CSV path builder
# ============================================================
def wta_csv_path(hemi: str) -> Path:
    if corr_type == "partial_corr":
        return wta_nl_folder / f"winning_seeds_by_subject_partial_corr_{hemi}_{variant}.csv"
    else:
        return wta_wb_folder / f"winning_seeds_by_subject_full_corr_{hemi}_{variant}_{mode}.csv"

# ============================================================
# Load WTA table — returns group row, consistency row, and n_subjects
# ============================================================
def load_wta_table(hemi: str):
    """
    Load the WTA CSV for a hemisphere and return:
        group_row        — Series: parcel → 1-based winner label (float)
        consistency_row  — Series: parcel → CONSISTENCY_% (float, 0–100)
        n_subjects       — int: number of subject rows (excludes GROUP/CONSISTENCY rows)
    """
    fpath = wta_csv_path(hemi)
    if not fpath.exists():
        raise FileNotFoundError(f"WTA CSV not found: {fpath}")

    df = pd.read_csv(fpath, index_col=0)

    required = {"GROUP", "CONSISTENCY_%"}
    missing  = required - set(df.index)
    if missing:
        raise ValueError(
            f"Expected rows {required} in {fpath.name} — missing: {missing}.\n"
            "Was this file produced by the current WTA pipeline?"
        )

    group_row       = df.loc["GROUP"].astype(float)
    consistency_row = df.loc["CONSISTENCY_%"].astype(float)

    # n_subjects = rows that look like subject IDs (start with "sub-")
    n_subjects = sum(1 for idx in df.index if str(idx).startswith("sub-"))

    return group_row, consistency_row, n_subjects

# ============================================================
# Build subplot data for one macro-region × hemisphere
# ============================================================
def build_subplot_data(
    target_macro: str,
    group_row: pd.Series,
    consistency_row: pd.Series,
    n_subjects: int,
) -> pd.DataFrame:
    """
    For a given target macro-region, collect parcel-level data:
        - winner_macro  : the GROUP winning seed
        - consistency   : CONSISTENCY_% (0–100)
        - n_winners     : number of subjects = consistency% × n_subjects / 100
        - color         : bar color (full opacity if in-scope, muted otherwise)

    Excludes self-seed parcels (GROUP winner == target macro-region).
    Sorted descending by n_winners (highest at top in plotly horizontal bars).
    """
    rows = []
    for parcel in seed_to_parcels[target_macro]:
        if parcel not in group_row.index:
            continue

        winner_label = group_row[parcel]
        if pd.isna(winner_label):
            continue

        winner_macro = number_to_seed[int(winner_label)]

        # Self-seed exclusion
        if winner_macro == target_macro:
            continue

        pct = consistency_row.get(parcel, np.nan)
        if pd.isna(pct):
            n_win = np.nan
        else:
            # Round to nearest integer subject count
            n_win = round(pct * n_subjects / 100.0)

        rows.append({
            "parcel":       parcel,
            "consistency":  pct,
            "n_winners":    n_win,
            "winner_macro": winner_macro,
            "color":        seed_bar_color(winner_macro),
        })

    if not rows:
        return pd.DataFrame(columns=["parcel", "consistency", "n_winners", "winner_macro", "color"])

    df = (
        pd.DataFrame(rows)
        .sort_values("n_winners", ascending=True)   # ascending → highest at top in plotly
    )
    return df

# ============================================================
# Build figure
# ============================================================
print("=" * 70)
print(
    f"WTA consistency figure — {corr_type} / {variant}"
    + (f" / {mode}" if mode != "none" else "")
)
print("=" * 70)

# Load both hemispheres
wta_data = {}
for hemi in ("lh", "rh"):
    print(f"  Loading {hemi}...")
    group_row, consistency_row, n_subjects = load_wta_table(hemi)
    wta_data[hemi] = (group_row, consistency_row, n_subjects)
    print(f"    n_subjects = {n_subjects}")

# Use the LH subject count for the x-axis range (both should match)
n_subjects_display = wta_data["lh"][2]

# Layout: 2 rows (LH, RH) × 5 cols (macro-regions)
N_ROWS = 2
N_COLS = len(TARGET_MACROS)

BAR_HEIGHT_PX = 22
TITLE_PAD_PX  = 40
max_parcels   = max(len(seed_to_parcels[m]) for m in TARGET_MACROS)
subplot_h_px  = max_parcels * BAR_HEIGHT_PX + TITLE_PAD_PX
fig_height    = N_ROWS * subplot_h_px + 140

subplot_titles = (
    [f"{m}  —  LH" for m in TARGET_MACROS]
    + [f"{m}  —  RH" for m in TARGET_MACROS]
)

fig = make_subplots(
    rows=N_ROWS,
    cols=N_COLS,
    subplot_titles=subplot_titles,
    horizontal_spacing=0.06,
    vertical_spacing=0.10,
    print_grid=False,
)

HEMIS_ROW = [("lh", 1), ("rh", 2)]

# Track which seeds have been added to the legend (one entry per seed globally)
legend_seeds_added = set()

for hemi, row_idx in HEMIS_ROW:
    group_row, consistency_row, n_subjects = wta_data[hemi]

    for col_idx, target_macro in enumerate(TARGET_MACROS, start=1):
        df_plot = build_subplot_data(
            target_macro, group_row, consistency_row, n_subjects
        )

        if df_plot.empty:
            fig.add_trace(go.Bar(x=[], y=[], showlegend=False), row=row_idx, col=col_idx)
            continue

        # One trace per winning seed (groups legend entries and keeps color consistent)
        for winner_macro, df_seed in df_plot.groupby("winner_macro", sort=False):
            show_legend = winner_macro not in legend_seeds_added
            legend_seeds_added.add(winner_macro)

            fig.add_trace(
                go.Bar(
                    x=df_seed["n_winners"],
                    y=df_seed["parcel"],
                    orientation="h",
                    name=winner_macro,
                    legendgroup=winner_macro,
                    # Bar fill color (muted if out-of-scope)
                    marker=dict(
                        color=df_seed["color"].tolist(),
                        line=dict(width=0),
                    ),
                    # Legend chip always full-opacity
                    legendgrouptitle=None,
                    width=0.7,
                    showlegend=show_legend,
                    hovertemplate=(
                        "<b>%{y}</b><br>"
                        f"Winner: {winner_macro}<br>"
                        "Subjects: %{x}<extra></extra>"
                    ),
                ),
                row=row_idx,
                col=col_idx,
            )

        # X-axis: integer subject counts, range 0–n_subjects
        x_tickvals = list(range(0, n_subjects_display + 1, 5))
        fig.update_xaxes(
            range=[0, n_subjects_display],
            tickvals=x_tickvals,
            ticktext=[str(v) for v in x_tickvals],
            title=dict(
                text="N subjects" if row_idx == N_ROWS else "",
                standoff=6,
            ),
            showline=True,
            linecolor="black",
            row=row_idx, col=col_idx,
        )
        fig.update_yaxes(
            showline=True,
            linecolor="black",
            tickfont=dict(size=11),
            row=row_idx, col=col_idx,
        )

# ============================================================
# Fix legend marker colors to full opacity
# (marker.color in the trace controls bars; legenditem color follows it,
#  but since bars per trace may be muted we override via update_traces)
# ============================================================
for seed in legend_seeds_added:
    fig.update_traces(
        selector=dict(name=seed, legendgroup=seed),
        legendrank=list(macro_regions).index(seed) if seed in macro_regions else 99,
    )

# ============================================================
# Global layout
# ============================================================
mode_label = f" [{mode}]" if mode != "none" else ""
fig_title  = (
    f"WTA consistency — {corr_type.replace('_', ' ')} / "
    f"{variant}{mode_label}  "
    f"(n={n_subjects_display} subjects, bars = N subjects with GROUP winner)"
)

fig.update_layout(
    template=fig_template,
    title=dict(text=fig_title, font=dict(size=15), x=0.5, xanchor="center"),
    height=int(fig_height),
    width=1600,
    barmode="stack",
    legend=dict(
        title=dict(text="Winning seed<br><sup>faded = outside target regions</sup>"),
        orientation="v",
        x=1.01,
        y=1.0,
        xanchor="left",
        font=dict(size=12),
        itemsizing="constant",
    ),
    margin=dict(l=80, r=200, t=100, b=60),
)

# ============================================================
# Export
# ============================================================
mode_token = f"_{mode}" if mode != "none" else ""
fig_fn     = fig_dir / f"wta_consistency_{corr_type}_{variant}{mode_token}.pdf"

print(f"\n  Saving: {fig_fn}")
fig.write_image(str(fig_fn))
print("  Done.")