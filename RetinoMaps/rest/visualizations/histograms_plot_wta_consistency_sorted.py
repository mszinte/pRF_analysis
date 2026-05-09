#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

histograms_plot_wta_consistency_sorted.py

-----------------------------------------------------------------------------------------
Goal:
    Plot WTA consistency figures from per-parcel WTA CSV outputs.

    Layout: 10 subplots — 5 macro-regions (columns) × 2 hemispheres (rows).
    Each subplot shows the parcels belonging to that macro-region × hemisphere:
        - bar length  = number of subjects with that GROUP winner
        - bar color   = GROUP winner seed color if winner is one of the 5 target
                        macro-regions (mPCS/sPCS/iPCS/sIPS/iIPS),
                        silver/grey otherwise
        - y-axis      = parcel names in a fixed order derived from LH data:
                          1. grouped by winner (in-scope winners first, then grey)
                          2. within each winner group, sorted by n_subjects descending
                        The same parcel order is reused for RH so hemispheres are comparable.

    Self-seed parcels (GROUP winner == target macro-region) are excluded.

    Primary output of interest is LH (top row); RH (bottom row) uses the same
    parcel order for direct comparison.

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
    SVG + PNG: {fig_dir}/wta_consistency_{corr_type}_{variant}{_mode}.svg / .png

To run:
    $ cd projects/pRF_analysis/RetinoMaps/rest/visualizations
    $ python wta_consistency_fig.py /scratch/mszinte/data RetinoMaps 327 partial_corr concat_clean none
-----------------------------------------------------------------------------------------
Written by Marco Bedini (marco.bedini@univ-amu.fr)
Feeding Uriel's code into Claude to get consistent style
-----------------------------------------------------------------------------------------
"""

import warnings
warnings.filterwarnings("ignore")

import os
import sys
import numpy as np
import pandas as pd
from pathlib import Path

import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.io as pio

# Disable MathJax to prevent the "Loading MathJax" watermark on static exports
pio.kaleido.scope.mathjax = None

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
    print(f"WARNING: mode='{mode}' ignored for partial_corr — set to 'none'.")
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

# The 5 target (in-scope) macro-regions — drive column order and color logic
TARGET_MACROS = ["mPCS", "sPCS", "iPCS", "sIPS", "iIPS"]

# ============================================================
# Colors
# In-scope winners: full seed color
# Out-of-scope winners: silver/grey (color identity removed intentionally)
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

OUT_OF_SCOPE_COLOR = "rgba(192, 192, 192, 0.2)"   # silver

def bar_color(winner_macro: str) -> str:
    """Full seed color for in-scope winners, silver for out-of-scope."""
    if winner_macro in TARGET_MACROS:
        r, g, b = SEED_COLORS_RGB[winner_macro]
        return f"rgb({r},{g},{b})"
    return OUT_OF_SCOPE_COLOR

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
# Load WTA table
# ============================================================
def load_wta_table(hemi: str):
    """
    Returns:
        group_row        — Series: parcel → 1-based winner label (float)
        consistency_row  — Series: parcel → CONSISTENCY_% (float, 0–100)
        n_subjects       — int: number of subject rows
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
    n_subjects      = sum(1 for idx in df.index if str(idx).startswith("sub-"))

    return group_row, consistency_row, n_subjects

# ============================================================
# Build per-parcel data for one macro-region × hemisphere (unsorted)
# ============================================================
def build_parcel_data(
    target_macro: str,
    group_row: pd.Series,
    consistency_row: pd.Series,
    n_subjects: int,
) -> pd.DataFrame:
    """
    Collect parcel-level winner data for a target macro-region.
    Excludes self-seed parcels.
    Returns an unsorted DataFrame with columns:
        parcel, n_winners, winner_macro, color, in_scope
    """
    rows = []
    for parcel in seed_to_parcels[target_macro]:
        if parcel not in group_row.index:
            continue

        winner_label = group_row[parcel]
        if pd.isna(winner_label):
            continue

        winner_macro = number_to_seed[int(winner_label)]
        if winner_macro == target_macro:
            continue   # self-seed exclusion

        pct  = consistency_row.get(parcel, np.nan)
        n_w  = round(pct * n_subjects / 100.0) if not pd.isna(pct) else np.nan
        in_s = winner_macro in TARGET_MACROS

        rows.append({
            "parcel":       parcel,
            "n_winners":    n_w,
            "winner_macro": winner_macro,
            "color":        bar_color(winner_macro),
            "in_scope":     in_s,
        })

    return pd.DataFrame(
        rows,
        columns=["parcel", "n_winners", "winner_macro", "color", "in_scope"],
    )

# ============================================================
# Derive canonical parcel order from LH data for one macro-region
#
# Ordering rules (applied in priority):
#   1. In-scope winners (TARGET_MACROS) come before out-of-scope (grey) winners
#   2. Within each group, sort by winner macro (in TARGET_MACROS order for in-scope,
#      alphabetically for out-of-scope) so same-colored bars are adjacent
#   3. Within each winner, sort by n_winners descending (most consistent at top)
#
# The resulting parcel list is then reused for RH unchanged.
# ============================================================
def derive_parcel_order(df: pd.DataFrame) -> list:
    """
    Given a parcel DataFrame (from LH), return a list of parcel names in the
    canonical display order (bottom-to-top in plotly, i.e. ascending in the list
    means top-of-chart appears last).

    We want top-of-chart = highest consistency, so we return the list in
    ascending order of display rank (plotly horizontal bar: first entry = bottom).
    """
    if df.empty:
        return []

    # Assign winner group sort key:
    #   in-scope: position in TARGET_MACROS (0-based, lower = earlier = bottom of chart)
    #   out-of-scope: large number so they appear at the bottom of the chart
    n_targets = len(TARGET_MACROS)

    def winner_sort_key(row):
        if row["in_scope"]:
            return TARGET_MACROS.index(row["winner_macro"])
        else:
            return n_targets + sorted(macro_regions).index(row["winner_macro"]) \
                   if row["winner_macro"] in macro_regions else n_targets + 99

    df = df.copy()
    df["_winner_key"] = df.apply(winner_sort_key, axis=1)

    # Sort: primary = winner group (ascending puts out-of-scope at bottom of chart),
    #       secondary = n_winners descending within each group (most consistent at top)
    df_sorted = df.sort_values(
        ["_winner_key", "n_winners"],
        ascending=[True, False],
    )

    # plotly horizontal bar renders first list entry at the BOTTOM of the y-axis,
    # last entry at the TOP. We want in-scope groups at TOP with highest consistency
    # at the very top. So reverse the sorted order before passing to plotly.
    return df_sorted["parcel"].tolist()[::-1]

# ============================================================
# Main — load data, build figure
# ============================================================
print("=" * 70)
print(
    f"WTA consistency figure — {corr_type} / {variant}"
    + (f" / {mode}" if mode != "none" else "")
)
print("=" * 70)

wta_data = {}
for hemi in ("lh", "rh"):
    print(f"  Loading {hemi}...")
    gr, cr, ns = load_wta_table(hemi)
    wta_data[hemi] = (gr, cr, ns)
    print(f"    n_subjects = {ns}")

n_subjects_display = wta_data["lh"][2]

# Pre-compute canonical parcel orders from LH data
# {target_macro: [parcel, ...]}  — used for both LH and RH subplots
parcel_orders = {}
for target_macro in TARGET_MACROS:
    gr_lh, cr_lh, ns_lh = wta_data["lh"]
    df_lh = build_parcel_data(target_macro, gr_lh, cr_lh, ns_lh)
    parcel_orders[target_macro] = derive_parcel_order(df_lh)

# ============================================================
# Figure layout
# ============================================================
N_ROWS = 2
N_COLS = len(TARGET_MACROS)

BAR_HEIGHT_PX = 22
TITLE_PAD_PX  = 45
max_parcels   = max(len(seed_to_parcels[m]) for m in TARGET_MACROS)
subplot_h_px  = max_parcels * BAR_HEIGHT_PX + TITLE_PAD_PX
fig_height    = N_ROWS * subplot_h_px + 160

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

# Track legend entries globally (one chip per seed across all subplots)
legend_added = set()

for hemi, row_idx in HEMIS_ROW:
    group_row, consistency_row, n_subjects = wta_data[hemi]

    for col_idx, target_macro in enumerate(TARGET_MACROS, start=1):
        canonical_order = parcel_orders[target_macro]   # from LH

        if not canonical_order:
            fig.add_trace(go.Bar(x=[], y=[], showlegend=False), row=row_idx, col=col_idx)
            continue

        df_all = build_parcel_data(target_macro, group_row, consistency_row, n_subjects)

        if df_all.empty:
            fig.add_trace(go.Bar(x=[], y=[], showlegend=False), row=row_idx, col=col_idx)
            continue

        # Reindex to canonical order; parcels missing in this hemi get n_winners=0
        df_all = df_all.set_index("parcel")
        ordered_parcels  = canonical_order  # bottom → top in plotly
        x_vals, colors   = [], []

        for parcel in ordered_parcels:
            if parcel in df_all.index:
                n_w = df_all.loc[parcel, "n_winners"]
                col = df_all.loc[parcel, "color"]
                x_vals.append(0 if pd.isna(n_w) else n_w)
                colors.append(col)
            else:
                x_vals.append(0)
                colors.append("rgba(0,0,0,0)")   # invisible placeholder

        # Single trace per subplot (per-bar colors via marker.color list)
        # We still need separate traces per winner for the legend — split accordingly
        df_ordered = pd.DataFrame({
            "parcel":       ordered_parcels,
            "n_winners":    x_vals,
            "color":        colors,
        })

        # Re-attach winner_macro for legend grouping
        df_ordered["winner_macro"] = df_ordered["parcel"].map(
            lambda p: df_all.loc[p, "winner_macro"] if p in df_all.index else "__missing__"
        )

        for winner_macro, df_grp in df_ordered.groupby("winner_macro", sort=False):
            if winner_macro == "__missing__":
                continue

            show_legend = winner_macro not in legend_added
            legend_added.add(winner_macro)

            # Legend color: always full seed color (even if out-of-scope bars are grey,
            # the legend chip uses the true color so the reader can identify the region)
            if winner_macro in TARGET_MACROS:
                r, g, b = SEED_COLORS_RGB[winner_macro]
                legend_color = f"rgb({r},{g},{b})"
            else:
                legend_color = OUT_OF_SCOPE_COLOR

            fig.add_trace(
                go.Bar(
                    x=df_grp["n_winners"],
                    y=df_grp["parcel"],
                    orientation="h",
                    name=winner_macro,
                    legendgroup=winner_macro,
                    marker=dict(
                        color=df_grp["color"].tolist(),
                        line=dict(width=0),
                    ),
                    # Override legend chip color to always be the true seed color
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

        # Axis formatting
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
            categoryorder="array",
            categoryarray=ordered_parcels,   # enforce canonical order
            showline=True,
            linecolor="black",
            tickfont=dict(size=11),
            row=row_idx, col=col_idx,
        )

# ============================================================
# Global layout
# ============================================================
mode_label = f" [{mode}]" if mode != "none" else ""
fig_title  = (
    f"WTA consistency — {corr_type.replace('_', ' ')} / "
    f"{variant}{mode_label}  "
    f"(n={n_subjects_display} subjects, bar = N with GROUP winner; "
    f"grey = winner outside the 5 target regions)"
)

fig.update_layout(
    template=fig_template,
    title=dict(text=fig_title, font=dict(size=14), x=0.5, xanchor="center"),
    height=int(fig_height),
    width=1600,
    barmode="overlay",
    legend=dict(
        title=dict(text="Winning seed"),
        orientation="v",
        x=1.01,
        y=1.0,
        xanchor="left",
        font=dict(size=12),
        itemsizing="constant",
        tracegroupgap=4,
    ),
    margin=dict(l=80, r=200, t=100, b=60),
)

# ============================================================
# Export SVG + PNG (no PDF, avoids MathJax watermark issues)
# ============================================================
mode_token = f"_{mode}" if mode != "none" else ""
stem = fig_dir / f"wta_consistency_{corr_type}_{variant}{mode_token}"

for fmt, scale in [("svg", 1), ("png", 2)]:
    out = Path(str(stem) + f".{fmt}")
    print(f"  Saving {fmt.upper()}: {out}")
    fig.write_image(str(out), format=fmt, scale=scale)

print("  Done.")