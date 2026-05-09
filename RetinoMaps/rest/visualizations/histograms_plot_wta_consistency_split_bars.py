
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
histograms_plot_wta_consistency_split_bars.py
-----------------------------------------------------------------------------------------
Goal:
    Plot WTA consistency figures from per-parcel WTA CSV outputs.

    Layout: 1 row × 5 columns (one column per target macro-region).
    Each subplot shows the parcels of that macro-region with split bars:
        - Upper half-bar : LH winner (colored by LH GROUP winner seed)
        - Lower half-bar : RH winner (colored by RH GROUP winner seed)
        - Bar length     : number of subjects with that GROUP winner
        - In-scope winners (mPCS/sPCS/iPCS/sIPS/iIPS): full seed color
        - Out-of-scope winners: semi-transparent silver (excluded from legend)

    Parcel order derived from LH data (same for RH):
        1. In-scope winners first (ordered by TARGET_MACROS), out-of-scope last
        2. Within each winner group, sorted by n_subjects descending
        → Highest consistency in-scope parcels appear at the TOP of each subplot

    Self-seed parcels (GROUP winner == target macro-region) are excluded.

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
    SVG + PNG: {fig_dir}/wta_consistency_{corr_type}_{variant}{_mode}.svg/.png

To run:
    $ cd projects/pRF_analysis/RetinoMaps/rest/visualizations
    $ python histograms_plot_wta_consistency_sorted.py \
        /scratch/mszinte/data RetinoMaps 327 partial_corr concat_clean none
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

# Disable MathJax watermark on static exports
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
    "Usage: python histograms_plot_wta_consistency_sorted.py "
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

# The 5 target (in-scope) macro-regions
TARGET_MACROS = ["mPCS", "sPCS", "iPCS", "sIPS", "iIPS"]

# ============================================================
# Colors
#   In-scope  : full seed color
#   Out-of-scope : semi-transparent silver (matches document version)
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

OUT_OF_SCOPE_COLOR = "rgba(192, 192, 192, 0.2)"

def bar_color(winner_macro: str) -> str:
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
# Build per-parcel data for one macro-region × hemisphere
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
    Returns unsorted DataFrame: parcel, n_winners, winner_macro, color, in_scope
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

        pct = consistency_row.get(parcel, np.nan)
        n_w = round(pct * n_subjects / 100.0) if not pd.isna(pct) else np.nan
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
# Derive canonical parcel order from LH data
#
# Sort priority:
#   1. in_scope first (True before False)
#   2. within in-scope: winner order follows TARGET_MACROS sequence
#   3. within each winner group: n_winners ascending
# Then REVERSE so plotly (which renders first list entry at BOTTOM) places
# the highest-priority parcel at the TOP of the chart.
# ============================================================
def derive_parcel_order(df: pd.DataFrame) -> list:
    if df.empty:
        return []

    n_targets = len(TARGET_MACROS)

    def winner_sort_key(row):
        if row["in_scope"]:
            return TARGET_MACROS.index(row["winner_macro"])
        return (
            n_targets + sorted(macro_regions).index(row["winner_macro"])
            if row["winner_macro"] in macro_regions
            else n_targets + 99
        )

    df = df.copy()
    df["_winner_key"] = df.apply(winner_sort_key, axis=1)

    # Sort: winner group asc, then n_winners ASCENDING within group.
    # After reversing: last entry (highest n_winners in first in-scope group)
    # becomes the TOP of the plotly chart.
    df_sorted = df.sort_values(
        ["_winner_key", "n_winners"],
        ascending=[True, True],   # ascending n_winners so reverse puts highest at top
    )

    return df_sorted["parcel"].tolist()[::-1]

# ============================================================
# Main
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

# Canonical parcel order from LH (reused for RH)
parcel_orders = {}
for target_macro in TARGET_MACROS:
    gr_lh, cr_lh, ns_lh = wta_data["lh"]
    df_lh = build_parcel_data(target_macro, gr_lh, cr_lh, ns_lh)
    parcel_orders[target_macro] = derive_parcel_order(df_lh)

# ============================================================
# Figure layout — 1 row × 5 columns
# Each subplot contains two half-height bar traces per parcel:
#   LH bar : offset=+0.2  (upper half of the parcel slot)
#   RH bar : offset=-0.2  (lower half of the parcel slot)
#   bar width = 0.4  (half of the default 1.0 slot height)
# ============================================================
N_ROWS = 1
N_COLS = len(TARGET_MACROS)

BAR_HEIGHT_PX = 38    # slightly taller to accommodate two half-bars per parcel
TITLE_PAD_PX  = 45
max_parcels   = max(len(seed_to_parcels[m]) for m in TARGET_MACROS)
fig_height    = max_parcels * BAR_HEIGHT_PX + TITLE_PAD_PX + 160

fig = make_subplots(
    rows=N_ROWS,
    cols=N_COLS,
    subplot_titles=TARGET_MACROS,
    horizontal_spacing=0.06,
    vertical_spacing=0.10,
    print_grid=False,
)

# Track in-scope legend entries (one per seed, shown once globally)
# Out-of-scope seeds are NOT added to the legend
legend_added = set()

# LH = upper half-bar, RH = lower half-bar
HEMI_CONFIG = [
    # (hemi,  offset,  bar_width,  label_suffix)
    ("lh",  0.20,   0.38,       " LH"),
    ("rh", -0.20,   0.38,       " RH"),
]

for col_idx, target_macro in enumerate(TARGET_MACROS, start=1):
    canonical_order = parcel_orders[target_macro]

    if not canonical_order:
        fig.add_trace(go.Bar(x=[], y=[], showlegend=False), row=1, col=col_idx)
        continue

    for hemi, offset, bar_width, label_suffix in HEMI_CONFIG:
        group_row, consistency_row, n_subjects = wta_data[hemi]
        df_hemi = build_parcel_data(target_macro, group_row, consistency_row, n_subjects)

        if df_hemi.empty:
            fig.add_trace(go.Bar(x=[], y=[], showlegend=False), row=1, col=col_idx)
            continue

        df_hemi = df_hemi.set_index("parcel")

        # Build value + color arrays in canonical order
        x_vals, colors, winner_macros = [], [], []
        for parcel in canonical_order:
            if parcel in df_hemi.index:
                n_w = df_hemi.loc[parcel, "n_winners"]
                x_vals.append(0 if pd.isna(n_w) else float(n_w))
                colors.append(df_hemi.loc[parcel, "color"])
                winner_macros.append(df_hemi.loc[parcel, "winner_macro"])
            else:
                x_vals.append(0.0)
                colors.append("rgba(0,0,0,0)")
                winner_macros.append("__missing__")

        # Split into one trace per in-scope winner (for legend) +
        # one trace for all out-of-scope (no legend entry)
        # We use a single multi-color trace per hemi per subplot for the bars,
        # plus invisible zero-length traces solely to drive legend entries.

        # --- Legend driver traces (in-scope only, shown once globally) ---
        for winner_macro in TARGET_MACROS:
            if winner_macro == target_macro:
                continue   # self-seed never appears
            legend_key = winner_macro   # same key for LH and RH — one entry
            if legend_key not in legend_added:
                r, g, b = SEED_COLORS_RGB[winner_macro]
                fig.add_trace(
                    go.Bar(
                        x=[None], y=[None],
                        orientation="h",
                        name=winner_macro,
                        legendgroup=winner_macro,
                        marker=dict(color=f"rgb({r},{g},{b})"),
                        showlegend=True,
                    ),
                    row=1, col=col_idx,
                )
                legend_added.add(legend_key)

        # --- Actual bar trace (all parcels, per-bar color, no legend) ---
        hover_texts = [
            (
                f"<b>{p}</b><br>"
                f"Hemi: {hemi.upper()}<br>"
                f"Winner: {wm}<br>"
                f"Subjects: {int(xv)}"
            )
            for p, wm, xv in zip(canonical_order, winner_macros, x_vals)
        ]

        fig.add_trace(
            go.Bar(
                x=x_vals,
                y=canonical_order,
                orientation="h",
                marker=dict(
                    color=colors,
                    line=dict(width=0),
                ),
                width=bar_width,
                offset=offset,
                showlegend=False,
                hovertemplate="%{customdata}<extra></extra>",
                customdata=hover_texts,
                # Small text label on bar for hemi identification
                text=[label_suffix.strip()] * len(canonical_order),
                textposition="inside",
                insidetextanchor="start",
                textfont=dict(size=8, color="rgba(0,0,0,0.4)"),
            ),
            row=1, col=col_idx,
        )

    # Axis formatting
    x_tickvals = list(range(0, n_subjects_display + 1, 5))
    fig.update_xaxes(
        range=[0, n_subjects_display],
        tickvals=x_tickvals,
        ticktext=[str(v) for v in x_tickvals],
        title=dict(text="N subjects", standoff=6),
        showline=True,
        linecolor="black",
        row=1, col=col_idx,
    )
    fig.update_yaxes(
        categoryorder="array",
        categoryarray=canonical_order,
        showline=True,
        linecolor="black",
        tickfont=dict(size=11),
        row=1, col=col_idx,
    )

# ============================================================
# Global layout
# ============================================================
mode_label = f" [{mode}]" if mode != "none" else ""
fig_title  = (
    f"WTA consistency — {corr_type.replace('_', ' ')} / {variant}{mode_label}  "
    f"(n={n_subjects_display} subjects | upper bar = LH, lower bar = RH | "
    f"grey = winner outside 5 target regions)"
)

fig.update_layout(
    template=fig_template,
    title=dict(text=fig_title, font=dict(size=13), x=0.5, xanchor="center"),
    height=int(fig_height),
    width=1400,
    barmode="overlay",
    legend=dict(
        title=dict(text="Winning seed<br><sup>(in-scope only)</sup>"),
        orientation="v",
        x=1.01,
        y=1.0,
        xanchor="left",
        font=dict(size=12),
        itemsizing="constant",
        tracegroupgap=2,
    ),
    margin=dict(l=80, r=200, t=100, b=60),
)

# ============================================================
# Export SVG + PNG
# ============================================================
mode_token = f"_{mode}" if mode != "none" else ""
stem = fig_dir / f"wta_consistency_{corr_type}_{variant}{mode_token}"

for fmt, scale in [("svg", 1), ("png", 2)]:
    out = Path(str(stem) + f".{fmt}")
    print(f"  Saving {fmt.upper()}: {out}")
    fig.write_image(str(out), format=fmt, scale=scale)

print("  Done.")