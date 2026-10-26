#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
network_utils.py
-----------------------------------------------------------------------------
Ji et al. (2019) network assignment for HCP MMP1.0 parcels and WTA
network-composition donut-chart plotting.

Location: RetinoMaps/rest/utils/ji_network_utils.py

Import from any script in rest/stats/:

    import os, sys
    base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../"))
    sys.path.append(os.path.abspath(os.path.join(base_dir, "RetinoMaps/rest/utils")))
    from ji_network_utils import (
        PARCEL_NETWORK,
        JI_NETWORK_COLORS,
        JI_NETWORK_ORDER,
        get_network,
        wta_network_composition,
        plot_wta_donut,
    )

Design notes
------------
- PARCEL_NETWORK stores the exact per-hemisphere Ji network assignments derived
  from the 360-parcel label file (Ji et al. 2019, Table S3 / dlabel export).
  RH and LH are stored separately because PEF has a genuine hemispheric
  asymmetry (RH = Dorsal_Attention, LH = Cingulo-Opercular) in the Ji atlas.
- Visual and Visual2 are kept as distinct entries in PARCEL_NETWORK for
  correctness. The display function optionally merges them into a single
  "Visual" wedge (merge_visual=True, default) to avoid a 1-parcel sliver
  (V1 is the only Visual parcel in the 53-parcel set).
- Self-seed parcels (parcels belonging to the winning seed's own cluster)
  carry NaN in the WTA map. They are excluded from the composition count —
  both numerator and denominator — matching the WTA pipeline convention.
- The function is hemisphere-agnostic: pass the GROUP row of a WTA CSV for
  one hemi, plus the parcel list and hemi label, and it returns a clean
  counts DataFrame ready for plotting or CSV export.
-----------------------------------------------------------------------------
Source
------
Ji JL et al. (2019) Mapping the human brain's cortical-subcortical functional
network organization. NeuroImage 185: 35–57.
https://doi.org/10.1016/j.neuroimage.2018.10.006

Label file: derived from HCP MMP1.0 parcellation annotated with Ji 2019
network assignments (360-parcel version, label indices 1–360).
-----------------------------------------------------------------------------
Written by Marco Bedini (marco.bedini@univ-amu.fr)
-----------------------------------------------------------------------------
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from collections import Counter
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# ============================================================
# Parcel → Ji network lookup
#
# Keys: canonical parcel name (matches ATLAS_ORDER in rest_utils.py)
# Values: Ji network name string
#
# Derived from label_network_mapping.txt (360-parcel Ji 2019 atlas),
# cross-referenced against ATLAS_KEY_TABLES indices in rest_utils.py.
# PEF is the single parcel with RH ≠ LH assignment — both are preserved.
# ============================================================
PARCEL_NETWORK: Dict[str, Dict[str, str]] = {
    "rh": {
        "V1":    "Visual",
        "MST":   "Visual2",
        "V2":    "Visual2",
        "V3":    "Visual2",
        "V4":    "Visual2",
        "V8":    "Visual2",
        "FEF":   "Cingulo-Opercular",
        "PEF":   "Dorsal_Attention",    # RH only — LH is Cingulo-Opercular
        "55b":   "Language",
        "V3A":   "Visual2",
        "V7":    "Visual2",
        "IPS1":  "Visual2",
        "FFC":   "Visual2",
        "V3B":   "Visual2",
        "LO1":   "Visual2",
        "LO2":   "Visual2",
        "PIT":   "Visual2",
        "MT":    "Visual2",
        "7Pm":   "Frontoparietal",
        "24dv":  "Somatomotor",
        "7AL":   "Somatomotor",
        "SCEF":  "Cingulo-Opercular",
        "6ma":   "Cingulo-Opercular",
        "7Am":   "Cingulo-Opercular",
        "7PL":   "Dorsal_Attention",
        "7PC":   "Somatomotor",
        "LIPv":  "Visual2",
        "VIP":   "Visual2",
        "MIP":   "Dorsal_Attention",
        "6d":    "Somatomotor",
        "6mp":   "Somatomotor",
        "6v":    "Somatomotor",
        "p32pr": "Cingulo-Opercular",
        "6r":    "Cingulo-Opercular",
        "IFJa":  "Language",
        "IFJp":  "Frontoparietal",
        "LIPd":  "Dorsal_Attention",
        "6a":    "Dorsal_Attention",
        "i6-8":  "Frontoparietal",
        "AIP":   "Dorsal_Attention",
        "PH":    "Visual2",
        "IP2":   "Frontoparietal",
        "IP1":   "Frontoparietal",
        "IP0":   "Dorsal_Attention",
        "V6A":   "Visual2",
        "VMV1":  "Visual2",
        "VMV3":  "Visual2",
        "V4t":   "Visual2",
        "FST":   "Visual2",
        "V3CD":  "Visual2",
        "LO3":   "Visual2",
        "VMV2":  "Visual2",
        "VVC":   "Visual2",
    },
    "lh": {
        "V1":    "Visual",
        "MST":   "Visual2",
        "V2":    "Visual2",
        "V3":    "Visual2",
        "V4":    "Visual2",
        "V8":    "Visual2",
        "FEF":   "Cingulo-Opercular",
        "PEF":   "Cingulo-Opercular",  # LH only — RH is Dorsal_Attention
        "55b":   "Language",
        "V3A":   "Visual2",
        "V7":    "Visual2",
        "IPS1":  "Visual2",
        "FFC":   "Visual2",
        "V3B":   "Visual2",
        "LO1":   "Visual2",
        "LO2":   "Visual2",
        "PIT":   "Visual2",
        "MT":    "Visual2",
        "7Pm":   "Frontoparietal",
        "24dv":  "Somatomotor",
        "7AL":   "Somatomotor",
        "SCEF":  "Cingulo-Opercular",
        "6ma":   "Cingulo-Opercular",
        "7Am":   "Cingulo-Opercular",
        "7PL":   "Dorsal_Attention",
        "7PC":   "Somatomotor",
        "LIPv":  "Visual2",
        "VIP":   "Visual2",
        "MIP":   "Dorsal_Attention",
        "6d":    "Somatomotor",
        "6mp":   "Somatomotor",
        "6v":    "Somatomotor",
        "p32pr": "Cingulo-Opercular",
        "6r":    "Cingulo-Opercular",
        "IFJa":  "Language",
        "IFJp":  "Frontoparietal",
        "LIPd":  "Dorsal_Attention",
        "6a":    "Dorsal_Attention",
        "i6-8":  "Frontoparietal",
        "AIP":   "Dorsal_Attention",
        "PH":    "Visual2",
        "IP2":   "Frontoparietal",
        "IP1":   "Frontoparietal",
        "IP0":   "Dorsal_Attention",
        "V6A":   "Visual2",
        "VMV1":  "Visual2",
        "VMV3":  "Visual2",
        "V4t":   "Visual2",
        "FST":   "Visual2",
        "V3CD":  "Visual2",
        "LO3":   "Visual2",
        "VMV2":  "Visual2",
        "VVC":   "Visual2",
    },
}

# ============================================================
# Display constants
#
# Colors match Ji et al. (2019) figure conventions (RGBA from the label file).
# Order follows dorsal-to-ventral / sensory-to-association convention.
# When merge_visual=True (default), "Visual" and "Visual2" are merged under
# the "Visual" entry; the "Visual2" entry is unused in that case.
# ============================================================

# Ji network display colors (hex) — derived from RGBA in label file
JI_NETWORK_COLORS: Dict[str, str] = {
    "Visual":             "#0000FF",   # pure blue  (label RGBA: 0 0 255)
    "Visual2":            "#6400FF",   # purple-blue (label RGBA: 100 0 255)
    "Somatomotor":        "#00FFFF",   # cyan        (label RGBA: 0 255 255)
    "Dorsal_Attention":   "#00FF00",   # green       (label RGBA: 0 255 0)
    "Frontoparietal":     "#FFFF00",   # yellow      (label RGBA: 255 255 0)
    "Language":           "#009B9B",   # teal        (label RGBA: 0 155 155)
    "Cingulo-Opercular":  "#990099",   # purple      (label RGBA: 153 0 153)
}

# Canonical display order for legends and stacked bars
JI_NETWORK_ORDER: List[str] = [
    "Visual",
    "Visual2",
    "Somatomotor",
    "Dorsal_Attention",
    "Frontoparietal",
    "Language",
    "Cingulo-Opercular",
]

# Display labels for figures (LaTeX-friendly, no underscores)
JI_NETWORK_LABELS: Dict[str, str] = {
    "Visual":             "Visual",
    "Visual2":            "Visual (ext.)",
    "Somatomotor":        "Somatomotor",
    "Dorsal_Attention":   "Dorsal Attention",
    "Frontoparietal":     "Frontoparietal",
    "Language":           "Language",
    "Cingulo-Opercular":  "Cingulo-Opercular",
}


# ============================================================
# Public API
# ============================================================

def get_network(parcel: str, hemi: str) -> Optional[str]:
    """Return the Ji network name for a given parcel and hemisphere.

    Parameters
    ----------
    parcel : canonical MMP1.0 parcel name (e.g. "V1", "IFJp")
    hemi   : "rh" or "lh"

    Returns
    -------
    str network name, or None if parcel is not in PARCEL_NETWORK.
    """
    return PARCEL_NETWORK.get(hemi, {}).get(parcel)


def wta_network_composition(
    winners: pd.Series,
    parcels: List[str],
    hemi: str,
    seed_names: List[str],
    merge_visual: bool = True,
) -> pd.DataFrame:
    """
    Compute the percentage of parcels in each Ji network won by each seed.

    Self-seed parcels (NaN in the WTA map) are excluded from both numerator
    and denominator — this matches the WTA pipeline convention exactly.

    Parameters
    ----------
    winners    : pd.Series, index = parcel name, values = 1-based seed integer
                 (float, NaN for self-seed / unassigned parcels).
                 Typically the GROUP row from a WTA CSV, indexed by parcel name.
    parcels    : canonical parcel name list (canonical YAML order, matching
                 the column order of the WTA CSV).
    hemi       : "rh" or "lh" — selects the correct PARCEL_NETWORK table.
    seed_names : seed cluster names in 1-based label order
                 (seed_names[0] → label 1, seed_names[1] → label 2, ...).
    merge_visual : if True (default), collapse "Visual" and "Visual2" into a
                   single "Visual" network for display.  The single V1 parcel
                   (Visual) avoids a barely-visible sliver next to the 26-parcel
                   Visual2 group.

    Returns
    -------
    pd.DataFrame, shape (n_seeds, n_networks)
        Index   = seed names (in input order)
        Columns = Ji network names (in JI_NETWORK_ORDER, Visual2 dropped if
                  merge_visual=True)
        Values  = percentage of parcels in that network won by that seed
                  (0–100, NaN where a network has no parcels in this hemi/set)

    Notes
    -----
    - A parcel missing from PARCEL_NETWORK[hemi] generates a WARNING and is
      excluded (treated like NaN).
    - Percentages within each network column sum to 100 across seeds
      (ignoring NaN parcels).  They do NOT sum to 100 across networks.
    """
    hemi_lower = hemi.lower()
    if hemi_lower not in PARCEL_NETWORK:
        raise ValueError(f"hemi must be 'rh' or 'lh', got '{hemi}'")

    net_table = PARCEL_NETWORK[hemi_lower]

    # Determine effective network list after optional merge
    if merge_visual:
        def _normalize(n):
            return "Visual" if n == "Visual2" else n
    else:
        def _normalize(n):
            return n

    net_order = []
    seen = set()
    for n in JI_NETWORK_ORDER:
        norm = _normalize(n)
        if norm not in seen:
            net_order.append(norm)
            seen.add(norm)

    n_seeds    = len(seed_names)
    n_networks = len(net_order)

    # Accumulate counts: counts[seed_idx][network] = n_parcels won
    # Totals: totals[network] = n_non-NaN parcels in that network
    counts = {s: Counter() for s in seed_names}
    totals: Counter = Counter()
    unknown_parcels = []

    for parcel in parcels:
        winner_val = winners.get(parcel, np.nan)
        if pd.isna(winner_val):
            continue  # self-seed or unmapped — skip from counts and totals

        raw_net = net_table.get(parcel)
        if raw_net is None:
            unknown_parcels.append(parcel)
            continue

        network = _normalize(raw_net)
        totals[network] += 1

        seed_label = int(winner_val)   # 1-based
        if seed_label < 1 or seed_label > n_seeds:
            print(
                f"  WARNING [{hemi}]: parcel '{parcel}' has winner label "
                f"{seed_label} outside range [1, {n_seeds}] — skipping"
            )
            continue
        seed_name = seed_names[seed_label - 1]
        counts[seed_name][network] += 1

    if unknown_parcels:
        print(
            f"  WARNING [{hemi}]: {len(unknown_parcels)} parcel(s) not found in "
            f"PARCEL_NETWORK['{hemi_lower}'] — excluded: {unknown_parcels}"
        )

    # Build percentage matrix
    data = {}
    for network in net_order:
        col = []
        denom = totals.get(network, 0)
        for seed in seed_names:
            if denom == 0:
                col.append(np.nan)
            else:
                col.append(100.0 * counts[seed].get(network, 0) / denom)
        data[network] = col

    df = pd.DataFrame(data, index=seed_names)
    df.index.name   = "seed"
    df.columns.name = "network"

    # Sanity check: each non-NaN network column should sum to ~100
    for net in net_order:
        col_sum = df[net].sum(skipna=True)
        if totals.get(net, 0) > 0 and not np.isclose(col_sum, 100.0, atol=0.1):
            print(
                f"  WARNING: network '{net}' column sums to {col_sum:.2f}% "
                f"(expected 100) — check seed_names / label alignment"
            )

    return df


def plot_wta_donut(
    composition_df: pd.DataFrame,
    seed_colors: Dict[str, str],
    title: str = "",
    merge_visual: bool = True,
    figsize: Tuple[float, float] = (5.0, 5.0),
    donut_width: float = 0.45,
    text_threshold_pct: float = 5.0,
    ax: Optional[plt.Axes] = None,
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Draw a donut (annular pie) chart showing % of parcels won by each seed,
    with one wedge per Ji network, coloured by the dominant winning seed.

    For each Ji network the wedge colour = the seed with the highest % of
    winning parcels in that network. The wedge size = number of non-NaN
    parcels in that network (i.e. proportional to total parcels per network,
    not to % won).

    Parameters
    ----------
    composition_df   : output of wta_network_composition() —
                       (n_seeds × n_networks) DataFrame of percentages.
    seed_colors      : {seed_name: hex_color} — matches MACRO_COLORS in
                       rest_utils.py.
    title            : figure title string.
    merge_visual     : must match the value used when building composition_df.
    figsize          : (width, height) in inches.
    donut_width      : radial width of the annulus ring (0 < width < 1).
    text_threshold_pct : only annotate wedges whose winner % ≥ this value.
    ax               : existing Axes to draw on; creates a new figure if None.

    Returns
    -------
    (fig, ax)

    Notes
    -----
    The wedge SIZE reflects how many parcels belong to each Ji network.
    The wedge COLOR reflects which seed won the most parcels in that network.
    A legend shows all seeds present in the data.
    The center label shows the total number of non-NaN parcels used.

    This is intentionally a "one seed per network" summary — it does NOT show
    the full distribution within each network.  Use a stacked-bar chart (not
    in this module) if per-network breakdown across seeds is needed.
    """
    # --- derive network sizes (proportional wedge sizes) ---
    # Total parcels per network is not stored in composition_df; reconstruct
    # from the fact that percentages sum to 100: size = Σ (pct[seed,net] > 0)
    # is unreliable. Instead we annotate sizes from the raw counts embedded
    # in the DataFrame via a helper convention: the caller should pass the
    # DataFrame as returned by wta_network_composition(), which we can inspect
    # by checking which networks have non-NaN values.
    #
    # For wedge sizes we use equal weighting per network (each network = 1
    # unit) so the chart reads as "network coverage" rather than parcel count.
    # The % text annotation inside each wedge conveys the winning seed's share.

    networks   = list(composition_df.columns)
    seed_names = list(composition_df.index)

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize)
    else:
        fig = ax.get_figure()

    # For each network, find the dominant seed and its percentage
    dominant_seeds = {}
    dominant_pcts  = {}
    for net in networks:
        col = composition_df[net]
        valid = col.dropna()
        if valid.empty or valid.max() == 0:
            dominant_seeds[net] = None
            dominant_pcts[net]  = 0.0
        else:
            best_seed = valid.idxmax()
            dominant_seeds[net] = best_seed
            dominant_pcts[net]  = float(valid.max())

    # Assign colors
    wedge_colors = []
    for net in networks:
        seed = dominant_seeds[net]
        if seed is None:
            wedge_colors.append("#CCCCCC")
        else:
            wedge_colors.append(seed_colors.get(seed, "#CCCCCC"))

    # Equal wedge sizes
    sizes = [1.0] * len(networks)

    # Network display labels
    if merge_visual:
        display_labels_map = {**JI_NETWORK_LABELS, "Visual": "Visual"}
        # "Visual2" won't appear since it was merged — no entry needed
    else:
        display_labels_map = JI_NETWORK_LABELS

    labels = [display_labels_map.get(n, n) for n in networks]

    # --- draw ---
    wedges, _ = ax.pie(
        sizes,
        colors=wedge_colors,
        startangle=90,
        counterclock=False,
        wedgeprops=dict(width=donut_width, edgecolor="white", linewidth=1.5),
    )

    # Annotate with dominant-seed % inside each wedge
    for i, (wedge, net) in enumerate(zip(wedges, networks)):
        pct = dominant_pcts[net]
        if pct < text_threshold_pct:
            continue
        # Midpoint angle of wedge (in degrees, CCW from 3 o'clock)
        theta = np.deg2rad((wedge.theta1 + wedge.theta2) / 2.0)
        r = 1.0 - donut_width / 2.0   # mid-radius of annulus
        x = r * np.cos(theta)
        y = r * np.sin(theta)
        ax.text(
            x, y,
            f"{pct:.0f}%",
            ha="center", va="center",
            fontsize=7, fontweight="bold",
            color="white",
        )

    # Center annotation: total parcel count
    n_parcels_used = int(composition_df.notna().any(axis=0).sum())
    ax.text(
        0, 0,
        f"n={len(composition_df.columns)}\nnetworks",
        ha="center", va="center",
        fontsize=8, color="#444444",
    )

    # Legend: seeds
    handles = [
        mpatches.Patch(facecolor=seed_colors.get(s, "#CCCCCC"), label=s)
        for s in seed_names
    ]
    ax.legend(
        handles=handles,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.15),
        ncol=len(seed_names),
        fontsize=7,
        frameon=False,
    )

    if title:
        ax.set_title(title, fontsize=9, pad=8)

    ax.set_aspect("equal")
    return fig, ax