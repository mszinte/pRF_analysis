#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
WTA full corr histogram plots 
Consistency of winner across subjects in PCEF-IPEF

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

fig_dir = os.path.join(output_folder, "figures/consistency_histograms")
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

# Load consistency outputs

