#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Group-level statistics for intra-hemispheric full correlations.

Operates on per-subject Fisher-z matrices produced by the computation script
Averaging is always performed in Fisher-z space; Pearson r is recovered only
at the final reporting stage via tanh()

This separation is intentional:
  - Fisher-z has approximately constant variance ~1/(n-3), making it the
    correct space for averaging and any subsequent parametric tests
  - Raw Pearson r values have r-dependent variance and should never be
    averaged directly

---------------------------------------------------
Written by Marco Bedini (marco.bedini@univ-amu.fr)
---------------------------------------------------
"""

import os
import sys
import numpy as np
import pandas as pd

# ============================================================
# Paths
# ============================================================

main_data             = "/scratch/mszinte/data/RetinoMaps/derivatives/pp_data"
partial_output_folder = (
    "/scratch/mszinte/data/RetinoMaps/derivatives/pp_data"
    "/group/91k/rest/partial_corr/by_hemi"
)
os.makedirs(partial_output_folder, exist_ok=True)

# Personal imports
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../"))
sys.path.append(os.path.abspath(os.path.join(base_dir, "analysis_code/utils")))
from settings_utils import load_settings

# ============================================================
# Settings
# ============================================================

project_dir       = "RetinoMaps"
settings_path     = os.path.join(base_dir, project_dir, "settings.yml")
prf_settings_path = os.path.join(base_dir, project_dir, "prf-analysis.yml")
settings          = load_settings([settings_path, prf_settings_path])
analysis_info     = settings[0]
subjects          = analysis_info["subjects"]

# ============================================================
# ROIs — needed for DataFrame index/columns and .npz metadata
# ============================================================

clusters = analysis_info["rois-drawn"]
seed_to_parcels = analysis_info["rois-group-mmp"]

clusters.reverse()   # mPCS first

parcels = []
for cl in clusters:
    parcels.extend(seed_to_parcels[cl])

# ============================================================
# Configuration
# ============================================================

HEMIS = [
    {"label": "LH", "tag": "lh"},
    {"label": "RH", "tag": "rh"},
]

RUNS = ["run-01", "run-02", ""]   # "" = full concatenated session

# ============================================================
# GROUP STATS
# ============================================================

print("=== GROUP STATS ===\n")