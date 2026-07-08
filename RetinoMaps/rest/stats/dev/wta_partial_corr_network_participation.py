#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
wta_partial_corr_network_participation.py
-----------------------------------------------------------------------------------------
Goal:
   Assign winner-take-all (WTA) partial correlation results
   To the resting-state network partition by Ji et al. (2019)

-----------------------------------------------------------------------------------------
Inputs (sys.argv):
    1: main project directory   (e.g. /scratch/mszinte/data)
    2: project name/directory   (e.g. RetinoMaps)
    3: server group             (e.g. 327)
    4: server project           (e.g. b327)

Output:


To run:
    $ cd projects/pRF_analysis/RetinoMaps/rest/stats
    $ python wta_partial_corr_by_subject_by_hemi.py /scratch/mszinte/data RetinoMaps 327 b327
-----------------------------------------------------------------------------------------
Written by Marco Bedini (marco.bedini@univ-amu.fr)
-----------------------------------------------------------------------------------------
"""

import warnings
warnings.filterwarnings("ignore")

import os
import sys
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, List, Optional

# ============================================================
# Personal imports
# ============================================================
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../"))

sys.path.append(os.path.abspath(os.path.join(base_dir, "analysis_code/utils")))
from settings_utils import load_settings

# Load rest-specific settings (runs, exclusions)
rest_settings_path = os.path.join(base_dir, project_dir, "rest-settings.yml")
rest_settings      = load_settings([rest_settings_path])[0]

RUNS          = rest_settings["runs"]["value"]
RUN02_EXCLUDED = frozenset(rest_settings["run02_excluded"]["value"])

# ============================================================
# rest_utils — shared pipeline constants and functions
# ============================================================
sys.path.append(os.path.abspath(os.path.join(base_dir, "RetinoMaps/rest/utils")))
from rest_utils import (
    RUN02_EXCLUDED,
    VARIANTS,
    compute_winners,
    append_group_and_consistency,
)

# ============================================================
# Parse and validate arguments
# ============================================================
USAGE = (
    "Usage: python wta_partial_corr_by_subject_by_hemi.py "
    "<main_dir> <project_dir> <group> <server>"
)

if len(sys.argv) != 5:
    print(f"ERROR: expected 4 arguments, got {len(sys.argv) - 1}.\n{USAGE}")
    sys.exit(1)

main_dir    = sys.argv[1]
project_dir = sys.argv[2]
group       = sys.argv[3]
server      = sys.argv[4]

print("=" * 80)
print("WTA — partial correlation (Nilearn parcellated fisher-z CSVs)")
print("=" * 80)
print(f"  main_dir    : {main_dir}")
print(f"  project_dir : {project_dir}")
print(f"  group       : {group}")
print(f"  server      : {server}")

# ============================================================
# Load settings
# ============================================================
settings_path     = os.path.join(base_dir, project_dir, "settings.yml")
prf_settings_path = os.path.join(base_dir, project_dir, "prf-analysis.yml")
settings          = load_settings([settings_path, prf_settings_path])
analysis_info     = settings[0]
subjects          = analysis_info["subjects"]

# ============================================================
# Paths
# ============================================================
main_data          = Path(main_dir) / project_dir / "derivatives/pp_data"
output_folder      = main_data / "group/91k/rest/wta/nilearn"
group_corr_folder  = main_data / "group/91k/rest/partial_corr/by_hemi"
output_folder.mkdir(parents=True, exist_ok=True)

# ============================================================
# ROI / cluster settings
# ============================================================
clusters: List[str]                   = list(analysis_info["rois-drawn"])
seed_to_parcels: Dict[str, List[str]] = analysis_info["rois-group-mmp"]
clusters.reverse()                    # mPCS first

seed_to_number: Dict[str, int] = {s: i + 1 for i, s in enumerate(clusters)}

n_clusters = len(clusters)

