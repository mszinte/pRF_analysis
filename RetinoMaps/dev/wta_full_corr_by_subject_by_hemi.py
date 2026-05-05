"""
Created on Feb 27, 2025

wta_full_corr_by_subject_by_hemi.py
-----------------------------------------------------------------------------------------
Goal of the script:
Compute winner take all on full correlation results from workbench outputs
Written to run the same pipeline as with Nilearn outputs
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: server group (e.g. 327)
sys.argv[4]: server project (eg b327)
sys.argv[5]: grab legacy outputs or default
-----------------------------------------------------------------------------------------
Output(s):
TSV to import into the generate workbench dlabel file scripts
-----------------------------------------------------------------------------------------
To run:
1. cd to function
2. run python command
-----------------------------------------------------------------------------------------
Examples:
cd 
-----------------------------------------------------------------------------------------
Written by Marco Bedini (marco.bedini@univ-amu.fr)
-----------------------------------------------------------------------------------------
"""
# Stop warnings
import warnings
warnings.filterwarnings("ignore")

# Debug
import ipdb
deb = ipdb.set_trace

import os
import sys
import numpy as np
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
output_folder = main_data / "group/91k/rest/wta/workbench"
output_folder.mkdir(parents=True, exist_ok=True)

# Personal imports
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../"))
sys.path.append(os.path.abspath(os.path.join(base_dir, "analysis_code/utils")))
from settings_utils import load_settings
from surface_utils import load_surface
from cifti_utils import from_91k_to_32k

# Load settings
project_dir = 'RetinoMaps'
settings_path = os.path.join(base_dir, project_dir, "settings.yml")
prf_settings_path = os.path.join(base_dir, project_dir, "prf-analysis.yml")
settings = load_settings([settings_path, prf_settings_path])
analysis_info = settings[0]
subjects = analysis_info["subjects"]

# ============================================================
# ROIs
# ============================================================

# Load seed clusters and parcel to cluster assignments
clusters = analysis_info['rois-drawn']
seed_to_parcels = analysis_info['rois-group-mmp']

# Have mPCS as the first cluster instead of V1
clusters.reverse()

parcels = []
for cl in clusters:
    parcels.extend(seed_to_parcels[cl])

seed_to_number = {s: i+1 for i,s in enumerate(clusters)}

# =============================================================================
# Safe WTA function (no NaN/negative values or filtering)
# =============================================================================

def compute_winners(df_corr, seed_to_parcels, seed_to_number):
    """
    df_corr : DataFrame (seeds x parcels)
    """
    df = df_corr.copy()

    # exclude self-seed parcels
    for seed, plist in seed_to_parcels.items():
        for p in plist:
            if seed in df.index and p in df.columns:
                df.loc[seed, p] = np.nan

    winners = []
    for parcel in df.columns:
        col = df[parcel]

        # if entire column is NaN → no winner
        if col.isna().all():
            winners.append(np.nan)
        else:
            winners.append(seed_to_number[col.idxmax()])

    return np.array(winners)

# =============================================================================
# Loop over hemispheres
# =============================================================================

hemispheres = ['lh', 'rh']

for hemi in hemispheres:
    print(f"\n{'='*80}")
    print(f"Processing hemisphere: {hemi.upper()}")
    print('='*80)
    
    all_winners = []
    subject_ids = []
    
    # =============================================================================
    # Subject-level results
    # =============================================================================
    
    print("Processing subjects...")
    
    for subject in subjects:
        subj_dir = main_data / subject / "91k/rest/corr/full_corr/workbench_full_corr/by_hemi"
        
        csv_file = subj_dir / f"sub-08_task-rest_space-fsLR_den-91k_desc-fisher-z_lh_iIPS_parcellated_legacy-mode.tsv"
        
        if not csv_file.exists():
            print(f"{subject}: WARNING - File not found: {csv_file}")
            continue
        
        df_full = pd.read_csv(csv_file, index_col=0)
        winners = compute_winners(df_full, seed_to_parcels, seed_to_number)
        
        all_winners.append(winners)
        subject_ids.append(subject)
        
        print(f"{subject}: OK")
    
    # =============================================================================
    # Group-level results
    # =============================================================================
    
    group_result_path = main_data / "group/91k/rest/full_corr/by_hemi"
    group_csv = group_result_path / f"group_median_cluster_by_mmp-parcel_partial_{hemi}_by_hemi.csv"
    
    if not group_csv.exists():
        print(f"WARNING: Group file not found: {group_csv}")
        continue
    
    df_group = pd.read_csv(group_csv, index_col=0)
    group_winners = compute_winners(df_group, seed_to_parcels, seed_to_number)
    all_winners.append(group_winners)
    subject_ids.append("GROUP")
    
    # =============================================================================
    # Save table
    # =============================================================================
    
    winners_df = pd.DataFrame(all_winners, columns=df_group.columns, index=subject_ids)
    out_csv = output_folder / f"winning_seeds_by_subject_full_corr_{hemi}.csv"
    winners_df.to_csv(out_csv)
    
    print(f"\nSaved: {out_csv}")
    
    # =============================================================================
    # Compute consistency by parcel across subjects
    # =============================================================================
    
    subject_df = winners_df.iloc[:-1]
    
    consistency = []
    for parcel in winners_df.columns:
        gw = winners_df.loc["GROUP", parcel]
        if np.isnan(gw):
            consistency.append(np.nan)
        else:
            consistency.append(100 * (subject_df[parcel] == gw).sum() / len(subject_df))
    
    consistency_df = pd.DataFrame({
        "Parcel": winners_df.columns,
        "Group_Winner": winners_df.loc["GROUP"].values,
        "Consistency_%": consistency
    }).sort_values("Consistency_%", ascending=False)
    
    consistency_file = output_folder / f"winner_consistency_by_parcel_full_corr_{hemi}.csv"
    consistency_df.to_csv(consistency_file, index=False)
    
    print(f"Saved: {consistency_file}")

print("\n" + "="*80)
print("ALL HEMISPHERES COMPLETE")
print("="*80)
print(f"\nOutput files:")
print(f"  LH winners: {output_folder}/winning_seeds_by_subject_full_corr_lh.csv")
print(f"  RH winners: {output_folder}/winning_seeds_by_subject_full_corr_rh.csv")
print(f"  LH consistency: {output_folder}/winner_consistency_by_parcel_full_corr_lh.csv")
print(f"  RH consistency: {output_folder}/winner_consistency_by_parcel_full_corr_rh.csv")
print("\nDone.")