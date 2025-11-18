#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
--------------------------------------------------
Written by Marco Bedini (marco.bedini@univ-amu.fr) 
--------------------------------------------------
"""

import nibabel as nib
import numpy as np
import os

# Data paths
base_path = "/scratch/mszinte/data/RetinoMaps/derivatives/pp_data/group/91k/rest/wta/surfaces"
atlas_path = "/home/${USER}/projects/pRF_analysis/RetinoMaps/rest/mmp1_clusters/"
output_path = "/scratch/mszinte/data/RetinoMaps/derivatives/pp_data/group/91k/rest/wta"
os.makedirs(output_path, exist_ok=True)

# Brain settings
hemi_list = ["lh", "rh"]
region_labels = ["mPCS", "sPCS", "iPCS", "sIPS", "iIPS", "hMT+", "VO", "LO", "V3AB", "V3", "V2", "V1"]

# Loop per hemisphere
for hemi in hemi_list:
    correlation_data = []
    template_img = None

    print(f"\nProcessing hemisphere: {hemi}")

    for region in region_labels:
        # Correlation file
        corr_path = os.path.join(base_path, f"group_ses-01_task-rest_space-fsLR_den-91k_desc-full_corr_{region}_hollow_seed_{hemi}.func.gii")

        # Load correlation data
        corr_img = nib.load(corr_path)
        corr_data = corr_img.darrays[0].data.copy()

        correlation_data.append(corr_data)

        if template_img is None:
            template_img = corr_img

    # Stack and compute WTA
    stacked_data = np.stack(correlation_data, axis=0)
    wta_indices = np.argmax(stacked_data, axis=0) + 1  # 1-based indexing
    max_values = np.max(stacked_data, axis=0)

    # Save WTA map
    wta_out_path = os.path.join(output_path, f"group_winner_take_all_{hemi}.shape.gii")
    wta_darray = nib.gifti.GiftiDataArray(wta_indices.astype(np.float32))
    nib.save(nib.GiftiImage(darrays=[wta_darray]), wta_out_path)

    # Save max value map
    max_val_out_path = os.path.join(output_path, f"group_winner_take_all_max_values_{hemi}.shape.gii")
    maxval_darray = nib.gifti.GiftiDataArray(max_values.astype(np.float32))
    nib.save(nib.GiftiImage(darrays=[maxval_darray]), max_val_out_path)

output_path
