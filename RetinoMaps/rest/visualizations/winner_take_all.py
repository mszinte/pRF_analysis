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
import sys

# Data paths
base_path = "/scratch/mszinte/data/RetinoMaps/derivatives/pp_data/group/91k/rest/wta"
atlas_path = "/home/${USER}/projects/pRF_analysis/RetinoMaps/rest/mmp1_clusters/"
output_path = "/scratch/mszinte/data/RetinoMaps/derivatives/pp_data/group/91k/rest/wta/results"
os.makedirs(output_path, exist_ok=True)

# Custom utils
main_codes = "~/disks/meso_H/projects/pRF_analysis/analysis_code/utils"
from surface_utils import load_surface
from cifti_utils import from_91k_to_32k, make_cifti_image

# Brain settings
hemi_list = ["lh", "rh"]
region_labels = ["mPCS", "sPCS", "iPCS", "sIPS", "iIPS", "hMT+", "VO", "LO", "V3AB", "V3", "V2", "V1"]

# Loop per clusters
correlation_data = []
template_img = None

for region in region_labels:
	# Correlation file
        corr_path = os.path.join(base_path, f"group_ses-01_task-rest_space-fsLR_den-91k_desc-full_corr_{region}_median.shape.gii")
        # Exclusion mask
        mask_path =os.path.join(atlas_path, f"atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_leaveout_{region}.shape.gii")

        # Load correlation data
        corr_img = load_surface(corr_path)
        corr_data = corr_img.darrays[0].data.copy()

        # Load mask and exclude those vertices by setting them to -inf
        mask_img = load_surface(mask_path)
        mask = mask_img.darrays[0].data.astype(bool)
        corr_data[mask] = -np.inf  # exclude from WTA

        correlation_data.append(corr_data)

        if template_img is None:
            template_img = corr_img

# Stack and compute WTA
stacked_data = np.stack(correlation_data, axis=0)
wta_indices = np.argmax(stacked_data, axis=0) + 1  # 1-based indexing
max_values = np.max(stacked_data, axis=0)

# Save WTA map
wta_out_path = os.path.join(output_path, f"group_winner_take_all.shape.nii")
wta_darray = nib.gifti.GiftiDataArray(wta_indices.astype(np.float32))
nib.save(nib.GiftiImage(darrays=[wta_darray]), wta_out_path)

# Save max value map
max_val_out_path = os.path.join(output_path, f"group_winner_take_all_max_values.shape.nii")
maxval_darray = nib.gifti.GiftiDataArray(max_values.astype(np.float32))
nib.save(nib.GiftiImage(darrays=[maxval_darray]), max_val_out_path)
    
## We can also save using the util
## make_cifti_image(cor_final_new, template_file, output_file)

output_path
