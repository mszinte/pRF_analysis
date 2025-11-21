#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

Winner-take-all routine for seed-based resting-state functional connectivity

--------------------------------------------------
Written by Marco Bedini (marco.bedini@univ-amu.fr) 
--------------------------------------------------

"""

import os
import sys
import numpy as np
import nibabel as nib

# ---------------------------------------
# Paths
# ---------------------------------------

USER = os.environ.get("USER", "unknown")

utils_path = os.path.expanduser(f"/home/{USER}/projects/pRF_analysis/analysis_code/utils")

sys.path.insert(0, utils_path)
from surface_utils import load_surface, make_cifti_image
from cifti_utils import from_91k_to_32k

base_path = "/scratch/mszinte/data/RetinoMaps/derivatives/pp_data/group/91k/rest/median_full_corr"
atlas_path = f"/home/{USER}/projects/pRF_analysis/RetinoMaps/rest/mmp1_clusters/leaveout"

output_path = ("/scratch/mszinte/data/RetinoMaps/derivatives/pp_data/group/91k/rest/wta/results")
os.makedirs(output_path, exist_ok=True)

region_labels = [
    "mPCS", "sPCS", "iPCS", "sIPS", "iIPS",
    "hMT+", "VO", "LO", "V3AB", "V3", "V2", "V1"
]

# -------------------------------------------------------
# Load data per region
# -------------------------------------------------------

correlation_data = []
template_32k = None

for region in region_labels:

    corr_path = os.path.join(base_path, f"group_ses-01_task-rest_space-fsLR_den-91k_desc-full_corr_{region}_median.dscalar.nii")

    mask_path = os.path.join(atlas_path, f"atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_leaveout_{region}.dlabel.nii")

    # ---- load 91k CIFTI ----
    cifti_img = nib.load(corr_path)
    corr_91k = np.asarray(cifti_img.get_fdata()).astype(np.float32)

    # ---- convert 91k → 32k ----
    conv = from_91k_to_32k(
        img=cifti_img,
        data=corr_91k,
        return_concat_hemis=True,
        return_32k_mask=False
    )
    corr_32k = conv["data_concat"].squeeze()   # shape ≈ 65k

    # ---- load 32k mask ----
    mask_img = nib.load(mask_path)
    mask_arr = np.asarray(mask_img.get_fdata()).squeeze().astype(int)
    mask = mask_arr > 0                         # shape = 65k

    if corr_32k.shape != mask.shape:
        raise RuntimeError(f"Size mismatch: corr={corr_32k.shape}, mask={mask.shape}")

    # ---- apply leave-one-out exclusion ----
    corr_32k[mask] = -np.inf

    correlation_data.append(corr_32k)

    if template_32k is None:
        template_32k = cifti_img   # keep original axis info for writing


# -------------------------------------------------------
# Winner-take-all computation
# -------------------------------------------------------

stacked = np.stack(correlation_data, axis=0)     # shape = (12 regions, 91k verts)
wta_indices = np.argmax(stacked, axis=0) + 1     # 1–12
max_values = np.max(stacked, axis=0)


# -------------------------------------------------------
# Save as CIFTI (proper fsLR dense file)
# -------------------------------------------------------

template_cifti = "/scratch/mszinte/data/RetinoMaps/derivatives//pp_data/group/91k/rest/median_full_corr/group_ses-01_task-rest_space-fsLR_den-91k_desc-full_corr_sPCS_median.dscalar.nii"

template_img = nib.load(template_cifti)

# WTA indices
wta_img = make_cifti_image(
    data=wta_indices[np.newaxis, :],   # shape = (1, vertices)
    source_img=template_img
)
wta_out = os.path.join(output_path, "group_winner_take_all.dscalar.nii")
nib.save(wta_img, wta_out)

# Max values
maxval_img = make_cifti_image(
    data=max_values[np.newaxis, :],
    source_img=template_img
)
maxval_out = os.path.join(output_path, "group_wta_max_values.dscalar.nii")
nib.save(maxval_img, maxval_out)


print("Saved:")
print(" -", wta_out)
print(" -", maxval_out)

