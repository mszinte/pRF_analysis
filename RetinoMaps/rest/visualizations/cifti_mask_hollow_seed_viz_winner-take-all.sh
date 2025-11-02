#!/bin/bash

## Here we are setting all values within the seed ROI to -0.999 so that they won't pop up in the winner take all visualization
## But this way the data keeps the same shape so we can stack it together and compute more easily
### WARNING: this script is not functional yet

# Define base paths
BASE_DIR="/scratch/mszinte/data/RetinoMaps/derivatives/pp_data/group/91k/rest/median_full_corr"
SEED_DIR="/home/${USER}/projects/pRF_analysis/RetinoMaps/rest/mmp1_clusters/leaveout"

# Create two output directories
mkdir -p "${BASE_DIR}/wta";
mkdir -p "${BASE_DIR}/thresholded";

# Loop over your regions
for ROI in mPCS sPCS iPCS sIPS iIPS hMT+ VO LO V3AB V3 V2 V1; do

# Threshold all full corr dscalar to include only non-zero vertices
wb_command -cifti-reduce "${BASE_DIR}/group_ses-01_task-rest_space-fsLR_den-91k_desc-full_corr_${ROI}_median.dscalar.nii" \
COUNT_NONZERO "${BASE_DIR}/thresholded/group_ses-01_task-rest_space-fsLR_den-91k_desc-full_corr_${ROI}_median_non-zero.dscalar.nii";
 
# Mask using the appropriate dscalar file
wb_command -cifti-math 'corr * (seed == 0) + (-0.999 * (seed > 0))' \
    "${BASE_DIR}/wta/group_ses-01_task-rest_space-fsLR_den-91k_desc-full_corr_${ROI}_median_wta.dscalar.nii" \
    -fixnan -0.999 \
    -var corr "${BASE_DIR}/thresholded/group_ses-01_task-rest_space-fsLR_den-91k_desc-full_corr_${ROI}_median_non-zero.dscalar.nii" \
    -var seed "${SEED_DIR}/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_discarded_${ROI}.dscalar.nii" \
    -select 2 INDEXMAX;

done

