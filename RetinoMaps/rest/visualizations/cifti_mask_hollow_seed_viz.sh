#!/bin/bash

## Important note: we are going to use only Pearson full correlations for viz not the Fisher-z ones (we are looking at the latter only to check the differences)
## Purpose: Mask correlation maps so that seed ROIs are hollowed out bilaterally

#####################################################
# Written by Marco Bedini (marco.bedini@univ-amu.fr)
#####################################################

# Define paths
DATA_PATH_CORR="/scratch/mszinte/data/RetinoMaps/derivatives/pp_data/group/91k/rest/median_full_corr"
DATA_PATH_FISHER="/scratch/mszinte/data/RetinoMaps/derivatives/pp_data/group/91k/rest/median_fisher-z"
ATLAS_PATH="/home/${USER}/projects/pRF_analysis/RetinoMaps/rest/mmp1_clusters/leaveout"
OUTPUT_PATH_CORR="/scratch/mszinte/data/RetinoMaps/derivatives/pp_data/group/91k/rest/hollow_seed_viz_full_corr"
OUTPUT_PATH_FISHER="/scratch/mszinte/data/RetinoMaps/derivatives/pp_data/group/91k/rest/hollow_seed_viz_fisher-z"

# Create output directories if needed
mkdir -p "$OUTPUT_PATH_CORR" "$OUTPUT_PATH_FISHER"

# Define ROIs (bilateral)
REGIONS=(mPCS sPCS iPCS sIPS iIPS hMT+ VO LO V3AB V3 V2 V1)

# Process correlation and Fisher-z maps
for region in "${REGIONS[@]}"; do
    echo "Masking: $region"

    # ROI masks
    roi_file_lh="${ATLAS_PATH}/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_leaveout_lh_${region}.shape.gii"
    roi_file_rh="${ATLAS_PATH}/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_leaveout_rh_${region}.shape.gii"

    # Full correlation
    input_corr="${DATA_PATH_CORR}/group_ses-01_task-rest_space-fsLR_den-91k_desc-full_corr_${region}_median.dscalar.nii"
    output_corr="${OUTPUT_PATH_CORR}/group_ses-01_task-rest_space-fsLR_den-91k_desc-full_corr_${region}_hollow_seed.dscalar.nii"
    wb_command -cifti-restrict-dense-map "$input_corr" COLUMN "$output_corr" \
        -left-roi "$roi_file_lh" \
        -right-roi "$roi_file_rh"

    # Fisher-z correlation
    input_fisher="${DATA_PATH_FISHER}/group_ses-01_task-rest_space-fsLR_den-91k_desc-fisher-z_${region}_median.dscalar.nii"
    output_fisher="${OUTPUT_PATH_FISHER}/group_ses-01_task-rest_space-fsLR_den-91k_desc-fisher-z_${region}_hollow_seed.dscalar.nii"
    wb_command -cifti-restrict-dense-map "$input_fisher" COLUMN "$output_fisher" \
        -left-roi "$roi_file_lh" \
        -right-roi "$roi_file_rh"
done

echo "All masking operations complete!"
