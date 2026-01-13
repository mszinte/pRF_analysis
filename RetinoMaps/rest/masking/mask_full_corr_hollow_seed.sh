#!/bin/bash

## Purpose: Mask correlation maps so that seed clusters (eg sPCS) are hollowed out bilaterally

#####################################################
# Written by Marco Bedini (marco.bedini@univ-amu.fr)
#####################################################

# Define paths
BASE_PATH="/scratch/mszinte/data/RetinoMaps/derivatives/pp_data"
ATLAS_PATH="/home/${USER}/projects/pRF_analysis/RetinoMaps/rest/mmp1_clusters/leaveout"

# Define ROIs (bilateral)
CLUSTERS=(mPCS sPCS iPCS sIPS iIPS hMT+ VO LO V3AB V3 V2 V1)

# Iterate through subjects
for sub in 01 02 03 04 05 06 07 08 09 11 12 13 14 17 20 21 22 23 24 25; do
    echo "Processing sub-${sub}..."

    FULL_CORR="${BASE_PATH}/sub-${sub}/91k/rest/corr/full_corr"

    # Create output directories if needed
    mkdir -p "${FULL_CORR}/hollow_seed" "${FULL_CORR}/fisher-z/hollow_seed"

    OUT_PATH_FULL="${FULL_CORR}/hollow_seed"
    OUT_PATH_FISHERZ="${FULL_CORR}/fisher-z/hollow_seed"

    # Process correlation and Fisher-z maps
    for region in "${CLUSTERS[@]}"; do
        echo "Masking: $region"

        # ROI masks
        roi_file_lh="${ATLAS_PATH}/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_leaveout_lh_${region}.shape.gii"
        roi_file_rh="${ATLAS_PATH}/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_leaveout_rh_${region}.shape.gii"

        # Full correlation
        input_corr="${FULL_CORR}/sub-${sub}_ses-01_task-rest_space-fsLR_den-91k_desc-full_corr_${region}_masked.dscalar.nii"
        output_corr="${OUT_PATH_FULL}/sub-${sub}_ses-01_task-rest_space-fsLR_den-91k_desc-full_corr_${region}_masked_hollow_seed.dscalar.nii"
        
        wb_command -cifti-restrict-dense-map "$input_corr" COLUMN "$output_corr" \
            -left-roi "$roi_file_lh" \
            -right-roi "$roi_file_rh"

        # Fisher-z correlation
        input_fisher="${FULL_CORR}/fisher-z/sub-${sub}_ses-01_task-rest_space-fsLR_den-91k_desc-fisher-z_${region}_masked.dscalar.nii"
        output_fisher="${OUT_PATH_FISHERZ}/sub-${sub}_ses-01_task-rest_space-fsLR_den-91k_desc-fisher-z_${region}_masked_hollow_seed.dscalar.nii"
        
        wb_command -cifti-restrict-dense-map "$input_fisher" COLUMN "$output_fisher" \
            -left-roi "$roi_file_lh" \
            -right-roi "$roi_file_rh"
    done
done

echo "All masking operations complete!"

