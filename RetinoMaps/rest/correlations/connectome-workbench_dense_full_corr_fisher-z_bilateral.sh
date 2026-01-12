#!/bin/bash

#####################################################
# Written by Marco Bedini (marco.bedini@univ-amu.fr)
#####################################################

# Remember to source your bashrc before running

# Define some paths
TASK_RESULTS="/scratch/mszinte/data/RetinoMaps/derivatives/pp_data"
ATLAS="/home/${USER}/projects/pRF_analysis/RetinoMaps/rest/mmp1_clusters"

# Iterate through subjects
for i in 01 02 03 04 05 06 07 08 09 11 12 13 14 17 20 21 22 23 24 25; 
do

SEED_DIR="$TASK_RESULTS/sub-${i}/91k/rest/seed"
OUT_DIR="$TASK_RESULTS/sub-${i}/91k/rest/corr/full_corr"

## Make sure all files are accessible
chmod -Rf 771 "$SEED_DIR"
chgrp -Rf 771 "$SEED_DIR"

mkdir "$TASK_RESULTS/sub-${i}/91k/rest/corr/full_corr/fisher-z"

    for ROI in mPCS sPCS iPCS sIPS iIPS hMT+ VO LO V3AB V3 V2 V1; 
    do

        # Compute correlation using both hemispheres simultaneously
        wb_command -cifti-correlation "$TASK_RESULTS/sub-${i}/91k/rest/timeseries/sub-${i}_ses-01_task-rest_space-fsLR_den-91k_desc-denoised_bold.dtseries.nii" \
            "$OUT_DIR/fisher-z/sub-${i}_ses-01_task-rest_space-fsLR_den-91k_desc-fisher-z_${ROI}.dconn.nii" -roi-override \
            -left-roi "$SEED_DIR/sub-${i}_91k_intertask_Sac_Pur_vision-pursuit-saccade_lh_${ROI}.shape.gii" \
            -right-roi "$SEED_DIR/sub-${i}_91k_intertask_Sac_Pur_vision-pursuit-saccade_rh_${ROI}.shape.gii" \
            -fisher-z -mem-limit 5;

        # Average correlation values within the seed ROI (both hemispheres)
        wb_command -cifti-average-dense-roi "$OUT_DIR/fisher-z/sub-${i}_ses-01_task-rest_space-fsLR_den-91k_desc-fisher-z_${ROI}.dscalar.nii" \
            -left-roi "$SEED_DIR/sub-${i}_91k_intertask_Sac_Pur_vision-pursuit-saccade_lh_${ROI}.shape.gii" \
            -right-roi "$SEED_DIR/sub-${i}_91k_intertask_Sac_Pur_vision-pursuit-saccade_rh_${ROI}.shape.gii" \
            -cifti "$OUT_DIR/fisher-z/sub-${i}_ses-01_task-rest_space-fsLR_den-91k_desc-fisher-z_${ROI}.dconn.nii";

        # Mask for visualization
        wb_command -cifti-restrict-dense-map "$OUT_DIR/fisher-z/sub-${i}_ses-01_task-rest_space-fsLR_den-91k_desc-fisher-z_${ROI}.dscalar.nii" COLUMN \
            "$OUT_DIR/fisher-z/sub-${i}_ses-01_task-rest_space-fsLR_den-91k_desc-fisher-z_${ROI}_masked.dscalar.nii" \
            -left-roi "$ATLAS/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_left_hemi.shape.gii" \
            -right-roi "$ATLAS/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_right_hemi.shape.gii";

        # Optional: Export to text for Python/Pandas use
        # wb_command -cifti-convert -to-text "$OUTPUT_PATH/sub-${i}/sub-${i}_ses-01_task-rest_space-fsLR_den-91k_desc-fisher-z_${ROI}.dconn.nii" \
        # "$OUTPUT_PATH/sub-${i}/sub-${i}_ses-01_task-rest_space-fsLR_den-91k_desc-fisher-z_${ROI}.dconn.txt"

    done

done
