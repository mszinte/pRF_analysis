#!/bin/bash

#####################################################
# Written by Marco Bedini (marco.bedini@univ-amu.fr)
#####################################################

# Remember to source your bashrc before running

# Define some paths
TASK_RESULTS="/scratch/mszinte/data/RetinoMaps/derivatives/pp_data"
ATLAS="/scratch/mszinte/data/RetinoMaps/derivatives/pp_data/atlas"

# Iterate through subjects
for i in 01 02 03 04 05 06 07 08 09 11 12 13 14 17 20 21 22 23 24 25;
do

SEED_DIR="$TASK_RESULTS/sub-${i}/91k/rest/seed"
OUT_DIR="$TASK_RESULTS/sub-${i}/91k/rest/corr/full_corr"

## Make sure all files are accessible
chmod -Rf 771 "$SEED_DIR"
chgrp -Rf 771 "$SEED_DIR"

mkdir -p "$OUT_DIR"

	for ROI in mPCS sPCS iPCS sIPS iIPS hMT+ VO LO V3AB V3 V2 V1;
	do
	  for HEMI in lh rh;
	  do

	    if [ "$HEMI" = "lh" ]; then
	      ROI_FLAG="-left-roi"
	    else
	      ROI_FLAG="-right-roi"
	    fi

	    ROI_FILE="$SEED_DIR/sub-${i}_91k_intertask_Sac_Pur_vision-pursuit-saccade_${HEMI}_${ROI}.shape.gii"

	    # Full correlation
	    wb_command -cifti-correlation \
	      "$TASK_RESULTS/sub-${i}/91k/rest/timeseries/sub-${i}_ses-01_task-rest_space-fsLR_den-91k_desc-denoised_bold.dtseries.nii" \
	      "$OUT_DIR/sub-${i}_task-rest_space-fsLR_den-91k_desc-full_corr_${HEMI}_${ROI}.dconn.nii" \
	      -roi-override $ROI_FLAG "$ROI_FILE" \
	      -mem-limit 20

	    # Fisher-z
	    wb_command -cifti-correlation \
	      "$TASK_RESULTS/sub-${i}/91k/rest/timeseries/sub-${i}_ses-01_task-rest_space-fsLR_den-91k_desc-denoised_bold.dtseries.nii" \
	      "$OUT_DIR/sub-${i}_task-rest_space-fsLR_den-91k_desc-fisher-z_${HEMI}_${ROI}.dconn.nii" \
	      -roi-override $ROI_FLAG "$ROI_FILE" \
	      -fisher-z -mem-limit 20

	    # Average within seed ROI
	    wb_command -cifti-average-dense-roi \
	      "$OUT_DIR/sub-${i}_task-rest_space-fsLR_den-91k_desc-full_corr_${HEMI}_${ROI}.dscalar.nii" \
	      $ROI_FLAG "$ROI_FILE" \
	      -cifti "$OUT_DIR/sub-${i}_task-rest_space-fsLR_den-91k_desc-full_corr_${HEMI}_${ROI}.dconn.nii"

	    wb_command -cifti-average-dense-roi \
	      "$OUT_DIR/sub-${i}_task-rest_space-fsLR_den-91k_desc-fisher-z_${HEMI}_${ROI}.dscalar.nii" \
	      $ROI_FLAG "$ROI_FILE" \
	      -cifti "$OUT_DIR/sub-${i}_task-rest_space-fsLR_den-91k_desc-fisher-z_${HEMI}_${ROI}.dconn.nii"

	    # Parcellate targets
	    wb_command -cifti-parcellate \
	      "$OUT_DIR/sub-${i}_task-rest_space-fsLR_den-91k_desc-full_corr_${HEMI}_${ROI}.dscalar.nii" \
	      "$ATLAS/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_dseg.dlabel.nii" COLUMN \
	      "$OUT_DIR/sub-${i}_task-rest_space-fsLR_den-91k_desc-full_corr_${HEMI}_${ROI}_parcellated.pscalar.nii" \
	      -method MEAN

	    wb_command -cifti-parcellate \
	      "$OUT_DIR/sub-${i}_task-rest_space-fsLR_den-91k_desc-fisher-z_${HEMI}_${ROI}.dscalar.nii" \
	      "$ATLAS/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_dseg.dlabel.nii" COLUMN \
	      "$OUT_DIR/sub-${i}_task-rest_space-fsLR_den-91k_desc-fisher-z_${HEMI}_${ROI}_parcellated.pscalar.nii" \
	      -method MEAN

	    # Parcellate excluding outliers
	    wb_command -cifti-parcellate \
	      "$OUT_DIR/sub-${i}_task-rest_space-fsLR_den-91k_desc-full_corr_${HEMI}_${ROI}.dscalar.nii" \
	      "$ATLAS/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_dseg.dlabel.nii" COLUMN \
	      "$OUT_DIR/sub-${i}_task-rest_space-fsLR_den-91k_desc-full_corr_${HEMI}_${ROI}_parcellated_no_outliers.pscalar.nii" \
	      -method MEAN -exclude-outliers 3 3

	    wb_command -cifti-parcellate \
	      "$OUT_DIR/sub-${i}_task-rest_space-fsLR_den-91k_desc-fisher-z_${HEMI}_${ROI}.dscalar.nii" \
	      "$ATLAS/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_dseg.dlabel.nii" COLUMN \
	      "$OUT_DIR/sub-${i}_task-rest_space-fsLR_den-91k_desc-fisher-z_${HEMI}_${ROI}_parcellated_no_outliers.pscalar.nii" \
	      -method MEAN -exclude-outliers 3 3
		
		# Mask vertex-wise results for visualizations in supplementary information
	    wb_command -cifti-restrict-dense-map "$OUT_DIR/sub-${i}_task-rest_space-fsLR_den-91k_desc-full_corr_${HEMI}_${ROI}.dscalar.nii" COLUMN \
	    	"$OUT_DIR/sub-${i}_task-rest_space-fsLR_den-91k_desc-full_corr_${HEMI}_${ROI}_masked.dscalar.nii" \
	    	-left-roi "$ATLAS/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_left_hemi.shape.gii" \
	    	-right-roi "$ATLAS/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_right_hemi.shape.gii";
		
	    # Same step as previous with Fisher-z outputs
	    wb_command -cifti-restrict-dense-map "$OUT_DIR/sub-${i}_task-rest_space-fsLR_den-91k_desc-fisher-z_${HEMI}_${ROI}.dscalar.nii" COLUMN \
	    	"$OUT_DIR/sub-${i}_task-rest_space-fsLR_den-91k_desc-fisher-z_${HEMI}_${ROI}_masked.dscalar.nii" \
	    	-left-roi "$ATLAS/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_left_hemi.shape.gii" \
	    	-right-roi "$ATLAS/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_right_hemi.shape.gii";
		
	  done
	done

# Grant permissions to output files
chmod -Rf 771 "$OUT_DIR"
chgrp -Rf 771 "$OUT_DIR"

done
