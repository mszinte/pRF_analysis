#!/bin/bash

# Debugging
# set -e 
# set -x

#####################################################
# Written by Marco Bedini (marco.bedini@univ-amu.fr)
#####################################################

# Get the intertask results path
TASK_RESULTS="/scratch/mszinte/data/RetinoMaps/derivatives/pp_data"
ATLAS="/home/${USER}/projects/pRF_analysis/RetinoMaps/rest/mmp1_clusters"

# Think about if we want to have a folder for templates
# TEMPLATE="/projects/pRF_analysis/RetinoMaps/?/fsLR"

	for i in 01 02 03 04 05 06 07 08 09 11 12 13 14 17 20 21 22 23 24 25 170k;
	do
	
	# Output dirs for every subject
	mkdir "/scratch/mszinte/data/RetinoMaps/derivatives/pp_data/sub-${i}/91k/rest/seed/source_170k"
	OUT_DIR1="$TASK_RESULTS/sub-${i}/91k/rest/seed/source_170k"
	mkdir "/scratch/mszinte/data/RetinoMaps/derivatives/pp_data/sub-${i}/91k/rest/seed/target_91k"
	OUT_DIR2="$TASK_RESULTS/sub-${i}/91k/rest/seed/target_91k"

		### 1. Separate the results by hemisphere

		wb_command -cifti-separate "$TASK_RESULTS/sub-${i}/170k/intertask/sub-${i}_intertask_Sac_Pur.dtseries.nii" COLUMN -metric CORTEX_LEFT \
		"$OUT_DIR1/sub-${i}_intertask_Sac_Pur.dtseries.left_hemi.shape.gii";

		wb_command -cifti-separate "$TASK_RESULTS/sub-${i}/170k/intertask/sub-${i}_intertask_Sac_Pur.dtseries.nii" COLUMN -metric CORTEX_RIGHT \
		"$OUT_DIR1/sub-${i}_intertask_Sac_Pur.dtseries.right_hemi.shape.gii";

		### 2. Resample from 170k to 91k - try this one for more advanced registration options wb_command -cifti-resample 
		### (in the step before maybe use -cifti-separate to get the dense labels dlabel.nii) 

		wb_command -metric-resample \
		    "$OUT_DIR1/sub-${i}_intertask_Sac_Pur.dtseries.left_hemi.shape.gii" \
		    fsLR/tpl-fsLR_hemi-L_den-59k_sphere.surf.gii \
		    fsLR/tpl-fsLR_hemi-L_den-32k_sphere.surf.gii \
		    BARYCENTRIC -largest \
		    "$OUT_DIR2/sub-${i}_91k_intertask_Sac_Pur.dtseries.left_hemi.bary_largest.shape.gii";
			
		wb_command -metric-resample \
		    "$OUT_DIR1/sub-${i}_intertask_Sac_Pur.dtseries.right_hemi.shape.gii" \
		    fsLR/tpl-fsLR_hemi-R_den-59k_sphere.surf.gii \
		    fsLR/tpl-fsLR_hemi-R_den-32k_sphere.surf.gii \
		    BARYCENTRIC -largest \
		    "$OUT_DIR2/sub-${i}_91k_intertask_Sac_Pur.dtseries.right_hemi.bary_largest.shape.gii";

		### 3. Import labels to rename them and change the color map (using Uriel's rgb convention)

		wb_command -metric-label-import "$OUT_DIR2/sub-${i}_91k_intertask_Sac_Pur.dtseries.left_hemi.bary_largest.shape.gii" \
		intertask_label_list.txt \
		"$OUT_DIR2/sub-${i}_91k_intertask_Sac_Pur_lh_renamed_cmap.label.gii" && \
		
		wb_command -metric-label-import "$OUT_DIR2/sub-${i}_91k_intertask_Sac_Pur.dtseries.right_hemi.bary_largest.shape.gii" \
		intertask_label_list.txt \
		"$OUT_DIR2/sub-${i}_91k_intertask_Sac_Pur_rh_renamed_cmap.label.gii";

		### 4. Mask the GLM results with the Glasser macro-regions (simply comment out these steps if you want to keep the results unmasked)

		wb_command -label-mask "$OUT_DIR2/sub-${i}_91k_intertask_Sac_Pur_lh_renamed_cmap.label.gii" \
		"$ATLAS/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_left_hemi.shape.gii" \
		"$OUT_DIR2/sub-${i}_91k_intertask_Sac_Pur_lh_renamed_masked_cmap.label.gii";

		wb_command -label-mask "$OUT_DIR2/sub-${i}_91k_intertask_Sac_Pur_rh_renamed_cmap.label.gii" \
		"$ATLAS/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_right_hemi.shape.gii" \
		"$OUT_DIR2/sub-${i}_91k_intertask_Sac_Pur_rh_renamed_masked_cmap.label.gii";

		### 5. Loop over each label to generate separate metric files (useful to mask the timeseries)

		# Pursuit
		wb_command -gifti-label-to-roi "$OUT_DIR2/sub-${i}_91k_intertask_Sac_Pur_lh_renamed_masked_cmap.label.gii" \
		"$OUT_DIR2/sub-${i}_91k_intertask_Sac_Pur_lh_pursuit.shape.gii" -name Pursuit -map 2 && \
		wb_command -gifti-label-to-roi "$OUT_DIR2/sub-${i}_91k_intertask_Sac_Pur_rh_renamed_masked_cmap.label.gii" \
		"$OUT_DIR2/sub-${i}_91k_intertask_Sac_Pur_rh_pursuit.shape.gii" -name Pursuit -map 2;

		# Saccade
		wb_command -gifti-label-to-roi "$OUT_DIR2/sub-${i}_91k_intertask_Sac_Pur_lh_renamed_masked_cmap.label.gii" \
		"$OUT_DIR2/sub-${i}_91k_intertask_Sac_Pur_lh_saccade.shape.gii" -name Saccade -map 3 && \
		wb_command -gifti-label-to-roi "$OUT_DIR2/sub-${i}_91k_intertask_Sac_Pur_rh_renamed_masked_cmap.label.gii" \
		"$OUT_DIR2/sub-${i}_91k_intertask_Sac_Pur_rh_saccade.shape.gii" -name Saccade -map 3;

		# Pursuit & Saccade
		wb_command -gifti-label-to-roi "$OUT_DIR2/sub-${i}_91k_intertask_Sac_Pur_lh_renamed_masked_cmap.label.gii" \
		"$OUT_DIR2/sub-${i}_91k_intertask_Sac_Pur_lh_pursuit-saccade.shape.gii" -name Pursuit_and_Saccade -map 4 && \
		wb_command -gifti-label-to-roi "$OUT_DIR2/sub-${i}_91k_intertask_Sac_Pur_rh_renamed_masked_cmap.label.gii" \
		"$OUT_DIR2/sub-${i}_91k_intertask_Sac_Pur_rh_pursuit-saccade.shape.gii" -name Pursuit_and_Saccade -map 4;

		# Vision
		wb_command -gifti-label-to-roi "$OUT_DIR2/sub-${i}_91k_intertask_Sac_Pur_lh_renamed_masked_cmap.label.gii" \
		"$OUT_DIR2/sub-${i}_91k_intertask_Sac_Pur_lh_vision.shape.gii" -name Vision -map 5 && \
		wb_command -gifti-label-to-roi "$OUT_DIR2/sub-${i}_91k_intertask_Sac_Pur_rh_renamed_masked_cmap.label.gii" \
		"$OUT_DIR2/sub-${i}_91k_intertask_Sac_Pur_rh_vision.shape.gii" -name Vision -map 5;

		# Vision & Pursuit
		wb_command -gifti-label-to-roi "$OUT_DIR2/sub-${i}_91k_intertask_Sac_Pur_lh_renamed_masked_cmap.label.gii" \
		"$OUT_DIR2/sub-${i}_91k_intertask_Sac_Pur_lh_vision-pursuit.shape.gii" -name Vision_and_Pursuit -map 6 && \
		wb_command -gifti-label-to-roi "$OUT_DIR2/sub-${i}_91k_intertask_Sac_Pur_rh_renamed_masked_cmap.label.gii" \
		"$OUT_DIR2/sub-${i}_91k_intertask_Sac_Pur_rh_vision-pursuit.shape.gii" -name Vision_and_Pursuit -map 6;

		# Vision & Saccade
		wb_command -gifti-label-to-roi "$OUT_DIR2/sub-${i}_91k_intertask_Sac_Pur_lh_renamed_masked_cmap.label.gii" \
		"$OUT_DIR2/sub-${i}_91k_intertask_Sac_Pur_lh_vision-saccade.shape.gii" -name Vision_and_Saccade -map 7 && \
		wb_command -gifti-label-to-roi "$OUT_DIR2/sub-${i}_91k_intertask_Sac_Pur_rh_renamed_masked_cmap.label.gii" \
		"$OUT_DIR2/sub-${i}_91k_intertask_Sac_Pur_rh_vision-saccade.shape.gii" -name Vision_and_Saccade -map 7;

		# Vision & Pursuit & Saccade
		wb_command -gifti-label-to-roi "$OUT_DIR2/sub-${i}_91k_intertask_Sac_Pur_lh_renamed_masked_cmap.label.gii" \
		"$OUT_DIR2/sub-${i}_91k_intertask_Sac_Pur_lh_vision-pursuit-saccade.shape.gii" -name Vision_and_Pursuit_and_Saccade -map 8 && \
		wb_command -gifti-label-to-roi "$OUT_DIR2/sub-${i}_91k_intertask_Sac_Pur_rh_renamed_masked_cmap.label.gii" \
		"$OUT_DIR2/sub-${i}_91k_intertask_Sac_Pur_rh_vision-pursuit-saccade.shape.gii" -name Vision_and_Pursuit_and_Saccade -map 8;

		### 6. Create two additional files: saccade minus pursuit, pursuit minus saccade

		# wb_command -cifti-merge out.dtseries.nii -cifti first.dtseries.nii -index 1 -cifti second.dtseries.nii;

		### 7. Convert back to label files (the empty string "" tells the command to take whatever is non-empty in the metric file)

		wb_command -metric-label-import "$OUT_DIR2/sub-${i}_91k_intertask_Sac_Pur_lh_vision-pursuit-saccade.shape.gii" "" \
		"$OUT_DIR2/sub-${i}_91k_intertask_Sac_Pur_lh_vision-pursuit-saccade.label.gii";
		wb_command -metric-label-import "$OUT_DIR2/sub-${i}_91k_intertask_Sac_Pur_rh_vision-pursuit-saccade.shape.gii" "" \
		"$OUT_DIR2/sub-${i}_91k_intertask_Sac_Pur_rh_vision-pursuit-saccade.label.gii";

		### Next we'll rename the label keys to make it nicer (this part didn't work as expected - probably need to use additional flags)

		wb_command -metric-label-import "$OUT_DIR2/sub-${i}_91k_intertask_Sac_Pur_lh_vision-pursuit-saccade.shape.gii" \
		7_vision-pursuit-saccade.txt "$OUT_DIR2/sub-${i}_91k_intertask_Sac_Pur_lh_vision-pursuit-saccade.label.gii";
		wb_command -metric-label-import "$OUT_DIR2/sub-${i}_91k_intertask_Sac_Pur_rh_vision-pursuit-saccade.shape.gii" \
		7_vision-pursuit-saccade.txt "$OUT_DIR2/sub-${i}_91k_intertask_Sac_Pur_rh_vision-pursuit-saccade.label.gii";

	done

# Change all files permissions
chmod -Rf 771 *
chgrp -Rf 771 *
