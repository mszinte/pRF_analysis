#!/bin/bash

# Get the intertask results path
TASK_RESULTS="/scratch/mszinte/data/RetinoMaps/derivatives/pp_data"

# Think about if we want to have a folder for templates
# TEMPLATE="/projects/pRF_analysis/RetinoMaps/?/fsLR"

# Create input and output folders
mkdir source_170k
mkdir target_91k

for i in 01 02 03 04 05 06 07 08 09 11 12 13 14 17 20 21 22 23 24 25 170k;
do

	### 1. Separate the results by hemisphere

	wb_command -cifti-separate "$TASK_RESULTS/sub-"$i"/170k/intertask/sub-"$i"_intertask_Sac_Pur.dtseries.nii COLUMN -metric CORTEX_LEFT source_170k/sub-"$i"_intertask_Sac_Pur.dtseries.left_hemi.shape.gii";

	wb_command -cifti-separate "$TASK_RESULTS/sub-"$i"/170k/intertask/sub-"$i"_intertask_Sac_Pur.dtseries.nii COLUMN -metric CORTEX_RIGHT source_170k/sub-"$i"_intertask_Sac_Pur.dtseries.right_hemi.shape.gii";

	### 2. Resample from 170k to 91k - try this one for more advanced registration options wb_command -cifti-resample (in the step before maybe use -cifti-separate to get the dense labels dlabel.nii) 

	wb_command -metric-resample \
	    source_170k/sub-"$i"_intertask_Sac_Pur.dtseries.left_hemi.shape.gii \
	    fsLR/tpl-fsLR_hemi-L_den-59k_sphere.surf.gii \
	    fsLR/tpl-fsLR_hemi-L_den-32k_sphere.surf.gii \
	    BARYCENTRIC -largest \
	    target_91k/sub-"$i"_91k_intertask_Sac_Pur.dtseries.left_hemi.bary_largest.shape.gii;
		
	wb_command -metric-resample \
	    source_170k/sub-"$i"_intertask_Sac_Pur.dtseries.right_hemi.shape.gii \
	    fsLR/tpl-fsLR_hemi-R_den-59k_sphere.surf.gii \
	    fsLR/tpl-fsLR_hemi-R_den-32k_sphere.surf.gii \
	    BARYCENTRIC -largest \
	    target_91k/sub-"$i"_91k_intertask_Sac_Pur.dtseries.right_hemi.bary_largest.shape.gii;

	### 3. Import labels to rename them and change the color map (using Uriel's rgb convention)

	wb_command -metric-label-import target_91k/sub-"$i"_91k_intertask_Sac_Pur.dtseries.left_hemi.bary_largest.shape.gii intertask_label_list.txt \
	target_91k/sub-"$i"_91k_intertask_Sac_Pur_lh_renamed_cmap.label.gii && \
	wb_command -metric-label-import target_91k/sub-"$i"_91k_intertask_Sac_Pur.dtseries.right_hemi.bary_largest.shape.gii intertask_label_list.txt \
	target_91k/sub-"$i"_91k_intertask_Sac_Pur_rh_renamed_cmap.label.gii;

	### 4. Mask the GLM results with the Glasser macro-regions (simply comment out these steps if you want to keep the results unmasked)

	wb_command -label-mask target_91k/sub-"$i"_91k_intertask_Sac_Pur_lh_renamed_cmap.label.gii "fsLR/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_left_hemi.shape.gii" \
	target_91k/sub-"$i"_91k_intertask_Sac_Pur_lh_renamed_cmap.label.gii;

	wb_command -label-mask target_91k/sub-"$i"_91k_intertask_Sac_Pur_lh_renamed_cmap.label.gii "fsLR/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_left_hemi.shape.gii" \
	target_91k/sub-"$i"_91k_intertask_Sac_Pur_lh_renamed_cmap.label.gii;

	wb_command -label-mask target_91k/sub-"$i"_91k_intertask_Sac_Pur_rh_renamed_cmap.label.gii "fsLR/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_right_hemi.shape.gii" \
	target_91k/sub-"$i"_91k_intertask_Sac_Pur_rh_renamed_cmap.label.gii;

	wb_command -label-mask target_91k/sub-"$i"_91k_intertask_Sac_Pur_rh_renamed_cmap.label.gii "fsLR/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_right_hemi.shape.gii" \
	target_91k/sub-"$i"_91k_intertask_Sac_Pur_rh_renamed_cmap.label.gii;

	### 5. Loop over each label to generate separate metric files (useful to mask the timeseries)

	# Pursuit
	wb_command -gifti-label-to-roi target_91k/sub-"$i"_91k_intertask_Sac_Pur_lh_renamed_cmap.label.gii \
	target_91k/sub-"$i"_91k_intertask_Sac_Pur_lh_pursuit.shape.gii -name Pursuit -map 2 && \
	wb_command -gifti-label-to-roi target_91k/sub-"$i"_91k_intertask_Sac_Pur_rh_renamed_cmap.label.gii \
	target_91k/sub-"$i"_91k_intertask_Sac_Pur_rh_pursuit.shape.gii -name Pursuit -map 2;

	# Saccade
	wb_command -gifti-label-to-roi target_91k/sub-"$i"_91k_intertask_Sac_Pur_lh_renamed_cmap.label.gii \
	target_91k/sub-"$i"_91k_intertask_Sac_Pur_lh_saccade.shape.gii -name Saccade -map 3 && \
	wb_command -gifti-label-to-roi target_91k/sub-"$i"_91k_intertask_Sac_Pur_rh_renamed_cmap.label.gii \
	target_91k/sub-"$i"_91k_intertask_Sac_Pur_rh_saccade.shape.gii -name Saccade -map 3;

	# Pursuit & Saccade
	wb_command -gifti-label-to-roi target_91k/sub-"$i"_91k_intertask_Sac_Pur_lh_renamed_cmap.label.gii \
	target_91k/sub-"$i"_91k_intertask_Sac_Pur_lh_pursuit-saccade.shape.gii -name Pursuit_and_Saccade -map 4 && \
	wb_command -gifti-label-to-roi target_91k/sub-"$i"_91k_intertask_Sac_Pur_rh_renamed_cmap.label.gii \
	target_91k/sub-"$i"_91k_intertask_Sac_Pur_rh_pursuit-saccade.shape.gii -name Pursuit_and_Saccade -map 4;

	# Vision
	wb_command -gifti-label-to-roi target_91k/sub-"$i"_91k_intertask_Sac_Pur_lh_renamed_cmap.label.gii \
	target_91k/sub-"$i"_91k_intertask_Sac_Pur_lh_vision.shape.gii -name Vision -map 5 && \
	wb_command -gifti-label-to-roi target_91k/sub-"$i"_91k_intertask_Sac_Pur_rh_renamed_cmap.label.gii \
	target_91k/sub-"$i"_91k_intertask_Sac_Pur_rh_vision.shape.gii -name Vision -map 5;

	# Vision & Pursuit
	wb_command -gifti-label-to-roi target_91k/sub-"$i"_91k_intertask_Sac_Pur_lh_renamed_cmap.label.gii \
	target_91k/sub-"$i"_91k_intertask_Sac_Pur_lh_vision-pursuit.shape.gii -name Vision_and_Pursuit -map 6 && \
	wb_command -gifti-label-to-roi target_91k/sub-"$i"_91k_intertask_Sac_Pur_rh_renamed_cmap.label.gii \
	target_91k/sub-"$i"_91k_intertask_Sac_Pur_rh_vision-pursuit.shape.gii -name Vision_and_Pursuit -map 6;

	# Vision & Saccade
	wb_command -gifti-label-to-roi target_91k/sub-"$i"_91k_intertask_Sac_Pur_lh_renamed_cmap.label.gii \
	target_91k/sub-"$i"_91k_intertask_Sac_Pur_lh_vision-saccade.shape.gii -name Vision_and_Saccade -map 7 && \
	wb_command -gifti-label-to-roi target_91k/sub-"$i"_91k_intertask_Sac_Pur_rh_renamed_cmap.label.gii \
	target_91k/sub-"$i"_91k_intertask_Sac_Pur_rh_vision-saccade.shape.gii -name Vision_and_Saccade -map 7;

	# Vision & Pursuit & Saccade
	wb_command -gifti-label-to-roi target_91k/sub-"$i"_91k_intertask_Sac_Pur_lh_renamed_cmap.label.gii \
	target_91k/sub-"$i"_91k_intertask_Sac_Pur_lh_vision-pursuit-saccade.shape.gii -name Vision_and_Pursuit_and_Saccade -map 8 && \
	wb_command -gifti-label-to-roi target_91k/sub-"$i"_91k_intertask_Sac_Pur_rh_renamed_cmap.label.gii \
	target_91k/sub-"$i"_91k_intertask_Sac_Pur_rh_vision-pursuit-saccade.shape.gii -name Vision_and_Pursuit_and_Saccade -map 8;

	### 6. Create two additional files: saccade minus pursuit, pursuit minus saccade

	# wb_command -cifti-merge out.dtseries.nii -cifti first.dtseries.nii -index 1 -cifti second.dtseries.nii;

	### 7. Convert back to label files (the empty string "" tells the command to take whatever is non-empty in the metric file)

	wb_command -metric-label-import target_91k/sub-"$i"_91k_intertask_Sac_Pur_lh_vision-pursuit-saccade.shape.gii "" \
	target_91k/sub-"$i"_91k_intertask_Sac_Pur_lh_vision-pursuit-saccade.label.gii;
	wb_command -metric-label-import target_91k/sub-"$i"_91k_intertask_Sac_Pur_rh_vision-pursuit-saccade.shape.gii "" \
	target_91k/sub-"$i"_91k_intertask_Sac_Pur_rh_vision-pursuit-saccade.label.gii;

	### Next we'll rename the label keys to make it nicer

	wb_command -metric-label-import target_91k/sub-"$i"_91k_intertask_Sac_Pur_lh_vision-pursuit-saccade.shape.gii \
	labels/7_vision-pursuit-saccade.txt sub-"$i"_91k_intertask_Sac_Pur_lh_vision-pursuit-saccade.label.gii;
	wb_command -metric-label-import target_91k/sub-"$i"_91k_intertask_Sac_Pur_rh_vision-pursuit-saccade.shape.gii \
	labels/7_vision-pursuit-saccade.txt sub-"$i"_91k_intertask_Sac_Pur_rh_vision-pursuit-saccade.label.gii;

done
