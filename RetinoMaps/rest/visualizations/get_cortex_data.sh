#!/bin/bash

### Get GIFTI files to mask with nibabel

DATA_PATH_CORR="/scratch/mszinte/data/RetinoMaps/derivatives/pp_data/group/91k/rest/hollow_seed_viz_full_corr"

for ROI in mPCS sPCS iPCS sIPS iIPS hMT+ VO LO V3AB V3 V2 V1; do

	wb_command -cifti-separate "group_ses-01_task-rest_space-fsLR_den-91k_desc-full_corr_${ROI}_median.dscalar.nii" COLUMN \
	-metric CORTEX_LEFT "group_ses-01_task-rest_space-fsLR_den-91k_desc-full_corr_lh_${ROI}_median.shape.gii" \
	-metric CORTEX_RIGHT "group_ses-01_task-rest_space-fsLR_den-91k_desc-full_corr_rh_${ROI}_median.shape.gii";
	
done	
