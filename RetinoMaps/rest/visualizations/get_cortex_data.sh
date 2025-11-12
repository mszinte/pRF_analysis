#!/bin/bash

### Get surfaces to use as input for nibabel winner take all function

DATA_PATH_CORR="/scratch/mszinte/data/RetinoMaps/derivatives/pp_data/group/91k/rest/hollow_seed_viz_full_corr"
OUT_PATH="/scratch/mszinte/data/RetinoMaps/derivatives/pp_data/group/91k/rest/wta/surfaces"

for ROI in mPCS sPCS iPCS sIPS iIPS hMT+ VO LO V3AB V3 V2 V1; do

	wb_command -cifti-separate "$DATA_PATH_CORR/group_ses-01_task-rest_space-fsLR_den-91k_desc-full_corr_${ROI}_hollow_seed.dscalar.nii" COLUMN \
	-metric CORTEX_LEFT "$OUT_PATH/group_ses-01_task-rest_space-fsLR_den-91k_desc-full_corr_lh_${ROI}_hollow_seed.shape.gii" \
	-metric CORTEX_RIGHT "$OUT_PATH/group_ses-01_task-rest_space-fsLR_den-91k_desc-full_corr_rh_${ROI}_hollow_seed.shape.gii";
	
done	
