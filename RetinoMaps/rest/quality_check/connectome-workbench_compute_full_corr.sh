#!/bin/bash

#####################################################
# Goal of the script:
# Compute parcellated full correlations using Connectome Workbench v. 2.1
# Written by Marco Bedini (marco.bedini@univ-amu.fr)
#####################################################

# Define some paths
TASK_RESULTS="/scratch/mszinte/data/RetinoMaps/derivatives/xcp-d_w_atlas"
ATLAS="/scratch/mszinte/data/RetinoMaps/derivatives/pp_data/atlas"
OUT_DIR="/scratch/mszinte/data/RetinoMaps/derivatives/pp_data/quality_check"

# Used to check subjects with odd denoising (ids 03, 04, 14, 21, 23) where in some cases XCP-D failed the parcellate steps
# Affects only run-02 but we look also at the concatenated timeseries to understand if z-scoring removes the issue 

for i in 01 02 03 04 05 06 07 08 09 11 12 13 14 17 20 21 22 23 24 25;
do

	# Full correlation (run-01)
	wb_command -cifti-correlation \
	    "$TASK_RESULTS/sub-${i}/ses-01/func/sub-${i}_ses-01_task-rest_run-01_space-fsLR_den-91k_desc-denoised_bold.dtseries.nii" \
	    "$OUT_DIR/sub-${i}_task-rest_run-01_space-fsLR_den-91k_desc-full_corr.dconn.nii" \
	    -mem-limit 12

	# Full correlation (run-02)
	wb_command -cifti-correlation \
	    "$TASK_RESULTS/sub-${i}/ses-01/func/sub-${i}_ses-01_task-rest_run-02_space-fsLR_den-91k_desc-denoised_bold.dtseries.nii" \
	    "$OUT_DIR/sub-${i}_task-rest_run-02_space-fsLR_den-91k_desc-full_corr.dconn.nii" \
	    -mem-limit 12
	
	# Full correlation (for concat timeseries)
	wb_command -cifti-correlation \
	    "$TASK_RESULTS/sub-${i}/ses-01/func/sub-${i}_ses-01_task-rest_space-fsLR_den-91k_desc-denoised_bold.dtseries.nii" \
	    "$OUT_DIR/sub-${i}_task-rest_space-fsLR_den-91k_desc-full_corr.dconn.nii" \
	    -mem-limit 12

done

# Grant permissions to output files
chmod -Rf 771 "$OUT_DIR"
chgrp -Rf 327 "$OUT_DIR"
