#!/bin/bash

#####################################################
# Goal of the script:
# Compute full correlation by vertex intra-hemispherically using Connectome Workbench v. 2.1
# Outputs are parcellated using the mean method and masked for supplementary info visualizations
# 
# Written by Marco Bedini (marco.bedini@univ-amu.fr)
#####################################################

# Define some paths
TASK_RESULTS="/scratch/mszinte/data/RetinoMaps/derivatives/xcp-d_w_atlas"
ATLAS="/scratch/mszinte/data/RetinoMaps/derivatives/pp_data/atlas"

# Iterate through subjects with odd denoising
# Affects only run-02 but we look also at the concatenated timeseries
# To understand if z-scoring removes it 

for i in 03 04 14 21 23;
do

OUT_DIR="/scratch/mszinte/data/RetinoMaps/derivatives/pp_data/quality_check"

	# Full correlation
	wb_command -cifti-correlation \
	    "$TASK_RESULTS/sub-${i}/ses-01/func/sub-${i}_ses-01_task-rest_run-02_space-fsLR_den-91k_desc-denoised_bold.dtseries.nii" \
	    "$OUT_DIR/sub-${i}_task-rest_run-02_space-fsLR_den-91k_desc-full_corr.dconn.nii" \
	    -mem-limit 30

	# Parcellate targets
	wb_command -cifti-parcellate \
	    "$OUT_DIR/sub-${i}_task-rest_run-02_space-fsLR_den-91k_desc-full_corr.dscalar.nii" \
	    "$ATLAS/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_dseg.dlabel.nii" COLUMN \
	    "$OUT_DIR/sub-${i}_task-rest_run-02_space-fsLR_den-91k_desc-full_corr_parcellated.pscalar.nii" \
	    -method MEAN -only-numeric;

	# Convert all parcellated outputs to TSV for visualizations
    wb_command -cifti-convert -to-text \
        "$OUT_DIR/sub-${i}_task-rest_run-02_space-fsLR_den-91k_desc-full_corr_parcellated.pscalar.nii" \
        "$OUT_DIR/sub-${i}_task-rest_run-02_space-fsLR_den-91k_desc-full_corr_parcellated.tsv"

	# Remove files that occupy excessive memory space
    rm "$OUT_DIR/sub-${i}_task-rest_run-02_space-fsLR_den-91k_desc-full_corr.dconn.nii"
	
	# Full correlation (for concat timeseries)
	wb_command -cifti-correlation \
	    "$TASK_RESULTS/sub-${i}/91k/rest/timeseries/sub-${i}_ses-01_task-rest_space-fsLR_den-91k_desc-denoised_bold.dtseries.nii" \
	    "$OUT_DIR/sub-${i}_task-rest_space-fsLR_den-91k_desc-full_corr.dconn.nii" \
	    -mem-limit 30

	# Parcellate targets
	wb_command -cifti-parcellate \
	    "$OUT_DIR/sub-${i}_task-rest_space-fsLR_den-91k_desc-full_corr.dscalar.nii" \
	    "$ATLAS/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_dseg.dlabel.nii" COLUMN \
	    "$OUT_DIR/sub-${i}_task-rest_space-fsLR_den-91k_desc-full_corr_parcellated.pscalar.nii" \
	    -method MEAN -only-numeric;

	# Convert all parcellated outputs to TSV for visualizations
    wb_command -cifti-convert -to-text \
        "$OUT_DIR/sub-${i}_task-rest_space-fsLR_den-91k_desc-full_corr_parcellated.pscalar.nii" \
        "$OUT_DIR/sub-${i}_task-rest_space-fsLR_den-91k_desc-full_corr_parcellated.tsv"

	# Remove files that occupy excessive memory space
    rm "$OUT_DIR/sub-${i}_task-rest_space-fsLR_den-91k_desc-full_corr.dconn.nii"

done

# Grant permissions to output files
chmod -Rf 771 "$OUT_DIR"
chgrp -Rf 327 "$OUT_DIR"
