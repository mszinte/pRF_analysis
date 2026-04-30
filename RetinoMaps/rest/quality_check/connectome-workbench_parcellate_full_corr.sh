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
	# Parcellate targets in two steps:
	# The first produces a dense by parcel matrix
	# The second parcel by parcel
	# See wb_command -cifti-help for a shorthand on the typical mappings
	wb_command -cifti-parcellate \
	    "$OUT_DIR/sub-${i}_task-rest_run-01_space-fsLR_den-91k_desc-full_corr.dconn.nii" \
	    "$ATLAS/atlas-Glasser_space-fsLR_den-32k_dseg.dlabel.nii" COLUMN \
	    "$OUT_DIR/sub-${i}_task-rest_run-01_space-fsLR_den-91k_desc-full_corr_parcellated.dpconn.nii" \
	    -method MEAN

	# Parcellate targets
	wb_command -cifti-parcellate \
	    "$OUT_DIR/sub-${i}_task-rest_run-01_space-fsLR_den-91k_desc-full_corr.dconn.nii" \
	    "$ATLAS/atlas-Glasser_space-fsLR_den-32k_dseg.dlabel.nii" ROW \
	    "$OUT_DIR/sub-${i}_task-rest_run-01_space-fsLR_den-91k_desc-full_corr_parcellated.pconn.nii" \
	    -method MEAN

	# Convert all parcellated outputs to TSV for visualizations
    wb_command -cifti-convert -to-text \
        "$OUT_DIR/sub-${i}_task-rest_run-01_space-fsLR_den-91k_desc-full_corr_parcellated.pconn.nii" \
        "$OUT_DIR/sub-${i}_task-rest_run-01_space-fsLR_den-91k_desc-full_corr_parcellated.tsv"

	# Remove files that occupy excessive memory space
    # rm "$OUT_DIR/sub-${i}_task-rest_run-01_space-fsLR_den-91k_desc-full_corr.dconn.nii"
	
	wb_command -cifti-parcellate \
	    "$OUT_DIR/sub-${i}_task-rest_run-02_space-fsLR_den-91k_desc-full_corr.dconn.nii" \
	    "$ATLAS/atlas-Glasser_space-fsLR_den-32k_dseg.dlabel.nii" COLUMN \
	    "$OUT_DIR/sub-${i}_task-rest_run-02_space-fsLR_den-91k_desc-full_corr_parcellated.dpconn.nii" \
	    -method MEAN

	# Parcellate targets
	wb_command -cifti-parcellate \
	    "$OUT_DIR/sub-${i}_task-rest_run-02_space-fsLR_den-91k_desc-full_corr.dconn.nii" \
	    "$ATLAS/atlas-Glasser_space-fsLR_den-32k_dseg.dlabel.nii" ROW \
	    "$OUT_DIR/sub-${i}_task-rest_run-02_space-fsLR_den-91k_desc-full_corr_parcellated.pconn.nii" \
	    -method MEAN

	# Convert all parcellated outputs to TSV for visualizations
    wb_command -cifti-convert -to-text \
        "$OUT_DIR/sub-${i}_task-rest_run-02_space-fsLR_den-91k_desc-full_corr_parcellated.pconn.nii" \
        "$OUT_DIR/sub-${i}_task-rest_run-02_space-fsLR_den-91k_desc-full_corr_parcellated.tsv"

	# Remove files that occupy excessive memory space
    # rm "$OUT_DIR/sub-${i}_task-rest_run-02_space-fsLR_den-91k_desc-full_corr.dconn.nii"

	# Parcellate targets
	wb_command -cifti-parcellate \
	    "$OUT_DIR/sub-${i}_task-rest_space-fsLR_den-91k_desc-full_corr.dconn.nii" \
	    "$ATLAS/atlas-Glasser_space-fsLR_den-32k_dseg.dlabel.nii" COLUMN \
	    "$OUT_DIR/sub-${i}_task-rest_space-fsLR_den-91k_desc-full_corr_parcellated.dpconn.nii" \
	    -method MEAN

	# Parcellate targets
	wb_command -cifti-parcellate \
	    "$OUT_DIR/sub-${i}_task-rest_space-fsLR_den-91k_desc-full_corr.dconn.nii" \
	    "$ATLAS/atlas-Glasser_space-fsLR_den-32k_dseg.dlabel.nii" ROW \
	    "$OUT_DIR/sub-${i}_task-rest_space-fsLR_den-91k_desc-full_corr_parcellated.pconn.nii" \
	    -method MEAN

	# Convert all parcellated outputs to TSV for visualizations
    wb_command -cifti-convert -to-text \
        "$OUT_DIR/sub-${i}_task-rest_space-fsLR_den-91k_desc-full_corr_parcellated.pconn.nii" \
        "$OUT_DIR/sub-${i}_task-rest_space-fsLR_den-91k_desc-full_corr_parcellated.tsv"

	# Remove files that occupy excessive memory space
    # rm "$OUT_DIR/sub-${i}_task-rest_space-fsLR_den-91k_desc-full_corr.dconn.nii"

done

# Grant permissions to output files
chmod -Rf 771 "$OUT_DIR"
chgrp -Rf 327 "$OUT_DIR"
