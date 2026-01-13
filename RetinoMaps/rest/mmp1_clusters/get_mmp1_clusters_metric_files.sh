#!/bin/bash

#####################################################
# Written by Marco Bedini (marco.bedini@univ-amu.fr)
#####################################################

# Define paths
ATLAS_DIR="/home/${USER}/projects/pRF_analysis/RetinoMaps/rest/mmp1_clusters"
LABELS_DIR="$ATLAS_DIR/label_files"

# Output dir
# mkdir -p "$ATLAS_DIR/clusters"

# Define regions and hemispheres
REGIONS=(mPCS sPCS iPCS sIPS iIPS hMT+ VO LO V3AB V3 V2 V1)
HEMISPHERES=(lh rh)

# 1. Merge labels for all macro-regions
for hemi in "${HEMISPHERES[@]}"; do
  for region in "${REGIONS[@]}"; do
    wb_command -cifti-label-import atlas-Glasser_space-fsLR_den-32k_dseg.dlabel.nii "$LABELS_DIR/${hemi}_${region}_labels.txt" \
      "$ATLAS_DIR/clusters/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_discarded_${hemi}_${region}.dlabel.nii" -discard-others
  
  # 2. Get the corresponding metric files and binarize outputs

		# Separate metric files
		wb_command -cifti-separate "$ATLAS_DIR/clusters/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_discarded_${hemi}_${region}.dlabel.nii" COLUMN \
		-metric CORTEX_LEFT "$ATLAS_DIR/clusters/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_discarded_${hemi}_${region}.shape.gii" \
		-metric CORTEX_RIGHT "$ATLAS_DIR/clusters/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_discarded_${hemi}_${region}.shape.gii"

	# Binarize outputs
	wb_command -cifti-reduce "$ATLAS_DIR/clusters/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_discarded_${hemi}_${region}.dlabel.nii" -direction ROW \
	COUNT_NONZERO "$ATLAS_DIR/clusters/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_discarded_${hemi}_${region}_bin.shape.gii"
	done

done

# Change permissions
chmod -Rf 771 *
chgrp -Rf 771 *
