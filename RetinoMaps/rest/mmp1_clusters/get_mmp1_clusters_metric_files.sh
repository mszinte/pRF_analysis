#!/bin/bash

#####################################################
# Written by Marco Bedini (marco.bedini@univ-amu.fr)
#####################################################

# Output dir
mkdir clusters

# Define base paths
ATLAS_FILE="atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_dseg.dlabel.nii"
LABELS_DIR="label_files"
CLUSTERS_DIR="clusters"

# Define regions
REGIONS=(mPCS sPCS iPCS sIPS iIPS hMT+ VO LO V3AB V3 V2 V1)
HEMISPHERES=(lh rh)

# 1. Merge labels for all macro-regions
for hemi in "${HEMISPHERES[@]}"; do
  for region in "${REGIONS[@]}"; do
    wb_command -cifti-label-import "$ATLAS_FILE" "$LABELS_DIR/${hemi}_${region}_labels.txt" \
      "$CLUSTERS_DIR/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_discarded_${hemi}_${region}.dlabel.nii" -discard-others;
  done
done

# 2. Get the corresponding metric files
for hemi in "${HEMISPHERES[@]}"; do
  CORTEX="CORTEX_LEFT"
  [ "$hemi" == "rh" ] && CORTEX="CORTEX_RIGHT"
  
  for region in "${REGIONS[@]}"; do
    wb_command -cifti-separate "$CLUSTERS_DIR/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_discarded_${hemi}_${region}.dlabel.nii" COLUMN \
      -metric "$CORTEX" "$CLUSTERS_DIR/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_discarded_${hemi}_${region}.shape.gii";
  done
done

# Change permissions
chmod -Rf 771 *
chgrp -Rf 771 *

