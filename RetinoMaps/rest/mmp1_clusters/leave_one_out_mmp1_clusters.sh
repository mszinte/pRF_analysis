#!/bin/bash

#####################################################
# Written by Marco Bedini (marco.bedini@univ-amu.fr)
#####################################################

# Define base paths
LABELS_DIR="label_files"
CLUSTERS_DIR="clusters"

# Input file
ATLAS_FILE="atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_dseg.dlabel.nii"

# Output dir
LEAVEOUT_DIR="leaveout"

# Define the clusters and hemisphere suffix
CLUSTERS=(mPCS sPCS iPCS sIPS iIPS hMT+ VO LO V3AB V3 V2 V1)
HEMISPHERES=(lh rh)

# 1a. Generate "leave-one-out" cluster label per hemisphere
mkdir -p "$LEAVEOUT_DIR"
for hemi in "${HEMISPHERES[@]}"; do
  for leave_out in "${CLUSTERS[@]}"; do
    COMBINED_LABEL_FILE=$(mktemp)
    
    echo "💡 Creating label file excluding $leave_out for hemisphere $hemi"
    for region in "${CLUSTERS[@]}"; do
      if [ "$region" != "$leave_out" ]; then
        cat "$LABELS_DIR/${hemi}_${region}_labels.txt" >> "$COMBINED_LABEL_FILE"
      fi
    done
    OUT_LABEL="${LEAVEOUT_DIR}/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_leaveout_${hemi}_${leave_out}.dlabel.nii"
    wb_command -cifti-label-import "$ATLAS_FILE" "$COMBINED_LABEL_FILE" "$OUT_LABEL" -discard-others
    rm "$COMBINED_LABEL_FILE"
  done
done

# 1b. Generate "leave-one-out" bilateral cluster
for leave_out in "${CLUSTERS[@]}"; do
  COMBINED_LABEL_FILE=$(mktemp)
  echo "💡 Creating BILATERAL label file excluding $leave_out"
  for hemi in "${HEMISPHERES[@]}"; do
    for region in "${CLUSTERS[@]}"; do
      if [ "$region" != "$leave_out" ]; then
        cat "$LABELS_DIR/${hemi}_${region}_labels.txt" >> "$COMBINED_LABEL_FILE"
      fi
    done
  done
  OUT_LABEL="${LEAVEOUT_DIR}/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_leaveout_${leave_out}.dlabel.nii"
  wb_command -cifti-label-import "$ATLAS_FILE" "$COMBINED_LABEL_FILE" "$OUT_LABEL" -discard-others
  rm "$COMBINED_LABEL_FILE"
done

# 2. Extract metric files from the bilateral leave-one-out labels
for leave_out in "${CLUSTERS[@]}"; do
  DLABEL_FILE="${LEAVEOUT_DIR}/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_leaveout_${leave_out}.dlabel.nii"
  
  wb_command -cifti-separate "$DLABEL_FILE" COLUMN \
    -metric CORTEX_LEFT "${LEAVEOUT_DIR}/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_leaveout_lh_${leave_out}.shape.gii" \
    -metric CORTEX_RIGHT "${LEAVEOUT_DIR}/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_leaveout_rh_${leave_out}.shape.gii"
done

# 3. Convert all labels to ROIs (dscalar format)
echo "💡 Converting labels to ROIs..."
for leave_out in "${CLUSTERS[@]}"; do
  DLABEL_FILE="${LEAVEOUT_DIR}/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_leaveout_${leave_out}.dlabel.nii"
  DSCALAR_FILE="${LEAVEOUT_DIR}/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_leaveout_${leave_out}.dscalar.nii"
  
  wb_command -cifti-all-labels-to-rois "$DLABEL_FILE" 1 "$DSCALAR_FILE"
done

# 4. Create binary masks from ROIs
echo "💡 Creating binary masks..."
for leave_out in "${CLUSTERS[@]}"; do
  DSCALAR_FILE="${LEAVEOUT_DIR}/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_leaveout_${leave_out}.dscalar.nii"
  BIN_FILE="${LEAVEOUT_DIR}/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_leaveout_${leave_out}_bin.dscalar.nii"
  
  wb_command -cifti-reduce "$DSCALAR_FILE" COUNT_NONZERO "$BIN_FILE"
done

echo "✅ All label, metric, ROI, and binary mask files created!"
# Change permissions
chmod -Rf 771 *
chgrp -Rf 771 *
