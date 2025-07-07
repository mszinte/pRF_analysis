#!/bin/bash

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

    OUT_LABEL="${LEAVEOUT_DIR}/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_leaveout_${leave_out}_${hemi}.dlabel.nii"

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
    -metric CORTEX_LEFT "${LEAVEOUT_DIR}/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_leaveout_${leave_out}_${hemi}.shape.gii" \
    -metric CORTEX_RIGHT "${LEAVEOUT_DIR}/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_leaveout_${leave_out}_${hemi}.shape.gii"
done

echo "✅ All label and metric files created!"

# Change permissions
chmod -Rf 771 *
chgrp -Rf 771 *
