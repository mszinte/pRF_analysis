#!/bin/bash

#####################################################
# Written by Marco Bedini (marco.bedini@univ-amu.fr)
# Goal: get seed vertices for each macro-region e.g., sPCS
#####################################################

# Get the paths
TASK_RESULTS="/scratch/mszinte/data/RetinoMaps/derivatives/pp_data"
ATLAS="/scratch/mszinte/data/RetinoMaps/derivatives/pp_data/atlas/macro_regions"

# Iterate through subjects
for i in 01 02 03 04 05 06 07 08 09 11 12 13 14 17 20 21 22 23 24 25; 
do

  # Mask all active vertices by MMP1 macro-regions
  OUT_DIR="$TASK_RESULTS/sub-${i}/91k/rest/seed"
  
  # Left Hemisphere masks
  for region in lh_mPCS lh_sPCS lh_iPCS lh_sIPS lh_iIPS lh_hMT+ lh_VO lh_LO lh_V3AB lh_V3 lh_V2 lh_V1; do
    wb_command -metric-mask "$OUT_DIR/target_91k/sub-${i}_91k_intertask_Sac-Pur-pRF_lh_largest.shape.gii" \
    "$ATLAS/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_discarded_${region}.shape.gii" \
    "$OUT_DIR/sub-${i}_91k_intertask_Sac-Pur-pRF_${region}.shape.gii"
  done

  # Right Hemisphere masks
  for region in rh_mPCS rh_sPCS rh_iPCS rh_sIPS rh_iIPS rh_hMT+ rh_VO rh_LO rh_V3AB rh_V3 rh_V2 rh_V1; do
    wb_command -metric-mask "$OUT_DIR/target_91k/sub-${i}_91k_intertask_Sac-Pur-pRF_rh_largest.shape.gii" \
    "$ATLAS/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_discarded_${region}.shape.gii" \
    "$OUT_DIR/sub-${i}_91k_intertask_Sac-Pur-pRF_${region}.shape.gii"
  done
  
done

# Change file permissions
chmod -Rf 771 "$OUT_DIR"
chgrp -Rf 327 "$OUT_DIR"
