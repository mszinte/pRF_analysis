#!/bin/bash

# Get the paths
TASK_RESULTS="/scratch/mszinte/data/RetinoMaps/derivatives/pp_data"
ATLAS="~/projects/pRF_analysis/RetinoMaps/rest/mmp1_clusters/clusters"

# Iterate through subjects
for i in 01 02 03 04 05 06 07 08 09 11 12 13 14 17 20 21 22 23 24 25 170k; 
do

  # Mask all active vertices by MMP1 clusters
  
  OUT_DIR="/scratch/mszinte/data/RetinoMaps/derivatives/pp_data/sub-${i}/91k/rest/seed"
  
  # Left Hemisphere masks
  for region in lh_mPCS lh_sPCS lh_iPCS lh_sIPS lh_iIPS lh_hMT+ lh_VO lh_LO lh_V3AB lh_V3 lh_V2 lh_V1; do
    wb_command -metric-mask "$TASK_RESULTS/sub-${i}/91k/rest/seed/target_91k/sub-${i}_91k_intertask_Sac_Pur_lh_vision-pursuit-saccade.shape.gii" \
    "$ATLAS/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_discarded_${region}.shape.gii" \
    "$OUT_DIR/sub-${i}_91k_intertask_Sac_Pur_vision-pursuit-saccade_${region}.shape.gii"
  done

  # Right Hemisphere masks
  for region in rh_mPCS rh_sPCS rh_iPCS rh_sIPS rh_iIPS rh_hMT+ rh_VO rh_LO rh_V3AB rh_V3 rh_V2 rh_V1; do
    wb_command -metric-mask "$TASK_RESULTS/sub-${i}_91k_intertask_Sac_Pur_rh_vision-pursuit-saccade.shape.gii" \
    "$ATLAS/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_discarded_${region}.shape.gii" \
    "$OUT_DIR/sub-${i}_91k_intertask_Sac_Pur_vision-pursuit-saccade_${region}.shape.gii"
  done
done

