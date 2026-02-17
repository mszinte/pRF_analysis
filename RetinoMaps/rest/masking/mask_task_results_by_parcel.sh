#!/bin/bash

#####################################################
# Written by Marco Bedini (marco.bedini@univ-amu.fr)
#####################################################

# Get the paths
TASK_RESULTS="/scratch/mszinte/data/RetinoMaps/derivatives/pp_data"
ATLAS="/scratch/mszinte/data/RetinoMaps/derivatives/pp_data/atlas/mmp1_clusters/parcels"

# Subjects
SUBJECTS=(01 02 03 04 05 06 07 08 09 11 12 13 14 17 20 21 22 23 24 25 170k)

# MMP1 parcels within the macro-regions
ROIS=(
  L_V1_ROI L_MST_ROI L_V2_ROI L_V3_ROI L_V4_ROI L_V8_ROI L_FEF_ROI L_PEF_ROI L_55b_ROI L_V3A_ROI
  L_V7_ROI L_IPS1_ROI L_FFC_ROI L_V3B_ROI L_LO1_ROI L_LO2_ROI L_PIT_ROI L_MT_ROI L_7Pm_ROI L_24dv_ROI
  L_7AL_ROI L_SCEF_ROI L_6ma_ROI L_7Am_ROI L_7PL_ROI L_7PC_ROI L_LIPv_ROI L_VIP_ROI L_MIP_ROI
  L_6d_ROI L_6mp_ROI L_6v_ROI L_p32pr_ROI L_6r_ROI L_IFJa_ROI L_IFJp_ROI L_LIPd_ROI L_6a_ROI
  L_i6-8_ROI L_AIP_ROI L_PH_ROI L_IP2_ROI L_IP1_ROI L_IP0_ROI L_V6A_ROI L_VMV1_ROI L_VMV3_ROI
  L_V4t_ROI L_FST_ROI L_V3CD_ROI L_LO3_ROI L_VMV2_ROI L_VVC_ROI
  R_V1_ROI R_MST_ROI R_V2_ROI R_V3_ROI R_V4_ROI R_V8_ROI R_FEF_ROI R_PEF_ROI R_55b_ROI R_V3A_ROI
  R_V7_ROI R_IPS1_ROI R_FFC_ROI R_V3B_ROI R_LO1_ROI R_LO2_ROI R_PIT_ROI R_MT_ROI R_7Pm_ROI R_24dv_ROI
  R_7AL_ROI R_SCEF_ROI R_6ma_ROI R_7Am_ROI R_7PL_ROI R_7PC_ROI R_LIPv_ROI R_VIP_ROI R_MIP_ROI
  R_6d_ROI R_6mp_ROI R_6v_ROI R_p32pr_ROI R_6r_ROI R_IFJa_ROI R_IFJp_ROI R_LIPd_ROI R_6a_ROI
  R_i6-8_ROI R_AIP_ROI R_PH_ROI R_IP2_ROI R_IP1_ROI R_IP0_ROI R_V6A_ROI R_VMV1_ROI R_VMV3_ROI
  R_V4t_ROI R_FST_ROI R_V3CD_ROI R_LO3_ROI R_VMV2_ROI R_VVC_ROI
)

# Iterate through subjects
for subj in "${SUBJECTS[@]}"; do

  OUT_DIR="$TASK_RESULTS/sub-${subj}/91k/rest/seed"
  mkdir -p "$OUT_DIR"

  # Loop through parcels
  for roi in "${ROIS[@]}"; do
    hemi=${roi:0:1}   # "L" or "R"
    
    if [ "$hemi" == "L" ]; then
      INPUT="$OUT_DIR/target_91k/sub-${subj}_91k_intertask_Sac_Pur_lh_vision-pursuit-saccade.shape.gii"
    else
      INPUT="$OUT_DIR/target_91k/sub-${subj}_91k_intertask_Sac_Pur_rh_vision-pursuit-saccade.shape.gii"
    fi

    wb_command -metric-mask \
      "$INPUT" \
      "$ATLAS/${roi}.shape.gii" \
      "$OUT_DIR/sub-${subj}_91k_intertask_Sac_Pur_vision-pursuit-saccade_${roi}.shape.gii"
  done
done

# Fix permissions
chmod -Rf 771 "$OUT_DIR"
chgrp -Rf 771 "$OUT_DIR"
