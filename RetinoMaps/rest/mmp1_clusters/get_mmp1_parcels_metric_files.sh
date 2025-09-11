#!/bin/bash

# Input atlas
ATLAS_FILE="atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_dseg.dlabel.nii"

# Output directory
OUTDIR="parcels"
mkdir -p "$OUTDIR"

# List of ROI names with corresponding numeric label IDs
# Format: ROI_NAME LABEL_ID
ROIS=(
  "R_V1_ROI 1"
  "R_MST_ROI 2"
  "R_V2_ROI 4"
  "R_V3_ROI 5"
  "R_V4_ROI 6"
  "R_V8_ROI 7"
  "R_FEF_ROI 10"
  "R_PEF_ROI 11"
  "R_55b_ROI 12"
  "R_V3A_ROI 13"
  "R_V7_ROI 16"
  "R_IPS1_ROI 17"
  "R_FFC_ROI 18"
  "R_V3B_ROI 19"
  "R_LO1_ROI 20"
  "R_LO2_ROI 21"
  "R_PIT_ROI 22"
  "R_MT_ROI 23"
  "R_7Pm_ROI 29"
  "R_24dv_ROI 41"
  "R_7AL_ROI 42"
  "R_SCEF_ROI 43"
  "R_6ma_ROI 44"
  "R_7Am_ROI 45"
  "R_7PL_ROI 46"
  "R_7PC_ROI 47"
  "R_LIPv_ROI 48"
  "R_VIP_ROI 49"
  "R_MIP_ROI 50"
  "R_6d_ROI 54"
  "R_6mp_ROI 55"
  "R_6v_ROI 56"
  "R_p32pr_ROI 60"
  "R_6r_ROI 78"
  "R_IFJa_ROI 79"
  "R_IFJp_ROI 80"
  "R_LIPd_ROI 95"
  "R_6a_ROI 96"
  "R_i6-8_ROI 97"
  "R_AIP_ROI 117"
  "R_PH_ROI 138"
  "R_IP2_ROI 144"
  "R_IP1_ROI 145"
  "R_IP0_ROI 146"
  "R_V6A_ROI 152"
  "R_VMV1_ROI 153"
  "R_VMV3_ROI 154"
  "R_V4t_ROI 156"
  "R_FST_ROI 157"
  "R_V3CD_ROI 158"
  "R_LO3_ROI 159"
  "R_VMV2_ROI 160"
  "R_VVC_ROI 163"
  "L_V1_ROI 181"
  "L_MST_ROI 182"
  "L_V2_ROI 184"
  "L_V3_ROI 185"
  "L_V4_ROI 186"
  "L_V8_ROI 187"
  "L_FEF_ROI 190"
  "L_PEF_ROI 191"
  "L_55b_ROI 192"
  "L_V3A_ROI 193"
  "L_V7_ROI 196"
  "L_IPS1_ROI 197"
  "L_FFC_ROI 198"
  "L_V3B_ROI 199"
  "L_LO1_ROI 200"
  "L_LO2_ROI 201"
  "L_PIT_ROI 202"
  "L_MT_ROI 203"
  "L_7Pm_ROI 209"
  "L_24dv_ROI 221"
  "L_7AL_ROI 222"
  "L_SCEF_ROI 223"
  "L_6ma_ROI 224"
  "L_7Am_ROI 225"
  "L_7PL_ROI 226"
  "L_7PC_ROI 227"
  "L_LIPv_ROI 228"
  "L_VIP_ROI 229"
  "L_MIP_ROI 230"
  "L_6d_ROI 234"
  "L_6mp_ROI 235"
  "L_6v_ROI 236"
  "L_p32pr_ROI 240"
  "L_6r_ROI 258"
  "L_IFJa_ROI 259"
  "L_IFJp_ROI 260"
  "L_LIPd_ROI 275"
  "L_6a_ROI 276"
  "L_i6-8_ROI 277"
  "L_AIP_ROI 297"
  "L_PH_ROI 318"
  "L_IP2_ROI 324"
  "L_IP1_ROI 325"
  "L_IP0_ROI 326"
  "L_V6A_ROI 332"
  "L_VMV1_ROI 333"
  "L_VMV3_ROI 334"
  "L_V4t_ROI 336"
  "L_FST_ROI 337"
  "L_V3CD_ROI 338"
  "L_LO3_ROI 339"
  "L_VMV2_ROI 340"
  "L_VVC_ROI 343"
)

# Loop over all ROIs and extract
for roi in "${ROIS[@]}"; do
  set -- $roi
  NAME=$1
  LABEL_ID=$2

  # Choose hemisphere from ROI name
  if [[ "$NAME" == L* ]]; then
    CORTEX="CORTEX_LEFT"
  elif [[ "$NAME" == R* ]]; then
    CORTEX="CORTEX_RIGHT"
  else
    echo "⚠️ Could not determine hemisphere for ROI $NAME"
    continue
  fi

  # Extract label to dscalar
  wb_command -cifti-label-to-roi "$ATLAS_FILE" \
    "$OUTDIR/${NAME}.dscalar.nii" -key "$LABEL_ID"

  # Separate to metric GIFTI
  wb_command -cifti-separate "$OUTDIR/${NAME}.dscalar.nii" COLUMN \
    -metric "$CORTEX" "$OUTDIR/${NAME}.shape.gii"
done

# Permissions
chmod -Rf 771 "$OUTDIR"
