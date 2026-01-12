#!/bin/bash
#####################################################
# Winner-Take-All analysis using Connectome-Workbench
# Note: INDEXMAX returns winning ROW (map) that corresponds to the seed cluster (e.g. sPCS)
# Written by Marco Bedini (marco.bedini@univ-amu.fr)
#####################################################

BASE_PATH="/scratch/mszinte/data/RetinoMaps/derivatives/pp_data"
OUTPUT_PATH="${BASE_PATH}/group/91k/rest/wta"
mkdir -p "${OUTPUT_PATH}"

# Define ROI order for merging (determines map index)
# Index starts at 1 (we add 1 to 0-based indices)
# Order: mPCS=1, sPCS=2, iPCS=3, sIPS=4, iIPS=5, hMT+=6, VO=7, LO=8, V3AB=9, V3=10, V2=11, V1=12
declare -a ROIS=("mPCS" "sPCS" "iPCS" "sIPS" "iIPS" "hMT+" "VO" "LO" "V3AB" "V3" "V2" "V1")

# Define all parcel names in order (these are the columns)
declare -a PARCELS=(
    "R_V1_ROI" "R_MST_ROI" "R_V2_ROI" "R_V3_ROI" "R_V4_ROI" "R_V8_ROI" 
    "R_FEF_ROI" "R_PEF_ROI" "R_55b_ROI" "R_V3A_ROI" "R_V7_ROI" "R_IPS1_ROI"
    "R_FFC_ROI" "R_V3B_ROI" "R_LO1_ROI" "R_LO2_ROI" "R_PIT_ROI" "R_MT_ROI"
    "R_7Pm_ROI" "R_24dv_ROI" "R_7AL_ROI" "R_SCEF_ROI" "R_6ma_ROI" "R_7Am_ROI"
    "R_7PL_ROI" "R_7PC_ROI" "R_LIPv_ROI" "R_VIP_ROI" "R_MIP_ROI" "R_6d_ROI"
    "R_6mp_ROI" "R_6v_ROI" "R_p32pr_ROI" "R_6r_ROI" "R_IFJa_ROI" "R_IFJp_ROI"
    "R_LIPd_ROI" "R_6a_ROI" "R_i6-8_ROI" "R_AIP_ROI" "R_PH_ROI" "R_IP2_ROI"
    "R_IP1_ROI" "R_IP0_ROI" "R_V6A_ROI" "R_VMV1_ROI" "R_VMV3_ROI" "R_V4t_ROI"
    "R_FST_ROI" "R_V3CD_ROI" "R_LO3_ROI" "R_VMV2_ROI" "R_VVC_ROI"
    "L_V1_ROI" "L_MST_ROI" "L_V2_ROI" "L_V3_ROI" "L_V4_ROI" "L_V8_ROI"
    "L_FEF_ROI" "L_PEF_ROI" "L_55b_ROI" "L_V3A_ROI" "L_V7_ROI" "L_IPS1_ROI"
    "L_FFC_ROI" "L_V3B_ROI" "L_LO1_ROI" "L_LO2_ROI" "L_PIT_ROI" "L_MT_ROI"
    "L_7Pm_ROI" "L_24dv_ROI" "L_7AL_ROI" "L_SCEF_ROI" "L_6ma_ROI" "L_7Am_ROI"
    "L_7PL_ROI" "L_7PC_ROI" "L_LIPv_ROI" "L_VIP_ROI" "L_MIP_ROI" "L_6d_ROI"
    "L_6mp_ROI" "L_6v_ROI" "L_p32pr_ROI" "L_6r_ROI" "L_IFJa_ROI" "L_IFJp_ROI"
    "L_LIPd_ROI" "L_6a_ROI" "L_i6-8_ROI" "L_AIP_ROI" "L_PH_ROI" "L_IP2_ROI"
    "L_IP1_ROI" "L_IP0_ROI" "L_V6A_ROI" "L_VMV1_ROI" "L_VMV3_ROI" "L_V4t_ROI"
    "L_FST_ROI" "L_V3CD_ROI" "L_LO3_ROI" "L_VMV2_ROI" "L_VVC_ROI"
)

# Iterate through subjects
for sub in 01 02 03 04 05 06 07 08 09 11 12 13 14 17 20 21 22 23 24 25; do
	    
echo "Processing sub-${sub}..."
	    
FULL_CORR="${BASE_PATH}/sub-${sub}/91k/rest/corr/full_corr"
	    
# Build the merge command - merge all ROI correlation maps
MERGE_CMD="wb_command -cifti-merge ${FULL_CORR}/sub-${sub}_task-rest_space-fsLR_den-91k_desc-full_corr_merged.pscalar.nii -direction ROW" \
	    
	    # Add each ROI in order (this creates rows 0-11 for mPCS-V1)
	    for roi in "${ROIS[@]}"; do
		MERGE_CMD="${MERGE_CMD} -cifti ${FULL_CORR}/sub-${sub}_task-rest_space-fsLR_den-91k_desc-full_corr_${roi}_parcellated.pscalar.nii"
	    done
	    
	    # Execute merge
	    eval $MERGE_CMD
	    
	    # Get the index of maximum correlation for each parcel (column)
	    # This returns which ROW (seed) has the maximum correlation with each parcel
	    # Output format: one line per parcel, value is the row index (0-11)
	    wb_command -cifti-reduce \
		${FULL_CORR}/sub-${sub}_task-rest_space-fsLR_den-91k_desc-full_corr_merged.pscalar.nii \
		INDEXMAX \
		${OUTPUT_PATH}/sub-${sub}_indexmax_raw.wta.nii \
		-direction COLUMN;

done
echo "" >> "${OUTPUT_PATH}/winning_seeds_full_corr.csv"
        
done
