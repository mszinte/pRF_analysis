#!/bin/bash
#####################################################
# Winner-Take-All analysis using Connectome-Workbench
# Note: INDEXMAX returns winning COLUMN (seed cluster, e.g. sPCS) by ROW (MMP  label) 
# Written by Marco Bedini (marco.bedini@univ-amu.fr)
#####################################################

BASE_PATH="/scratch/mszinte/data/RetinoMaps/derivatives/pp_data"
ATLAS_DIR="/home/${USER}/projects/pRF_analysis/RetinoMaps/rest/mmp1_clusters/leaveout"
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

echo "=========================================="
echo "STEP 1: Processing individual subjects"
echo "=========================================="

# Iterate through subjects to compute individual WTA maps
for sub in 01 02 03 04 05 06 07 08 09 11 12 13 14 17 20 21 22 23 24 25; do
    
    echo "Processing sub-${sub}..."
    
    FULL_CORR="${BASE_PATH}/sub-${sub}/91k/rest/corr/full_corr"
    
    # Build the merge command - merge all ROI correlation maps
    MERGE_CMD="wb_command -cifti-merge ${FULL_CORR}/sub-${sub}_task-rest_space-fsLR_den-91k_desc-full_corr_merged.pscalar.nii -direction ROW"
    
    # Add each ROI in order (this creates rows 0-11 for mPCS-V1)
    for roi in "${ROIS[@]}"; do
        MERGE_CMD="${MERGE_CMD} -cifti ${FULL_CORR}/sub-${sub}_task-rest_space-fsLR_den-91k_desc-full_corr_${roi}_parcellated.pscalar.nii"
        
        # Hollow seed cluster to 0s so it doesn't pop up in the winner-take-all map
    	roi_file_lh="${ATLAS_PATH}/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_leaveout_lh_${roi}.shape.gii"
    	roi_file_rh="${ATLAS_PATH}/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_leaveout_rh_${roi}.shape.gii"
    
		wb_command -cifti-restrict-dense-map ${FULL_CORR}/sub-${sub}_task-rest_space-fsLR_den-91k_desc-full_corr_${roi}_parcellated.pscalar.nii COLUMN \
		${FULL_CORR}/sub-${sub}_task-rest_space-fsLR_den-91k_desc-full_corr_merged_hollow_seed.pscalar.nii \
		-left-roi "$roi_file_lh" \
    	-right-roi "$roi_file_rh"
        
    done
    
    # Execute merge
    eval $MERGE_CMD
    
    # Transpose so clusters are columns and parcels are rows
    wb_command -cifti-transpose \
        ${FULL_CORR}/sub-${sub}_task-rest_space-fsLR_den-91k_desc-full_corr_merged_hollow_seed.pscalar.nii \
        ${FULL_CORR}/sub-${sub}_task-rest_space-fsLR_den-91k_desc-full_corr_merged.transpose.nii \
        -mem-limit 8
    
    # Get the index of maximum correlation for each parcel
    # Output: 106 parcels × 1 column with winner index (1-12)
    wb_command -cifti-reduce \
        ${FULL_CORR}/sub-${sub}_task-rest_space-fsLR_den-91k_desc-full_corr_merged.transpose.nii \
        INDEXMAX \
        ${OUTPUT_PATH}/sub-${sub}_indexmax_wta.raw.nii \
        -direction COLUMN
    
    # Transpose back to standard orientation
    wb_command -cifti-transpose \
        ${OUTPUT_PATH}/sub-${sub}_indexmax_wta.raw.nii \
        ${OUTPUT_PATH}/sub-${sub}_indexmax_wta.pscalar.nii \
        -mem-limit 8
    
    echo "  Completed sub-${sub}"
done

echo ""
echo "=========================================="
echo "STEP 2: Computing group-level WTA (MODE)"
echo "=========================================="

# Build merge command for all subjects
# This creates a file with 106 parcels (rows) × 20 subjects (columns)
MERGE_GROUP_CMD="wb_command -cifti-merge ${OUTPUT_PATH}/all_subjects_wta_merged.pscalar.nii -direction ROW"

for sub in 01 02 03 04 05 06 07 08 09 11 12 13 14 17 20 21 22 23 24 25; do
    MERGE_GROUP_CMD="${MERGE_GROUP_CMD} -cifti ${OUTPUT_PATH}/sub-${sub}_indexmax_wta.pscalar.nii"
done

# Execute group merge
eval $MERGE_GROUP_CMD

echo "Merged all subjects into single file"

# Transpose so subjects are rows and parcels are columns
wb_command -cifti-transpose \
    ${OUTPUT_PATH}/all_subjects_wta_merged.pscalar.nii \
    ${OUTPUT_PATH}/all_subjects_wta_merged.transpose.nii \
    -mem-limit 8

echo "Transposed for MODE calculation"

# Calculate MODE across subjects (rows) for each parcel (column)
# This gives the most frequent winner across subjects for each parcel
wb_command -cifti-reduce \
    ${OUTPUT_PATH}/all_subjects_wta_merged.transpose.nii \
    MODE \
    ${OUTPUT_PATH}/group_wta_mode_raw.pscalar.nii \
    -direction COLUMN

echo "Calculated MODE across subjects"

# Transpose back to standard orientation (parcels as rows)
wb_command -cifti-transpose \
    ${OUTPUT_PATH}/group_wta_mode.raw.nii \
    ${OUTPUT_PATH}/group_wta_mode.pscalar.nii \
    -mem-limit 8

echo "Group-level WTA map created!"

echo ""
echo "=========================================="
echo "ANALYSIS COMPLETE"
echo "=========================================="
echo ""
echo "Output files:"
echo "  Individual subjects: ${OUTPUT_PATH}/sub-*_indexmax_wta.pscalar.nii"
echo "  Group MODE:          ${OUTPUT_PATH}/group_wta_mode.pscalar.nii"
echo ""
echo "Winner seed mapping:"
for i in "${!ROIS[@]}"; do
    echo "  $((i+1)) = ${ROIS[$i]}"
done
echo ""
echo "Next steps:"
echo "1. Convert group_wta_mode.pscalar.nii to label file for visualization"
echo "2. Compare with Python-derived results"
echo "3. Calculate consistency metrics"
