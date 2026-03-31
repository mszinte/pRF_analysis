#!/bin/bash
#####################################################
# Winner-Take-All analysis using Connectome-Workbench
# Flexible for different correlation types and hemisphere modes
# Note: INDEXMAX returns winning COLUMN (seed cluster, e.g. sPCS) by ROW (MMP  label)
#
# Written by Marco Bedini (marco.bedini@univ-amu.fr) with Claude's help
# 
# Usage examples:
#   ./workbench_compute_wta_full_corr_dev.sh full_corr bilateral default
#   ./workbench_compute_wta_full_corr_dev.sh fisher-z by_hemi legacy
#####################################################

# Parse command-line arguments
CORR_TYPE="${1:-full_corr}"      # 'full_corr' or 'fisher-z'
HEMI_TYPE="${2:-bilateral}"      # 'bilateral' or 'by_hemi'
LEGACY_MODE="${3:-default}"      # 'default' or 'legacy'

# Validate inputs
if [[ "$CORR_TYPE" != "full_corr" ]] && [[ "$CORR_TYPE" != "fisher-z" ]]; then
    echo "Error: CORR_TYPE must be 'full_corr' or 'fisher_z'"
    echo "Usage: $0 [full_corr|fisher-z] [bilateral|by_hemi] [default|legacy]"
    exit 1
fi

if [[ "$HEMI_TYPE" != "bilateral" ]] && [[ "$HEMI_TYPE" != "by_hemi" ]]; then
    echo "Error: HEMI_TYPE must be 'bilateral' or 'by_hemi'"
    echo "Usage: $0 [full_corr|fisher-z] [bilateral|by_hemi] [default|legacy]"
    exit 1
fi

# Set legacy flag for workbench
if [[ "$LEGACY_MODE" == "legacy" ]]; then
    LEGACY_FLAG="-legacy-mode"
else
    LEGACY_FLAG=""
fi

echo "=========================================="
echo "WTA ANALYSIS CONFIGURATION"
echo "=========================================="
echo "Correlation type: ${CORR_TYPE}"
echo "Hemisphere mode: ${HEMI_TYPE}"
echo "Legacy mode: ${LEGACY_MODE}"
echo "=========================================="
echo ""

# Main folders
BASE_PATH="/scratch/mszinte/data/RetinoMaps/derivatives/pp_data"
ATLAS_DIR="${BASE_PATH}/atlas"
OUTPUT_PATH="${BASE_PATH}/group/91k/rest/wta/workbench"
mkdir -p "${OUTPUT_PATH}"

# Define ROI order
declare -a ROIS=("mPCS" "sPCS" "iPCS" "sIPS" "iIPS" "hMT+" "VO" "LO" "V3AB" "V3" "V2" "V1")

# Construct file suffix based on options
FILE_SUFFIX="${CORR_TYPE}"
[[ "$HEMI_TYPE" == "by_hemi" ]] && FILE_SUFFIX="${FILE_SUFFIX}_by_hemi"

# Determine which hemispheres to process
if [[ "$HEMI_TYPE" == "by_hemi" ]]; then
    HEMIS=("lh" "rh")
else
    HEMIS=("bilateral")
fi

# ============================================
# HEMISPHERE LOOP
# ============================================

for HEMI in "${HEMIS[@]}"; do
    
    if [[ "$HEMI" != "by_hemi" ]]; then
        echo ""
        echo "=========================================="
        echo "PROCESSING HEMISPHERE: ${HEMI^^}"
        echo "=========================================="
        HEMI_SUFFIX="_${HEMI}"
    else
        HEMI_SUFFIX=""
    fi
    
    echo ""
    echo "=========================================="
    echo "STEP 1: Processing individual subjects"
    echo "=========================================="
    
    # Iterate through subjects
    for sub in 01 02 03 04 05 06 07 08 09 11 12 13 14 17 20 21 22 23 24 25; do
        
        echo "Processing sub-${sub}${HEMI_SUFFIX}..."
        
        # Create output directory
        FULL_CORR="${BASE_PATH}/sub-${sub}/91k/rest/corr/full_corr/wta"
        mkdir -p "${FULL_CORR}"
        
        # Process each ROI
        for roi in "${ROIS[@]}"; do
            
            # Process combined hemispheres
            INPUT_FILE="${BASE_PATH}/sub-${sub}/91k/rest/corr/full_corr/sub-${sub}_task-rest_space-fsLR_den-91k_desc-${CORR_TYPE}_${roi}.dscalar.nii"
                
             # Separate hemispheres for masking
             wb_command -cifti-separate "${INPUT_FILE}" COLUMN \
                -metric CORTEX_LEFT "${FULL_CORR}/sub-${sub}_${CORR_TYPE}_lh_${roi}.shape.gii" \
                -metric CORTEX_RIGHT "${FULL_CORR}/sub-${sub}_${CORR_TYPE}_rh_${roi}.shape.gii"
                
            # Mask left hemisphere
            wb_command -metric-math \
                '(corr * (hollow_seed > 0)) + ((hollow_seed == 0) * -1)' \
                "${FULL_CORR}/sub-${sub}_${CORR_TYPE}_lh_${roi}_hollow.shape.gii" \
                -var corr "${FULL_CORR}/sub-${sub}_${CORR_TYPE}_lh_${roi}.shape.gii" \
                -var hollow_seed "${ATLAS_DIR}/mmp1_clusters/leaveout/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_leaveout_lh_${roi}_bin.shape.gii"
                
            # Mask right hemisphere
            wb_command -metric-math \
                '(corr * (hollow_seed > 0)) + ((hollow_seed == 0) * -1)' \
                "${FULL_CORR}/sub-${sub}_${CORR_TYPE}_rh_${roi}_hollow.shape.gii" \
                -var corr "${FULL_CORR}/sub-${sub}_${CORR_TYPE}_rh_${roi}.shape.gii" \
                -var hollow_seed "${ATLAS_DIR}/mmp1_clusters/leaveout/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_leaveout_rh_${roi}_bin.shape.gii"
                
            # Recombine hemispheres
            wb_command -cifti-create-dense-scalar \
                "${FULL_CORR}/sub-${sub}_${CORR_TYPE}_${roi}_hollow.dscalar.nii" \
                -left-metric "${FULL_CORR}/sub-${sub}_${CORR_TYPE}_lh_${roi}_hollow.shape.gii" \
                -right-metric "${FULL_CORR}/sub-${sub}_${CORR_TYPE}_rh_${roi}_hollow.shape.gii"
                
                
            # Mask the hemisphere
            wb_command -metric-math \
                '(corr * (hollow_seed > 0)) + ((hollow_seed == 0) * -1)' \
                "${FULL_CORR}/sub-${sub}_${CORR_TYPE}_${HEMI}_${roi}_hollow.shape.gii" \
                -var corr "${FULL_CORR}/sub-${sub}_${CORR_TYPE}_${HEMI}_${roi}.shape.gii" \
                -var hollow_seed "${ATLAS_DIR}/mmp1_clusters/leaveout/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_leaveout_${HEMI}_${roi}_bin.shape.gii"
            
            # Parcellate
            OUTPUT_PARCELLATED="${FULL_CORR}/sub-${sub}_${CORR_TYPE}${HEMI_SUFFIX}_${roi}.pscalar.nii"
            
            INPUT_TO_PARCELLATE="${FULL_CORR}/sub-${sub}_${CORR_TYPE}_${roi}_hollow.dscalar.nii"
            
            wb_command -cifti-parcellate \
                "${INPUT_TO_PARCELLATE}" \
                "${ATLAS_DIR}/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_dseg.dlabel.nii" COLUMN \
                "${OUTPUT_PARCELLATED}" -method MEAN \
                -only-numeric ${LEGACY_FLAG}
        done
        
        # Merge all ROIs
        MERGE_OUTPUT="${FULL_CORR}/sub-${sub}_${FILE_SUFFIX}${HEMI_SUFFIX}_merged.pscalar.nii"
        MERGE_CMD="wb_command -cifti-merge ${MERGE_OUTPUT} -direction ROW"
        
        for roi in "${ROIS[@]}"; do
            MERGE_CMD="${MERGE_CMD} -cifti ${FULL_CORR}/sub-${sub}_${CORR_TYPE}${HEMI_SUFFIX}_${roi}.pscalar.nii"
        done
        
        eval $MERGE_CMD
        
        # Transpose
        wb_command -cifti-transpose \
            "${MERGE_OUTPUT}" \
            "${FULL_CORR}/sub-${sub}_${FILE_SUFFIX}${HEMI_SUFFIX}_merged.transpose.nii" \
            -mem-limit 20
        
        # INDEXMAX
        wb_command -cifti-reduce \
            "${FULL_CORR}/sub-${sub}_${FILE_SUFFIX}${HEMI_SUFFIX}_merged.transpose.nii" \
            INDEXMAX \
            "${OUTPUT_PATH}/sub-${sub}_indexmax_wta_${FILE_SUFFIX}${HEMI_SUFFIX}_raw.transpose.nii" \
            -direction COLUMN
        
        # Transpose back
        wb_command -cifti-transpose \
            "${OUTPUT_PATH}/sub-${sub}_indexmax_wta_${FILE_SUFFIX}${HEMI_SUFFIX}_raw.transpose.nii" \
            "${OUTPUT_PATH}/sub-${sub}_indexmax_wta_${FILE_SUFFIX}${HEMI_SUFFIX}.pscalar.nii" \
            -mem-limit 20
        
        echo "  Completed sub-${sub}${HEMI_SUFFIX}"
    done
    
    # ============================================
    # GROUP-LEVEL MODE
    # ============================================
    
    echo ""
    echo "=========================================="
    echo "STEP 2: Computing group-level WTA (MODE)"
    echo "=========================================="
    
    # Merge all subjects
    MERGE_GROUP_CMD="wb_command -cifti-merge ${OUTPUT_PATH}/all_subjects_wta_${FILE_SUFFIX}${HEMI_SUFFIX}_merged.pscalar.nii -direction ROW"
    
    for sub in 01 02 03 04 05 06 07 08 09 11 12 13 14 17 20 21 22 23 24 25; do
        MERGE_GROUP_CMD="${MERGE_GROUP_CMD} -cifti ${OUTPUT_PATH}/sub-${sub}_indexmax_wta_${FILE_SUFFIX}${HEMI_SUFFIX}.pscalar.nii"
    done
    
    eval $MERGE_GROUP_CMD
    echo "Merged all subjects"
    
    # Transpose
    wb_command -cifti-transpose \
        "${OUTPUT_PATH}/all_subjects_wta_${FILE_SUFFIX}${HEMI_SUFFIX}_merged.pscalar.nii" \
        "${OUTPUT_PATH}/all_subjects_wta_${FILE_SUFFIX}${HEMI_SUFFIX}_merged_transpose.nii" \
        -mem-limit 20
    echo "Transposed for MODE calculation"
    
    # Calculate MODE
    wb_command -cifti-reduce \
        "${OUTPUT_PATH}/all_subjects_wta_${FILE_SUFFIX}${HEMI_SUFFIX}_merged_transpose.nii" \
        MODE \
        "${OUTPUT_PATH}/group_wta_${FILE_SUFFIX}${HEMI_SUFFIX}_mode_raw.nii" \
        -direction COLUMN
    echo "Calculated MODE"
    
    # Transpose back
    wb_command -cifti-transpose \
        "${OUTPUT_PATH}/group_wta_${FILE_SUFFIX}${HEMI_SUFFIX}_mode_raw.nii" \
        "${OUTPUT_PATH}/group_wta_${FILE_SUFFIX}${HEMI_SUFFIX}.pscalar.nii" \
        -mem-limit 20
    
    # Convert to TSV
    wb_command -cifti-convert -to-text \
        "${OUTPUT_PATH}/group_wta_${FILE_SUFFIX}${HEMI_SUFFIX}.pscalar.nii" \
        "${OUTPUT_PATH}/group_wta_${FILE_SUFFIX}${HEMI_SUFFIX}.tsv"
    
    echo "Group-level WTA map created for ${HEMI}${HEMI_SUFFIX}!"
    
done  # End hemisphere loop

# ============================================
# FINAL SUMMARY
# ============================================

echo ""
echo "=========================================="
echo "ANALYSIS COMPLETE"
echo "=========================================="
echo ""
echo "Configuration:"
echo "  Correlation: ${CORR_TYPE}"
echo "  Hemisphere: ${HEMI_TYPE}"
echo "  Legacy mode: ${LEGACY_MODE}"
echo ""
echo "Output files:"
if [[ "$HEMI_TYPE" == "bilateral" ]]; then
    echo "  Subjects: ${OUTPUT_PATH}/sub-*_indexmax_wta_${FILE_SUFFIX}.pscalar.nii"
    echo "  Group: ${OUTPUT_PATH}/group_wta_${FILE_SUFFIX}.pscalar.nii"
    echo "  TSV: ${OUTPUT_PATH}/group_wta_${FILE_SUFFIX}.tsv"
else
    echo "  LH Subjects: ${OUTPUT_PATH}/sub-*_indexmax_wta_${FILE_SUFFIX}_lh.pscalar.nii"
    echo "  RH Subjects: ${OUTPUT_PATH}/sub-*_indexmax_wta_${FILE_SUFFIX}_rh.pscalar.nii"
    echo "  LH Group: ${OUTPUT_PATH}/group_wta_${FILE_SUFFIX}_lh.pscalar.nii"
    echo "  RH Group: ${OUTPUT_PATH}/group_wta_${FILE_SUFFIX}_rh.pscalar.nii"
    echo "  LH TSV: ${OUTPUT_PATH}/group_wta_${FILE_SUFFIX}_lh.tsv"
    echo "  RH TSV: ${OUTPUT_PATH}/group_wta_${FILE_SUFFIX}_rh.tsv"
fi
echo ""
echo "Winner seed mapping:"
for i in "${!ROIS[@]}"; do
    echo "  $((i+1)) = ${ROIS[$i]}"
done
echo ""

# Change permissions
chmod -Rf 771 ${BASE_PATH}
chgrp -Rf 327 ${BASE_PATH}