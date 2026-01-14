#!/bin/bash
#####################################################
# Winner-Take-All analysis using Connectome-Workbench
# Note: INDEXMAX returns winning COLUMN (seed cluster, e.g. sPCS) by ROW (MMP  label) 
# Written by Marco Bedini (marco.bedini@univ-amu.fr)
#####################################################

BASE_PATH="/scratch/mszinte/data/RetinoMaps/derivatives/pp_data"
ATLAS_DIR="/home/${USER}/projects/pRF_analysis/RetinoMaps/rest/mmp1_clusters"
OUTPUT_PATH="${BASE_PATH}/group/91k/rest/wta"
mkdir -p "${OUTPUT_PATH}"

# Define ROI order for merging (determines map index)
# Index starts at 1 (we add 1 to 0-based indices)
# Order: mPCS=1, sPCS=2, iPCS=3, sIPS=4, iIPS=5, hMT+=6, VO=7, LO=8, V3AB=9, V3=10, V2=11, V1=12
declare -a ROIS=("mPCS" "sPCS" "iPCS" "sIPS" "iIPS" "hMT+" "VO" "LO" "V3AB" "V3" "V2" "V1")

echo "=========================================="
echo "STEP 1: Processing individual subjects"
echo "=========================================="

# Iterate through subjects to compute individual WTA maps
for sub in 01 02 03 04 05 06 07 08 09 11 12 13 14 17 20 21 22 23 24 25; do
    
    echo "Processing sub-${sub}..."
    
    FULL_CORR="${BASE_PATH}/sub-${sub}/91k/rest/corr/full_corr"
    
    ################# Developing this part ###############################################
    
    # Remove volume info so we end up with approximately 59k vertices (use the fsLR file not the atlas it's cleaner conceptually)
    for roi in "${ROIS[@]}"; do
    	wb_command -cifti-restrict-dense-map \
    	"${FULL_CORR}/sub-${sub}_ses-01_task-rest_space-fsLR_den-91k_desc-full_corr_${roi}.dscalar.nii" \
    	COLUMN \
    	-left-roi "$ATLAS_DIR/tpl-fsLR_hemi-L_den-59k_desc-vaavg_midthickness.shape.gii" \
    	-right-roi "$ATLAS_DIR/tpl-fsLR_hemi-R_den-59k_desc-vaavg_midthickness.shape.gii" \
    	"${FULL_CORR}/sub-${sub}_task-rest_space-fsLR_den-91k_desc-full_corr_${roi}_cortex.dscalar.nii" 
    	
    # Zero out values outside the atlas
    	wb_command -cifti-math \
    	'corr * (atlas > 0)' \
        "${FULL_CORR}/sub-${sub}_task-rest_space-fsLR_den-91k_desc-full_corr_${roi}_cortex.dscalar.nii" \
        -var corr "${FULL_CORR}/sub-${sub}_task-rest_space-fsLR_den-91k_desc-full_corr_${roi}_cortex.dscalar.nii" \
        -var atlas "${ATLAS_DIR}/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_dseg_bin.dscalar.nii"
    
    	# Copy file for debugging
    	cp "${FULL_CORR}/sub-${sub}_task-rest_space-fsLR_den-91k_desc-full_corr_${roi}_cortex.dscalar.nii" \
    	"${FULL_CORR}/sub-${sub}_task-rest_space-fsLR_den-91k_desc-full_corr_${roi}_cortex_b.dscalar.nii"
    
    # Zero out values for the seed cluster (HERE USE THE SHAPE BIN FILE) 	
    	wb_command -cifti-math \
    	'(corr + (-0.999 * (cluster_lh > 0)) + (-0.999 * (cluster_rh > 0)))' \
        "${FULL_CORR}/sub-${sub}_task-rest_space-fsLR_den-91k_desc-full_corr_${roi}_cortex_b.dscalar.nii" \
        -var corr "${FULL_CORR}/sub-${sub}_task-rest_space-fsLR_den-91k_desc-full_corr_${roi}_cortex_b.dscalar.nii" \
        -var cluster_lh "${ATLAS_DIR}/clusters/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_discarded_lh_${roi}_bin.shape.gii" \
        -var cluster_rh "${ATLAS_DIR}/clusters/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_discarded_rh_${roi}_bin.shape.gii"
        
    # Parcellate the outputs now
	wb_command -cifti-parcellate \
	"$FULL_CORR/sub-${sub}_task-rest_space-fsLR_den-91k_desc-full_corr_${roi}_cortex_b.dscalar.nii" \
        "$ATLAS_DIR/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_dseg.dlabel.nii" COLUMN \
        "$FULL_CORR/sub-${sub}_task-rest_space-fsLR_den-91k_desc-full_corr_${roi}_hollow_seed.pscalar.nii" -method MEAN; 
        ## to exclude outliers, add this flag at the end -exclude-outliers 3 3	
        
	done
    #######################################################################################
    
    # Build the merge command - merge all ROI correlation maps
    MERGE_CMD="wb_command -cifti-merge ${FULL_CORR}/sub-${sub}_task-rest_space-fsLR_den-91k_desc-full_corr_merged.pscalar.nii -direction ROW"
    
    # Add each ROI in order (this creates rows 0-11 for mPCS-V1)
    for roi in "${ROIS[@]}"; do
        MERGE_CMD="${MERGE_CMD} -cifti ${FULL_CORR}/sub-${sub}_task-rest_space-fsLR_den-91k_desc-full_corr_${roi}_hollow_seed.pscalar.nii"
    done
    
    # Execute merge
    eval $MERGE_CMD
    
    # Transpose so clusters are columns and parcels are rows
    wb_command -cifti-transpose \
        ${FULL_CORR}/sub-${sub}_task-rest_space-fsLR_den-91k_desc-full_corr_merged.pscalar.nii \
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
