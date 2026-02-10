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
    
    # Create separate output dir
    mkdir -p "${BASE_PATH}/sub-${sub}/91k/rest/corr/full_corr/wta"
    FULL_CORR="${BASE_PATH}/sub-${sub}/91k/rest/corr/full_corr/wta"
    
    ################# Debugging this part ###############################################
    ## add print statements with wb_command -file-information to inspect the outputs
    for roi in "${ROIS[@]}"; do
        
        # Get metric files for both hemi (was used to remove medial wall but I'm not sure is neeeded)
    	# wb_command -cifti-separate "${BASE_PATH}/sub-${sub}/91k/rest/corr/full_corr/sub-${sub}_task-rest_space-fsLR_den-91k_desc-full_corr_${roi}.dscalar.nii" COLUMN \
    	# -metric CORTEX_LEFT "$FULL_CORR/sub-${sub}_task-rest_space-fsLR_den-91k_desc-full_corr_lh_${roi}.shape.gii" \
    	# -metric CORTEX_RIGHT "$FULL_CORR/sub-${sub}_task-rest_space-fsLR_den-91k_desc-full_corr_rh_${roi}.shape.gii"
    
    	# Replace seed cluster so it doesn't show up in the wta (hollow seed)
    	# wb_command -metric-mask "$FULL_CORR/sub-${sub}_task-rest_space-fsLR_den-91k_desc-full_corr_lh_${roi}.shape.gii" \
    	# "$ATLAS_DIR/leaveout/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_leaveout_lh_${roi}.shape.gii" \
    	# "$FULL_CORR/sub-${sub}_task-rest_space-fsLR_den-91k_desc-full_corr_lh_${roi}_hollow_seed.shape.gii"

        # wb_command -metric-mask "$FULL_CORR/sub-${sub}_task-rest_space-fsLR_den-91k_desc-full_corr_rh_${roi}.shape.gii" \
    	# "$ATLAS_DIR/leaveout/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_leaveout_rh_${roi}.shape.gii" \
    	# "$FULL_CORR/sub-${sub}_task-rest_space-fsLR_den-91k_desc-full_corr_rh_${roi}_hollow_seed.shape.gii"

        # It's safer to use cifti-mask
        wb_command -cifti-math \
    	'corr * (hollow_seed = 0)' \
        "${FULL_CORR}/sub-${sub}_task-rest_space-fsLR_den-91k_desc-full_corr_${roi}.dscalar.nii" \
        -var corr "${FULL_CORR}/sub-${sub}_task-rest_space-fsLR_den-91k_desc-full_corr_${roi}.dscalar.nii" \
        -var hollow_seed "${ATLAS_DIR}/leaveout/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_leaveout_${roi}_bin.dscalar.nii"

    	# Go back to the combined hemi version in cifti format
    	# wb_command -cifti-create-dense-scalar "${FULL_CORR}/sub-${sub}_task-rest_space-fsLR_den-91k_desc-full_corr_${roi}_hollow_seed.dscalar.nii" \
    	# -left-metric "$FULL_CORR/sub-${sub}_task-rest_space-fsLR_den-91k_desc-full_corr_lh_${roi}_hollow_seed.shape.gii" \
        # -right-metric "$FULL_CORR/sub-${sub}_task-rest_space-fsLR_den-91k_desc-full_corr_rh_${roi}_hollow_seed.shape.gii"
    	
    	# Parcellate the full corr outputs
        wb_command -cifti-parcellate \
        "$FULL_CORR/sub-${sub}_task-rest_space-fsLR_den-91k_desc-full_corr_${roi}_hollow_seed.dscalar.nii" \
        "$ATLAS_DIR/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_dseg.dlabel.nii" COLUMN \
        "$FULL_CORR/sub-${sub}_task-rest_space-fsLR_den-91k_desc-full_corr_${roi}.pscalar.nii" -method MEAN; 
        ## to exclude outliers, add this flag at the end -exclude-outliers 3 3
        
        # Double check this
        # NOTE: the parcels in the output file are sorted by the numeric label keys, in ascending order.
        
	done
    #######################################################################################
    
    # Build the merge command - merge all ROI correlation maps
    # Probably here the output has to be a .pconn not .pscalar file
    MERGE_CMD="wb_command -cifti-merge ${FULL_CORR}/sub-${sub}_task-rest_space-fsLR_den-91k_desc-full_corr_merged.pscalar.nii -direction ROW"
    
    # Add each ROI in order (this creates rows 0-11 for mPCS-V1)
    for roi in "${ROIS[@]}"; do
        MERGE_CMD="${MERGE_CMD} -cifti ${FULL_CORR}/sub-${sub}_task-rest_space-fsLR_den-91k_desc-full_corr_${roi}.pscalar.nii"
    done
    
    # Execute merge
    eval $MERGE_CMD
    
    # Transpose so parcels are rows and clusters are columns 
    wb_command -cifti-transpose \
        "${FULL_CORR}/sub-${sub}_task-rest_space-fsLR_den-91k_desc-full_corr_merged.pscalar.nii" \
        "${FULL_CORR}/sub-${sub}_task-rest_space-fsLR_den-91k_desc-full_corr_merged.transpose.nii" \
        -mem-limit 20
    
    # Get the index of maximum correlation for each parcel
    # Output: 106 parcels × 1 column with winner index (1-12)
    wb_command -cifti-reduce \
        "${FULL_CORR}/sub-${sub}_task-rest_space-fsLR_den-91k_desc-full_corr_merged.transpose.nii" \
        INDEXMAX \
        "${OUTPUT_PATH}/sub-${sub}_indexmax_wta_full_corr.raw.nii" \
        -direction COLUMN
    
    # Transpose back to standard orientation
    wb_command -cifti-transpose \
        "${OUTPUT_PATH}/sub-${sub}_indexmax_wta_full_corr.raw.nii" \
        "${OUTPUT_PATH}/sub-${sub}_indexmax_wta_full_corr.pscalar.nii" \
        -mem-limit 20
    
    echo "  Completed sub-${sub}"
done

echo ""
echo "=========================================="
echo "STEP 2: Computing group-level WTA (MODE)"
echo "=========================================="

# Build merge command for all subjects
# This creates a file with 106 parcels (rows) × 20 subjects (columns)
MERGE_GROUP_CMD="wb_command -cifti-merge ${OUTPUT_PATH}/all_subjects_wta_full_corr_merged.pscalar.nii -direction ROW"

for sub in 01 02 03 04 05 06 07 08 09 11 12 13 14 17 20 21 22 23 24 25; do
    MERGE_GROUP_CMD="${MERGE_GROUP_CMD} -cifti ${OUTPUT_PATH}/sub-${sub}_indexmax_wta_full_corr.pscalar.nii"
done

# Inspect the shape of the output
wb_command -file-information ${OUTPUT_PATH}/sub-${sub}_indexmax_wta_full_corr.pscalar.nii

# Execute group merge
eval $MERGE_GROUP_CMD

echo "Merged all subjects into single file"

# Transpose so subjects are rows and parcels are columns
wb_command -cifti-transpose \
    "${OUTPUT_PATH}/all_subjects_wta_full_corr_merged.pscalar.nii" \
    "${OUTPUT_PATH}/all_subjects_wta_full_corr_merged.transpose.nii" \
    -mem-limit 20

echo "Transposed for MODE calculation"

# Calculate MODE across subjects (rows) for each parcel (column)
# This gives the most frequent winner across subjects for each parcel
wb_command -cifti-reduce \
    "${OUTPUT_PATH}/all_subjects_wta_full_corr_merged.transpose.nii" \
    MODE \
    "${OUTPUT_PATH}/group_wta_full_corr_mode.raw.nii" \
    -direction COLUMN

echo "Calculated MODE across subjects"

# Transpose back to standard orientation (parcels as rows)
wb_command -cifti-transpose \
    "${OUTPUT_PATH}/group_wta_full_corr_mode.raw.nii" \
    "${OUTPUT_PATH}/group_wta_full_corr.pscalar.nii" \
    -mem-limit 20

# Dump the main output into another format
wb_command -cifti-convert -to-text "${OUTPUT_PATH}/group_wta_full_corr.pscalar.nii" "${OUTPUT_PATH}/group_wta_full_corr.tsv"

# Verify that the label order matches the original key ascending order


# Change file permissions
chmod -Rf 771 ${BASE_PATH}
chgrp -Rf 327 ${BASE_PATH}

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
