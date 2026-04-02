#!/bin/bash
#####################################################
# Winner-Take-All analysis using Connectome-Workbench
# Flexible for different correlation types and hemisphere modes
# Note: INDEXMAX returns winning COLUMN (seed cluster, e.g. sPCS) by ROW (MMP label)
#
# Written by Marco Bedini (marco.bedini@univ-amu.fr) with Claude's help
#
# Usage examples:
#   ./workbench_compute_wta.sh full_corr bilateral default
#   ./workbench_compute_wta.sh fisher-z by_hemi legacy
#####################################################

# Parse command-line arguments
CORR_TYPE="${1:-full_corr}"      # 'full_corr' or 'fisher-z'
HEMI_TYPE="${2:-bilateral}"      # 'bilateral' or 'by_hemi'
LEGACY_MODE="${3:-default}"      # 'default' or 'legacy'

# Validate inputs
if [[ "$CORR_TYPE" != "full_corr" ]] && [[ "$CORR_TYPE" != "fisher-z" ]]; then
    echo "Error: CORR_TYPE must be 'full_corr' or 'fisher-z'"
    echo "Usage: $0 [full_corr|fisher-z] [bilateral|by_hemi] [default|legacy]"
    exit 1
fi

if [[ "$HEMI_TYPE" != "bilateral" ]] && [[ "$HEMI_TYPE" != "by_hemi" ]]; then
    echo "Error: HEMI_TYPE must be 'bilateral' or 'by_hemi'"
    echo "Usage: $0 [full_corr|fisher-z] [bilateral|by_hemi] [default|legacy]"
    exit 1
fi

if [[ "$LEGACY_MODE" != "default" ]] && [[ "$LEGACY_MODE" != "legacy" ]]; then
    echo "Error: LEGACY_MODE must be 'default' or 'legacy'"
    echo "Usage: $0 [full_corr|fisher-z] [bilateral|by_hemi] [default|legacy]"
    exit 1
fi

# Set legacy flag for workbench
[[ "$LEGACY_MODE" == "legacy" ]] && LEGACY_FLAG="-legacy-mode" || LEGACY_FLAG=""

echo "=========================================="
echo "WTA ANALYSIS CONFIGURATION"
echo "=========================================="
echo "Correlation type: ${CORR_TYPE}"
echo "Hemisphere mode:  ${HEMI_TYPE}"
echo "Legacy mode:      ${LEGACY_MODE}"
echo "=========================================="
echo ""

# Main folders
BASE_PATH="/scratch/mszinte/data/RetinoMaps/derivatives/pp_data"
ATLAS_DIR="${BASE_PATH}/atlas"
OUTPUT_PATH="${BASE_PATH}/group/91k/rest/wta/workbench"
mkdir -p "${OUTPUT_PATH}"

# Define ROI order (seed clusters; INDEXMAX winner values map to 1-based position here)
declare -a ROIS=("mPCS" "sPCS" "iPCS" "sIPS" "iIPS" "hMT+" "VO" "LO" "V3AB" "V3" "V2" "V1")

# Construct file suffix based on options
FILE_SUFFIX="${CORR_TYPE}"
[[ "$HEMI_TYPE" == "by_hemi" ]] && FILE_SUFFIX="${FILE_SUFFIX}_by_hemi"
[[ "$LEGACY_MODE" == "legacy" ]] && FILE_SUFFIX="${FILE_SUFFIX}_legacy"

# Determine which hemispheres to process
if [[ "$HEMI_TYPE" == "by_hemi" ]]; then
    HEMIS=("lh" "rh")
else
    HEMIS=("bilateral")
fi

declare -a subjects=("01" "02" "03" "04" "05" "06" "07" "08" "09" "11" "12" "13" "14" "17" "20" "21" "22" "23" "24" "25")

# ============================================
# HEMISPHERE LOOP
# ============================================

for HEMI in "${HEMIS[@]}"; do

    if [[ "$HEMI" != "bilateral" ]]; then
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

    for sub in "${subjects[@]}"; do

        echo "Processing sub-${sub}${HEMI_SUFFIX}..."

        FULL_CORR="${BASE_PATH}/sub-${sub}/91k/rest/corr/full_corr/wta"
        mkdir -p "${FULL_CORR}"

        for roi in "${ROIS[@]}"; do

            INPUT_FILE="${BASE_PATH}/sub-${sub}/91k/rest/corr/full_corr/sub-${sub}_task-rest_space-fsLR_den-91k_desc-${CORR_TYPE}_${roi}.dscalar.nii"

            # Separate hemispheres for masking
            wb_command -cifti-separate "${INPUT_FILE}" COLUMN \
                -metric CORTEX_LEFT  "${FULL_CORR}/sub-${sub}_${CORR_TYPE}_lh_${roi}.shape.gii" \
                -metric CORTEX_RIGHT "${FULL_CORR}/sub-${sub}_${CORR_TYPE}_rh_${roi}.shape.gii"

            # Mask left hemisphere: set non-seed-cluster vertices to -1
            # Safer to use this rather than masking with cifti-mask (this produces 0s that we don't want when parcellating the output)
            wb_command -metric-math \
                '(corr * (hollow_seed > 0)) + ((hollow_seed == 0) * -1)' \
                "${FULL_CORR}/sub-${sub}_${CORR_TYPE}_lh_${roi}_hollow.shape.gii" \
                -var corr        "${FULL_CORR}/sub-${sub}_${CORR_TYPE}_lh_${roi}.shape.gii" \
                -var hollow_seed "${ATLAS_DIR}/mmp1_clusters/leaveout/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_leaveout_lh_${roi}_bin.shape.gii"

            # Mask right hemisphere
            wb_command -metric-math \
                '(corr * (hollow_seed > 0)) + ((hollow_seed == 0) * -1)' \
                "${FULL_CORR}/sub-${sub}_${CORR_TYPE}_rh_${roi}_hollow.shape.gii" \
                -var corr        "${FULL_CORR}/sub-${sub}_${CORR_TYPE}_rh_${roi}.shape.gii" \
                -var hollow_seed "${ATLAS_DIR}/mmp1_clusters/leaveout/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_leaveout_rh_${roi}_bin.shape.gii"

            # Recombine hemispheres into dense scalar
            wb_command -cifti-create-dense-scalar \
                "${FULL_CORR}/sub-${sub}_${CORR_TYPE}_${roi}_hollow.dscalar.nii" \
                -left-metric  "${FULL_CORR}/sub-${sub}_${CORR_TYPE}_lh_${roi}_hollow.shape.gii" \
                -right-metric "${FULL_CORR}/sub-${sub}_${CORR_TYPE}_rh_${roi}_hollow.shape.gii"

            # Parcellate: mean correlation per MMP parcel
            wb_command -cifti-parcellate \
                "${FULL_CORR}/sub-${sub}_${CORR_TYPE}_${roi}_hollow.dscalar.nii" \
                "${ATLAS_DIR}/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_dseg.dlabel.nii" COLUMN \
                "${FULL_CORR}/sub-${sub}_${CORR_TYPE}${HEMI_SUFFIX}_${roi}.pscalar.nii" \
                -method MEAN -only-numeric ${LEGACY_FLAG}

            # Convert to TSV for violin plots
            wb_command -cifti-convert -to-text \
                "${FULL_CORR}/sub-${sub}_${CORR_TYPE}${HEMI_SUFFIX}_${roi}.pscalar.nii" \
                "${FULL_CORR}/sub-${sub}_${CORR_TYPE}${HEMI_SUFFIX}_${roi}.tsv"

        done

        # Merge all ROIs into a single pscalar (rows = ROIs, cols = parcels)
        MERGE_INPUT="${FULL_CORR}/sub-${sub}_${FILE_SUFFIX}${HEMI_SUFFIX}_merged.pscalar.nii"
        MERGE_CMD="wb_command -cifti-merge ${MERGE_INPUT} -direction ROW"
        for roi in "${ROIS[@]}"; do
            MERGE_CMD="${MERGE_CMD} -cifti ${FULL_CORR}/sub-${sub}_${CORR_TYPE}${HEMI_SUFFIX}_${roi}.pscalar.nii"
        done
        eval $MERGE_CMD

        # Transpose so INDEXMAX operates along the ROI dimension per parcel
        wb_command -cifti-transpose \
            "${MERGE_INPUT}" \
            "${FULL_CORR}/sub-${sub}_${FILE_SUFFIX}${HEMI_SUFFIX}_merged_transpose.nii" \
            -mem-limit 20

        # INDEXMAX: for each parcel, which ROI had the highest correlation?
        # Output: 106 parcels × 1 column with winner index (1-12)
        wb_command -cifti-reduce \
            "${FULL_CORR}/sub-${sub}_${FILE_SUFFIX}${HEMI_SUFFIX}_merged_transpose.nii" \
            INDEXMAX \
            "${OUTPUT_PATH}/sub-${sub}_indexmax_wta_${FILE_SUFFIX}${HEMI_SUFFIX}_raw.nii" \
            -direction COLUMN

        # Transpose back to standard pscalar orientation
        wb_command -cifti-transpose \
            "${OUTPUT_PATH}/sub-${sub}_indexmax_wta_${FILE_SUFFIX}${HEMI_SUFFIX}_raw.nii" \
            "${OUTPUT_PATH}/sub-${sub}_indexmax_wta_${FILE_SUFFIX}${HEMI_SUFFIX}.pscalar.nii" \
            -mem-limit 20

        echo "  ✓ Completed sub-${sub}${HEMI_SUFFIX}"
    done

    # ============================================
    # GROUP-LEVEL: compute MODE across subjects
    # ============================================

    echo ""
    echo "=========================================="
    echo "STEP 2: Computing group-level WTA (MODE)"
    echo "=========================================="

    # Merge all subjects' winner maps (rows = subjects, cols = parcels)
    GROUP_MERGED="${OUTPUT_PATH}/all_subjects_wta_${FILE_SUFFIX}${HEMI_SUFFIX}_merged.pscalar.nii"
    MERGE_GROUP_CMD="wb_command -cifti-merge ${GROUP_MERGED} -direction ROW"
    for sub in "${subjects[@]}"; do
        MERGE_GROUP_CMD="${MERGE_GROUP_CMD} -cifti ${OUTPUT_PATH}/sub-${sub}_indexmax_wta_${FILE_SUFFIX}${HEMI_SUFFIX}.pscalar.nii"
    done
    eval $MERGE_GROUP_CMD
    echo "Merged all subjects"

    # Transpose for column-wise MODE
    wb_command -cifti-transpose \
        "${GROUP_MERGED}" \
        "${OUTPUT_PATH}/all_subjects_wta_${FILE_SUFFIX}${HEMI_SUFFIX}_merged_transpose.nii" \
        -mem-limit 20
    echo "Transposed for MODE calculation"

    # MODE: most common winner per parcel across subjects
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

    # Convert to TSV for downstream colormap script
    wb_command -cifti-convert -to-text \
        "${OUTPUT_PATH}/group_wta_${FILE_SUFFIX}${HEMI_SUFFIX}.pscalar.nii" \
        "${OUTPUT_PATH}/group_wta_${FILE_SUFFIX}${HEMI_SUFFIX}.tsv"

    echo "Group-level WTA map created for hemi=${HEMI}!"

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
echo "  Hemisphere:  ${HEMI_TYPE}"
echo "  Legacy mode: ${LEGACY_MODE}"
echo ""
echo "Output files:"
if [[ "$HEMI_TYPE" == "bilateral" ]]; then
    echo "  Subjects: ${OUTPUT_PATH}/sub-*_indexmax_wta_${FILE_SUFFIX}.pscalar.nii"
    echo "  Group:    ${OUTPUT_PATH}/group_wta_${FILE_SUFFIX}.pscalar.nii"
    echo "  TSV:      ${OUTPUT_PATH}/group_wta_${FILE_SUFFIX}.tsv"
else
    for hemi in lh rh; do
        echo "  ${hemi^^} Subjects: ${OUTPUT_PATH}/sub-*_indexmax_wta_${FILE_SUFFIX}_${hemi}.pscalar.nii"
        echo "  ${hemi^^} Group:    ${OUTPUT_PATH}/group_wta_${FILE_SUFFIX}_${hemi}.pscalar.nii"
        echo "  ${hemi^^} TSV:      ${OUTPUT_PATH}/group_wta_${FILE_SUFFIX}_${hemi}.tsv"
    done
fi
echo ""
echo "Winner seed mapping:"
for i in "${!ROIS[@]}"; do
    echo "  $((i+1)) = ${ROIS[$i]}"
done
echo ""

# Set permissions
chmod -Rf 771 "${BASE_PATH}"
chgrp -Rf 327 "${BASE_PATH}"