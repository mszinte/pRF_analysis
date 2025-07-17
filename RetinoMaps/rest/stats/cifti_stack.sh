#!/bin/bash

## H Stack subject-specific correlation matrices into one group file (rows = subjects)
## Useful so we can use -cifti-reduce later to compute the median values

# Define base paths
TASK_RESULTS="/scratch/mszinte/data/RetinoMaps/derivatives/pp_data"
OUTPUT_PATH="$TASK_RESULTS/group/91k/rest/"

# Define subjects, ROIs, data types
SUBJECTS=(01 02 03 04 05 06 07 08 09 11 12 13 14 17 20 21 22 23 24 25)
ROIS=(mPCS sPCS iPCS sIPS iIPS hMT+ VO LO V3AB V3 V2 V1)
# HEMIS=(lh rh both)
TYPES=("full_corr" "fisher-z")

# Loop over ROIs, hemispheres, and data types
for ROI in "${ROIS[@]}"; do
    for TYPE in "${TYPES[@]}"; do

    # Choose subfolder and descriptor prefix based on type
    if [[ "$TYPE" == "full_corr" ]]; then
      SUBFOLDER="full_corr"
      DESC="desc-full_corr"
    elif [[ "$TYPE" == "fisher-z" ]]; then
      SUBFOLDER="fisher-z"
      DESC="desc-fisher-z"
    fi

    # Final output filename
    OUTPUT_FILE="$OUTPUT_PATH/group_ses-01_task-rest_space-fsLR_den-91k_${DESC}_stacked.dscalar.nii"

    # Build merge command
    MERGE_CMD=("wb_command" "-cifti-merge" "$OUTPUT_FILE")

    for SUB in "${SUBJECTS[@]}"; do
      INPUT_FILE="$TASK_RESULTS/sub-${SUB}/91k/rest/corr/${SUBFOLDER}/sub-${SUB}_ses-01_task-rest_space-fsLR_den-91k_${DESC}_${ROI}_masked.dscalar.nii"
      if [[ -f "$INPUT_FILE" ]]; then
        MERGE_CMD+=("-cifti" "$INPUT_FILE")
      else
        echo "⚠️  Missing: $INPUT_FILE"
      fi
    done

    echo "🔄 Merging TYPE=${TYPE} | ROI=${ROI}..."
      "${MERGE_CMD[@]}"

    done
done

