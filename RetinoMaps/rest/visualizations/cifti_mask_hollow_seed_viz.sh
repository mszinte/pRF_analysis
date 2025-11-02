#!/bin/bash

## Important note: we are going to use only Pearson full correlations for viz not the Fisher-z ones (we are looking at the latter only to check the differences)

# Define paths
DATA_PATH_FISHER="/scratch/mszinte/data/RetinoMaps/derivatives/pp_data/group/91k/rest/median_fisher-z"
DATA_PATH_CORR="/scratch/mszinte/data/RetinoMaps/derivatives/pp_data/group/91k/rest/median_full_corr"
ATLAS_PATH="/home/${USER}/projects/pRF_analysis/RetinoMaps/rest/mmp1_clusters"
OUTPUT_PATH_CORR="${DATA_PATH_CORR}/hollow_seed_viz_full_corr"
OUTPUT_PATH_FISHER="${DATA_PATH_FISHER}/hollow_seed_viz_fisher-z"

# Create output directories if they don't exist
mkdir -p "$OUTPUT_PATH_CORR"
mkdir -p "$OUTPUT_PATH_FISHER"

# Process full correlation files
echo "Processing full correlation files..."
for region in lh_mPCS lh_sPCS lh_iPCS lh_sIPS lh_iIPS lh_hMT+ lh_VO lh_LO lh_V3AB lh_V3 lh_V2 lh_V1 rh_mPCS rh_sPCS rh_iPCS rh_sIPS rh_iIPS rh_hMT+ rh_VO rh_LO rh_V3AB rh_V3 rh_V2 rh_V1 mPCS sPCS iPCS sIPS iIPS hMT+ VO LO V3AB V3 V2 V1; do
    if [[ "$region" == lh_* ]]; then
        specific_roi_flag="-left-roi"
        specific_roi_file="${ATLAS_PATH}/leaveout/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_leaveout_${region}.shape.gii"
        opposite_roi_flag="-right-roi"
        opposite_roi_file="${ATLAS_PATH}/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_right_hemi.shape.gii"
    elif [[ "$region" == rh_* ]]; then
        specific_roi_flag="-right-roi"
        specific_roi_file="${ATLAS_PATH}/leaveout/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_leaveout_${region}.shape.gii"
        opposite_roi_flag="-left-roi"
        opposite_roi_file="${ATLAS_PATH}/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_left_hemi.shape.gii"
    else
        specific_roi_flag="-left-roi"
        specific_roi_file="${ATLAS_PATH}/leaveout/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_leaveout_lh_${region}.shape.gii"
        opposite_roi_flag="-right-roi"
        opposite_roi_file="${ATLAS_PATH}/leaveout/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_leaveout_rh_${region}.shape.gii"
    fi

    input_file="${DATA_PATH_CORR}/group_ses-01_task-rest_space-fsLR_den-91k_desc-full_corr_${region}_median.dscalar.nii"
    output_file="${OUTPUT_PATH_CORR}/group_ses-01_task-rest_space-fsLR_den-91k_desc-full_corr_${region}_hollow_seed.dscalar.nii"

    if [[ -f "$input_file" && -f "$specific_roi_file" && -f "$opposite_roi_file" ]]; then
        echo "Masking: $input_file"
        wb_command -cifti-restrict-dense-map "$input_file" COLUMN "$output_file" \
            "${specific_roi_flag}" "$specific_roi_file" \
            "${opposite_roi_flag}" "$opposite_roi_file"
    else
        echo "Warning: Missing file for ${region} - input: $([ -f "$input_file" ] && echo "exists" || echo "missing"), specific mask: $([ -f "$specific_roi_file" ] && echo "exists" || echo "missing"), opposite mask: $([ -f "$opposite_roi_file" ] && echo "exists" || echo "missing")"
    fi

    input_masked_file="${DATA_PATH_CORR}/group_ses-01_task-rest_space-fsLR_den-91k_desc-full_corr_${region}_masked_median.dscalar.nii"
    output_masked_file="${OUTPUT_PATH_CORR}/group_ses-01_task-rest_space-fsLR_den-91k_desc-full_corr_${region}_masked_hollow_seed.dscalar.nii"

    if [[ -f "$input_masked_file" && -f "$specific_roi_file" && -f "$opposite_roi_file" ]]; then
        echo "Masking: $input_masked_file"
        wb_command -cifti-restrict-dense-map "$input_masked_file" COLUMN "$output_masked_file" \
            "${specific_roi_flag}" "$specific_roi_file" \
            "${opposite_roi_flag}" "$opposite_roi_file"
    fi
done

# Process fisher-z correlation files
echo "Processing fisher-z correlation files..."
for region in lh_mPCS lh_sPCS lh_iPCS lh_sIPS lh_iIPS lh_hMT+ lh_VO lh_LO lh_V3AB lh_V3 lh_V2 lh_V1 rh_mPCS rh_sPCS rh_iPCS rh_sIPS rh_iIPS rh_hMT+ rh_VO rh_LO rh_V3AB rh_V3 rh_V2 rh_V1 mPCS sPCS iPCS sIPS iIPS hMT+ VO LO V3AB V3 V2 V1; do
    if [[ "$region" == lh_* ]]; then
        specific_roi_flag="-left-roi"
        specific_roi_file="${ATLAS_PATH}/leaveout/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_leaveout_${region}.shape.gii"
        opposite_roi_flag="-right-roi"
        opposite_roi_file="${ATLAS_PATH}/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_right_hemi.shape.gii"
    elif [[ "$region" == rh_* ]]; then
        specific_roi_flag="-right-roi"
        specific_roi_file="${ATLAS_PATH}/leaveout/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_leaveout_${region}.shape.gii"
        opposite_roi_flag="-left-roi"
        opposite_roi_file="${ATLAS_PATH}/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_left_hemi.shape.gii"
    else
        specific_roi_flag="-left-roi"
        specific_roi_file="${ATLAS_PATH}/leaveout/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_leaveout_lh_${region}.shape.gii"
        opposite_roi_flag="-right-roi"
        opposite_roi_file="${ATLAS_PATH}/leaveout/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_leaveout_rh_${region}.shape.gii"
    fi

    input_file="${DATA_PATH_FISHER}/group_ses-01_task-rest_space-fsLR_den-91k_desc-full_corr_fisher-z_${region}_median.dscalar.nii"
    output_file="${OUTPUT_PATH_FISHER}/group_ses-01_task-rest_space-fsLR_den-91k_desc-full_corr_fisher-z_${region}_hollow_seed.dscalar.nii"

    if [[ -f "$input_file" && -f "$specific_roi_file" && -f "$opposite_roi_file" ]]; then
        echo "Masking: $input_file"
        wb_command -cifti-restrict-dense-map "$input_file" COLUMN "$output_file" \
            "${specific_roi_flag}" "$specific_roi_file" \
            "${opposite_roi_flag}" "$opposite_roi_file"
    else
        echo "Warning: Missing file for ${region} - input: $([ -f "$input_file" ] && echo "exists" || echo "missing"), specific mask: $([ -f "$specific_roi_file" ] && echo "exists" || echo "missing"), opposite mask: $([ -f "$opposite_roi_file" ] && echo "exists" || echo "missing")"
    fi

    input_masked_file="${DATA_PATH_FISHER}/group_ses-01_task-rest_space-fsLR_den-91k_desc-full_corr_fisher-z_${region}_masked_median.dscalar.nii"
    output_masked_file="${OUTPUT_PATH_FISHER}/group_ses-01_task-rest_space-fsLR_den-91k_desc-full_corr_fisher-z_${region}_masked_hollow_seed.dscalar.nii"

    if [[ -f "$input_masked_file" && -f "$specific_roi_file" && -f "$opposite_roi_file" ]]; then
        echo "Masking: $input_masked_file"
        wb_command -cifti-restrict-dense-map "$input_masked_file" COLUMN "$output_masked_file" \
            "${specific_roi_flag}" "$specific_roi_file" \
            "${opposite_roi_flag}" "$opposite_roi_file"
    fi
done

echo "All masking operations complete!"

