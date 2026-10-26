#!/bin/bash

## Important note: we are going to use only Fisher-z maps for stats not Pearson correlations (uncomment the lines if you need them)

#####################################################
# Written by Marco Bedini (marco.bedini@univ-amu.fr)
#####################################################

# Define paths
DATA_PATH_FISHER="/media/marc_be/marc_be_vault1/RetinoMaps/derivatives/func_connectivity/Retinomaps_ROIs/data/correlations/group/stacked_fisher-z"
DATA_PATH_CORR="/media/marc_be/marc_be_vault1/RetinoMaps/derivatives/func_connectivity/Retinomaps_ROIs/data/correlations/group/stacked_full_corr"
ATLAS_PATH="/media/marc_be/marc_be_vault1/RetinoMaps/derivatives/func_connectivity/Retinomaps_ROIs/atlas/macro_regions/leaveout"
HEMI_ATLAS_PATH="/media/marc_be/marc_be_vault1/RetinoMaps/derivatives/func_connectivity/Retinomaps_ROIs/atlas"
OUTPUT_PATH_FISHER="${DATA_PATH_FISHER}/hollow_seed_for_stats"
# OUTPUT_PATH_CORR="${DATA_PATH_CORR}/hollow_seed_for_stats"

# Create output directories if they don't exist
mkdir -p "$OUTPUT_PATH_FISHER"
# mkdir -p "$OUTPUT_PATH_CORR"

# Define ROI list (hemi-specific and bilateral)
ALL_REGIONS=(lh_mPCS lh_sPCS lh_iPCS lh_sIPS lh_iIPS lh_hMT+ lh_VO lh_LO lh_V3AB lh_V3 lh_V2 lh_V1 rh_mPCS rh_sPCS rh_iPCS rh_sIPS rh_iIPS rh_hMT+ rh_VO rh_LO rh_V3AB rh_V3 rh_V2 rh_V1 mPCS sPCS iPCS sIPS iIPS hMT+ VO LO V3AB V3 V2 V1)

# Helper function to apply masks
apply_mask() {
    local input="$1"
    local output="$2"
    local left_mask="$3"
    local right_mask="$4"
    local apply_left="$5"
    local apply_right="$6"

    if [[ -f "$input" ]]; then
        CMD=("wb_command" "-cifti-restrict-dense-map" "$input" "COLUMN" "$output")
        if [[ "$apply_left" == "1" && -f "$left_mask" ]]; then
            CMD+=("-left-roi" "$left_mask")
        fi
        if [[ "$apply_right" == "1" && -f "$right_mask" ]]; then
            CMD+=("-right-roi" "$right_mask")
        fi
        echo "Masking: $input → $output"
        "${CMD[@]}"
    else
        echo "⚠️  Missing input: $input"
    fi
}

# Loop over all ROIs
for region in "${ALL_REGIONS[@]}"; do
    if [[ "$region" == lh_* ]]; then
        base_roi="${region#lh_}"
        left_mask="${ATLAS_PATH}/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_discarded_leaveout_${region}.shape.gii"
        right_mask="${HEMI_ATLAS_PATH}/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_discarded_right_hemi.shape.gii"
        apply_left=1
        apply_right=1  # still apply full mask in both hemispheres
    elif [[ "$region" == rh_* ]]; then
        base_roi="${region#rh_}"
        left_mask="${HEMI_ATLAS_PATH}/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_discarded_left_hemi.shape.gii"
        right_mask="${ATLAS_PATH}/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_discarded_leaveout_${region}.shape.gii"
        apply_left=1
        apply_right=1
    else
        base_roi="$region"
        left_mask="${ATLAS_PATH}/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_discarded_leaveout_lh_${region}.shape.gii"
        right_mask="${ATLAS_PATH}/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_discarded_leaveout_rh_${region}.shape.gii"
        apply_left=1
        apply_right=1
    fi

    # === FISHER-Z ===
    # Unmasked
    in_file_fz="${DATA_PATH_FISHER}/group_ses-01_task-rest_space-fsLR_den-91k_desc-full_corr_fisher-z_${region}_stacked.dscalar.nii"
    out_file_fz="${OUTPUT_PATH_FISHER}/group_ses-01_task-rest_space-fsLR_den-91k_desc-full_corr_fisher-z_${region}_stacked_hollow_seed.dscalar.nii"
    apply_mask "$in_file_fz" "$out_file_fz" "$left_mask" "$right_mask" "$apply_left" "$apply_right"

    # Masked
    in_masked_fz="${DATA_PATH_FISHER}/group_ses-01_task-rest_space-fsLR_den-91k_desc-full_corr_fisher-z_${region}_masked_stacked.dscalar.nii"
    out_masked_fz="${OUTPUT_PATH_FISHER}/group_ses-01_task-rest_space-fsLR_den-91k_desc-full_corr_fisher-z_${region}_masked_stacked_hollow_seed.dscalar.nii"
    apply_mask "$in_masked_fz" "$out_masked_fz" "$left_mask" "$right_mask" "$apply_left" "$apply_right"

    # === FULL CORR ===
    # Unmasked
#   in_file_corr="${DATA_PATH_CORR}/group_ses-01_task-rest_space-fsLR_den-91k_desc-full_corr_${region}_stacked.dscalar.nii"
#   out_file_corr="${OUTPUT_PATH_CORR}/group_ses-01_task-rest_space-fsLR_den-91k_desc-full_corr_${region}_stacked_hollow_seed.dscalar.nii"
#   apply_mask "$in_file_corr" "$out_file_corr" "$left_mask" "$right_mask" "$apply_left" "$apply_right"

    # Masked
#   in_masked_corr="${DATA_PATH_CORR}/group_ses-01_task-rest_space-fsLR_den-91k_desc-full_corr_${region}_masked_stacked.dscalar.nii"
#   out_masked_corr="${OUTPUT_PATH_CORR}/group_ses-01_task-rest_space-fsLR_den-91k_desc-full_corr_${region}_masked_stacked_hollow_seed.dscalar.nii"
#   apply_mask "$in_masked_corr" "$out_masked_corr" "$left_mask" "$right_mask" "$apply_left" "$apply_right"
done

echo "✅ All masking operations complete!"

