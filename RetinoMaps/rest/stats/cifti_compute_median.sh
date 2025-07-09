#!/bin/bash

# With this script we go from the stacked files to a median ROI by all target vertices connectivity matrix

# Define base paths
BASE_DIR="/scratch/mszinte/data/RetinoMaps/derivatives/pp_data/group/91k/rest"

# Define source and destination directories
CORR_SOURCE="${BASE_DIR}/stacked_full_corr"
CORR_DEST="${BASE_DIR}/median_full_corr"

FISHER_SOURCE="${BASE_DIR}/stacked_fisher-z"
FISHER_DEST="${BASE_DIR}/median_fisher-z"

# Create destination directories if they don't exist
mkdir -p "$CORR_DEST"
mkdir -p "$FISHER_DEST"

# Process full correlation files
echo "Processing full correlation files..."
cd "$CORR_SOURCE"

for stacked_file in group_ses-01_task-rest_space-fsLR_den-91k_desc-full_corr_*_stacked.dscalar.nii; do
    if [[ -f "$stacked_file" ]]; then
        # Create name for median file by replacing "stacked" with "median"
        median_file="${stacked_file/stacked/median}"
        
        echo "Processing: $stacked_file"
        wb_command -cifti-reduce "$stacked_file" MEDIAN "${CORR_DEST}/${median_file}"
    fi
done

# Process fisher-z files
echo "Processing fisher-z files..."
cd "$FISHER_SOURCE"

for stacked_file in group_ses-01_task-rest_space-fsLR_den-91k_desc-full_corr_fisher-z_*_stacked.dscalar.nii; do
    if [[ -f "$stacked_file" ]]; then
        # Create name for median file by replacing "stacked" with "median"
        median_file="${stacked_file/stacked/median}"
        
        echo "Processing: $stacked_file"
        wb_command -cifti-reduce "$stacked_file" MEDIAN "${FISHER_DEST}/${median_file}"
    fi
done

echo "All processing complete!"
