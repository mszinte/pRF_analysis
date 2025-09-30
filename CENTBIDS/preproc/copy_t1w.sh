#!/bin/bash

# =========================================================================
# Copy T1w files from presurfer derivatives of MP2RAGE to BIDS anat folder
# =========================================================================

# Path to the settings file
SETTINGS_JSON="settings.json"

# Directories (adjust if different)
PROJECT_DIR="./"                # path to your BIDS root
DERIVATIVES_DIR="${PROJECT_DIR}/derivatives"

# Parse subjects and anatomical sessions from JSON
subjects=$(jq -r '.subjects[]' "${SETTINGS_JSON}")
anat_sessions=$(jq -r '.anat_session[]' "${SETTINGS_JSON}")

# Loop over subjects and sessions
for sub in $subjects; do
    sub_num=$(echo $sub | sed 's/sub-//')
    for ses in $anat_sessions; do
        echo "Processing $sub $ses ..."

        t1w_dest="${PROJECT_DIR}/${sub}/${ses}/anat/${sub}_${ses}_T1w.nii.gz"
        
        # Skip if destination already exists
        if [ -f "${t1w_dest}" ]; then
            echo "Skipping: T1w already exists for ${sub} ${ses}"
            continue
        fi

        # Source file
        mprage_source="${DERIVATIVES_DIR}/presurfer/${sub}/${ses}/presurf_MPRAGEise/${sub}_${ses}_acq-UNI_run-01_MP2RAGE_MPRAGEised.nii"

        if [ -f "${mprage_source}" ]; then
            mkdir -p "$(dirname "${t1w_dest}")"
            gzip -c "${mprage_source}" > "${t1w_dest}"
            echo "Copied: ${mprage_source} -> ${t1w_dest}"
        else
            echo "Warning: MPRAGE file not found for ${sub} ${ses}"
        fi
    done
done
