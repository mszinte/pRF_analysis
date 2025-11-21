#!/bin/bash

# Goal: create GIFTI surfaces to import data in Python (nibabel-compatible)
# Outputs go into the winner-take-all script

#####################################################
# Written by Marco Bedini (marco.bedini@univ-amu.fr)
#####################################################

# Base paths
FISHERZ_PATH="/scratch/mszinte/data/RetinoMaps/derivatives/pp_data/group/91k/rest/hollow_seed_viz_fisher-z"
FULLCORR_PATH="/scratch/mszinte/data/RetinoMaps/derivatives/pp_data/group/91k/rest/hollow_seed_viz_full_corr"

# Output folder
OUT_GIFTI="/scratch/mszinte/data/RetinoMaps/derivatives/pp_data/group/91k/rest/wta/surfaces"
mkdir -p "$OUT_GIFTI"

# === Loop 1: Fisher-z files ===
for FILE in "$FISHERZ_PATH"/group_ses-01_task-rest_space-fsLR_den-91k_desc-full_corr_fisher-z_*_hollow_seed.dscalar.nii; do
  [ -e "$FILE" ] || continue
  BASENAME=$(basename "$FILE" .dscalar.nii)
  echo "🧠 Splitting Fisher-z: $BASENAME"

  OUT_LH="$OUT_GIFTI/${BASENAME}_lh.func.gii"
  OUT_RH="$OUT_GIFTI/${BASENAME}_rh.func.gii"

  wb_command -cifti-separate "$FILE" COLUMN \
    -metric CORTEX_LEFT "$OUT_LH" \
    -metric CORTEX_RIGHT "$OUT_RH"
done

# === Loop 2: Full correlation files ===
for FILE in "$FULLCORR_PATH"/group_ses-01_task-rest_space-fsLR_den-91k_desc-full_corr_*_hollow_seed.dscalar.nii; do
  [ -e "$FILE" ] || continue
  BASENAME=$(basename "$FILE" .dscalar.nii)
  echo "🧠 Splitting Full correlation: $BASENAME"

  OUT_LH="$OUT_GIFTI/${BASENAME}_lh.func.gii"
  OUT_RH="$OUT_GIFTI/${BASENAME}_rh.func.gii"

  wb_command -cifti-separate "$FILE" COLUMN \
    -metric CORTEX_LEFT "$OUT_LH" \
    -metric CORTEX_RIGHT "$OUT_RH"
done

