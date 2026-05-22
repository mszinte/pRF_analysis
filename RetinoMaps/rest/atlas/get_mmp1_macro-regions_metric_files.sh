#!/bin/bash

#####################################################
# Written by Marco Bedini (marco.bedini@univ-amu.fr)
#####################################################

# Define paths
ATLAS_DIR="/scratch/mszinte/data/RetinoMaps/derivatives/pp_data/atlas/mmp1"
LABELS_DIR="$ATLAS_DIR/label_files"

# Output dir
mkdir -p "$ATLAS_DIR/macro_regions"

# Define regions and hemispheres
REGIONS=(mPCS sPCS iPCS sIPS iIPS hMT+ VO LO V3AB V3 V2 V1)

# 1. Merge labels for all macro-regions
for region in "${REGIONS[@]}"; do

    wb_command -cifti-label-import "$ATLAS_DIR/atlas-Glasser_space-fsLR_den-32k_dseg.dlabel.nii" "$LABELS_DIR/lh_${region}_labels.txt" \
    	"$ATLAS_DIR/macro_regions/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_discarded_lh_${region}.dlabel.nii" -discard-others

	wb_command -cifti-label-import "$ATLAS_DIR/atlas-Glasser_space-fsLR_den-32k_dseg.dlabel.nii" "$LABELS_DIR/rh_${region}_labels.txt" \
    	"$ATLAS_DIR/macro_regions/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_discarded_rh_${region}.dlabel.nii" -discard-others

# 2. Get the dscalar files (useful later to make new dlabel files for macro-regions)

	wb_command -cifti-reduce "$ATLAS_DIR/macro_regions/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_discarded_lh_${region}.dlabel.nii" \
	COUNT_NONZERO "$ATLAS_DIR/macro_regions/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_discarded_lh_${region}.dscalar.nii"

	wb_command -cifti-reduce "$ATLAS_DIR/macro_regions/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_discarded_rh_${region}.dlabel.nii" \
	COUNT_NONZERO "$ATLAS_DIR/macro_regions/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_discarded_rh_${region}.dscalar.nii"

# 3. Get the corresponding metric files and binarize outputs

	# Separate metric files
	wb_command -cifti-separate "$ATLAS_DIR/macro_regions/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_discarded_lh_${region}.dlabel.nii" COLUMN \
		-metric CORTEX_LEFT "$ATLAS_DIR/macro_regions/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_discarded_lh_${region}.shape.gii"
		
	wb_command -cifti-separate "$ATLAS_DIR/macro_regions/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_discarded_rh_${region}.dlabel.nii" COLUMN \
		-metric CORTEX_RIGHT "$ATLAS_DIR/macro_regions/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_discarded_rh_${region}.shape.gii"

	# Binarize macro-regions by hemi metric files
	wb_command -metric-reduce "$ATLAS_DIR/macro_regions/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_discarded_lh_${region}.shape.gii" COUNT_NONZERO \
		"$ATLAS_DIR/macro_regions/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_discarded_lh_${region}_bin.shape.gii" -only-numeric

	# Binarize macro-regions by hemi metric files
	wb_command -metric-reduce "$ATLAS_DIR/macro_regions/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_discarded_rh_${region}.shape.gii" COUNT_NONZERO \
		"$ATLAS_DIR/macro_regions/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_discarded_rh_${region}_bin.shape.gii" -only-numeric
done

# Change permissions
chmod -Rf 771 *
chgrp -Rf 327 *
