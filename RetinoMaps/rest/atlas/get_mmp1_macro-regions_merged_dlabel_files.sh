#!/bin/bash
#####################################################
# Written by Marco Bedini (marco.bedini@univ-amu.fr)
#####################################################

# Define paths
ATLAS_DIR="/scratch/mszinte/data/RetinoMaps/derivatives/pp_data/atlas/mmp1"
LABELS_DIR="$ATLAS_DIR/label_files"
MACRO_DIR="$ATLAS_DIR/macro_regions"
TEMPLATE="$ATLAS_DIR/atlas-Glasser_space-fsLR_den-32k_dseg.dlabel.nii"

PARCEL="$MACRO_DIR/atlas-Glasser_space-fsLR_den-32k_macro-regions_parcellation.dscalar.nii"
OUTPUT="$ATLAS_DIR/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_discarded_macro-regions.dlabel.nii"
COMBINED_LABELS="$LABELS_DIR/all_macro-regions_label.txt"

# Output dir
mkdir -p "$MACRO_DIR"

# Define regions: key = index+1 (mPCS=1 ... V1=12)
REGIONS=(mPCS sPCS iPCS sIPS iIPS "hMT+" VO LO V3AB V3 V2 V1)

# ---------------------------------------------------------------------------
# 1. Merge per-hemisphere binary scalars into one bilateral dscalar per region
# ---------------------------------------------------------------------------
for region in "${REGIONS[@]}"; do
    wb_command -cifti-create-dense-from-template \
        "$TEMPLATE" \
        "$MACRO_DIR/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_discarded_${region}.dscalar.nii" \
        -metric CORTEX_LEFT  "$MACRO_DIR/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_discarded_lh_${region}_bin.shape.gii" \
        -metric CORTEX_RIGHT "$MACRO_DIR/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_discarded_rh_${region}_bin.shape.gii"
done

# ---------------------------------------------------------------------------
# 2. Multiply each binary dscalar by its key value (1-12)
# ---------------------------------------------------------------------------
key=1
for region in "${REGIONS[@]}"; do
    wb_command -cifti-math "x * $key" \
        "$MACRO_DIR/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_discarded_${region}_key.dscalar.nii" \
        -var x "$MACRO_DIR/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_discarded_${region}.dscalar.nii"
    (( key++ ))
done

# ---------------------------------------------------------------------------
# 3. Sum all 12 keyed dscalars into one parcellation scalar
# ---------------------------------------------------------------------------
merge_cmd=(wb_command -cifti-merge "$PARCEL")
for region in "${REGIONS[@]}"; do
    merge_cmd+=(
        -cifti "$MACRO_DIR/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_discarded_${region}_key.dscalar.nii"
    )
done
"${merge_cmd[@]}"

wb_command -cifti-reduce "$PARCEL" SUM "$PARCEL"

# ---------------------------------------------------------------------------
# 4. Write combined label file (all 12 regions, keys 1-12)
# ---------------------------------------------------------------------------
declare -A roi_r=( [mPCS]=255 [sPCS]=255 [iPCS]=151 [sIPS]=44  [iIPS]=0   ["hMT+"]=0   [VO]=0   [LO]=150 [V3AB]=235 [V3]=248 [V2]=250 [V1]=243 )
declare -A roi_g=( [mPCS]=111 [sPCS]=234 [iPCS]=255 [sIPS]=255 [iIPS]=152 ["hMT+"]=25  [VO]=0   [LO]=0   [V3AB]=127 [V3]=160 [V2]=196 [V1]=231 )
declare -A roi_b=( [mPCS]=0   [sPCS]=0   [iPCS]=0   [sIPS]=150 [iIPS]=255 ["hMT+"]=255 [VO]=200 [LO]=90  [V3AB]=134 [V3]=126 [V2]=132 [V1]=155 )

> "$COMBINED_LABELS"   # truncate/create
key=1
for region in "${REGIONS[@]}"; do
    r="${roi_r[$region]}"
    g="${roi_g[$region]}"
    b="${roi_b[$region]}"
    printf '%s\n%d %d %d %d %d\n' "$region" "$key" "$r" "$g" "$b" 255 >> "$COMBINED_LABELS"
    (( key++ ))
done

# ---------------------------------------------------------------------------
# 5. Import all labels at once → single dlabel with 12 named color regions
# ---------------------------------------------------------------------------
wb_command -cifti-label-import \
    "$PARCEL" \
    "$COMBINED_LABELS" \
    "$OUTPUT" \
    -discard-others

# Fix permissions
chmod -Rf 771 "$ATLAS_DIR"
chgrp -Rf 327 "$ATLAS_DIR"

echo "Done. Output: $OUTPUT"