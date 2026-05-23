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
OUTPUT="$ATLAS_DIR/atlas-Glasser_space-fsLR_den-32k_macro-regions.dlabel.nii"
COMBINED_LABELS="$LABELS_DIR/all_macro-regions_label.txt"

# Output dir
mkdir -p "$MACRO_DIR"

# Define regions and hemispheres
# Keys run 1-24: lh_mPCS=1 ... lh_V1=12, rh_mPCS=13 ... rh_V1=24
REGIONS=(mPCS sPCS iPCS sIPS iIPS "hMT+" VO LO V3AB V3 V2 V1)
HEMIS=(lh rh)

declare -A HEMI_STRUCTURE=( [lh]=CORTEX_LEFT [rh]=CORTEX_RIGHT )
declare -A HEMI_PREFIX=(    [lh]=LH          [rh]=RH           )

declare -A roi_r=( [mPCS]=255 [sPCS]=255 [iPCS]=151 [sIPS]=44  [iIPS]=0   ["hMT+"]=0   [VO]=0   [LO]=150 [V3AB]=235 [V3]=248 [V2]=250 [V1]=243 )
declare -A roi_g=( [mPCS]=111 [sPCS]=234 [iPCS]=255 [sIPS]=255 [iIPS]=152 ["hMT+"]=25  [VO]=0   [LO]=0   [V3AB]=127 [V3]=160 [V2]=196 [V1]=231 )
declare -A roi_b=( [mPCS]=0   [sPCS]=0   [iPCS]=0   [sIPS]=150 [iIPS]=255 ["hMT+"]=255 [VO]=200 [LO]=90  [V3AB]=134 [V3]=126 [V2]=132 [V1]=155 )

# ---------------------------------------------------------------------------
# 1. Create one dscalar per hemisphere/region using a single cortical structure
# ---------------------------------------------------------------------------
for hemi in "${HEMIS[@]}"; do
    structure="${HEMI_STRUCTURE[$hemi]}"
    for region in "${REGIONS[@]}"; do
        wb_command -cifti-create-dense-from-template \
            "$TEMPLATE" \
            "$MACRO_DIR/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_discarded_${hemi}_${region}.dscalar.nii" \
            -metric "$structure" \
            "$MACRO_DIR/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_discarded_${hemi}_${region}_bin.shape.gii"
    done
done

# ---------------------------------------------------------------------------
# 2. Multiply each dscalar by its unique key (lh: 1-12, rh: 13-24)
# ---------------------------------------------------------------------------
key=1
for hemi in "${HEMIS[@]}"; do
    for region in "${REGIONS[@]}"; do
        wb_command -cifti-math "x * $key" \
            "$MACRO_DIR/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_discarded_${hemi}_${region}_key.dscalar.nii" \
            -var x "$MACRO_DIR/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_discarded_${hemi}_${region}.dscalar.nii"
        (( key++ ))
    done
done

# ---------------------------------------------------------------------------
# 3. Sum all 24 keyed dscalars into one parcellation scalar
# ---------------------------------------------------------------------------
merge_cmd=(wb_command -cifti-merge "$PARCEL")
for hemi in "${HEMIS[@]}"; do
    for region in "${REGIONS[@]}"; do
        merge_cmd+=(
            -cifti "$MACRO_DIR/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_discarded_${hemi}_${region}_key.dscalar.nii"
        )
    done
done
"${merge_cmd[@]}"

wb_command -cifti-reduce "$PARCEL" SUM "$PARCEL"

# ---------------------------------------------------------------------------
# 4. Write combined label file (24 entries, keys 1-24, prefixed LH_/RH_)
# ---------------------------------------------------------------------------
> "$COMBINED_LABELS"
key=1
for hemi in "${HEMIS[@]}"; do
    prefix="${HEMI_PREFIX[$hemi]}"
    for region in "${REGIONS[@]}"; do
        r="${roi_r[$region]}"
        g="${roi_g[$region]}"
        b="${roi_b[$region]}"
        printf '%s\n%d %d %d %d %d\n' "${prefix}_${region}" "$key" "$r" "$g" "$b" 255 >> "$COMBINED_LABELS"
        (( key++ ))
    done
done

# ---------------------------------------------------------------------------
# 5. Import all 24 labels → single dlabel with hemisphere-specific regions
# ---------------------------------------------------------------------------
wb_command -cifti-label-import \
    "$PARCEL" \
    "$COMBINED_LABELS" \
    "$OUTPUT"

# Fix permissions
chmod -Rf 771 "$ATLAS_DIR"
chgrp -Rf 327 "$ATLAS_DIR"

echo "Done. Output: $OUTPUT"