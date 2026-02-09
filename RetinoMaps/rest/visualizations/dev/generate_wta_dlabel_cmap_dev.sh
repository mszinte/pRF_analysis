#!/bin/bash
#####################################################
# Generate label file with winner-take-all colors from TSV
# Written by Marco Bedini (marco.bedini@univ-amu.fr)
#####################################################

BASE_PATH="/scratch/mszinte/data/RetinoMaps/derivatives/pp_data"
ATLAS_DIR="/home/${USER}/projects/pRF_analysis/RetinoMaps/rest/mmp1_clusters"
OUTPUT_PATH="${BASE_PATH}/group/91k/rest/wta"
TSV_FILE="${OUTPUT_PATH}/group_wta_full_corr.tsv"
OUTPUT_FILE="${OUTPUT_PATH}/wta_full_corr_labels.txt"

# Define winner colors (seed number -> RGB + alpha)
declare -A WINNER_COLORS
WINNER_COLORS[1]="255 111 0 255"     # mPCS
WINNER_COLORS[2]="255 234 0 255"     # sPCS
WINNER_COLORS[3]="151 255 0 255"     # iPCS
WINNER_COLORS[4]="44 255 150 255"    # sIPS
WINNER_COLORS[5]="0 152 255 255"     # iIPS
WINNER_COLORS[6]="0 25 255 255"      # hMT+
WINNER_COLORS[7]="0 0 200 255"       # VO
WINNER_COLORS[8]="150 0 90 255"      # LO
WINNER_COLORS[9]="235 127 134 255"   # V3AB
WINNER_COLORS[10]="248 160 126 255"  # V3
WINNER_COLORS[11]="250 196 132 255"  # V2
WINNER_COLORS[12]="243 231 155 255"  # V1

# Define parcel mapping (in order as they appear in TSV - should match row order)
declare -a PARCELS=(
    "V1" "MST" "V2" "V3" "V4" "V8" 
    "FEF" "PEF" "55b" "V3A" "V7" "IPS1"
    "FFC" "V3B" "LO1" "LO2" "PIT" "MT"
    "7Pm" "24dv" "7AL" "SCEF" "6ma" "7Am"
    "7PL" "7PC" "LIPv" "VIP" "MIP" "6d"
    "6mp" "6v" "p32pr" "6r" "IFJa" "IFJp"
    "LIPd" "6a" "i6-8" "AIP" "PH" "IP2"
    "IP1" "IP0" "V6A" "VMV1" "VMV3" "V4t"
    "FST" "V3CD" "LO3" "VMV2" "VVC"
)

# Atlas keys - Right hemisphere
declare -A R_KEYS
R_KEYS["V1"]=1; R_KEYS["MST"]=2; R_KEYS["V2"]=4; R_KEYS["V3"]=5; R_KEYS["V4"]=6
R_KEYS["V8"]=7; R_KEYS["FEF"]=10; R_KEYS["PEF"]=11; R_KEYS["55b"]=12; R_KEYS["V3A"]=13
R_KEYS["V7"]=16; R_KEYS["IPS1"]=17; R_KEYS["FFC"]=18; R_KEYS["V3B"]=19; R_KEYS["LO1"]=20
R_KEYS["LO2"]=21; R_KEYS["PIT"]=22; R_KEYS["MT"]=23; R_KEYS["7Pm"]=29; R_KEYS["24dv"]=41
R_KEYS["7AL"]=42; R_KEYS["SCEF"]=43; R_KEYS["6ma"]=44; R_KEYS["7Am"]=45; R_KEYS["7PL"]=46
R_KEYS["7PC"]=47; R_KEYS["LIPv"]=48; R_KEYS["VIP"]=49; R_KEYS["MIP"]=50; R_KEYS["6d"]=54
R_KEYS["6mp"]=55; R_KEYS["6v"]=56; R_KEYS["p32pr"]=60; R_KEYS["6r"]=78; R_KEYS["IFJa"]=79
R_KEYS["IFJp"]=80; R_KEYS["LIPd"]=95; R_KEYS["6a"]=96; R_KEYS["i6-8"]=97; R_KEYS["AIP"]=117
R_KEYS["PH"]=138; R_KEYS["IP2"]=144; R_KEYS["IP1"]=145; R_KEYS["IP0"]=146; R_KEYS["V6A"]=152
R_KEYS["VMV1"]=153; R_KEYS["VMV3"]=154; R_KEYS["V4t"]=156; R_KEYS["FST"]=157; R_KEYS["V3CD"]=158
R_KEYS["LO3"]=159; R_KEYS["VMV2"]=160; R_KEYS["VVC"]=163

# Atlas keys - Left hemisphere
declare -A L_KEYS
L_KEYS["V1"]=181; L_KEYS["MST"]=182; L_KEYS["V2"]=184; L_KEYS["V3"]=185; L_KEYS["V4"]=186
L_KEYS["V8"]=187; L_KEYS["FEF"]=190; L_KEYS["PEF"]=191; L_KEYS["55b"]=192; L_KEYS["V3A"]=193
L_KEYS["V7"]=196; L_KEYS["IPS1"]=197; L_KEYS["FFC"]=198; L_KEYS["V3B"]=199; L_KEYS["LO1"]=200
L_KEYS["LO2"]=201; L_KEYS["PIT"]=202; L_KEYS["MT"]=203; L_KEYS["7Pm"]=209; L_KEYS["24dv"]=221
L_KEYS["7AL"]=222; L_KEYS["SCEF"]=223; L_KEYS["6ma"]=224; L_KEYS["7Am"]=225; L_KEYS["7PL"]=226
L_KEYS["7PC"]=227; L_KEYS["LIPv"]=228; L_KEYS["VIP"]=229; L_KEYS["MIP"]=230; L_KEYS["6d"]=234
L_KEYS["6mp"]=235; L_KEYS["6v"]=236; L_KEYS["p32pr"]=240; L_KEYS["6r"]=258; L_KEYS["IFJa"]=259
L_KEYS["IFJp"]=260; L_KEYS["LIPd"]=275; L_KEYS["6a"]=276; L_KEYS["i6-8"]=277; L_KEYS["AIP"]=297
L_KEYS["PH"]=318; L_KEYS["IP2"]=324; L_KEYS["IP1"]=325; L_KEYS["IP0"]=326; L_KEYS["V6A"]=332
L_KEYS["VMV1"]=333; L_KEYS["VMV3"]=334; L_KEYS["V4t"]=336; L_KEYS["FST"]=337; L_KEYS["V3CD"]=338
L_KEYS["LO3"]=339; L_KEYS["VMV2"]=340; L_KEYS["VVC"]=343

echo "Reading TSV file: ${TSV_FILE}"

# Read winners from TSV file (one value per row)
# Each row corresponds to a parcel in the order defined in PARCELS array
mapfile -t WINNERS < "${TSV_FILE}"

echo "Found ${#WINNERS[@]} winner values"
echo "Expected ${#PARCELS[@]} parcels"

# Clear output file
> "${OUTPUT_FILE}"

# Generate labels for each parcel
parcel_count=0
for i in "${!PARCELS[@]}"; do
    parcel="${PARCELS[$i]}"
    winner="${WINNERS[$i]}"
    
    # Skip if no winner
    if [[ -z "$winner" ]] || [[ "$winner" == "nan" ]]; then
        echo "Warning: No winner for ${parcel} (row $((i+1)))"
        continue
    fi
    
    # Remove any whitespace and round winner to integer
    winner=$(echo "$winner" | tr -d '[:space:]')
    winner=$(printf "%.0f" "$winner" 2>/dev/null)
    
    # Check if winner is valid (1-12)
    if [[ ! "$winner" =~ ^[0-9]+$ ]] || [[ "$winner" -lt 1 ]] || [[ "$winner" -gt 12 ]]; then
        echo "Warning: Invalid winner value '$winner' for ${parcel} (row $((i+1)))"
        continue
    fi
    
    # Get the color for this winner
    color="${WINNER_COLORS[$winner]}"
    
    # Get keys for this parcel
    r_key="${R_KEYS[$parcel]}"
    l_key="${L_KEYS[$parcel]}"
    
    # Write right hemisphere entry
    if [[ -n "$r_key" ]]; then
        echo "R_${parcel}_ROI" >> "${OUTPUT_FILE}"
        echo "$r_key $color" >> "${OUTPUT_FILE}"
        ((parcel_count++))
    fi
    
    # Write left hemisphere entry
    if [[ -n "$l_key" ]]; then
        echo "L_${parcel}_ROI" >> "${OUTPUT_FILE}"
        echo "$l_key $color" >> "${OUTPUT_FILE}"
        ((parcel_count++))
    fi
done

echo ""
echo "Label file created: ${OUTPUT_FILE}"
echo "Total label entries: $parcel_count"
echo ""
echo "First 20 lines:"
head -n 20 "${OUTPUT_FILE}"
echo ""

# Import labels into CIFTI file
echo "Importing labels into CIFTI format..."
wb_command -cifti-label-import \
    "${ATLAS_DIR}/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_dseg.dlabel.nii" \
    "${OUTPUT_FILE}" \
    "${OUTPUT_PATH}/group_wta_full_corr.dlabel.nii" \
    -discard-others

echo ""
echo "CIFTI label file created: ${OUTPUT_PATH}/group_wta_full_corr.dlabel.nii"
echo "Done!"
