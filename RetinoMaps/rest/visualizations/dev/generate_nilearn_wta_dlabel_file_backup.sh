#!/bin/bash
#####################################################
# Generate label file with winner-take-all colors from Nilearn CSV
# Written by Marco Bedini (marco.bedini@univ-amu.fr)
# Example use for main output of interest
# $ ./generate_nilearn_wta_dlabel_file.sh fisher-z by_hemi
#####################################################

# Parse command-line arguments
CORR_TYPE="${1:-full_corr}"      # full_corr or fisher-z
HEMI_TYPE="${2:-bilateral}"      # bilateral or by_hemi
PARC_MODE="${3:-default}"        # default or legacy

# Validate inputs
if [[ "$CORR_TYPE" != "full_corr" ]] && [[ "$CORR_TYPE" != "fisher-z" ]]; then
    echo "Error: CORR_TYPE must be 'full_corr' or 'fisher-z'"
    echo "Usage: $0 [full_corr|fisher-z] [bilateral|by_hemi] [default|legacy]"
    exit 1
fi

if [[ "$HEMI_TYPE" != "bilateral" ]] && [[ "$HEMI_TYPE" != "by_hemi" ]]; then
    echo "Error: HEMI_TYPE must be 'bilateral' or 'by_hemi'"
    echo "Usage: $0 [full_corr|fisher-z] [bilateral|by_hemi] [default|legacy]"
    exit 1
fi

if [[ "$PARC_MODE" != "default" ]] && [[ "$PARC_MODE" != "legacy" ]]; then
    echo "Error: PARC_MODE must be 'default' or 'legacy'"
    echo "Usage: $0 [full_corr|fisher-z] [bilateral|by_hemi] [default|legacy]"
    exit 1
fi

echo "============================================"
echo "WTA LABEL GENERATION (Nilearn)"
echo "============================================"
echo "Correlation type:  ${CORR_TYPE}"
echo "Hemisphere type:   ${HEMI_TYPE}"
echo "Parcellation mode: ${PARC_MODE}"
echo "============================================"
echo ""

# Construct file suffix based on options
FILE_SUFFIX="${CORR_TYPE}"
[[ "$HEMI_TYPE" == "by_hemi" ]] && FILE_SUFFIX="${FILE_SUFFIX}_by_hemi"
[[ "$PARC_MODE" == "legacy"  ]] && FILE_SUFFIX="${FILE_SUFFIX}_legacy"

BASE_PATH="/scratch/mszinte/data/RetinoMaps/derivatives/pp_data"
ATLAS_DIR="${BASE_PATH}/atlas/mmp1_clusters"
OUTPUT_PATH="${BASE_PATH}/group/91k/rest/wta/nilearn"
OUTPUT_FILE="${ATLAS_DIR}/wta_mmp1_labels_nilearn.txt"
OUTPUT_FILE_PCEF_IPEF="${ATLAS_DIR}/wta_mmp1_labels_nilearn_pcef-ipef.txt"

# Winner colors (seed number -> RGB + alpha)
# Seed numbers correspond to ROI order: 1=mPCS ... 12=V1
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

# Cluster membership (for self-win validation)
declare -A CLUSTER_PARCELS
CLUSTER_PARCELS["mPCS"]="SCEF p32pr 24dv"
CLUSTER_PARCELS["sPCS"]="FEF i6-8 6a 6d 6mp 6ma"
CLUSTER_PARCELS["iPCS"]="PEF IFJp 6v 6r IFJa 55b"
CLUSTER_PARCELS["sIPS"]="VIP LIPv LIPd IP2 7PC AIP 7AL 7Am 7Pm"
CLUSTER_PARCELS["iIPS"]="IP0 IPS1 V7 MIP IP1 V6A 7PL"
CLUSTER_PARCELS["hMT+"]="V4t MST MT FST"
CLUSTER_PARCELS["VO"]="V8 PIT PH FFC VMV1 VMV2 VMV3 VVC"
CLUSTER_PARCELS["LO"]="LO1 LO2 LO3"
CLUSTER_PARCELS["V3AB"]="V3CD V3A V3B"
CLUSTER_PARCELS["V3"]="V3 V4"
CLUSTER_PARCELS["V2"]="V2"
CLUSTER_PARCELS["V1"]="V1"

# Map cluster names to winner seed numbers
declare -A CLUSTER_TO_WINNER
CLUSTER_TO_WINNER["mPCS"]=1
CLUSTER_TO_WINNER["sPCS"]=2
CLUSTER_TO_WINNER["iPCS"]=3
CLUSTER_TO_WINNER["sIPS"]=4
CLUSTER_TO_WINNER["iIPS"]=5
CLUSTER_TO_WINNER["hMT+"]=6
CLUSTER_TO_WINNER["VO"]=7
CLUSTER_TO_WINNER["LO"]=8
CLUSTER_TO_WINNER["V3AB"]=9
CLUSTER_TO_WINNER["V3"]=10
CLUSTER_TO_WINNER["V2"]=11
CLUSTER_TO_WINNER["V1"]=12

# Parcel order expected by write_hemi_labels (must match atlas key order)
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

# Column order as they appear left-to-right in the Nilearn CSV header
# Verified against the image pasted in the manuscript notes
declare -a CSV_COL_ORDER=(
    "SCEF" "p32pr" "24dv" "FEF" "i6-8" "6a" "6d" "6mp" "6ma" "PEF"
    "IFJp" "6v" "6r" "IFJa" "55b" "VIP" "LIPv" "LIPd" "IP2" "7PC"
    "AIP" "7AL" "7Am" "7Pm" "IP0" "IPS1" "V7" "MIP" "IP1" "V6A"
    "7PL" "V4t" "MST" "MT" "FST" "V8" "PIT" "PH" "FFC" "VMV1"
    "VMV2" "VMV3" "VVC" "LO1" "LO2" "LO3" "V3A" "V3B" "V3CD" "V3"
    "V4" "V2" "V1"
)

# Pre-compute reindex map: PARCEL_TO_CSV_COL[i] = column index in CSV for PARCELS[i]
# Built at startup so the CSV row extraction is a simple indexed lookup
declare -A _CSV_COL_IDX   # parcel name -> 0-based CSV column index
for _ci in "${!CSV_COL_ORDER[@]}"; do
    _CSV_COL_IDX["${CSV_COL_ORDER[$_ci]}"]="$_ci"
done

declare -a PARCEL_TO_CSV_COL=()
for _pi in "${!PARCELS[@]}"; do
    _parcel="${PARCELS[$_pi]}"
    if [[ -z "${_CSV_COL_IDX[$_parcel]+x}" ]]; then
        echo "FATAL: PARCELS entry '${_parcel}' not found in CSV_COL_ORDER. Check column definitions." >&2
        exit 1
    fi
    PARCEL_TO_CSV_COL[$_pi]="${_CSV_COL_IDX[$_parcel]}"
done

# Atlas keys - Right hemisphere
declare -A R_KEYS
R_KEYS["V1"]=1;   R_KEYS["MST"]=2;   R_KEYS["V2"]=4;    R_KEYS["V3"]=5;    R_KEYS["V4"]=6
R_KEYS["V8"]=7;   R_KEYS["FEF"]=10;  R_KEYS["PEF"]=11;  R_KEYS["55b"]=12;  R_KEYS["V3A"]=13
R_KEYS["V7"]=16;  R_KEYS["IPS1"]=17; R_KEYS["FFC"]=18;  R_KEYS["V3B"]=19;  R_KEYS["LO1"]=20
R_KEYS["LO2"]=21; R_KEYS["PIT"]=22;  R_KEYS["MT"]=23;   R_KEYS["7Pm"]=29;  R_KEYS["24dv"]=41
R_KEYS["7AL"]=42; R_KEYS["SCEF"]=43; R_KEYS["6ma"]=44;  R_KEYS["7Am"]=45;  R_KEYS["7PL"]=46
R_KEYS["7PC"]=47; R_KEYS["LIPv"]=48; R_KEYS["VIP"]=49;  R_KEYS["MIP"]=50;  R_KEYS["6d"]=54
R_KEYS["6mp"]=55; R_KEYS["6v"]=56;   R_KEYS["p32pr"]=60; R_KEYS["6r"]=78;  R_KEYS["IFJa"]=79
R_KEYS["IFJp"]=80; R_KEYS["LIPd"]=95; R_KEYS["6a"]=96;  R_KEYS["i6-8"]=97; R_KEYS["AIP"]=117
R_KEYS["PH"]=138; R_KEYS["IP2"]=144; R_KEYS["IP1"]=145; R_KEYS["IP0"]=146; R_KEYS["V6A"]=152
R_KEYS["VMV1"]=153; R_KEYS["VMV3"]=154; R_KEYS["V4t"]=156; R_KEYS["FST"]=157; R_KEYS["V3CD"]=158
R_KEYS["LO3"]=159; R_KEYS["VMV2"]=160; R_KEYS["VVC"]=163

# Atlas keys - Left hemisphere
declare -A L_KEYS
L_KEYS["V1"]=181;  L_KEYS["MST"]=182; L_KEYS["V2"]=184;  L_KEYS["V3"]=185;  L_KEYS["V4"]=186
L_KEYS["V8"]=187;  L_KEYS["FEF"]=190; L_KEYS["PEF"]=191; L_KEYS["55b"]=192; L_KEYS["V3A"]=193
L_KEYS["V7"]=196;  L_KEYS["IPS1"]=197; L_KEYS["FFC"]=198; L_KEYS["V3B"]=199; L_KEYS["LO1"]=200
L_KEYS["LO2"]=201; L_KEYS["PIT"]=202; L_KEYS["MT"]=203;  L_KEYS["7Pm"]=209; L_KEYS["24dv"]=221
L_KEYS["7AL"]=222; L_KEYS["SCEF"]=223; L_KEYS["6ma"]=224; L_KEYS["7Am"]=225; L_KEYS["7PL"]=226
L_KEYS["7PC"]=227; L_KEYS["LIPv"]=228; L_KEYS["VIP"]=229; L_KEYS["MIP"]=230; L_KEYS["6d"]=234
L_KEYS["6mp"]=235; L_KEYS["6v"]=236;  L_KEYS["p32pr"]=240; L_KEYS["6r"]=258; L_KEYS["IFJa"]=259
L_KEYS["IFJp"]=260; L_KEYS["LIPd"]=275; L_KEYS["6a"]=276; L_KEYS["i6-8"]=277; L_KEYS["AIP"]=297
L_KEYS["PH"]=318;  L_KEYS["IP2"]=324; L_KEYS["IP1"]=325; L_KEYS["IP0"]=326; L_KEYS["V6A"]=332
L_KEYS["VMV1"]=333; L_KEYS["VMV3"]=334; L_KEYS["V4t"]=336; L_KEYS["FST"]=337; L_KEYS["V3CD"]=338
L_KEYS["LO3"]=339; L_KEYS["VMV2"]=340; L_KEYS["VVC"]=343

# ============================================
# HELPER: validate winner value is integer in [1,12]
# Returns 0 (true) if valid
# ============================================
validate_winner() {
    local winner="$1"
    [[ "$winner" =~ ^[0-9]+$ ]] && [[ "$winner" -ge 1 ]] && [[ "$winner" -le 12 ]]
}

# ============================================
# HELPER: check if parcel won its own cluster (self-win warning)
# Returns 0 if self-win detected
# ============================================
check_self_win() {
    local parcel="$1"
    local winner="$2"
    local cluster expected
    for cluster in mPCS sPCS iPCS sIPS iIPS "hMT+" VO LO V3AB V3 V2 V1; do
        if echo "${CLUSTER_PARCELS[$cluster]}" | grep -wq "$parcel"; then
            expected="${CLUSTER_TO_WINNER[$cluster]}"
            [[ "$winner" -eq "$expected" ]] && return 0
            return 1
        fi
    done
    return 1
}

# ============================================
# HELPER: write label entries for one hemisphere side
# Uses eval-based indirection instead of local -n for bash 3.2 compatibility
# Arguments:
#   $1  output_file
#   $2  hemi_char ("R" or "L")
#   $3  name of keys associative array ("R_KEYS" or "L_KEYS")
#   $4  name of winners indexed array
#   $5  name of integer parcel counter variable (incremented in place)
#   $6  name of integer warning counter variable (incremented in place)
# ============================================
write_hemi_labels() {
    local output_file="$1"
    local hemi_char="$2"
    local keys_arr_name="$3"
    local winners_arr_name="$4"
    local count_var="$5"
    local warn_var="$6"

    local i parcel winner color key

    for i in "${!PARCELS[@]}"; do
        parcel="${PARCELS[$i]}"

        eval "winner=\"\${${winners_arr_name}[$i]:-}\""

        [[ -z "$winner" || "$winner" == "nan" ]] && continue

        winner=$(echo "$winner" | tr -d '[:space:]')
        winner=$(printf "%.0f" "$winner" 2>/dev/null)

        if ! validate_winner "$winner"; then
            echo "  WARNING: Parcel '${hemi_char}_${parcel}' has invalid winner '$winner' — skipping" >&2
            continue
        fi

        if check_self_win "$parcel" "$winner"; then
            echo "  WARNING: Parcel '${hemi_char}_${parcel}' won its own cluster (winner=$winner)!" >&2
            eval "${warn_var}=$(( ${!warn_var} + 1 ))"
        fi

        color="${WINNER_COLORS[$winner]}"
        eval "key=\"\${${keys_arr_name}[$parcel]:-}\""

        if [[ -z "$key" ]]; then
            echo "  WARNING: No atlas key for ${hemi_char}_${parcel} — skipping" >&2
            continue
        fi

        printf '%s\n%s %s\n' "${hemi_char}_${parcel}_ROI" "${key} ${color}" >> "$output_file"
        eval "${count_var}=$(( ${!count_var} + 1 ))"
    done
}

# ============================================
# HELPER: build label file from rh and lh winner arrays
# Arguments:
#   $1  output label file path
#   $2  name of rh winners array
#   $3  name of lh winners array
#   $4  name of parcel counter variable
#   $5  name of warning counter variable
# ============================================
build_label_file() {
    local output_file="$1"
    local rh_arr_name="$2"
    local lh_arr_name="$3"
    local count_var="$4"
    local warn_var="$5"

    > "$output_file"

    write_hemi_labels "$output_file" "R" "R_KEYS" "$rh_arr_name" "$count_var" "$warn_var"
    write_hemi_labels "$output_file" "L" "L_KEYS" "$lh_arr_name" "$count_var" "$warn_var"
}

# ============================================
# HELPER: write label entries for one hemisphere side — PCEF/IPEF filter
# Identical to write_hemi_labels but only emits parcels whose winner
# is in [1,5] (mPCS=1, sPCS=2, iPCS=3, sIPS=4, iIPS=5).
# Arguments: same as write_hemi_labels
# ============================================
write_hemi_labels_pcef_ipef() {
    local output_file="$1"
    local hemi_char="$2"
    local keys_arr_name="$3"
    local winners_arr_name="$4"
    local count_var="$5"
    local warn_var="$6"

    local i parcel winner color key

    for i in "${!PARCELS[@]}"; do
        parcel="${PARCELS[$i]}"

        eval "winner=\"\${${winners_arr_name}[$i]:-}\""

        [[ -z "$winner" || "$winner" == "nan" ]] && continue

        winner=$(echo "$winner" | tr -d '[:space:]')
        winner=$(printf "%.0f" "$winner" 2>/dev/null)

        if ! validate_winner "$winner"; then
            continue   # already warned in the full-map pass
        fi

        # PCEF/IPEF filter: keep only winners 1–5
        [[ "$winner" -lt 1 || "$winner" -gt 5 ]] && continue

        if check_self_win "$parcel" "$winner"; then
            eval "${warn_var}=$(( ${!warn_var} + 1 ))"
        fi

        color="${WINNER_COLORS[$winner]}"
        eval "key=\"\${${keys_arr_name}[$parcel]:-}\""

        [[ -z "$key" ]] && continue   # already warned in the full-map pass

        printf '%s\n%s %s\n' "${hemi_char}_${parcel}_ROI" "${key} ${color}" >> "$output_file"
        eval "${count_var}=$(( ${!count_var} + 1 ))"
    done
}

# ============================================
# HELPER: build PCEF/IPEF-filtered label file
# Arguments: same as build_label_file
# ============================================
build_label_file_pcef_ipef() {
    local output_file="$1"
    local rh_arr_name="$2"
    local lh_arr_name="$3"
    local count_var="$4"
    local warn_var="$5"

    > "$output_file"

    write_hemi_labels_pcef_ipef "$output_file" "R" "R_KEYS" "$rh_arr_name" "$count_var" "$warn_var"
    write_hemi_labels_pcef_ipef "$output_file" "L" "L_KEYS" "$lh_arr_name" "$count_var" "$warn_var"
}

# ============================================
# HELPER: read one row from a Nilearn CSV into a winners array
# The CSV has a header row then one row per subject; the last data row
# is the group result. Columns are in CSV_COL_ORDER; values are
# reindexed into PARCELS order using PARCEL_TO_CSV_COL.
#
# Arguments:
#   $1  path to CSV file
#   $2  row to extract: "group" (last data row) or a subject ID e.g. "sub-01"
#       Subject rows are matched on the first column (row label).
#   $3  name of winners array to populate
# Returns 1 on error
# ============================================
read_csv_row_into_array() {
    local csv_path="$1"
    local row_id="$2"
    local arr_name="$3"

    if [[ ! -f "$csv_path" ]]; then
        echo "Error: CSV file not found: $csv_path" >&2
        return 1
    fi

    # Load all lines (skip header row 0)
    local -a _all_lines=()
    mapfile -t _all_lines < "$csv_path"
    local n_lines="${#_all_lines[@]}"

    # Validate: need at least header + 1 data row
    if [[ "$n_lines" -lt 2 ]]; then
        echo "Error: CSV has fewer than 2 lines (no data rows): $csv_path" >&2
        return 1
    fi

    # Locate the target data row
    local _target_line=""
    if [[ "$row_id" == "group" ]]; then
        # Last line is always the group result
        _target_line="${_all_lines[$(( n_lines - 1 ))]}"
    else
        # Match on first field (subject ID)
        local _line
        for _line in "${_all_lines[@]:1}"; do   # skip header
            local _first_field
            _first_field=$(echo "$_line" | cut -d',' -f1)
            if [[ "$_first_field" == "$row_id" ]]; then
                _target_line="$_line"
                break
            fi
        done
        if [[ -z "$_target_line" ]]; then
            echo "Error: Row '${row_id}' not found in ${csv_path}" >&2
            return 1
        fi
    fi

    # Split the target line into fields (comma-separated)
    local -a _fields=()
    IFS=',' read -r -a _fields <<< "$_target_line"

    # Validate field count: expect 1 label column + 53 parcel columns
    local expected_fields=$(( 1 + ${#PARCELS[@]} ))
    if [[ "${#_fields[@]}" -ne "$expected_fields" ]]; then
        echo "Error: Expected $expected_fields fields (1 label + ${#PARCELS[@]} parcels) but got ${#_fields[@]} in row '${row_id}': $csv_path" >&2
        return 1
    fi

    # Reindex: map each PARCELS[i] position to its CSV column
    # CSV data starts at field index 1 (field 0 is the row label)
    eval "${arr_name}=()"
    local _pi
    for _pi in "${!PARCELS[@]}"; do
        local _csv_col="${PARCEL_TO_CSV_COL[$_pi]}"
        local _field_idx=$(( _csv_col + 1 ))   # +1 to skip row-label field
        eval "${arr_name}[$_pi]=\"\${_fields[$_field_idx]}\""
    done
}

# ============================================
# HELPER: read winners for group or a subject from Nilearn CSVs
#
# For bilateral: reads a single CSV, copies winners to both rh and lh arrays.
# For by_hemi:   reads two CSVs (one per hemisphere).
#
# CSV naming convention (current):
#   bilateral: winning_seeds_by_subject_partial_corr_bilateral.csv
#   by_hemi:   winning_seeds_by_subject_partial_corr_lh.csv / _rh.csv
#
# Arguments:
#   $1  row_id: "group" or subject ID (e.g. "sub-01")
#   $2  name of rh winners array to populate
#   $3  name of lh winners array to populate
# Returns 1 on error
# ============================================
read_all_winners() {
    local row_id="$1"
    local rh_arr_name="$2"
    local lh_arr_name="$3"

    local csv_path hemi idx

    if [[ "$HEMI_TYPE" == "by_hemi" ]]; then
        for hemi in rh lh; do
            csv_path="${OUTPUT_PATH}/winning_seeds_by_subject_partial_corr_${hemi}.csv"

            if [[ ! -f "$csv_path" ]]; then
                echo "Error: CSV not found: $csv_path" >&2
                return 1
            fi

            if [[ "$hemi" == "rh" ]]; then
                read_csv_row_into_array "$csv_path" "$row_id" "$rh_arr_name" || return 1
            else
                read_csv_row_into_array "$csv_path" "$row_id" "$lh_arr_name" || return 1
            fi
        done

    else
        # Bilateral: single CSV, same winners for both hemispheres
        csv_path="${OUTPUT_PATH}/winning_seeds_by_subject_partial_corr_bilateral.csv"

        if [[ ! -f "$csv_path" ]]; then
            echo "Error: CSV not found: $csv_path" >&2
            return 1
        fi

        read_csv_row_into_array "$csv_path" "$row_id" "$rh_arr_name" || return 1

        # Copy rh into lh
        local rh_len
        eval "rh_len=\${#${rh_arr_name}[@]}"
        eval "${lh_arr_name}=()"
        for idx in $(seq 0 $(( rh_len - 1 ))); do
            eval "${lh_arr_name}[$idx]=\"\${${rh_arr_name}[$idx]}\""
        done
    fi
}

# ============================================
# GROUP-LEVEL PROCESSING
# ============================================

echo "============================================"
echo "GROUP-LEVEL PROCESSING"
echo "============================================"

declare -a group_rh_winners=()
declare -a group_lh_winners=()

# Group result is always the last row of the CSV
if ! read_all_winners "group" "group_rh_winners" "group_lh_winners"; then
    echo "Error: Failed to read group winners. Aborting." >&2
    exit 1
fi

echo "RH winners loaded: ${#group_rh_winners[@]}"
echo "LH winners loaded: ${#group_lh_winners[@]}"
echo "Expected: ${#PARCELS[@]} parcels per hemisphere"
echo ""

group_parcel_count=0
group_warning_count=0

build_label_file "$OUTPUT_FILE" \
    "group_rh_winners" "group_lh_winners" \
    "group_parcel_count" "group_warning_count"

echo "============================================"
echo "GROUP-LEVEL SUMMARY"
echo "============================================"
echo "Label file:   ${OUTPUT_FILE}"
echo "Total labels: $group_parcel_count"
echo "Warnings:     $group_warning_count"
echo ""

wb_command -cifti-label-import \
    "${ATLAS_DIR}/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_dseg.dlabel.nii" \
    "${OUTPUT_FILE}" \
    "${OUTPUT_PATH}/group_wta_${FILE_SUFFIX}.dlabel.nii" \
    -discard-others

echo "CIFTI file: ${OUTPUT_PATH}/group_wta_${FILE_SUFFIX}.dlabel.nii"
echo ""

# --- PCEF/IPEF filtered map (winners 1–5 only) ---
group_pcef_ipef_count=0
group_pcef_ipef_warn=0

build_label_file_pcef_ipef "$OUTPUT_FILE_PCEF_IPEF" \
    "group_rh_winners" "group_lh_winners" \
    "group_pcef_ipef_count" "group_pcef_ipef_warn"

wb_command -cifti-label-import \
    "${ATLAS_DIR}/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_dseg.dlabel.nii" \
    "${OUTPUT_FILE_PCEF_IPEF}" \
    "${OUTPUT_PATH}/group_wta_${FILE_SUFFIX}_pcef-ipef.dlabel.nii" \
    -discard-others

echo "PCEF/IPEF labels: $group_pcef_ipef_count"
echo "CIFTI file: ${OUTPUT_PATH}/group_wta_${FILE_SUFFIX}_pcef-ipef.dlabel.nii"
echo ""

# ============================================
# SUBJECT-LEVEL PROCESSING
# ============================================

echo "============================================"
echo "SUBJECT-LEVEL PROCESSING"
echo "============================================"
echo ""

subjects=("01" "02" "03" "04" "05" "06" "07" "08" "09" "11" "12" "13" "14" "17" "20" "21" "22" "23" "24" "25")

for sub in "${subjects[@]}"; do
    echo "Processing sub-${sub}..."

    declare -a sub_rh_winners=()
    declare -a sub_lh_winners=()

    # Row label in the CSV matches the subject folder name
    if ! read_all_winners "sub-${sub}" "sub_rh_winners" "sub_lh_winners"; then
        echo "  ERROR: Skipping sub-${sub} (failed to read winners)" >&2
        unset sub_rh_winners sub_lh_winners
        continue
    fi

    SUB_LABELS="${OUTPUT_PATH}/sub-${sub}_wta_${FILE_SUFFIX}_labels.txt"
    sub_parcel_count=0
    sub_warning_count=0

    build_label_file "$SUB_LABELS" \
        "sub_rh_winners" "sub_lh_winners" \
        "sub_parcel_count" "sub_warning_count"

    wb_command -cifti-label-import \
        "${ATLAS_DIR}/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_dseg.dlabel.nii" \
        "${SUB_LABELS}" \
        "${OUTPUT_PATH}/sub-${sub}_wta_${FILE_SUFFIX}.dlabel.nii" \
        -discard-others

    # --- PCEF/IPEF filtered map (winners 1–5 only) ---
    SUB_LABELS_PCEF_IPEF="${OUTPUT_PATH}/sub-${sub}_wta_${FILE_SUFFIX}_pcef-ipef_labels.txt"
    sub_pcef_ipef_count=0
    sub_pcef_ipef_warn=0

    build_label_file_pcef_ipef "$SUB_LABELS_PCEF_IPEF" \
        "sub_rh_winners" "sub_lh_winners" \
        "sub_pcef_ipef_count" "sub_pcef_ipef_warn"

    wb_command -cifti-label-import \
        "${ATLAS_DIR}/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_dseg.dlabel.nii" \
        "${SUB_LABELS_PCEF_IPEF}" \
        "${OUTPUT_PATH}/sub-${sub}_wta_${FILE_SUFFIX}_pcef-ipef.dlabel.nii" \
        -discard-others

    if [[ $sub_warning_count -gt 0 ]]; then
        echo "  ⚠️  Labels: $sub_parcel_count | Warnings: $sub_warning_count"
    else
        echo "  ✓ Labels: $sub_parcel_count"
    fi

    unset sub_rh_winners sub_lh_winners
done

chmod -Rf 771 /scratch/mszinte/data/RetinoMaps
chgrp -Rf 327 /scratch/mszinte/data/RetinoMaps


echo ""
echo "============================================"
echo "ALL PROCESSING COMPLETE"
echo "============================================"
echo "Group (full):         ${OUTPUT_PATH}/group_wta_${FILE_SUFFIX}.dlabel.nii"
echo "Group (PCEF/IPEF):    ${OUTPUT_PATH}/group_wta_${FILE_SUFFIX}_pcef-ipef.dlabel.nii"
echo "Subjects (full):      ${OUTPUT_PATH}/sub-*_wta_${FILE_SUFFIX}.dlabel.nii"
echo "Subjects (PCEF/IPEF): ${OUTPUT_PATH}/sub-*_wta_${FILE_SUFFIX}_pcef-ipef.dlabel.nii"
echo ""