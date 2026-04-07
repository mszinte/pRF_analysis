#!/bin/bash
#####################################################
# Generate label file with winner-take-all colors from TSV
# Written by Marco Bedini (marco.bedini@univ-amu.fr)
# Example use for main output of interest
# $ ./generate_workbench_wta_dlabel_file.sh fisher-z by_hemi
#####################################################

# Parse command-line arguments
CORR_TYPE="${1:-full_corr}"
HEMI_TYPE="${2:-bilateral}"
PARC_MODE="${3:-default}"

# Validate inputs
if [[ "$CORR_TYPE" != "full_corr" ]] && [[ "$CORR_TYPE" != "fisher-z" ]]; then
    echo "Error: CORR_TYPE must be 'full_corr' or 'fisher-z'"
    exit 1
fi

if [[ "$HEMI_TYPE" != "bilateral" ]] && [[ "$HEMI_TYPE" != "by_hemi" ]]; then
    echo "Error: HEMI_TYPE must be 'bilateral' or 'by_hemi'"
    exit 1
fi

# ============================================
# FILE SUFFIX
# ============================================
FILE_SUFFIX="${CORR_TYPE}"
[[ "$HEMI_TYPE" == "by_hemi" ]] && FILE_SUFFIX="${FILE_SUFFIX}_by_hemi"
[[ "$PARC_MODE" == "legacy" ]] && FILE_SUFFIX="${FILE_SUFFIX}_legacy"   ### ADDED

BASE_PATH="/scratch/mszinte/data/RetinoMaps/derivatives/pp_data"
ATLAS_DIR="${BASE_PATH}/atlas/mmp1_clusters"
OUTPUT_PATH="${BASE_PATH}/group/91k/rest/wta/workbench"
OUTPUT_FILE="${ATLAS_DIR}/wta_mmp1_labels.txt"

# ============================================
# COLORS
# ============================================
declare -A WINNER_COLORS
WINNER_COLORS[1]="255 111 0 255"
WINNER_COLORS[2]="255 234 0 255"
WINNER_COLORS[3]="151 255 0 255"
WINNER_COLORS[4]="44 255 150 255"
WINNER_COLORS[5]="0 152 255 255"
WINNER_COLORS[6]="0 25 255 255"
WINNER_COLORS[7]="0 0 200 255"
WINNER_COLORS[8]="150 0 90 255"
WINNER_COLORS[9]="235 127 134 255"
WINNER_COLORS[10]="248 160 126 255"
WINNER_COLORS[11]="250 196 132 255"
WINNER_COLORS[12]="243 231 155 255"

# ============================================
# PARCELS
# ============================================
declare -a PARCELS=( "V1" "MST" "V2" "V3" "V4" "V8" "FEF" "PEF" "55b" "V3A" "V7" "IPS1"
"FFC" "V3B" "LO1" "LO2" "PIT" "MT" "7Pm" "24dv" "7AL" "SCEF" "6ma" "7Am"
"7PL" "7PC" "LIPv" "VIP" "MIP" "6d" "6mp" "6v" "p32pr" "6r" "IFJa" "IFJp"
"LIPd" "6a" "i6-8" "AIP" "PH" "IP2" "IP1" "IP0" "V6A" "VMV1" "VMV3" "V4t"
"FST" "V3CD" "LO3" "VMV2" "VVC" )

# ============================================
# VALIDATION
# ============================================
validate_winner() {
    local winner="$1"
    [[ "$winner" =~ ^[0-9]+$ ]] && [[ "$winner" -ge 1 ]] && [[ "$winner" -le 12 ]]
}

# ============================================
# WRITE LABELS
# ============================================
write_hemi_labels() {
    local output_file="$1"
    local hemi_char="$2"
    local keys_arr_name="$3"
    local winners_arr_name="$4"
    local count_var="$5"
    local warn_var="$6"
    local filter_mode="$7"   ### ADDED

    local i parcel winner color key

    for i in "${!PARCELS[@]}"; do
        parcel="${PARCELS[$i]}"
        eval "winner=\"\${${winners_arr_name}[$i]:-}\""

        [[ -z "$winner" || "$winner" == "nan" ]] && continue

        # Sanitize
        winner=$(echo "$winner" | tr -d '[:space:]')

        # Skip non-numeric entries (fixes "55b" crash)
        if ! [[ "$winner" =~ ^[0-9]+([.][0-9]+)?$ ]]; then
            continue
        fi

        # Convert to integer
        winner=$(printf "%.0f" "$winner" 2>/dev/null)

        # Validate range
        if ! validate_winner "$winner"; then
            continue
        fi

        ### ADDED: filter winners 1–5
        if [[ "$filter_mode" == "pcef_ipef" ]]; then
            if [[ "$winner" -lt 1 || "$winner" -gt 5 ]]; then
                continue
            fi
        fi

        color="${WINNER_COLORS[$winner]}"
        eval "key=\"\${${keys_arr_name}[$parcel]:-}\""
        [[ -z "$key" ]] && continue

        printf '%s\n%s %s\n' "${hemi_char}_${parcel}_ROI" "${key} ${color}" >> "$output_file"
        eval "${count_var}=$(( ${!count_var} + 1 ))"
    done
}

# ============================================
# BUILD LABEL FILE
# ============================================
build_label_file() {
    local output_file="$1"
    local rh_arr_name="$2"
    local lh_arr_name="$3"
    local count_var="$4"
    local warn_var="$5"
    local filter_mode="$6"   ### ADDED

    > "$output_file"

    write_hemi_labels "$output_file" "R" "R_KEYS" "$rh_arr_name" "$count_var" "$warn_var" "$filter_mode"
    write_hemi_labels "$output_file" "L" "L_KEYS" "$lh_arr_name" "$count_var" "$warn_var" "$filter_mode"
}

# ============================================
# HELPER: read a TSV file into a named array variable
# Uses eval for bash 3.2 compatibility
# Arguments:
#   $1  path to TSV file
#   $2  name of array variable to populate
# Returns 1 on error
# ============================================
read_tsv_into_array() {
    local tsv_path="$1"
    local arr_name="$2"

    if [[ ! -f "$tsv_path" ]]; then
        echo "Error: TSV file not found: $tsv_path" >&2
        return 1
    fi

    # Read lines into a temp array then copy via eval
    local -a _tmp_arr=()
    mapfile -t _tmp_arr < "$tsv_path"

    # Serialize into named array in caller scope
    eval "${arr_name}=()"
    local idx
    for idx in "${!_tmp_arr[@]}"; do
        eval "${arr_name}[$idx]=\"\${_tmp_arr[$idx]}\""
    done
}

# ============================================
# READ WINNERS (patched for 106 rows)
# ============================================
read_all_winners() {
    local entity_prefix="$1"
    local rh_arr_name="$2"
    local lh_arr_name="$3"

    local tsv_path hemi

    if [[ "$HEMI_TYPE" == "by_hemi" ]]; then
        for hemi in rh lh; do
            tsv_path="${OUTPUT_PATH}/${entity_prefix}_${FILE_SUFFIX}_${hemi}.tsv"
            if [[ "$hemi" == "rh" ]]; then
                read_tsv_into_array "$tsv_path" "$rh_arr_name" || return 1
            else
                read_tsv_into_array "$tsv_path" "$lh_arr_name" || return 1
            fi
        done

        eval "${rh_arr_name}=(\"\${rh_tmp[@]}\")"
        eval "${lh_arr_name}=(\"\${lh_tmp[@]}\")"

        ### ADDED: split 106 rows
        local rh_len
        eval "rh_len=\${#${rh_arr_name}[@]}"

        if [[ "$rh_len" -eq 106 ]]; then
            local idx
            for idx in $(seq 0 52); do
                eval "tmp_rh[$idx]=\${${rh_arr_name}[$idx]}"
                eval "tmp_lh[$idx]=\${${lh_arr_name}[$((idx+53))]}"
            done
            eval "${rh_arr_name}=(\"\${tmp_rh[@]}\")"
            eval "${lh_arr_name}=(\"\${tmp_lh[@]}\")"
        fi
    fi
}

# ============================================
# GROUP PROCESSING
# ============================================
declare -a group_rh_winners=()
declare -a group_lh_winners=()

read_all_winners "group_wta" "group_rh_winners" "group_lh_winners"

group_parcel_count=0
group_warning_count=0

build_label_file "$OUTPUT_FILE" \
    "group_rh_winners" "group_lh_winners" \
    "group_parcel_count" "group_warning_count" \
    ""   ### MODIFIED

wb_command -cifti-label-import \
    "${ATLAS_DIR}/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_dseg.dlabel.nii" \
    "${OUTPUT_FILE}" \
    "${OUTPUT_PATH}/group_wta_${FILE_SUFFIX}.dlabel.nii" \
    -discard-others

### ADDED: filtered output
FILTERED_OUTPUT_FILE="${ATLAS_DIR}/wta_mmp1_labels_pcef-ipef.txt"

build_label_file "$FILTERED_OUTPUT_FILE" \
    "group_rh_winners" "group_lh_winners" \
    group_parcel_count group_warning_count \
    "pcef_ipef"

wb_command -cifti-label-import \
    "${ATLAS_DIR}/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_dseg.dlabel.nii" \
    "${FILTERED_OUTPUT_FILE}" \
    "${OUTPUT_PATH}/group_wta_${FILE_SUFFIX}_pcef-ipef.dlabel.nii" \
    -discard-others