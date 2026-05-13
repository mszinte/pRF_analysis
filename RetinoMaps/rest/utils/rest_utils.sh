#!/bin/bash
# rest_utils.sh
# -----------------------------------------------------------------------------
# Shared constants and helper functions for WTA dlabel generation.
#
# Location: RetinoMaps/rest/utils/rest_utils.sh
#
# Source from any script in rest/stats/ (or any sibling of utils/):
#
#   SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
#   source "${SCRIPT_DIR}/../utils/rest_utils.sh"
#
# Provides:
#   Constants : WINNER_COLORS, CLUSTER_PARCELS, CLUSTER_TO_WINNER,
#               R_KEYS, L_KEYS, SUBJECTS
#   Functions : validate_winner, check_self_win,
#               read_csv_row_into_array,
#               write_hemi_labels, write_hemi_labels_pcef_ipef,
#               build_label_file, build_label_file_pcef_ipef
#
# Not provided (script-specific):
#   Paths (BASE_PATH, INPUT_PATH, OUTPUT_PATH …)
#   CSV filename construction  (read_all_winners)
#   process_variant()
# -----------------------------------------------------------------------------
# Written by Marco Bedini (marco.bedini@univ-amu.fr)
# -----------------------------------------------------------------------------

# ============================================================
# Winner colors  (seed number → "R G B A")
# 1=mPCS  2=sPCS  3=iPCS  4=sIPS  5=iIPS  6=hMT+
# 7=VO    8=LO    9=V3AB  10=V3   11=V2   12=V1
# ============================================================
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

# ============================================================
# Cluster membership  (for self-win validation)
# ============================================================
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

# ============================================================
# Atlas keys
# Used only to look up the integer label value for wb_command.
# NOT used to determine column order in CSVs (that is name-based).
# ============================================================
declare -A R_KEYS
R_KEYS["V1"]=1;    R_KEYS["MST"]=2;   R_KEYS["V2"]=4;    R_KEYS["V3"]=5;    R_KEYS["V4"]=6
R_KEYS["V8"]=7;    R_KEYS["FEF"]=10;  R_KEYS["PEF"]=11;  R_KEYS["55b"]=12;  R_KEYS["V3A"]=13
R_KEYS["V7"]=16;   R_KEYS["IPS1"]=17; R_KEYS["FFC"]=18;  R_KEYS["V3B"]=19;  R_KEYS["LO1"]=20
R_KEYS["LO2"]=21;  R_KEYS["PIT"]=22;  R_KEYS["MT"]=23;   R_KEYS["7Pm"]=29;  R_KEYS["24dv"]=41
R_KEYS["7AL"]=42;  R_KEYS["SCEF"]=43; R_KEYS["6ma"]=44;  R_KEYS["7Am"]=45;  R_KEYS["7PL"]=46
R_KEYS["7PC"]=47;  R_KEYS["LIPv"]=48; R_KEYS["VIP"]=49;  R_KEYS["MIP"]=50;  R_KEYS["6d"]=54
R_KEYS["6mp"]=55;  R_KEYS["6v"]=56;   R_KEYS["p32pr"]=60; R_KEYS["6r"]=78;  R_KEYS["IFJa"]=79
R_KEYS["IFJp"]=80; R_KEYS["LIPd"]=95; R_KEYS["6a"]=96;   R_KEYS["i6-8"]=97; R_KEYS["AIP"]=117
R_KEYS["PH"]=138;  R_KEYS["IP2"]=144; R_KEYS["IP1"]=145; R_KEYS["IP0"]=146; R_KEYS["V6A"]=152
R_KEYS["VMV1"]=153; R_KEYS["VMV3"]=154; R_KEYS["V4t"]=156; R_KEYS["FST"]=157; R_KEYS["V3CD"]=158
R_KEYS["LO3"]=159; R_KEYS["VMV2"]=160; R_KEYS["VVC"]=163

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

# ============================================================
# Canonical subject list  (shared across all scripts)
# ============================================================
SUBJECTS=(
    "sub-01" "sub-02" "sub-03" "sub-04" "sub-05"
    "sub-06" "sub-07" "sub-08" "sub-09" "sub-11"
    "sub-12" "sub-13" "sub-14" "sub-17" "sub-20"
    "sub-21" "sub-22" "sub-23" "sub-24" "sub-25"
)

# ============================================================
# validate_winner
# Returns 0 (true) if $1 is an integer in [1, 12].
# ============================================================
validate_winner() {
    local winner="$1"
    [[ "$winner" =~ ^[0-9]+$ ]] && [[ "$winner" -ge 1 ]] && [[ "$winner" -le 12 ]]
}

# ============================================================
# check_self_win
# Returns 0 if parcel $1 was won by its own cluster seed (winner $2).
# ============================================================
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

# ============================================================
# read_csv_row_into_array
#
# Read one named row from a WTA CSV into an associative array keyed
# by parcel name:  arr["V1"]="3"  arr["MT"]="6"  etc.
#
# CSV format (produced by both wta_full_corr_by_subject_by_hemi.py
# and wta_partial_corr_by_subject_by_hemi.py):
#   Row 0      : header — ,<parcel1>,<parcel2>,...  (YAML canonical order)
#   Rows 1..N  : one subject per row  (label = sub-XX)
#   Row N+1    : GROUP
#   Last row   : CONSISTENCY_%  (never use as winner data)
#
# Column order is YAML canonical for BOTH pipelines. Values are looked up
# by parcel name from the header — completely order-independent.
#
# Arguments:
#   $1  path to CSV file
#   $2  row label  ("GROUP", "sub-01", …)
#   $3  name of associative array to populate
# Returns 1 on error.
# ============================================================
read_csv_row_into_array() {
    local csv_path="$1"
    local row_label="$2"
    local arr_name="$3"

    if [[ ! -f "$csv_path" ]]; then
        echo "ERROR: CSV file not found: $csv_path" >&2
        return 1
    fi

    local -a _all_lines=()
    mapfile -t _all_lines < "$csv_path"
    local n_lines="${#_all_lines[@]}"

    if [[ "$n_lines" -lt 3 ]]; then
        echo "ERROR: CSV has fewer than 3 lines (need header + data + GROUP): $csv_path" >&2
        return 1
    fi

    # Build name → column-index map from header
    # Header: ,<parcel1>,<parcel2>,...  (field 0 is the empty row-label slot)
    local -a _header_fields=()
    IFS=',' read -r -a _header_fields <<< "${_all_lines[0]}"

    declare -A _col_index=()
    local _hi
    for _hi in "${!_header_fields[@]}"; do
        local _pname="${_header_fields[$_hi]}"
        [[ -n "$_pname" ]] && _col_index["$_pname"]="$_hi"
    done

    if [[ "${#_col_index[@]}" -eq 0 ]]; then
        echo "ERROR: No parcel columns found in header of $csv_path" >&2
        return 1
    fi

    # Find target row by first-field label
    local _target_line="" _line _first_field
    for _line in "${_all_lines[@]:1}"; do
        _first_field="${_line%%,*}"
        if [[ "$_first_field" == "$row_label" ]]; then
            _target_line="$_line"
            break
        fi
    done

    if [[ -z "$_target_line" ]]; then
        echo "ERROR: Row label '${row_label}' not found in ${csv_path}" >&2
        return 1
    fi

    # Split target row into fields
    local -a _fields=()
    IFS=',' read -r -a _fields <<< "$_target_line"

    local expected_fields=$(( 1 + ${#_col_index[@]} ))
    if [[ "${#_fields[@]}" -ne "$expected_fields" ]]; then
        echo "ERROR: Expected $expected_fields fields but got ${#_fields[@]}" \
             "in row '${row_label}': $csv_path" >&2
        return 1
    fi

    # Populate output associative array: parcel_name → winner value
    eval "${arr_name}=()"
    local _pname _cidx
    for _pname in "${!_col_index[@]}"; do
        _cidx="${_col_index[$_pname]}"
        eval "${arr_name}[\"${_pname}\"]=\"\${_fields[${_cidx}]}\""
    done
}

# ============================================================
# write_hemi_labels
#
# Iterate over all parcels in the atlas-key table for one hemisphere,
# look up each winner by parcel name, and append one label entry per
# valid parcel to the output file.
#
# Arguments:
#   $1  output label file path
#   $2  hemi char  ("R" or "L")
#   $3  name of atlas-key associative array  ("R_KEYS" or "L_KEYS")
#   $4  name of winners associative array    (parcel_name → winner_value)
#   $5  name of integer parcel-count variable  (incremented in place)
#   $6  name of integer warning-count variable (incremented in place)
# ============================================================
write_hemi_labels() {
    local output_file="$1"
    local hemi_char="$2"
    local keys_arr_name="$3"
    local winners_arr_name="$4"
    local count_var="$5"
    local warn_var="$6"

    local parcel winner color key
    local -a _parcel_list=()
    eval "_parcel_list=(\"\${!${keys_arr_name}[@]}\")"

    for parcel in "${_parcel_list[@]}"; do

        eval "winner=\"\${${winners_arr_name}[\"${parcel}\"]:-}\""

        [[ -z "$winner" || "$winner" == "nan" ]] && continue

        winner=$(echo "$winner" | tr -d '[:space:]')
        winner=$(printf "%.0f" "$winner" 2>/dev/null)

        if ! validate_winner "$winner"; then
            echo "  WARNING: ${hemi_char}_${parcel} has invalid winner '$winner' — skipping" >&2
            continue
        fi

        if check_self_win "$parcel" "$winner"; then
            echo "  WARNING: ${hemi_char}_${parcel} won its own cluster (winner=$winner)" >&2
            eval "${warn_var}=$(( ${!warn_var} + 1 ))"
        fi

        color="${WINNER_COLORS[$winner]}"
        eval "key=\"\${${keys_arr_name}[\"${parcel}\"]:-}\""

        if [[ -z "$key" ]]; then
            echo "  WARNING: No atlas key for ${hemi_char}_${parcel} — skipping" >&2
            continue
        fi

        printf '%s\n%s %s\n' "${hemi_char}_${parcel}_ROI" "${key} ${color}" >> "$output_file"
        eval "${count_var}=$(( ${!count_var} + 1 ))"
    done
}

# ============================================================
# write_hemi_labels_pcef_ipef
#
# Identical to write_hemi_labels but only emits parcels whose winner
# is in [1,5]  (mPCS=1, sPCS=2, iPCS=3, sIPS=4, iIPS=5).
# Arguments: identical to write_hemi_labels.
# ============================================================
write_hemi_labels_pcef_ipef() {
    local output_file="$1"
    local hemi_char="$2"
    local keys_arr_name="$3"
    local winners_arr_name="$4"
    local count_var="$5"
    local warn_var="$6"

    local parcel winner color key
    local -a _parcel_list=()
    eval "_parcel_list=(\"\${!${keys_arr_name}[@]}\")"

    for parcel in "${_parcel_list[@]}"; do

        eval "winner=\"\${${winners_arr_name}[\"${parcel}\"]:-}\""

        [[ -z "$winner" || "$winner" == "nan" ]] && continue

        winner=$(echo "$winner" | tr -d '[:space:]')
        winner=$(printf "%.0f" "$winner" 2>/dev/null)

        if ! validate_winner "$winner"; then
            continue   # already warned in the full-map pass
        fi

        [[ "$winner" -lt 1 || "$winner" -gt 5 ]] && continue

        if check_self_win "$parcel" "$winner"; then
            eval "${warn_var}=$(( ${!warn_var} + 1 ))"
        fi

        color="${WINNER_COLORS[$winner]}"
        eval "key=\"\${${keys_arr_name}[\"${parcel}\"]:-}\""
        [[ -z "$key" ]] && continue

        printf '%s\n%s %s\n' "${hemi_char}_${parcel}_ROI" "${key} ${color}" >> "$output_file"
        eval "${count_var}=$(( ${!count_var} + 1 ))"
    done
}

# ============================================================
# build_label_file
# Write a full (all winners) label file for both hemispheres.
# ============================================================
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

# ============================================================
# build_label_file_pcef_ipef
# Write a PCEF/IPEF-filtered (winners 1–5 only) label file.
# ============================================================
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