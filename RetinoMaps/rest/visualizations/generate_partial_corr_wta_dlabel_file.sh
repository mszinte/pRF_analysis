#!/bin/bash
#####################################################
# Generate WTA dlabel files from partial-correlation CSVs
#
# Reads winning_seeds_by_subject_partial_corr_{hemi}_{variant}.csv
# produced by wta_partial_corr_by_subject_by_hemi.py.
#
# Usage:
#   $ ./generate_partial_corr_wta_dlabel_file.sh
#
# Written by Marco Bedini (marco.bedini@univ-amu.fr)
#####################################################

echo "============================================"
echo "WTA DLABEL GENERATION — partial correlation"
echo "============================================"
echo ""

# ============================================================
# Source shared utilities
# (WINNER_COLORS, R_KEYS, L_KEYS, SUBJECTS, validate_winner,
#  check_self_win, read_csv_row_into_array,
#  write_hemi_labels, write_hemi_labels_pcef_ipef,
#  build_label_file, build_label_file_pcef_ipef)
# ============================================================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../utils/rest_utils.sh"

# ============================================================
# Paths
# ============================================================
BASE_PATH="/scratch/mszinte/data/RetinoMaps/derivatives/pp_data"
ATLAS_DIR="${BASE_PATH}/atlas"
INPUT_PATH="${BASE_PATH}/group/91k/rest/wta/nilearn"
OUTPUT_PATH="${BASE_PATH}/group/91k/rest/wta/nilearn/dlabel"
ATLAS_DLABEL="${ATLAS_DIR}/mmp1/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_dseg.dlabel.nii"

# Temporary label text files (overwritten on each wb_command call)
_LABEL_FILE="${ATLAS_DIR}/mmp1/wta_mmp1_labels_partial_corr.txt"
_LABEL_FILE_PCEF="${ATLAS_DIR}/mmp1/wta_mmp1_labels_partial_corr_pcef-ipef.txt"

mkdir -p "${OUTPUT_PATH}"

# ============================================================
# read_all_winners
#
# Populate rh and lh winner associative arrays for one row label
# and variant from the per-hemisphere partial-corr CSVs.
#
# CSV naming convention (from wta_partial_corr_by_subject_by_hemi.py):
#   winning_seeds_by_subject_partial_corr_{hemi}_{variant}.csv
#
# Arguments:
#   $1  row_label  ("GROUP", "sub-01", …)
#   $2  variant    ("concat", "concat_clean", "run-01", "run-02")
#   $3  name of rh winners associative array to populate
#   $4  name of lh winners associative array to populate
# Returns 1 on error.
# ============================================================
read_all_winners() {
    local row_label="$1"
    local variant="$2"
    local rh_arr_name="$3"
    local lh_arr_name="$4"

    local hemi csv_path

    for hemi in rh lh; do
        csv_path="${INPUT_PATH}/winning_seeds_by_subject_partial_corr_${hemi}_${variant}.csv"
        if [[ "$hemi" == "rh" ]]; then
            read_csv_row_into_array "$csv_path" "$row_label" "$rh_arr_name" || return 1
        else
            read_csv_row_into_array "$csv_path" "$row_label" "$lh_arr_name" || return 1
        fi
    done
}

# ============================================================
# process_variant
#
# Generate group + per-subject dlabel files (full map and
# PCEF/IPEF-filtered) for one variant.
# ============================================================
process_variant() {
    local variant="$1"

    echo ""
    echo "============================================"
    echo "Variant: ${variant}"
    echo "============================================"

    # ── GROUP ────────────────────────────────────────────────
    echo "  Processing GROUP..."

    declare -A grp_rh=()
    declare -A grp_lh=()

    if ! read_all_winners "GROUP" "$variant" "grp_rh" "grp_lh"; then
        echo "  ERROR: Failed to read GROUP winners for '${variant}' — skipping." >&2
        return 1
    fi

    local grp_count=0 grp_warn=0
    build_label_file "$_LABEL_FILE" "grp_rh" "grp_lh" "grp_count" "grp_warn"

    wb_command -cifti-label-import \
        "$ATLAS_DLABEL" "$_LABEL_FILE" \
        "${OUTPUT_PATH}/group_wta_partial_corr_${variant}.dlabel.nii" \
        -discard-others

    local grp_pcef_count=0 grp_pcef_warn=0
    build_label_file_pcef_ipef "$_LABEL_FILE_PCEF" "grp_rh" "grp_lh" \
        "grp_pcef_count" "grp_pcef_warn"

    wb_command -cifti-label-import \
        "$ATLAS_DLABEL" "$_LABEL_FILE_PCEF" \
        "${OUTPUT_PATH}/group_wta_partial_corr_${variant}_pcef-ipef.dlabel.nii" \
        -discard-others

    echo "  GROUP: labels=${grp_count} warnings=${grp_warn}" \
         "| PCEF/IPEF labels=${grp_pcef_count}"

    # ── SUBJECTS ─────────────────────────────────────────────
    local sub
    for sub in "${SUBJECTS[@]}"; do

        declare -A sub_rh=()
        declare -A sub_lh=()

        if ! read_all_winners "$sub" "$variant" "sub_rh" "sub_lh"; then
            echo "  ${sub}: ERROR — skipping" >&2
            unset sub_rh sub_lh
            continue
        fi

        local sub_labels="${OUTPUT_PATH}/${sub}_wta_partial_corr_${variant}_labels.txt"
        local sub_count=0 sub_warn=0
        build_label_file "$sub_labels" "sub_rh" "sub_lh" "sub_count" "sub_warn"

        wb_command -cifti-label-import \
            "$ATLAS_DLABEL" "$sub_labels" \
            "${OUTPUT_PATH}/${sub}_wta_partial_corr_${variant}.dlabel.nii" \
            -discard-others

        local sub_labels_pcef="${OUTPUT_PATH}/${sub}_wta_partial_corr_${variant}_pcef-ipef_labels.txt"
        local sub_pcef_count=0 sub_pcef_warn=0
        build_label_file_pcef_ipef "$sub_labels_pcef" "sub_rh" "sub_lh" \
            "sub_pcef_count" "sub_pcef_warn"

        wb_command -cifti-label-import \
            "$ATLAS_DLABEL" "$sub_labels_pcef" \
            "${OUTPUT_PATH}/${sub}_wta_partial_corr_${variant}_pcef-ipef.dlabel.nii" \
            -discard-others

        if [[ $sub_warn -gt 0 ]]; then
            echo "  ${sub}: labels=${sub_count} warnings=${sub_warn}"
        else
            echo "  ${sub}: OK (labels=${sub_count})"
        fi

        unset sub_rh sub_lh
    done
}

# ============================================================
# Main — loop over all variants
# ============================================================
for variant in concat concat_clean run-01 run-02; do
    process_variant "$variant"
done

echo ""
echo "============================================"
echo "ALL VARIANTS COMPLETE"
echo "============================================"
echo "Outputs: ${OUTPUT_PATH}/"