#!/bin/bash
#####################################################
# Written by Marco Bedini (marco.bedini@univ-amu.fr)
#####################################################

# Output dir
out_dir="/scratch/mszinte/data/RetinoMaps/derivatives/pp_data/atlas/mmp1/label_files"
mkdir -p "$out_dir"

# ---------------------------------------------------------------------------
# Write macro-region label files for wb_command -metric-label-import
# Format:
#     <region_name>
#     <key> <R> <G> <B> <A>
# Keys: mPCS=1 ... V1=12
# Colors sourced from figure-settings.yml > rois_colors > value
# TODO: replace hardcoded values with dynamic parsing of figure-settings.yml
# ---------------------------------------------------------------------------

write_label() {
    local hemi=$1 region=$2 key=$3 r=$4 g=$5 b=$6
    printf '%s\n%d %d %d %d %d\n' "$region" "$key" "$r" "$g" "$b" 255 \
        > "$out_dir/${hemi}_${region}_macro-region_label.txt"
    echo "  wrote [key=$key] ${hemi}_${region}_macro-region_label.txt  RGB($r,$g,$b)"
}

#          hemi  region   key  R    G    B
write_label lh   mPCS      1  255  111    0
write_label rh   mPCS      1  255  111    0
write_label lh   sPCS      1  255  234    0
write_label rh   sPCS      1  255  234    0
write_label lh   iPCS      1  151  255    0
write_label rh   iPCS      1  151  255    0
write_label lh   sIPS      1   44  255  150
write_label rh   sIPS      1   44  255  150
write_label lh   iIPS      1    0  152  255
write_label rh   iIPS      1    0  152  255
write_label lh   "hMT+"    1    0   25  255
write_label rh   "hMT+"    1    0   25  255
write_label lh   VO        1    0    0  200
write_label rh   VO        1    0    0  200
write_label lh   LO        1  150    0   90
write_label rh   LO        1  150    0   90
write_label lh   V3AB      1  235  127  134
write_label rh   V3AB      1  235  127  134
write_label lh   V3        1  248  160  126
write_label rh   V3        1  248  160  126
write_label lh   V2        1  250  196  132
write_label rh   V2        1  250  196  132
write_label lh   V1        1  243  231  155
write_label rh   V1        1  243  231  155

# Fix permissions
chmod -Rf 771 "$out_dir"
chgrp -Rf 327 "$out_dir"

echo ""
echo "Done. 24 label files written to: $out_dir"