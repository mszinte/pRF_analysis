#!/bin/bash

# Input file
input_file="Glasser_filtered_labels.txt"

# Output dir
mkdir label_files

# Declare associative array of clusters
declare -A clusters=(
  ["mPCS"]="SCEF 24dv p32pr"
  ["sPCS"]="FEF i6-8 6a 6d 6mp 6ma"
  ["iPCS"]="PEF IFJp 6v 6r IFJa 55b"
  ["sIPS"]="VIP LIPv LIPd IP2 7PC AIP 7AL 7Am 7Pm"
  ["iIPS"]="IP0 IPS1 V7 MIP IP1 V6A 7PL"
  ["hMT+"]="V4t MST MT FST"
  ["VO"]="V8 PIT PH FFC VMV1 VMV2 VMV3 VVC"
  ["LO"]="LO1 LO2 LO3"
  ["V3AB"]="V3A V3B V3CD"
  ["V3"]="V3 V4"
  ["V2"]="V2"
  ["V1"]="V1"
)

# Read the input file line-by-line
while read -r label && read -r values; do
  # Extract hemisphere (lh/rh) and base name of label
  if [[ $label =~ ^(R|L)_(.+)_ROI$ ]]; then
    hemi=${BASH_REMATCH[1]}
    roi=${BASH_REMATCH[2]}
    prefix="rh"  # Default hemisphere label
    [[ $hemi == "L" ]] && prefix="lh"

    # Check which cluster this ROI belongs to
    for cluster in "${!clusters[@]}"; do
      for roi_match in ${clusters[$cluster]}; do
        if [[ "$roi" == "$roi_match" ]]; then
          # Append to appropriate output file
          output_file="label_files/${prefix}_${cluster}_labels.txt"
          echo "$label" >> "$output_file"
          echo "$values" >> "$output_file"
        fi
      done
    done
  fi
done < "$input_file"

# Change permissions
chmod -Rf 771 *
chgrp -Rf 771 *

