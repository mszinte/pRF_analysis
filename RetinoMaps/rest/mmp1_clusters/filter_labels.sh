#!/bin/bash

# Define the file containing the label table
label_file="Glasser_labels.txt"

# Specify the labels you want to filter (without L_ or R_)
labels_to_keep=("SCEF" "24dv" "p32pr" "FEF" "i6-8" "6a" "6d" "6mp" "6ma" "PEF" "IFJp" 
"6v" "6r" "IFJa" "55b" "VIP" "LIPv" "LIPd" "IP2" "7PC" "AIP" "7AL" "7Am" "7Pm" "IP0" "IPS1" 
"V7" "MIP" "IP1" "V6A" "7PL" "V4t" "MST" "MT" "FST" "V8" "PIT" "PH" "FFC" "VMV1" "VMV2" 
"VMV3" "VVC" "LO1" "LO2" "LO3" "V3A" "V3B" "V3CD" "V3" "V4" "V2" "V1")

# Redirect the filtered output to Glasser_filtered_labels.txt
output_file="Glasser_filtered_labels.txt"
> "$output_file"  # Clear the output file before appending new lines

# Read the file line by line, and process pairs of label and RGB line
while IFS= read -r line; do
  for label in "${labels_to_keep[@]}"; do
    # Match only exact labels with the hemisphere prefix (L_ or R_)
    if [[ $line == "L_${label}_ROI" ]] || [[ $line == "R_${label}_ROI" ]]; then
      echo "$line" >> "$output_file"   # Write the label
      IFS= read -r next_line           # Safely read the next line containing index and RGB
      echo "$next_line" >> "$output_file"  # Write the next line (index and RGB)
    fi
  done
done < "$label_file"

echo "Filtered labels and RGB values have been written to $output_file"

# Change permissions
chmod -Rf 771 *
chgrp -Rf 771 *
