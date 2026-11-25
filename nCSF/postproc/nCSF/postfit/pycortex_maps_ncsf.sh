#!/bin/bash
#-----------------------------------------------------------------------------------------
# pycortex_maps_ncsf.sh
# -----------------------------------------------------------------------------------------
# Goal of the script:
# Launch across subjects the function pycortex_maps_ncsf.py
# -----------------------------------------------------------------------------------------
# Input(s):
# input[1]: project code directory
# input[2]: project name (correspond to directory)
# input[3]: main data directory (correspond to directory)
# input[4]: Save maps in the overlay (y/n)
# (analysis folder read by pycortex_maps_ncsf.py from output_folder of nCSF-analysis.yml)
# -----------------------------------------------------------------------------------------
# Output(s):
# All pycortex maps for nCSF
# -----------------------------------------------------------------------------------------
# To run:
# 0. TO RUN ON INVIBE SERVER (with Inkscape)
# 1. cd to function
# >> cd ~/disks/meso_H/projects/pRF_analysis/nCSF/postproc/nCSF/postfit
# 2. run python command
# >> sh pycortex_maps_ncsf.sh [code directory] [project name] 
#                             [main directory] [save_in_overlay]
# -----------------------------------------------------------------------------------------
# Exemple:
# cd ~/disks/meso_H/projects/pRF_analysis/nCSF/postproc/nCSF/postfit
# sh pycortex_maps_ncsf.sh ~/disks/meso_H/projects nCSF ~/disks/meso_shared n
# -----------------------------------------------------------------------------------------
# Written by Martin Szinte (martin.szinte@gmail.com)
# Edited by Uriel Lascombes (uriel.lascombes@laposte.net)
# -----------------------------------------------------------------------------------------

# Check if the base path, project name, and data path are provided as arguments
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <base_path> <project_name> <data_path> <save_in_overlay>"
    exit 1
fi

# Define the base path, project name, and data path
base_path="$1"
project_name="$2"
data_path="$3"
save_in_overlay="$4"

# Define the path to the settings.yml file
settings_file="${base_path}/pRF_analysis/${project_name}/settings.yml"

# Define current directory
# [NCSF-PRFMODEL] AVANT : cd "${base_path}/pRF_analysis/analysis_code/postproc/prf/postfit"
cd "${base_path}/pRF_analysis/${project_name}/postproc/nCSF/postfit"

# Read the subjects from settings.yml using Python
subjects=$(python -c "
import yaml
with open('$settings_file', 'r') as file:
    data = yaml.safe_load(file)
    print('\n'.join(data['subjects']['value']))
")

# Loop through each subject and run the Python code
for subject in $subjects
do
    echo "Processing pycortex_maps_ncsf.py for: $subject"
    # [NCSF-PRFMODEL] AVANT : python pycortex_maps_css.py "$data_path" "$project_name" "$subject" "$save_in_overlay" "$output_folder"
    python pycortex_maps_ncsf.py "$data_path" "$project_name" "$subject" "$save_in_overlay"
done