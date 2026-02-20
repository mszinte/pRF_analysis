#!/bin/bash
#-----------------------------------------------------------------------------------------
# pycortex_maps_gauss.sh
# -----------------------------------------------------------------------------------------
# Goal of the script:
# Launch across subjects the function pycortex_maps_gauss.py
# -----------------------------------------------------------------------------------------
# Input(s):
# input[1]: project code directory
# input[2]: project name (correspond to directory)
# input[3]: main data directory (correspond to directory)
# input[4]: Save maps in the overlay (y/n)
# -----------------------------------------------------------------------------------------
# Output(s):
# All pycortex maps for grid fit
# -----------------------------------------------------------------------------------------
# To run:
# 0. TO RUN ON INVIBE SERVER (with Inkscape)
# 1. cd to function
# >> cd ~/disks/meso_H/projects/pRF_analysis/analysis_code/postproc/prf/postfit
# 2. run python command
# >> sh pycortex_maps_gauss.sh [code directory] [project name] [main directory] 
#                                [save_in_overlay] [analysis folder]
# -----------------------------------------------------------------------------------------
# Exemple:
# cd ~/disks/meso_H/projects/pRF_analysis/analysis_code/postproc/prf/postfit
# sh pycortex_maps_gauss.sh ~/disks/meso_H/projects RetinoMaps ~/disks/meso_S/data n
# -----------------------------------------------------------------------------------------
# Written by Martin Szinte (martin.szinte@gmail.com)
# Edited by Uriel Lascombes (uriel.lascombes@laposte.net)
# -----------------------------------------------------------------------------------------

# Define the base path, project name, and data path
base_path="$1"
project_name="$2"
data_path="$3"
save_in_overlay="$4"

# Define the path to the settings.yml file
settings_file="${base_path}/pRF_analysis/${project_name}/settings.yml"

# Define current directory
cd "${base_path}/pRF_analysis/analysis_code/postproc/prf/postfit"

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
    echo "Processing pycortex_maps_gauss.py for: $subject"
    python pycortex_maps_gauss.py "$data_path" "$project_name" "$subject" "$save_in_overlay" 
done