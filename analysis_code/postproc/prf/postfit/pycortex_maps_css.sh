#!/bin/bash
#-----------------------------------------------------------------------------------------
# pycortex_maps_css.sh
# -----------------------------------------------------------------------------------------
# Goal of the script:
# Launch across subjects the function pycortex_maps_rois.py
# -----------------------------------------------------------------------------------------
# Input(s):
# input[1]: project code directory
# input[2]: project name (correspond to directory)
# input[3]: main data directory (correspond to directory)
# -----------------------------------------------------------------------------------------
# Output(s):
# All pycortex maps for CSS
# -----------------------------------------------------------------------------------------
# To run:
# 0. TO RUN ON INVIBE SERVER (with Inkscape)
# 1. cd to function
# >> cd ~/disks/meso_H/projects/[PROJECT]/analysis_code/postproc/prf/postfit
# 2. run python command
# >> sh pycortex_maps_css.sh [code directory] [project name] [main directory]
# -----------------------------------------------------------------------------------------
# Exemple:
# cd ~/disks/meso_H/projects/pRF_analysis/analysis_code/postproc/prf/postfit
# sh pycortex_maps_css.sh ~/disks/meso_H/projects MotConf ~/disks/meso_S/data
# -----------------------------------------------------------------------------------------
# Written by Martin Szinte (martin.szinte@gmail.com)
# Edited by Uriel Lascombes (uriel.lascombes@laposte.net)
# -----------------------------------------------------------------------------------------

# Check if the base path, project name, and data path are provided as arguments
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <base_path> <project_name> <data_path>"
    exit 1
fi

# Define the base path, project name, and data path
base_path="$1"
project_name="$2"
data_path="$3"

# Define the path to the settings.json file
settings_file="${base_path}/${project_name}/analysis_code/settings.json"

# Define current directory
cd "${base_path}/${project_name}/analysis_code/postproc/prf/postfit"

# Read the subjects from settings.json using Python
subjects=$(python -c "import json; data = json.load(open('$settings_file')); print('\n'.join(data['subjects']))")

# Loop through each subject and run the Python code
for subject in $subjects
do
    echo "Processing pycortex_maps_css.py for: $subject"
    python pycortex_maps_css.py "$data_path" "$project_name" "$subject" n
done
