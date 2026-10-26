#!/bin/bash
#-----------------------------------------------------------------------------------------
# pycortex_webgl_retinomaps.sh
# -----------------------------------------------------------------------------------------
# Goal of the script:
# Launch across subjects the function pycortex_webgl_css.py
# -----------------------------------------------------------------------------------------
# Input(s):
# input[1]: project code directory
# input[2]: project name (correspond to directory)
# input[3]: main data directory (correspond to directory)
# -----------------------------------------------------------------------------------------
# Output(s):
# All pycortex maps for ROIs
# -----------------------------------------------------------------------------------------
# To run:
# 1. cd to function
# >> cd ~/projects/pRF_analysis/RetinoMaps/webgl/
# 2. run python command
# >> sh pycortex_webgl_retinomaps.sh [code directory] [project name] [main directory]
# -----------------------------------------------------------------------------------------
# Exemple:
# cd ~/projects/pRF_analysis/RetinoMaps/webgl/
# sh pycortex_webgl_retinomaps.sh ~/projects RetinoMaps /scratch/mszinte/data
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

# Define the path to the settings.yml file
settings_file="${base_path}/pRF_analysis/${project_name}/settings.yml"

# Define current directory
cd "${base_path}/pRF_analysis/${project_name}/webgl"

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
    echo "Processing pycortex_webgl_css.py for: $subject"
    python pycortex_webgl_retinomaps.py "$data_path" "$project_name" "$subject" 327 1
done


