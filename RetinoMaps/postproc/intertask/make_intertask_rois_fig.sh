#!/bin/bash
#-----------------------------------------------------------------------------------------
# make_intertask_rois_fig.sh
# -----------------------------------------------------------------------------------------
# Goal of the script:
# Launch across subjects the function make_intertask_rois_fig.py
# -----------------------------------------------------------------------------------------
# Input(s):
# input[1]: project code directory
# input[2]: project name (correspond to directory)
# input[3]: main data directory (correspond to directory)
# -----------------------------------------------------------------------------------------
# Output(s):
# All ROI based figures 
# -----------------------------------------------------------------------------------------
# To run:
# 1. cd to function
# >> cd ~/projects/pRF_analysis/RetinoMaps/postproc/intertask/
# 2. run shell command
# >> sh make_rois_fig.sh [code directory] [project name] [main directory]
# -----------------------------------------------------------------------------------------
# Exemple:
# cd ~/projects/pRF_analysis/RetinoMaps/postproc/intertask/
# sh make_intertask_rois_fig.sh ~/projects RetinoMaps /scratch/mszinte/data
# -----------------------------------------------------------------------------------------
# Written by Martin Szinte (martin.szinte@gmail.com)
# Edited by Uriel Lascombes (uriel.lascombes@laposte.net)
# -----------------------------------------------------------------------------------------

# Define the base path, project name, and data path
base_path="$1"
project_name="$2"
data_path="$3"

# Define the path to the settings.yml file
settings_file="${base_path}/pRF_analysis/${project_name}/settings.yml"

# Define current directory
cd "${base_path}/pRF_analysis/${project_name}/postproc/intertask"

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
    echo "Processing make_intertask_rois_fig.py for: $subject"
    python make_intertask_rois_fig.py "$data_path" "$project_name" "$subject" 327
done