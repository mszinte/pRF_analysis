#!/bin/bash
#-----------------------------------------------------------------------------------------
# make_od_img.sh
#-----------------------------------------------------------------------------------------
# Goal of the script:
# Launch across subjects the function make_od_img.py
#-----------------------------------------------------------------------------------------
# Input(s):
# input[1]: project code directory
# input[2]: project name (correspond to directory)
# input[3]: main data directory
# input[4]: analysis name (e.g prf)
# input[5]: server group (e.g. 327)
#-----------------------------------------------------------------------------------------
# Output(s):
# Ocular dominance surface images per subject
#-----------------------------------------------------------------------------------------
# To run:
# 1. cd to function
# >> cd ~/projects/pRF_analysis/amblyo7T_prf/postproc/prf/postfit
# 2. run shell command
# >> sh make_od_img.sh [code directory] [project name] [main directory] [analysis name] [group]
#-----------------------------------------------------------------------------------------
# Example:
# cd ~/projects/pRF_analysis/amblyo7T_prf/postproc/prf/postfit
# sh make_od_img.sh ~/projects amblyo7T_prf /scratch/mszinte/data prf 327
# sh make_od_img.sh ~/disks/meso_H/projects amblyo7T_prf ~/disks/meso_shared prf 327
#-----------------------------------------------------------------------------------------
# Written by Uriel Lascombes (uriel.lascombes@laposte.net)
#-----------------------------------------------------------------------------------------

# Check inputs
if [ "$#" -ne 5 ]; then
    echo "Usage: $0 <base_path> <project_name> <data_path> <analysis_name> <group>"
    exit 1
fi

# Define inputs
base_path="$1"
project_name="$2"
data_path="$3"
analysis_name="$4"
group="$5"

# Define settings file path
settings_file="${base_path}/pRF_analysis/${project_name}/settings.yml"

# cd to script location
cd "${base_path}/pRF_analysis/${project_name}/postproc/prf/postfit"

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
    echo "Processing make_od_img.py for: $subject"
    python make_od_img.py "$data_path" "$project_name" "$subject" "$analysis_name" "$group"
done
