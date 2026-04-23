#!/bin/bash
#-----------------------------------------------------------------------------------------
# make_corr_fig.sh
#-----------------------------------------------------------------------------------------
# Goal of the script:
# Launch across subjects the function make_corr_fig.py
#-----------------------------------------------------------------------------------------
# Input(s):
# input[1]: project code directory
# input[2]: project name (correspond to directory)
# input[3]: main data directory
# input[4]: server group (e.g. 327)
#-----------------------------------------------------------------------------------------
# Output(s):
# Correlation figures across eyes per subject and group
#-----------------------------------------------------------------------------------------
# To run:
# 1. cd to function
# >> cd ~/projects/pRF_analysis/amblyo7T_prf/postproc/prf/postfit
# 2. run shell command
# >> sh make_corr_fig.sh [code directory] [project name] [main directory] [group]
#-----------------------------------------------------------------------------------------
# Example:
# cd ~/projects/pRF_analysis/amblyo7T_prf/postproc/prf/postfit
# sh make_corr_fig.sh /home/mszinte/projects amblyo7T_prf /scratch/mszinte/data 327
#-----------------------------------------------------------------------------------------
# Written by Martin Szinte (martin.szinte@gmail.com)
#-----------------------------------------------------------------------------------------

# Check inputs
if [ "$#" -ne 4 ]; then
    echo "Usage: $0 <base_path> <project_name> <data_path> <group>"
    exit 1
fi

# Define inputs
base_path="$1"
project_name="$2"
data_path="$3"
group="$4"

# Define settings file path
settings_file="${base_path}/pRF_analysis/${project_name}/settings.yml"

# cd to script location
cd "${base_path}/pRF_analysis/${project_name}/postproc/prf/postfit"

# Read individual subjects from settings.yml
subjects=$(python -c "
import yaml
with open('$settings_file', 'r') as file:
    data = yaml.safe_load(file)
    print('\n'.join(data['subjects']['value']))
")

# Run per individual subject
for subject in $subjects
do
    echo "Processing make_corr_fig.py for: $subject"
    python make_corr_fig.py "$data_path" "$project_name" "$subject" "$group"
done

# Run for groups
for group_label in group-patient group-control
do
    echo "Processing make_corr_fig.py for: $group_label"
    python make_corr_fig.py "$data_path" "$project_name" "$group_label" "$group"
done