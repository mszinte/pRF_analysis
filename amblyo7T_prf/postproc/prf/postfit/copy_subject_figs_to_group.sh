#!/bin/bash
#-----------------------------------------------------------------------------------------
# copy_subject_figs_to_group.sh
#-----------------------------------------------------------------------------------------
# Goal of the script:
# Launch the copy_subject_figs_to_group.py function for group-patient and group-control
#-----------------------------------------------------------------------------------------
# Input(s):
# input[1]: main data directory
# input[2]: project name (correspond to directory)
#-----------------------------------------------------------------------------------------
# Output(s):
# Figures copied to group folders organized by type
#-----------------------------------------------------------------------------------------
# To run:
# 1. cd to function
# >> cd ~/projects/pRF_analysis/amblyo7T_prf/postproc/prf/postfit
# 2. run shell command
# >> sh copy_subject_figs_to_group.sh [main directory] [project name]
#-----------------------------------------------------------------------------------------
# Example:
# cd ~/projects/pRF_analysis/amblyo7T_prf/postproc/prf/postfit
# sh copy_subject_figs_to_group.sh /scratch/mszinte/data amblyo7T_prf
#-----------------------------------------------------------------------------------------
# Written by Martin Szinte (martin.szinte@gmail.com)
#-----------------------------------------------------------------------------------------
# Check inputs
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <data_path> <project_name>"
    exit 1
fi
# Define inputs
data_path="$1"
project_name="$2"
# cd to script location
cd ~/projects/pRF_analysis/${project_name}/postproc/prf/postfit
# Run for both groups
for group_label in group-patient group-control
do
    echo "Processing copy_subject_figs_to_group.py for: $group_label"
    python copy_subject_figs_to_group.py "$data_path" "$project_name" "$group_label"
done