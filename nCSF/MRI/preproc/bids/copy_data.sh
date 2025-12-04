#!/bin/bash
#-----------------------------------------------------------------------------------------
# copy_data.sh
# -----------------------------------------------------------------------------------------
# Goal of the script:
# Copy data from CRMBM server and create directories
# -----------------------------------------------------------------------------------------
# Input(s):
# input[1]: origin data directory
# input[2]: desired data directory
# input[3]: code directory
# -----------------------------------------------------------------------------------------
# Output(s):
# Transfer data into sourcedata
# -----------------------------------------------------------------------------------------
# To run:
# 1. cd to function
# >> cd ~/projects/[PROJECT]/analysis_code/preproc/bids
# 2. run shell command
# >> sh copy_data.sh [origin dir] [destination dir] [code directory]
# -----------------------------------------------------------------------------------------
# Exemple:
# cd ~/projects/nCSF/analysis_code/preproc/bids
# sh copy_data.sh /Volumes/DATA_CEMEREM/data/protocols/terra/brain/terrastart /Users/uriel/disks/meso_shared/nCSF /Users/uriel/disks/meso_H/projects/nCSFexp
# -----------------------------------------------------------------------------------------
# Written by Uriel Lascombes (uriel.lascombes@laposte.net) 
# Edited by Martin Szinte (martin.szinte@gmail.com)
# -----------------------------------------------------------------------------------------

# Check if the origin dir and destination dir are provided as arguments
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <origin_dir> <destination_dir> <code_dir>"
    exit 1
fi

# Define the origin dir and destination dir
origin_dir="$1"
destination_dir="$2"
code_dir="$3"

# Create additional BIDS folders
mkdir "${destination_dir}/code"
mkdir "${destination_dir}/derivatives"
mkdir "${destination_dir}/sourcedata"

# Define the path to the settings_MRI.json file
settings_file="${code_dir}/analysis_code/settings_MRI.json"

# Read the subjects from settings.json using Python
subjects=$(python -c "import json; data = json.load(open('$settings_file')); print('\n'.join(data['source_subject_name']))")

for subject in subjects
do
    # Copy raw data
    rsync -avuz --progress "" /scratch/mszinte/data/RetinoMaps/code
    



