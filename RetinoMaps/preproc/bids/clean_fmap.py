"""
-----------------------------------------------------------------------------------------
clean_fmap.py
-----------------------------------------------------------------------------------------
Goal of the script:
Remove tasks not in RetinoMaps for open neuro in the intended for
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject ID (e.g. sub-01)
sys.argv[4]: server group (e.g. 327)

-----------------------------------------------------------------------------------------
Output(s):
clean fmap 
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/pRF_analysis/RetinoMaps/preproc/bids
2. run python command
python clean_fmap.py [main directory] [project name] [subject] [group]
-----------------------------------------------------------------------------------------
Example:
cd ~/projects/pRF_analysis/RetinoMaps/preproc/bids
python clean_fmap.py /scratch/mszinte/data RetinoMapsOpenNeuro sub-01 327
-----------------------------------------------------------------------------------------
Written by Uriel Lascombes (uriel.lascombes@laposte.net)
Edited by Martin Szinte (mail@martinszinte.net)
-----------------------------------------------------------------------------------------
"""

# Stop warnings
import warnings
warnings.filterwarnings("ignore")

# Debug
import ipdb
deb = ipdb.set_trace

# General imports
import os
import sys
import glob
import json

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
group = sys.argv[4]

tasks_to_keep = ['pRF','SacLoc','PurLoc','rest']

fmap_dir = '{}/{}'.format(main_dir, project_dir)

epi_files = glob.glob(os.path.join(fmap_dir, subject, 'ses-*', 'fmap', '*_epi.json'))
fieldmap_files = glob.glob(os.path.join(fmap_dir, subject, 'ses-*', 'fmap', '*_fieldmap.json'))

for json_file in epi_files + fieldmap_files:
    with open(json_file, 'r') as f:
        data = json.load(f)

    if 'IntendedFor' in data:
        data['IntendedFor'] = [
            item for item in data['IntendedFor']
            if any(task in item for task in tasks_to_keep)
        ]

    with open(json_file, 'w') as f:
        json.dump(data, f, indent=4)



print("BIDS correction done.")

# # Change permissions (kept exactly as requested)
# print("Changing files permissions in {}/{}".format(main_dir, project_dir))
# os.system("chmod -Rf 771 {}/{}".format(main_dir, project_dir))
# os.system("chgrp -Rf {} {}/{}".format(group, main_dir, project_dir))