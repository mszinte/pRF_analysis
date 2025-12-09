"""
-----------------------------------------------------------------------------------------
bids_conversion.py
-----------------------------------------------------------------------------------------
Goal of the script:
Rename subject and session
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: server group (e.g. 327)
-----------------------------------------------------------------------------------------
Output(s):
rename files and folders
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/pRF_analysis/nCSF/MRI/preproc/bids
2. run python command
python bids_conversion.py [main directory] [project name] [group]
-----------------------------------------------------------------------------------------
Example:
cd ~/projects/pRF_analysis/nCSF/MRI/preproc/bids
python bids_conversion.py /scratch/mszinte/data nCSF 327
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
import json

# Personal import
sys.path.append("{}/../../../utils".format(os.getcwd()))
from bids_utils import bidsify_fmap, bidsify_anat

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
group = sys.argv[3]

# Load settings and analysis parameters
with open('../../../settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)
subjects = analysis_info['subjects']
sessions = analysis_info['sessions']

bids_data_dir = "{}/{}".format(main_dir, project_dir)

for n_subject, subject in enumerate(subjects):
    for n_session, session in enumerate(sessions):

        bidsify_fmap(
            sub=subject,
            ses=session,
            base_dir=bids_data_dir
        )

        bidsify_anat(
            sub=subject,
            ses=session,
            base_dir=bids_data_dir
        )

# Define permission cmd
print('Changing files permissions in {}/{}'.format(main_dir, project_dir))
os.system("chmod -Rf 771 {}/{}".format(main_dir, project_dir))
os.system("chgrp -Rf {} {}/{}".format(group, main_dir, project_dir))