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
>> cd ~/projects/pRF_analysis/RetinoMaps/preproc/bids
2. run python command
python bids_conversion.py [main directory] [project name] [group]
-----------------------------------------------------------------------------------------
Example:
cd ~/projects/pRF_analysis/RetinoMaps/preproc/bids
python bids_conversion.py /scratch/mszinte/data RetinoMaps 327
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


# Personal imports
sys.path.append("{}/../../../analysis_code/utils".format(os.getcwd()))
from settings_utils import load_settings

sys.path.append("{}/../../utils".format(os.getcwd()))
from bids_utils import bidsify_fmap, bidsify_anat, bidsify_func, correct_task_columns_event

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
group = sys.argv[3]

# Load settings
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../"))
settings_path = os.path.join(base_dir, project_dir, "settings.yml")
settings = load_settings([settings_path])
analysis_info = settings[0]

subjects = analysis_info["subjects"]
sessions = analysis_info["sessions"]
tasks = analysis_info["task_names"]

base_dir = "{}/{}".format(main_dir, project_dir)

for subject in subjects:
    print(subject)
    for session in sessions:

        # Apply BIDS fixes
        bidsify_anat(subject, session, base_dir)
        bidsify_fmap(subject, session, base_dir)
        bidsify_func(subject, session, base_dir)
        
        # Change task column in DeepMReye events for trial_type
        correct_task_columns_event(base_dir, subject, session)
        


print("BIDS correction done.")

# Change permissions (kept exactly as requested)
print("Changing files permissions in {}/{}".format(main_dir, project_dir))
os.system("chmod -Rf 771 {}/{}".format(main_dir, project_dir))
os.system("chgrp -Rf {} {}/{}".format(group, main_dir, project_dir))