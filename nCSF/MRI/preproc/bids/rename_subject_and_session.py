"""
-----------------------------------------------------------------------------------------
rename_subject_and_session.py
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
python rename_subject_and_session.py [main directory] [project name] [group]
-----------------------------------------------------------------------------------------
Example:
cd ~/projects/pRF_analysis/nCSF/MRI/preproc/bids
python rename_subject_and_session.py /scratch/mszinte/data nCSF 327
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
sys.path.append("{}/../../../utils".format(os.getcwd()))
from settings_utils import load_settings


# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
group = sys.argv[3]

# Load settings
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../../"))
settings_path = os.path.join(base_dir, "nCSF", "settings.yml")
settings = load_settings([settings_path])
analysis_info = settings[0]

subject_names = analysis_info['subject_names']
session_names = analysis_info['session_names']

# Define directories
sourcedata_dir = "{}/{}/sourcedata".format(main_dir, project_dir)

# Rename subjects directly in sourcedata
renamed_folders = []

for old_name, new_name in subject_names.items():
    old_path = os.path.join(sourcedata_dir, old_name)
    new_path = os.path.join(sourcedata_dir, new_name)

    if not os.path.exists(old_path):
        print("Source not found: {}".format(old_path))
        continue

    if os.path.exists(new_path):
        print("Destination already exists: {} — skipping folder rename.".format(new_path))
        renamed_folders.append(new_path)
        continue

    print("Renaming folder {} to {}".format(old_path, new_path))
    os.rename(old_path, new_path)
    renamed_folders.append(new_path)

# Rename all files and folders inside the renamed folders
full_mapping = {}
full_mapping.update(subject_names)
full_mapping.update(session_names)

for folder in renamed_folders:
    print("Renaming inside folder: {}".format(folder))

    for dirpath, dirnames, filenames in os.walk(folder, topdown=False):

        # Rename files
        for fname in filenames:
            new_fname = fname
            for old, new in full_mapping.items():
                if old in new_fname:
                    new_fname = new_fname.replace(old, new)

            if new_fname != fname:
                old_path = os.path.join(dirpath, fname)
                new_path = os.path.join(dirpath, new_fname)
                os.rename(old_path, new_path)
                print("File renamed: {} to {}".format(fname, new_fname))

        # Rename directories
        for dname in dirnames:
            new_dname = dname
            for old, new in full_mapping.items():
                if old in new_dname:
                    new_dname = new_dname.replace(old, new)

            if new_dname != dname:
                old_path = os.path.join(dirpath, dname)
                new_path = os.path.join(dirpath, new_dname)
                os.rename(old_path, new_path)
                print("Directory renamed: {} to {}".format(dname, new_dname))

print("Finished: all subject folders and internal files renamed directly in sourcedata.")

# Permissions
print("Changing files permissions in {}/{}".format(main_dir, project_dir))
os.system("chmod -Rf 771 {}/{}".format(main_dir, project_dir))
os.system("chgrp -Rf {} {}/{}".format(group, main_dir, project_dir))