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
import json
import subprocess

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
group = sys.argv[3]

# Load settings and analysis parameters
with open('../../../settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)

subject_names = analysis_info['subject_names']
session_names = analysis_info['session_names']

# Define directories
sourcedata_dir = "{}/{}/sourcedata".format(main_dir, project_dir)
bids_dir = "{}/{}".format(main_dir, project_dir)

# Function: fast copy using rsync -azuv
def fast_copy(src, dst):
    """
    Copy a folder tree quickly using rsync with archive, compression, update, verbose.
    """
    if not os.path.exists(dst):
        os.makedirs(dst)
    cmd = ["rsync", "-azuv", src + "/", dst + "/"]
    subprocess.run(cmd, check=True)

# Copy data from sourcedata_dir to bids_dir
created_folders = []

for old_name, new_name in subject_names.items():
    src = os.path.join(sourcedata_dir, old_name)
    dst = os.path.join(bids_dir, new_name)

    if not os.path.exists(src):
        print("Source not found: {}".format(src))
        continue

    if os.path.exists(dst):
        print("Destination already exists: {} — skipping copy.".format(dst))
        continue

    print("Copying {} → {}".format(src, dst))
    fast_copy(src, dst)
    created_folders.append(dst)

# Rename all files and folders inside the copied folders (subject + session)
full_mapping = {**subject_names, **session_names}

for folder in created_folders:
    print("Renaming inside folder: {}".format(folder))

    # Walk through the folder tree from bottom to top
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
                print("File renamed: {} → {}".format(fname, new_fname))

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
                print("Directory renamed: {} → {}".format(dname, new_dname))

print("\nFinished: all folders and files copied and renamed (subject + session).")

# # Define permission cmd
# print('Changing files permissions in {}/{}'.format(main_dir, project_dir))
# os.system("chmod -Rf 771 {}/{}".format(main_dir, project_dir))
# os.system("chgrp -Rf {} {}/{}".format(group, main_dir, project_dir))