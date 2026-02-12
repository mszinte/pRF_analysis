"""
-----------------------------------------------------------------------------------------
copy_events_and_matfiles.py
-----------------------------------------------------------------------------------------
Goal of the script:
Copy events files from the experiment data folders (nCSF, pRF, DeepMReyeCalib) to the BIDS dataset.
Matches runs correctly and handles differences in task names across datasets.
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: nCSF experiment data directory
sys.argv[2]: pRF experiment data directory
sys.argv[3]: DeepMReyeCalib experiment data directory
sys.argv[4]: BIDS dataset directory (meso_shared/nCSF)
sys.argv[5]: server group (e.g. 327)
-----------------------------------------------------------------------------------------
Output(s):
Copies and fills BIDS _events.tsv files for all subjects and sessions
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/pRF_analysis/nCSF/MRI/preproc/bids
2. run python command
python copy_events_and_matfiles.py [nCSF data dir] [pRF data dir] [DeepMReyeCalib data dir] [BIDS dataset dir] [group]
-----------------------------------------------------------------------------------------
Example:
cd ~/projects/pRF_analysis/nCSF/MRI/preproc/bids
python copy_events_and_matfiles.py ~/projects/nCSFexp/experiment_code/data ~/projects/prfexp7t/data/nCSF ~/projects/deepmreyecalib/experiment_code/data /scratch/mszinte/data/nCSF 327
-----------------------------------------------------------------------------------------
Written by Uriel Lascombes (uriel.lascombes@laposte.net)
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
import subprocess

# Personal imports
sys.path.append("{}/../../../utils".format(os.getcwd()))
from settings_utils import load_settings

# Inputs
ncsf_dir = sys.argv[1]
prf_dir = sys.argv[2]
deepm_dir = sys.argv[3]
bids_dir = sys.argv[4]
group = sys.argv[5]

# Load settings
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../../"))
settings_path = os.path.join(base_dir, "nCSF", "settings.yml")
settings = load_settings([settings_path])
analysis_info = settings[0]

subjects = analysis_info['subjects']
sessions = analysis_info['sessions']

for subject in subjects:
    for session in sessions:

        prf_func_dir = '{}/{}/{}/func'.format(prf_dir, subject, session)
        deepmreye_func_dir = '{}/{}/{}/func'.format(deepm_dir, subject, session)
        ncsf_func_dir = '{}/{}/{}/func'.format(ncsf_dir, subject, session)

        prf_add_dir = '{}/{}/{}/add'.format(prf_dir, subject, session)
        deepmreye_add_dir = '{}/{}/{}/add'.format(deepm_dir, subject, session)
        ncsf_add_dir = '{}/{}/{}/add'.format(ncsf_dir, subject, session)

        funcs_dir = [prf_func_dir, deepmreye_func_dir, ncsf_func_dir]
        adds_dir = [prf_add_dir, deepmreye_add_dir, ncsf_add_dir]

        events_fn = [
            f for func_dir in funcs_dir
            for f in glob.glob('{}/*events.tsv'.format(func_dir))
        ]

        adds_fn = (
            [f for add_dir in adds_dir for f in glob.glob('{}/*mat*.mat'.format(add_dir))] +
            [f for func_dir in funcs_dir for f in glob.glob('{}/*mat*.mat'.format(func_dir))]
        )

        copied_files = events_fn + adds_fn

        task_mapping = {
            "DeepMReyeCalib7t": {
                "experiment_task": "DeepMReyeCalib7t",
                "bids_task": "DeepMReye",
                "scan_direction": "PA"
            },
            "pRF": {
                "experiment_task": "pRF",
                "bids_task": "pRF",
                "scan_direction": "PA"
            },
            "nCSF": {
                "experiment_task": "nCSF",
                "bids_task": "nCSF",
                "scan_direction": "PA"
            }
        }

        files_to_copy = []

        for src_file in copied_files:

            if not os.path.exists(src_file):
                print("Warning: source file does not exist : {}".format(src_file))
                continue

            dst_filename = os.path.basename(src_file)

            for task_info in task_mapping.values():
                dst_filename = dst_filename.replace(
                    "task-{}".format(task_info['experiment_task']),
                    "task-{}_dir-{}".format(
                        task_info['bids_task'],
                        task_info['scan_direction']
                    )
                )

            dst_filename = dst_filename.replace(
                "_matlab.mat",
                "_matFile.mat"
            )

            dst_file = "{}/sourcedata/{}/{}/func/{}".format(
                bids_dir, subject, session, dst_filename
            )

            files_to_copy.append((src_file, dst_file))

        for src, dst in files_to_copy:
            print("copying {} -> {}".format(src, dst))
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            subprocess.run(
                ["rsync", "-azv", src, dst],
                check=True
            )

# Change permissions and group
print('Changing files permissions in {}'.format(bids_dir))
os.system("chmod -Rf 771 {}".format(bids_dir))
os.system("chgrp -Rf {} {}".format(group, bids_dir))