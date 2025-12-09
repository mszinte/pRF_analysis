"""
-----------------------------------------------------------------------------------------
copy_events.py
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
python copy_events.py [nCSF data dir] [pRF data dir] [DeepMReyeCalib data dir] [BIDS dataset dir] [group]
-----------------------------------------------------------------------------------------
Example:
cd ~/projects/pRF_analysis/nCSF/MRI/preproc/bids
python copy_events.py ~/projects/nCSFexp/experiment_code/data ~/projects/prfexp7t/data/nCSF ~/projects/deepmreyecalib/experiment_code/data /scratch/mszinte/data/nCSF 327
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
import json
import shutil

# Inputs
ncsf_dir = sys.argv[1]
prf_dir = sys.argv[2]
deepm_dir = sys.argv[3]
bids_dir = sys.argv[4]
group = sys.argv[5]

# Load settings
with open('../../../settings.json') as f:
    analysis_info = json.load(f)

subjects = analysis_info['subjects']
sessions = analysis_info['sessions']

# Mapping BIDS task -> source task folder
task_mapping = {
    "nCSF": {"folder": ncsf_dir, "source_task": "nCSF"},
    "pRF": {"folder": prf_dir, "source_task": "pRF"},
    "DeepMReye": {"folder": deepm_dir, "source_task": "DeepMReyeCalib7t"}
}

def find_source_file(src_dir, source_task, sub, ses, run_number):
    """Return the source events file path matching subject/session/run"""
    func_dir = os.path.join(src_dir, sub, ses, 'func')
    files = [f for f in os.listdir(func_dir) if source_task in f and f.endswith('_events.tsv') and 'Ams' not in f]
    files.sort()
    for f in files:
        if '_run-{0}_'.format(run_number) in f:
            return os.path.join(func_dir, f)
    return None

for sub in subjects:
    for ses in sessions:
        func_dir = os.path.join(bids_dir, sub, ses, 'func')
        bids_files = [f for f in os.listdir(func_dir) if f.endswith('_events.tsv')]
        bids_files.sort()
        for bids_file in bids_files:
            for task_dir in task_mapping:
                if task_dir in bids_file:
                    run_number = bids_file.split('_run-')[1].split('_')[0]
                    src_info = task_mapping[task_dir]
                    src_file = find_source_file(src_info["folder"], src_info["source_task"], sub, ses, run_number)
                    if src_file:
                        shutil.copy2(src_file, os.path.join(func_dir, bids_file))
                        print("Filled {0} with {1}".format(bids_file, os.path.basename(src_file)))
                    else:
                        print("WARNING: source file not found for {0} run {1}".format(task_dir, run_number))

# Change permissions and group
print('Changing files permissions in {}'.format(bids_dir))
os.system("chmod -Rf 771 {}".format(bids_dir))
os.system("chgrp -Rf {} {}".format(group, bids_dir))