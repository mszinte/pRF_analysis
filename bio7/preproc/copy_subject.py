"""
-----------------------------------------------------------------------------------------
copy_subject.py
-----------------------------------------------------------------------------------------
Goal of the script:
Copy subject from RetinoMaps to RetinoMapsOpenNeuro. Copy only the usefull task. 
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: origin directory 
sys.argv[2]: destination directory
sys.argv[3]: subject name
sys.argv[4]: tasks to copy
sys.argv[5]: group of shared data (e.g. 327)
-----------------------------------------------------------------------------------------
Output(s):
Copy subject from RetinoMaps to RetinoMapsOpenNeuro. Copy only the usefull task. 
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/pRF_analysis/RetinoMaps/preproc
2. run python command
python copy_events_and_matfiles.py [origin dir] [destination dir] [subject to copy] [list of task to copy] [group]
-----------------------------------------------------------------------------------------
Example:
cd ~/projects/pRF_analysis/RetinoMaps/preproc
python copy_subject.py /scratch/mszinte/data/RetinoMaps /scratch/mszinte/data/RetinoMapsOpenNeuro sub-01 [pRF,SacLoc,PurLoc,rest] 327
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

# Inputs
ori_dir = sys.argv[1]
dest_dir = sys.argv[2]
subject = sys.argv[3]
tasks_arg = sys.argv[4]
group = sys.argv[5]


# Parse tasks argument
tasks_arg = tasks_arg.strip("[]")
tasks = [t.strip() for t in tasks_arg.split(",")]

print("Tasks to copy: {}".format(tasks))

# Paths
src_sub = os.path.join(ori_dir, subject)
dst_sub = os.path.join(dest_dir, subject)

os.makedirs(dst_sub, exist_ok=True)

# Loop over sessions
sessions = sorted(glob.glob(os.path.join(src_sub, "ses-*")))

if len(sessions) == 0:
    raise RuntimeError("No ses-* found for {}".format(subject))

for ses_path in sessions:

    ses = os.path.basename(ses_path)
    print("Processing {}".format(ses))

    src_ses = os.path.join(src_sub, ses)
    dst_ses = os.path.join(dst_sub, ses)
    os.makedirs(dst_ses, exist_ok=True)

    # Copy everything except func
    for item in os.listdir(src_ses):

        if item == "func" or item == "dwi" or item.startswith("."):
            continue

        src = os.path.join(src_ses, item)
        dst = os.path.join(dst_ses, item)

        if os.path.isdir(src):
            subprocess.run(
                ["rsync", "-azv", src + "/", dst],
                check=True
            )
        else:
            subprocess.run(
                ["rsync", "-azv", src, dst],
                check=True
            )

    # Copy selected tasks in func
    src_func = os.path.join(src_ses, "func")
    dst_func = os.path.join(dst_ses, "func")
    os.makedirs(dst_func, exist_ok=True)

    for task in tasks:
        pattern = os.path.join(
            src_func,
            "{}_{}_task-{}_*".format(subject, ses, task)
        )

        files = glob.glob(pattern)

        if len(files) == 0:
            print("  ⚠️  No files for {} in {}".format(task, ses))
            continue

        for fpath in files:
            subprocess.run(
                ["rsync", "-azv", fpath, dst_func],
                check=True
            )

# Change permissions and group
print('Changing files permissions in {}'.format(dest_dir))
os.system("chmod -Rf 771 {}".format(dest_dir))
os.system("chgrp -Rf {} {}".format(group, dest_dir))
