"""
-----------------------------------------------------------------------------------------
deface_sbatch.py
-----------------------------------------------------------------------------------------
Goal of the script:
Deface anatomical images
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject (e.g. sub-001)
sys.argv[4]: group (e.g. 327)
sys.argv[5]: server job or not (1 = server, 0 = terminal)
sys.argv[6]: server project (e.g. b327)
-----------------------------------------------------------------------------------------
Output(s):
Defaced images, originals are copied
-----------------------------------------------------------------------------------------
To run: run python commands
>> cd ~/projects/pRF_analysis/analysis_code/preproc/bids/
>> python deface_sbatch.py [main directory] [project name] [subject num] [group] 
                           [server_in] [server_project]
--------------------------------------------------------------------------------------------
Exemple:
cd ~/projects/pRF_analysis/analysis_code/preproc/bids/
python deface_sbatch.py /scratch/mszinte/data centbids sub-t043 327 1 b327
-----------------------------------------------------------------------------------------
Written by Martin Szinte (martin.szinte@gmail.com)
adapted by Sina Kling (sina.kling@outlook.de)
-----------------------------------------------------------------------------------------
"""
import os
import sys
import shutil
import time
from bids import BIDSLayout
import pdb
deb = pdb.set_trace


# -------------------- INPUTS --------------------
if len(sys.argv) < 5:
    print("Usage: python deface_sbatch_flex.py <main_dir> <project_name> <subject> <server_in>")
    sys.exit(1)

main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3] if sys.argv[3].startswith('sub-') else f"sub-{sys.argv[3]}"
group = sys.argv[4]
server_in = int(sys.argv[5])
server_project = sys.argv[6]

# Paths
bids_root = os.path.join(main_dir, project_dir)
log_dir = os.path.join(bids_root, "derivatives", "pp_data", "logs")
sh_folder = os.path.join(bids_root, "derivatives", "pp_data", "jobs")
sourcedata_dir = os.path.join(bids_root, "sourcedata")

os.makedirs(log_dir, exist_ok=True)
os.makedirs(sh_folder, exist_ok=True)
os.makedirs(sourcedata_dir, exist_ok=True)

# -------------------- FIND ANATOMY FILES --------------------
anat_files = []

layout = BIDSLayout(bids_root, validate=False, derivatives=False)
files = layout.get(subject=subject.replace("sub-", ""), datatype="anat", return_type="filename")
anat_files = [f for f in files if f.endswith((".nii", ".nii.gz"))]


# fallback to glob if pybids not available or failed
if not anat_files:
    for root, dirs, files in os.walk(os.path.join(bids_root, subject)):
        for f in files:
            if f.endswith((".nii", ".nii.gz")) and "anat" in root:
                anat_files.append(os.path.join(root, f))

if not anat_files:
    print(f"No anatomical files found for {subject}. Exiting.")
    sys.exit(0)

print(f"Found {len(anat_files)} anatomical file(s).")

#  -------------------- COPY TO SOURCEDATA --------------------
for f in anat_files:
    rel_path = os.path.relpath(f, bids_root)
    dst_path = os.path.join(sourcedata_dir, rel_path)
    os.makedirs(os.path.dirname(dst_path), exist_ok=True)

    if os.path.exists(dst_path):
        print(f"File already exists in sourcedata, skipping: {dst_path}")
    else:
        shutil.copy2(f, dst_path)
        print(f"Copied {f} -> {dst_path}")

# -------------------- CREATE SH FILE --------------------
sh_file = os.path.join(sh_folder, f"{subject}_deface.sh")

# define FSL comand
fsl_cmd = """\
export FSLDIR='{main_dir}/{project_dir}/code/fsl' export PATH=$PATH:$FSLDIR/bin source $FSLDIR/etc/fslconf/fsl.sh\n""".format(main_dir=main_dir, project_dir=project_dir)

# define change mode change group comand
chmod_cmd = """chmod -Rf 771 {main_dir}/{project_dir}\n""".format(main_dir=main_dir, project_dir=project_dir) 
chgrp_cmd = """chgrp -Rf {group} {main_dir}/{project_dir}\n""".format(main_dir=main_dir, project_dir=project_dir, group=group)

with open(sh_file, "w") as f:
    if server_in:
        f.write("#!/bin/bash\n")
        f.write(f"#SBATCH -p skylake\n")
        f.write(f"#SBATCH -A {server_project}\n")
        f.write(f"#SBATCH --nodes=1\n")
        f.write(f"#SBATCH --cpus-per-task=8\n")
        f.write(f"#SBATCH --time=04:00:00\n")
        f.write(f"#SBATCH -e {log_dir}/{subject}_deface_%N_%j.err\n")
        f.write(f"#SBATCH -o {log_dir}/{subject}_deface_%N_%j.out\n")
        f.write(f"#SBATCH -J {subject}_deface\n\n")
    else:
        f.write("#!/bin/bash\n\n")

    # FSL environment
    f.write(fsl_cmd)

    # pydeface
    for anat in anat_files:
        f.write(f"pydeface {anat} --verbose\n")

    # permissions
    f.write(chmod_cmd)
    f.write(chgrp_cmd)


os.chmod(sh_file, 0o755)
print(f"Deface shell script written to {sh_file}")


# -------------------- RUN OR SUBMIT --------------------
if server_in:
    print(f"Submitting {sh_file} via sbatch...")
    os.system(f"sbatch {sh_file}")
else:
    print(f"Running {sh_file} locally...")
    os.system(f"sh {sh_file}")
