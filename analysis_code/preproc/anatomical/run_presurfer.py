#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 27 11:48:39 2023
-----------------------------------------------------------------------------------------
run_presurfer.py
-----------------------------------------------------------------------------------------
Goal of the script:
Create bids directory from tarfiles on mesocentre
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: mesocentre project data directory 
sys.argv[2]: project name
sys.argv[3]: individual subject name
-----------------------------------------------------------------------------------------
To run:
On mesocentre
>> conda activate snakemake3.9
>> cd ~/projects/bio7/analysis_code/preproc/bids
>> python run_presurfer.py [meso_proj_dir] [subject_number] 
python run_presurfer.py /scratch/mszinte/data/ centbids sub-t044
----------------------------------------------------------------------------------------
@author: penny
adapted by Sina Kling (sina.kling@outlook.de) 
----------------------------------------------------------------------------------------
"""

import os
import sys
import shutil

#conda config --add channels conda-forge
#conda config --add channels bioconda
#conda config --set channel_priority strict
#conda create -c conda-forge -c bioconda -n snakemake3.9 snakemake python=3.9
#conda activate snakemake3.9 
#pip install git+https://github.com/khanlab/snakebids.git



#ADAPT config (snakebids.yml) before: add correct bids folder and output folder

# Define directories 
main_dir = sys.argv[1]
project_name = sys.argv[2]
project_dir = os.path.join(main_dir,  project_name)
subject = sys.argv[3]
sub_num = subject[4:] 
print(f"processing subject: {sub_num}")

# Make all directories
presurfer_dir = os.path.join(project_dir, "derivatives/presurfer")
os.makedirs(presurfer_dir, exist_ok=True)

output_dir = presurfer_dir
print(f"saving data to: {output_dir}")

sh_folder = os.path.join(presurfer_dir, "jobs")
os.makedirs(sh_folder, exist_ok=True)

sh_file = os.path.join(sh_folder, f"{subject}_presurf.sh")

# Get directory of the current script
script_dir = os.path.dirname(os.path.realpath(__file__))


# Write shell script
presurfer_outputs = os.path.join(presurfer_dir, "log_outputs")
os.makedirs(presurfer_outputs, exist_ok=True)

org_script = os.path.join(script_dir, "organise_presurfer.py")
org_cmd = f"python {org_script} {project_dir} {subject}\n"

with open(sh_file, 'w') as f:
    pres_cmd = (
        f"python {project_dir}/code/presurfer-smk-master/presurfer/run.py "
        f"{project_dir} {output_dir} participant --participant_label {sub_num} "
        f"--cores 1 --use-singularity --force-output\n"  #--force-output
    )
    f.write(pres_cmd)
    f.write(org_cmd)

print(f"Created {sh_file} with the following commands:\n")
print(open(sh_file).read())


# Run the script directly
os.system(f"sh {sh_file}")

cache_dir = os.path.expanduser("~/.cache/snakemake/")
if os.path.exists(cache_dir):
    shutil.rmtree(cache_dir)
