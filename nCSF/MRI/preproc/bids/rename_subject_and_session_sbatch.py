"""
-----------------------------------------------------------------------------------------
rename_subject_and_session_sbatch.py
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
python rename_subject_and_session_sbatch.py /scratch/mszinte/data nCSF 327 b327
-----------------------------------------------------------------------------------------
Written by Uriel Lascombes (uriel.lascombes@laposte.net)
Edited by Martin Szinte (mail@martinszinte.net)
-----------------------------------------------------------------------------------------
"""

# Stop warnings
import warnings
warnings.filterwarnings("ignore")

# debug 
import ipdb 
deb = ipdb.set_trace

# general imports
import os
import sys
import json

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
group = sys.argv[3]
server_project = sys.argv[4]

# Load settings
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../../"))
settings_path = os.path.join(base_dir, "nCSF", "settings.yml")
settings = load_settings([settings_path])
analysis_info = settings[0]
cluster_name  = analysis_info['cluster_name']

nb_procs = 1
memory_val = 48
hour_proc = 2

# set folders
log_dir = "{}/{}/derivatives/log_outputs".format(main_dir, project_dir)
job_dir = "{}/{}/derivatives/jobs".format(main_dir, project_dir)
os.makedirs(log_dir, exist_ok=True)
os.makedirs(job_dir, exist_ok=True)

slurm_cmd = """\
#!/bin/bash
#SBATCH -p {cluster_name}
#SBATCH -A {server_project}
#SBATCH --nodes=1
#SBATCH --mem={memory_val}gb
#SBATCH --cpus-per-task={nb_procs}
#SBATCH --time={hour_proc}:00:00
#SBATCH -e {log_dir}/rename_subj_ses.err
#SBATCH -o {log_dir}/rename_subj_ses.out
#SBATCH -J rename_subj_ses
""".format(
server_project=server_project, cluster_name=cluster_name, nb_procs=nb_procs, 
hour_proc=hour_proc, memory_val=memory_val, log_dir=log_dir)

rename_subj_ses_cmd = "python rename_subject_and_session.py {} {} {}".format(main_dir, project_dir, group)

# create sh fn
sh_fn = "{}/rename_subj_ses.sh".format(job_dir)
of = open(sh_fn, 'w')

of.write("{} \n{}".format(slurm_cmd, rename_subj_ses_cmd))
of.close()

# Submit jobs
print("Submitting {} to queue".format(sh_fn))
os.system("sbatch {}".format(sh_fn))