"""
-----------------------------------------------------------------------------------------
preproc_end_sbatch.py
-----------------------------------------------------------------------------------------
Goal of the script:
Run preproc end on mesocenter 
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject name (e.g. sub-01)
sys.argv[4]: group (e.g. 327)
sys.argv[5]: server project (e.g. b327)
-----------------------------------------------------------------------------------------
Output(s):
sh file for running batch command
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/[PROJECT]/analysis_code/preproc/functional
2. run python command
>> python preproc_end_sbatch.py [main directory] [project] [subject] [group [server num]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/projects/pRF_analysis/analysis_code/preproc/functional
python preproc_end_sbatch.py /scratch/mszinte/data RetinoMaps sub-01 327 b327
python preproc_end_sbatch.py /scratch/mszinte/data centbids sub-2100247523 327 b327
-----------------------------------------------------------------------------------------
Written by Martin Szinte (martin.szinte@gmail.com)
and Uriel Lascombes (uriel.lascombes@laposte.net)
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
sys.path.append("{}/../../utils".format(os.getcwd()))
from settings_utils import load_settings

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
group = sys.argv[4]
server_project = sys.argv[5]
memory_val = 48
nb_procs = 8
hour_proc = 10

# Load settings
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../"))
settings_path = os.path.join(base_dir, project_dir, "settings.yml")
prf_settings_path = os.path.join(base_dir, project_dir, "prf-analysis.yml")
settings = load_settings([settings_path, prf_settings_path])
analysis_info = settings[0]

cluster_name  = analysis_info['cluster_name']

# Set folders
log_dir = "{}/{}/derivatives/pp_data/{}/log_outputs".format(main_dir, project_dir, subject)
job_dir = "{}/{}/derivatives/pp_data/{}/jobs".format(main_dir, project_dir, subject)
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
#SBATCH -e {log_dir}/{subject}_preproc_end_%N_%j_%a.err
#SBATCH -o {log_dir}/{subject}_preproc_end_%N_%j_%a.out
#SBATCH -J {subject}_preproc_end
""".format(server_project=server_project, cluster_name=cluster_name,
           nb_procs=nb_procs, hour_proc=hour_proc, 
           subject=subject, memory_val=memory_val, log_dir=log_dir)
    
preproc_end_surf_cmd = "python preproc_end.py {} {} {} {}".format(main_dir, project_dir, subject, group)

# Create sh fn
sh_fn = "{}/{}_preproc_end.sh".format(job_dir, subject)

of = open(sh_fn, 'w')
of.write("{} \n{}".format(slurm_cmd, preproc_end_surf_cmd))
of.close()

# Submit jobs
print("Submitting {} to queue".format(sh_fn))
os.system("sbatch {}".format(sh_fn))