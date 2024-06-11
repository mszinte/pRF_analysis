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
python preproc_end_sbatch.py /scratch/mszinte/data MotConf sub-01 327 b327
-----------------------------------------------------------------------------------------
Written by Martin Szinte (martin.szinte@gmail.com)
Edited by Uriel Lascombes (uriel.lascombes@laposte.net)
-----------------------------------------------------------------------------------------
"""

# Stop warnings
import warnings
warnings.filterwarnings("ignore")

# General imports
import json
import os
import sys
import ipdb
deb = ipdb.set_trace

# Define analysis parameters
with open('../../settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
group = sys.argv[4]
server_project = sys.argv[5]

# Define cluster/server specific parameters
cluster_name  = analysis_info['cluster_name']
proj_name = analysis_info['project_name']
nb_procs = 8
memory_val = 48
hour_proc = 10

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