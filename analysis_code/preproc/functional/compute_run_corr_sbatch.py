"""
-----------------------------------------------------------------------------------------
compute_run_corr_sbatch.py
-----------------------------------------------------------------------------------------
Goal of the script:
Run computation of correlation inter-run 
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
>> python compute_run_corr_sbatch.py [main directory] [project] [subject] [group] 
                                                                            [server num]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/projects/pRF_analysis/analysis_code/preproc/functional

python compute_run_corr_sbatch.py /scratch/mszinte/data MotConf sub-01 327 b327

python compute_run_corr_sbatch.py /scratch/mszinte/data RetinoMaps sub-01 327 b327

python compute_run_corr_sbatch.py /scratch/mszinte/data amblyo_prf sub-01 327 b327

python compute_run_corr_sbatch.py /scratch/mszinte/data centbids sub-2100247523 327 b327
-----------------------------------------------------------------------------------------
Written by Martin Szinte (martin.szinte@gmail.com)
Edited by Uriel Lascombes (uriel.lascombes@laposte.net)
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

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
group = sys.argv[4]
server_project = sys.argv[5]

# Define analysis parameters
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../"))
settings_path = os.path.join(base_dir, project_dir, "settings.json")

with open(settings_path) as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)

# Define cluster/server specific parameters
cluster_name  = analysis_info['cluster_name']
proj_name = analysis_info['project_name']
nb_procs = 1
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
#SBATCH -e {log_dir}/{subject}_run_corr_%N_%j_%a.err
#SBATCH -o {log_dir}/{subject}_run_corr_%N_%j_%a.out
#SBATCH -J {subject}_run_corr
""".format(server_project=server_project, cluster_name=cluster_name,
           nb_procs=nb_procs, hour_proc=hour_proc, 
           subject=subject, memory_val=memory_val, log_dir=log_dir)
    
run_corr_surf_cmd = "python compute_run_corr.py {} {} {} {}".format(main_dir, project_dir, subject, group)

# Create sh fn
sh_fn = "{}/{}_run_corr.sh".format(job_dir, subject)

of = open(sh_fn, 'w')
of.write("{} \n{}".format(slurm_cmd, run_corr_surf_cmd))
of.close()

# Submit jobs
print("Submitting {} to queue".format(sh_fn))
os.system("sbatch {}".format(sh_fn))