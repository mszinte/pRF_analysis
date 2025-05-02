"""
-----------------------------------------------------------------------------------------
submit_stats_corr_jobs.py
-----------------------------------------------------------------------------------------
Goal of the script:
Create and submit jobscript to make a gaussian grid fit for pRF analysis
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject name (e.g. sub-01)
sys.argv[4]: group (e.g. 327)
sys.argv[5]: server project (e.g. b327)
-----------------------------------------------------------------------------------------
Output(s):
.sh file to execute in server
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/[PROJECT]/analysis_code/postproc/rest/
2. run python command
python submit_stats_corr_jobs.py [main directory] [project name] [subject] 
                                 [group] [server project]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/projects/pRF_analysis/analysis_code/postproc/rest/
python submit_stats_corr_jobs.py /scratch/mszinte/data RetinoMaps sub-01 327 b327
-----------------------------------------------------------------------------------------
Written by Martin Szinte (martin.szinte@gmail.com)
Edited by Marco Bedini (marco.bedini@univ-amu.fr)
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
import json

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
group = sys.argv[4]
server_project = sys.argv[5]
memory_val = 30
hour_proc = 1
nb_procs = 8

# Cluster settings
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../"))
settings_path = os.path.join(base_dir, project_dir, "settings.json")
with open(settings_path) as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)
cluster_name  = analysis_info['cluster_name']
task_names = analysis_info['task_names']
task_name = task_names[0]

# Define permission cmd
chmod_cmd = "chmod -Rf 771 {}/{}".format(main_dir, project_dir)
chgrp_cmd = "chgrp -Rf {} {}/{}".format(group, main_dir, project_dir)

# Make directories 
logs_dir =  f'{main_dir}/{project_dir}/derivatives/pp_data/{subject}/91k/{task_name}/log_outputs'
os.makedirs(logs_dir, exist_ok=True)
job_dir =  f'{main_dir}/{project_dir}/derivatives/pp_data/{subject}/91k/{task_name}/jobs'
os.makedirs(job_dir, exist_ok=True)

slurm_cmd = """\
#!/bin/bash
#SBATCH -p {cluster_name}
#SBATCH -A {server_project}
#SBATCH --nodes=1
#SBATCH --mem={memory_val}gb
#SBATCH --cpus-per-task={nb_procs}
#SBATCH --time={hour_proc}:00:00
#SBATCH -e {logs_dir}/{subject}_stats_corr_%N_%j_%a.err
#SBATCH -o {logs_dir}/{subject}_stats_corr_%N_%j_%a.out
#SBATCH -J {subject}_rest_stats_corr
""".format(server_project=server_project, 
           cluster_name=cluster_name,
           nb_procs=nb_procs, 
           hour_proc=hour_proc, 
           subject=subject, 
           memory_val=memory_val, 
           logs_dir=logs_dir)

# Define stats_corr cmd
stats_corr_cmd = "python stats_corr.py {} {} {}".format(main_dir, project_dir, subject)

# Create sh
sh_fn = "{}/{}_stats_corr.sh".format(job_dir,subject)
of = open(sh_fn, 'w')
of.write("{} \n{} \n{} \n{}".format(slurm_cmd, stats_corr_cmd, chmod_cmd, chgrp_cmd))
of.close()

# Submit jobs
print("Submitting {} to queue".format(sh_fn))
os.system("sbatch {}".format(sh_fn))