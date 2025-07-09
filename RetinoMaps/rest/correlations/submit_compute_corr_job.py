"""
-----------------------------------------------------------------------------------------
submit_compute_corr_job.py
-----------------------------------------------------------------------------------------
Goal of the script:
Submit the compute dtseries correlation scripts to SLURM queue
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (corresponds to directory)
sys.argv[3]: subject name (e.g. sub-01)
sys.argv[4]: group (e.g. 327)
sys.argv[5]: server project (e.g. b327)
sys.argv[6]: correlation type ('full' or 'fisher-z')
-----------------------------------------------------------------------------------------
Output(s):
.sh file to execute in server
-----------------------------------------------------------------------------------------
Example:
cd ~/projects/pRF_analysis/analysis_code/postproc/rest/correlations/
python submit_compute_corr_job.py /scratch/mszinte/data RetinoMaps sub-01 327 b327 full
python submit_compute_corr_job.py /scratch/mszinte/data RetinoMaps sub-01 327 b327 fisher-z
-----------------------------------------------------------------------------------------
Written by Marco Bedini (marco.bedini@univ-amu.fr) adapting previous examples in the RetinoMaps project
"""

# Stop warnings
import warnings
warnings.filterwarnings("ignore")

# General imports
import os
import sys

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
group = sys.argv[4]
server_project = sys.argv[5]
corr_type = sys.argv[6]  # 'full' or 'fisher-z'

# Validate correlation type
if corr_type not in ['full', 'fisher-z']:
    print("Error: correlation type must be 'full' or 'fisher-z'")
    sys.exit(1)

# Job parameters
memory_val = 16        # GB
hour_proc = 2          # hours
nb_procs = 32          # number of CPUs
cluster_name = 'skylake'

# Define your actual command based on correlation type
if corr_type == 'full':
    script_name = "compute_dtseries_corr_bilateral.sh"
    job_suffix = "full"
elif corr_type == 'fisher-z':
    script_name = "compute_dtseries_corr_fisher-z_bilateral.sh"
    job_suffix = "fisher_z"

# Define permission commands
chmod_cmd = "chmod -Rf 771 {}/{}".format(main_dir, project_dir)
chgrp_cmd = "chgrp -Rf {} {}/{}".format(group, main_dir, project_dir)

# Make directories 
logs_dir =  f'{main_dir}/{project_dir}/derivatives/pp_data/{subject}/91k/rest/log_outputs'
os.makedirs(logs_dir, exist_ok=True)
job_dir =  f'{main_dir}/{project_dir}/derivatives/pp_data/{subject}/91k/rest/jobs'
os.makedirs(job_dir, exist_ok=True)

# SLURM header
slurm_cmd = """\
#!/bin/bash
#SBATCH -p {cluster_name}
#SBATCH -A {server_project}
#SBATCH --nodes=1
#SBATCH --mem={memory_val}gb
#SBATCH --cpus-per-task={nb_procs}
#SBATCH --time={hour_proc}:00:00
#SBATCH -e {logs_dir}/{subject}_job_%N_%j_%a.err
#SBATCH -o {logs_dir}/{subject}_job_%N_%j_%a.out
#SBATCH -J {subject}_{job_suffix}
""".format(server_project=server_project, 
           cluster_name=cluster_name,
           nb_procs=nb_procs, 
           hour_proc=hour_proc, 
           subject=subject, 
           memory_val=memory_val, 
           logs_dir=logs_dir,
           job_suffix=job_suffix)

# Assuming scripts are in the same directory as this Python script
script_dir = os.path.dirname(os.path.abspath(__file__))
your_script_path = os.path.join(script_dir, script_name)

# Check if script exists
if not os.path.exists(your_script_path):
    print(f"Error: Script {your_script_path} not found!")
    sys.exit(1)

main_cmd = "bash {}".format(your_script_path)

# Create job script
sh_fn = "{}/{}_{}.sh".format(job_dir, subject, job_suffix)
of = open(sh_fn, 'w')
of.write("{} \n{} \n{} \n{}".format(slurm_cmd, main_cmd, chmod_cmd, chgrp_cmd))
of.close()

# Submit job
print("Submitting {} ({}) to queue".format(sh_fn, corr_type))
os.system("sbatch {}".format(sh_fn))
