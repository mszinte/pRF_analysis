"""
-----------------------------------------------------------------------------------------
submit_workbench_compute_full_corr_job.py
-----------------------------------------------------------------------------------------
Goal of the script:
Submit the compute dtseries correlation scripts to SLURM queue
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (corresponds to directory)
sys.argv[3]: group (e.g. 327)
sys.argv[4]: server project (e.g. b327)
sys.argv[5]: runtime (e.g. 1:00:00)
-----------------------------------------------------------------------------------------
Output(s):
.sh file to execute in server
-----------------------------------------------------------------------------------------
Example:
conda activate pRF_env
cd projects/pRF_analysis/RetinoMaps/rest/correlations/
python submit_workbench_compute_full_corr_job.py /scratch/mszinte/data RetinoMaps 327 b327 1:00:00
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
group = sys.argv[3]
server_project = sys.argv[4]
proc_time = sys.argv[5] # usually an hour is enough

memory_val = 20 # GB
nb_procs = 16 # number of CPUs
cluster_name = 'skylake'

# Define your command
script_name = "connectome-workbench_dense_full_corr_bilateral.sh"
job_suffix = "workbench"
corr_type = "full_corr"

# Define permission commands
chmod_cmd = "chmod -Rf 771 {}/{}".format(main_dir, project_dir)
chgrp_cmd = "chgrp -Rf {} {}/{}".format(group, main_dir, project_dir)

# Make directories
logs_dir = f'{main_dir}/{project_dir}/derivatives/pp_data/log_outputs'
os.makedirs(logs_dir, exist_ok=True)

job_dir = f'{main_dir}/{project_dir}/derivatives/pp_data/jobs'
os.makedirs(job_dir, exist_ok=True)

# SLURM header
slurm_cmd = f"""#!/bin/bash
#SBATCH -p {cluster_name}
#SBATCH -A {server_project}
#SBATCH --nodes=1
#SBATCH --mem={memory_val}gb
#SBATCH --cpus-per-task={nb_procs}
#SBATCH --time={proc_time}
#SBATCH -e {logs_dir}/job_%N_%j_%a.err
#SBATCH -o {logs_dir}/job_%N_%j_%a.out
#SBATCH -J all_{job_suffix}
"""

# Assuming scripts are in the same directory as this Python script
script_dir = os.path.dirname(os.path.abspath(__file__))
your_script_path = os.path.join(script_dir, script_name)

# Check if script exists
if not os.path.exists(your_script_path):
    print(f"Error: Script {your_script_path} not found!")
    sys.exit(1)

main_cmd = "bash {}".format(your_script_path)
chmod_cmd = f"chmod -Rf 771 {main_dir}/{project_dir}"
chgrp_cmd = f"chgrp -Rf {group} {main_dir}/{project_dir}"

# Create job script
sh_fn = f"{job_dir}/all_{job_suffix}.sh"
with open(sh_fn, 'w') as of:
    of.write(f"{slurm_cmd}\n{main_cmd}\n{chmod_cmd}\n{chgrp_cmd}\n")

# Submit job
print("Submitting {} ({}) to queue".format(sh_fn, corr_type))
os.system("sbatch {}".format(sh_fn))
