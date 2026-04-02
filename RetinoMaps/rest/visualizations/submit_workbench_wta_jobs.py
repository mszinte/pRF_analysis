"""
-----------------------------------------------------------------------------------------
submit_workbench_wta_jobs.py
-----------------------------------------------------------------------------------------
Goal of the script:
Submit the compute winner take all scripts to SLURM queue
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (corresponds to directory)
sys.argv[3]: group (e.g. 327)
sys.argv[4]: server project (e.g. b327)
sys.argv[5]: script name to run
sys.argv[6]: runtime (e.g. 0:30:00)
sys.argv[7:] flexible options to compute wta on different outputs
-----------------------------------------------------------------------------------------
Output(s):
.sh file to execute in server
-----------------------------------------------------------------------------------------
Example:
source .bashrc
conda activate pRF_env
cd projects/pRF_analysis/RetinoMaps/rest/visualizations
python submit_workbench_wta_jobs.py /scratch/mszinte/data RetinoMaps 327 b327 workbench_compute_wta_full_corr_dev.sh 0:30:00 fisher-z by_hemi legacy
-----------------------------------------------------------------------------------------
Written by Marco Bedini (marco.bedini@univ-amu.fr)
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
script_name = sys.argv[5]
proc_time = sys.argv[6]
extra_args = sys.argv[7:]

memory_val = 20 # GB
nb_procs = 16 # number of CPUs
cluster_name = 'skylake'
job_suffix = "wta_workbench"

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
#SBATCH -J {job_suffix}
"""

# Assuming scripts are in the same directory as this Python script
script_dir = os.path.dirname(os.path.abspath(__file__))
your_script_path = os.path.join(script_dir, script_name)

# Check if script exists
if not os.path.exists(your_script_path):
    print(f"Error: Script {your_script_path} not found!")
    sys.exit(1)

# Add flexible options for wta (fisher-z by hemi legacy)
args_str = " ".join(extra_args)
main_cmd = f"bash {your_script_path} {args_str}"

chmod_cmd = f"chmod -Rf 771 {main_dir}/{project_dir}"
chgrp_cmd = f"chgrp -Rf {group} {main_dir}/{project_dir}"

# Create job script
sh_fn = f"{job_dir}/{job_suffix}.sh"
with open(sh_fn, 'w') as of:
    of.write(f"{slurm_cmd}\n{main_cmd}\n{chmod_cmd}\n{chgrp_cmd}\n")

# Submit job
print("Submitting ({}) to queue".format(sh_fn))
os.system("sbatch {}".format(sh_fn))

