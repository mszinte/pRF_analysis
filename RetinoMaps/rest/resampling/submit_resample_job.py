"""
-----------------------------------------------------------------------------------------
submit_resample_job.py
-----------------------------------------------------------------------------------------
Goal of the script:
Submit the resample from 170k to 91k scripts to SLURM queue
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (corresponds to directory)
sys.argv[3]: group (e.g. 327)
sys.argv[4]: server project (e.g. b327)
sys.argv[5]: script name to run:
sys.argv[6]: runtime (e.g. 1:00:00)
-----------------------------------------------------------------------------------------
Output(s):
.sh file to execute in server
-----------------------------------------------------------------------------------------
Example:
cd projects/pRF_analysis/RetinoMaps/rest/resampling
python submit_resample_job.py /scratch/mszinte/data RetinoMaps 327 b327 adap-bary 2:00:00
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
mode = sys.argv[5]   # "adap-bary" or "bary"
proc_time = sys.argv[6] #

# Define which way you want to resample
# Adap-bary is recommended for this scenario
if mode == "adap-bary":
    script_name = "resample_to_fsLR91k_adap-bary.sh -largest"
elif mode == "bary":
    script_name = "resample_to_fsLR91k_bary.sh"
else:
    raise ValueError(
        "Invalid mode. Use either 'adap-bary' or 'bary'. "
        f"Got: {mode}"
    )

memory_val = 20 # GB
nb_procs = 16 # number of CPUs
cluster_name = 'skylake'
job_suffix = "workbench"

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

main_cmd = "bash {}".format(your_script_path)
chmod_cmd = f"chmod -Rf 771 {main_dir}/{project_dir}"
chgrp_cmd = f"chgrp -Rf {group} {main_dir}/{project_dir}"

# Create job script
sh_fn = f"{job_dir}/resample_to_91k_all_{job_suffix}.sh"
with open(sh_fn, 'w') as of:
    of.write(f"{slurm_cmd}\n{main_cmd}\n{chmod_cmd}\n{chgrp_cmd}\n")

# Submit job
print("Submitting {} to queue".format(sh_fn))
os.system("sbatch {}".format(sh_fn))
