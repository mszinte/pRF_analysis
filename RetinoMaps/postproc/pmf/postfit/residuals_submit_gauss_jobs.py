"""
-----------------------------------------------------------------------------------------
residuals_submit_gauss_jobs.py
-----------------------------------------------------------------------------------------
Goal of the script:
Create and submit jobscript to make a prf gaussian fit on PMF residuals
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
>> cd ~/projects/pRF_analysis/RetinoMaps/postproc/pmf/postfit
2. run python command
python residuals_submit_gauss_jobs.py [main directory] [project name] [subject]
                                [group] [server project]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/projects/pRF_analysis/RetinoMaps/postproc/pmf/fit
python residuals_submit_gauss_jobs.py /scratch/mszinte/data RetinoMaps sub-01 327 b327
-----------------------------------------------------------------------------------------
Written by Sina Kling (sina.kling@outlook.de)
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

# Personal imports
script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.abspath(os.path.join(script_dir, "../../../../analysis_code/utils")))
from pycortex_utils import set_pycortex_config_file
from settings_utils import load_settings

# Inputs
main_dir       = sys.argv[1]
project_dir    = sys.argv[2]
subject        = sys.argv[3]
group          = sys.argv[4]
server_project = sys.argv[5]

memory_val = 30
hour_proc  = 6
nb_procs   = 8

# Load settings
base_dir          = os.path.abspath(os.path.join(script_dir, "../../../../"))
settings_path     = os.path.join(base_dir, project_dir, "pmf-settings.yml")
prf_settings_path = os.path.join(base_dir, project_dir, "prf-analysis.yml")
settings          = load_settings([settings_path, prf_settings_path])
analysis_info     = settings[0]
task = 'SacLoc'

cluster_name = analysis_info['cluster_name']

# Set pycortex db and colormaps
cortex_dir = "{}/{}/derivatives/pp_data/cortex".format(main_dir, project_dir)
set_pycortex_config_file(cortex_dir)

# Directories
pp_dir = "{}/{}/derivatives/pp_data".format(main_dir, project_dir)

# Permission commands
chmod_cmd = "chmod -Rf 771 {}/{}".format(main_dir, project_dir)
chgrp_cmd = "chgrp -Rf {} {}/{}".format(group, main_dir, project_dir)

# Find residual files (fsnative .gii and 170k .nii)
residual_fns = []
residual_fns.extend(glob.glob(
    "{}/{}/fsnative/pmf/pmf_residuals/*_task-{}*.func.gii".format(pp_dir, subject, task)))
residual_fns.extend(glob.glob(
    "{}/{}/170k/pmf/pmf_residuals/*_task-{}*.dtseries.nii".format(pp_dir, subject, task)))

if not residual_fns:
    raise FileNotFoundError(
        f"No residual files found for {subject}.\n"
        f"  Checked fsnative: {pp_dir}/{subject}/fsnative/pmf/pmf_residuals/\n"
        f"  Checked 170k:     {pp_dir}/{subject}/170k/pmf/pmf_residuals/"
    )

print(f"Found {len(residual_fns)} residual file(s):")
for fn in residual_fns:
    print(f"  {fn}")

for fit_num, residual_fn in enumerate(residual_fns):

    if residual_fn.endswith('.nii'):
        fit_dir      = "{}/{}/170k/pmf".format(pp_dir, subject)
        jobs_dir     = "{}/{}/170k/pmf/jobs".format(pp_dir, subject)
        logs_dir     = "{}/{}/170k/pmf/log_outputs".format(pp_dir, subject)

    elif residual_fn.endswith('.gii'):
        fit_dir      = "{}/{}/fsnative/pmf".format(pp_dir, subject)
        jobs_dir     = "{}/{}/fsnative/pmf/jobs".format(pp_dir, subject)
        logs_dir     = "{}/{}/fsnative/pmf/log_outputs".format(pp_dir, subject)

    for d in [fit_dir, jobs_dir, logs_dir]:
        os.makedirs(d, exist_ok=True)

    slurm_cmd = """\
#!/bin/bash
#SBATCH -p {cluster_name}
#SBATCH -A {server_project}
#SBATCH --nodes=1
#SBATCH --mem={memory_val}gb
#SBATCH --cpus-per-task={nb_procs}
#SBATCH --time={hour_proc}:00:00
#SBATCH -e {logs_dir}/{subject}_residuals_gauss_%N_%j_%a.err
#SBATCH -o {logs_dir}/{subject}_residuals_gauss_%N_%j_%a.out
#SBATCH -J {subject}_res_gaussfit
""".format(cluster_name=cluster_name, server_project=server_project,
           nb_procs=nb_procs, hour_proc=hour_proc,
           subject=subject, memory_val=memory_val, logs_dir=logs_dir)

    fit_cmd = "python residuals_gaussfit.py {} {} {} {} {}".format(
        main_dir, project_dir, subject, residual_fn, nb_procs)

    sh_fn = "{}/jobs/{}_residuals_gaussfit-{}.sh".format(fit_dir, subject, fit_num)
    with open(sh_fn, 'w') as f:
        f.write("{} \n{} \n{} \n{}".format(slurm_cmd, fit_cmd, chmod_cmd, chgrp_cmd))

    print("Submitting {} to queue".format(sh_fn))
    os.system("sbatch {}".format(sh_fn))