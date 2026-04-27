"""
-----------------------------------------------------------------------------------------
residuals_submit_gauss_jobs.py
-----------------------------------------------------------------------------------------
Goal of the script:
Create and submit jobscript to make a prf gaussian fit on PMF residuals or original BOLD
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject name(s) (e.g. sub-01 or sub-01,sub-02,sub-03)
sys.argv[4]: group (e.g. 327)
sys.argv[5]: server project (e.g. b327)
sys.argv[6]: input type - 'residuals' or 'bold' (default: residuals)
-----------------------------------------------------------------------------------------
Output(s):
.sh file to execute in server
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/pRF_analysis/RetinoMaps/postproc/pmf/postfit
2. run python command
python residuals_submit_gauss_jobs.py [main directory] [project name] [subject(s)]
                                [group] [server project] [input type]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/projects/pRF_analysis/RetinoMaps/postproc/pmf/postfit

# single subject, residuals (default)
python residuals_submit_gauss_jobs.py /scratch/mszinte/data RetinoMaps sub-01 327 b327

# multiple subjects, residuals
python residuals_submit_gauss_jobs.py /scratch/mszinte/data RetinoMaps "sub-01,sub-02,sub-03" 327 b327

# multiple subjects, original bold
python residuals_submit_gauss_jobs.py /scratch/mszinte/data RetinoMaps "sub-01,sub-02" 327 b327 bold
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
subjects = [s.strip() for s in sys.argv[3].split(',')] 
group          = sys.argv[4]
server_project = sys.argv[5]
input_type     = sys.argv[6] if len(sys.argv) > 6 else 'residuals'  # 'residuals' or 'bold'

if input_type not in ('residuals', 'bold'):
    raise ValueError(f"input_type must be 'residuals' or 'bold', got '{input_type}'")

memory_val = 30
hour_proc  = 6
nb_procs   = 8

# Load settings
base_dir          = os.path.abspath(os.path.join(script_dir, "../../../../"))
settings_path     = os.path.join(base_dir, project_dir, "pmf-settings.yml")
prf_settings_path = os.path.join(base_dir, project_dir, "prf-analysis.yml")
settings          = load_settings([settings_path, prf_settings_path])
analysis_info     = settings[0]
task              = 'SacLoc'
cluster_name      = analysis_info['cluster_name']

# Set pycortex db and colormaps
cortex_dir = "{}/{}/derivatives/pp_data/cortex".format(main_dir, project_dir)
set_pycortex_config_file(cortex_dir)

pp_dir = "{}/{}/derivatives/pp_data".format(main_dir, project_dir)

# Permission commands
chmod_cmd = "chmod -Rf 771 {}/{}".format(main_dir, project_dir)
chgrp_cmd = "chgrp -Rf {} {}/{}".format(group, main_dir, project_dir)

# ── main loop over subjects ───────────────────────────────────────────────────
for subject in subjects:

    print(f"\n{'='*60}")
    print(f"Subject: {subject}  |  input_type: {input_type}")
    print(f"{'='*60}")

    # find input files depending on input_type
    input_fns = []

    if input_type == 'residuals':
        input_fns.extend(glob.glob(
            "{}/{}/fsnative/pmf/pmf_residuals/*_task-{}*.func.gii".format(
                pp_dir, subject, task)))
        input_fns.extend(glob.glob(
            "{}/{}/170k/pmf/pmf_residuals/*_task-{}*.dtseries.nii".format(
                pp_dir, subject, task)))

    elif input_type == 'bold':
        input_fns.extend(glob.glob(
            "{}/{}/fsnative/func/fmriprep_dct_z-score_concat/*_task-{}*.func.gii".format(
                pp_dir, subject, task)))
        input_fns.extend(glob.glob(
            "{}/{}/170k/func/fmriprep_dct_z-score_concat/*_task-{}*.dtseries.nii".format(
                pp_dir, subject, task)))

    if not input_fns:
        print(f"WARNING: No {input_type} files found for {subject}, skipping.")
        print(f"  Checked fsnative: {pp_dir}/{subject}/fsnative/pmf/pmf_residuals/")
        print(f"  Checked 170k:     {pp_dir}/{subject}/170k/pmf/pmf_residuals/")
        continue

    print(f"Found {len(input_fns)} file(s):")
    for fn in input_fns:
        print(f"  {fn}")

    # ── loop over files ───────────────────────────────────────────────────────
    for fit_num, input_fn in enumerate(input_fns):

        if input_fn.endswith('.nii'):
            fit_dir  = "{}/{}/170k/pmf".format(pp_dir, subject)
            jobs_dir = "{}/{}/170k/pmf/jobs".format(pp_dir, subject)
            logs_dir = "{}/{}/170k/pmf/log_outputs".format(pp_dir, subject)

        elif input_fn.endswith('.gii'):
            fit_dir  = "{}/{}/fsnative/pmf".format(pp_dir, subject)
            jobs_dir = "{}/{}/fsnative/pmf/jobs".format(pp_dir, subject)
            logs_dir = "{}/{}/fsnative/pmf/log_outputs".format(pp_dir, subject)

        for d in [fit_dir, jobs_dir, logs_dir]:
            os.makedirs(d, exist_ok=True)

        # job name and script name reflect input_type
        job_label = 'res' if input_type == 'residuals' else 'bold'

        slurm_cmd = """\
#!/bin/bash
#SBATCH -p {cluster_name}
#SBATCH -A {server_project}
#SBATCH --nodes=1
#SBATCH --mem={memory_val}gb
#SBATCH --cpus-per-task={nb_procs}
#SBATCH --time={hour_proc}:00:00
#SBATCH -e {logs_dir}/{subject}_{job_label}_gauss_%N_%j_%a.err
#SBATCH -o {logs_dir}/{subject}_{job_label}_gauss_%N_%j_%a.out
#SBATCH -J {subject}_{job_label}_gaussfit
""".format(cluster_name=cluster_name, server_project=server_project,
           nb_procs=nb_procs, hour_proc=hour_proc,
           subject=subject, memory_val=memory_val,
           logs_dir=logs_dir, job_label=job_label)

        fit_cmd = "python residuals_gaussfit.py {} {} {} {} {}".format(
            main_dir, project_dir, subject, input_fn, nb_procs)

        sh_fn = "{}/jobs/{}_{}_gaussfit-{}.sh".format(
            fit_dir, subject, job_label, fit_num)

        with open(sh_fn, 'w') as f:
            f.write("{} \n{} \n{} \n{}".format(
                slurm_cmd, fit_cmd, chmod_cmd, chgrp_cmd))

        print(f"Submitting {sh_fn} to queue")
        os.system("sbatch {}".format(sh_fn))