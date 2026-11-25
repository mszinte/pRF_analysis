"""
-----------------------------------------------------------------------------------------
submit_ncsf_jobs_gpu.py
-----------------------------------------------------------------------------------------
Goal of the script:
Create and submit jobscripts to make the nCSF fit on GPU (prfmodel, ncsf_fit_gpu.py)
Same structure as analysis_code/postproc/prf/fit/prf_submit_css_jobs_gpu_v2.py
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject name (e.g. sub-01)
sys.argv[4]: analysis name (e.g. nCSF, loads nCSF-analysis.yml)
sys.argv[5]: group (e.g. 327)
sys.argv[6]: server project (e.g. b327)
-----------------------------------------------------------------------------------------
Output(s):
.sh file to execute in server
-----------------------------------------------------------------------------------------
To run:
0. use a python environment where prfmodel is installed (https://github.com/popylar-org/prfmodel);
   the jobs run with the same python interpreter as this script
1. cd to function
>> cd ~/projects/pRF_analysis/nCSF/postproc/nCSF/fit
2. run python command
python submit_ncsf_jobs_gpu.py [main directory] [project name] [subject]
                               [analysis name] [group] [server project]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/projects/pRF_analysis/nCSF/postproc/nCSF/fit
python submit_ncsf_jobs_gpu.py /scratch/mszinte/data nCSF sub-01 nCSF 327 b327
-----------------------------------------------------------------------------------------
Written by Martin Szinte (martin.szinte@gmail.com)
and Uriel Lascombes (uriel.lascombes@laposte.net)
-----------------------------------------------------------------------------------------
"""
# Stop warnings
import warnings
warnings.filterwarnings("ignore")

# General imports
import os
import sys
import glob
from importlib.util import find_spec

# stop here if prfmodel or its dependencies are missing (checked without importing them)
PRFMODEL_URL = "https://github.com/popylar-org/prfmodel"
missing_packages = [package for package in ['prfmodel', 'keras', 'torch'] if find_spec(package) is None]
if missing_packages:
    raise ImportError(
        "Missing python package(s) for the GPU nCSF fit: {}\n"
        "Please install prfmodel and its dependencies in the python environment used to run this script,\n"
        "see {}\n"
        "(current python: {})".format(', '.join(missing_packages), PRFMODEL_URL, sys.executable))

# Personal imports
script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.abspath(os.path.join(script_dir, "../../../../analysis_code/utils")))
from settings_utils import load_settings

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
analysis_name = sys.argv[4]
group = sys.argv[5]
server_project = sys.argv[6]

# GPU job settings
cluster_name = 'hopper'     # GPU
nb_gpus = 1                 # prfmodel uses a single GPU
nb_procs = 4                # CPUs for data loading / saving
memory_val = 40             # RAM (GB)
hour_proc = 2
fit_script = 'ncsf_fit_gpu.py'

# Load settings (same order as ncsf_fit_gpu.py)
base_dir = os.path.abspath(os.path.join(script_dir, "../../../../"))
settings_fns = [os.path.join(base_dir, project_dir, "settings.yml"),
                os.path.join(base_dir, project_dir, "prf-analysis.yml"),
                os.path.join(base_dir, project_dir, f"{analysis_name}-analysis.yml")]
analysis_info = load_settings([fn for fn in settings_fns if os.path.isfile(fn)])[0]

task_names = analysis_info['analysis_task_names']
preproc_prep = analysis_info['preproc_prep']
filtering = analysis_info['filtering']
normalizations = analysis_info['normalization']
avg_methods = analysis_info['avg_methods']
output_folder = analysis_info['output_folder']
formats_to_fit = analysis_info.get('ncsf_formats_to_fit', ['fsnative', '170k'])

# Define directories
pp_dir = "{}/{}/derivatives/pp_data".format(main_dir, project_dir)

# define permission cmd
chmod_cmd = "chmod -Rf 771 {}/{}".format(main_dir, project_dir)
chgrp_cmd = "chgrp -Rf {} {}/{}".format(group, main_dir, project_dir)

env_cmd = """\
cd {script_dir}
echo "python: {python}"
nvidia-smi""".format(script_dir=script_dir, python=sys.executable)

# Define fns (filenames)
pp_fns = []
for normalization in normalizations:
    for avg_method in avg_methods:
        for task_name in task_names:
            print("=========================")
            print(f"nCSF GPU fit: normalization: {normalization}, data: {avg_method}, task: {task_name}, analysis: {analysis_name}")
            print("=========================\n")
            if 'fsnative' in formats_to_fit:
                pp_fns.extend(glob.glob("{}/{}/fsnative/func/{}_{}_{}_{}/*_task-{}_*{}*.func.gii".format(
                    pp_dir, subject, preproc_prep, filtering, normalization, avg_method, task_name, avg_method)))
            if '170k' in formats_to_fit:
                pp_fns.extend(glob.glob("{}/{}/170k/func/{}_{}_{}_{}/*_task-{}_*{}*.dtseries.nii".format(
                    pp_dir, subject, preproc_prep, filtering, normalization, avg_method, task_name, avg_method)))

for fit_num, pp_fn in enumerate(pp_fns):
    fmt = '170k' if pp_fn.endswith('.nii') else 'fsnative'
    ncsf_dir = "{}/{}/{}/{}".format(pp_dir, subject, fmt, output_folder)
    ncsf_logs_dir = "{}/log_outputs".format(ncsf_dir)
    os.makedirs("{}/jobs".format(ncsf_dir), exist_ok=True)
    os.makedirs(ncsf_logs_dir, exist_ok=True)

    # averaging method of this file, from its folder name (e.g. fmriprep_dct_z-score_loo-avg -> loo-avg)
    avg_method = os.path.basename(os.path.dirname(pp_fn)).split('_')[-1]

    slurm_cmd = """\
#!/bin/bash
#SBATCH -p {cluster_name}
#SBATCH -A {server_project}
#SBATCH --nodes=1
#SBATCH --gres=gpu:{nb_gpus}
#SBATCH --mem={memory_val}gb
#SBATCH --cpus-per-task={nb_procs}
#SBATCH --time={hour_proc}:00:00
#SBATCH -e {log_dir}/{subject}_{avg_method}-{analysis_name}-ncsf_gpu_fit_%N_%j.err
#SBATCH -o {log_dir}/{subject}_{avg_method}-{analysis_name}-ncsf_gpu_fit_%N_%j.out
#SBATCH -J {subject}_{avg_method}-{analysis_name}-ncsf_gpu_fit
""".format(server_project=server_project, cluster_name=cluster_name, nb_gpus=nb_gpus, nb_procs=nb_procs,
           hour_proc=hour_proc, subject=subject, memory_val=memory_val, log_dir=ncsf_logs_dir,
           avg_method=avg_method, analysis_name=analysis_name)

    # same python interpreter as the submitting shell (where prfmodel was found)
    fit_cmd = "{} {} {} {} {} {} {} {}".format(
        sys.executable, fit_script, main_dir, project_dir, subject, pp_fn, analysis_name, nb_procs)

    sh_fn = "{}/jobs/{}_{}-{}-ncsf_gpu_fit-{}.sh".format(ncsf_dir, subject, avg_method, analysis_name, fit_num)
    with open(sh_fn, 'w') as of:
        of.write("{} \n{} \n{} \n{} \n{}".format(slurm_cmd, env_cmd, fit_cmd, chmod_cmd, chgrp_cmd))

    # Submit jobs
    print("Submitting {} to queue".format(sh_fn))
    os.system("sbatch {}".format(sh_fn))
