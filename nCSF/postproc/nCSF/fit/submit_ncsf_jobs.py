"""
-----------------------------------------------------------------------------------------
submit_ncsf_jobs.py
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
>> cd ~/projects/pRF_analysis/nCSF/postproc/nCSF/fit
2. run python command
python submit_ncsf_jobs.py [main directory] [project name] [subject] 
                              [group] [server project]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/projects/pRF_analysis/nCSF/postproc/nCSF/fit
python submit_ncsf_jobs.py /scratch/mszinte/data nCSF sub-01 327 b327
-----------------------------------------------------------------------------------------
Written by Martin Szinte (martin.szinte@gmail.com)
and Uriel Lascombes (uriel.lascombes@laposte.net)
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
sys.path.append("{}/../../../../analysis_code/utils".format(os.getcwd()))
from settings_utils import load_settings
# from pycortex_utils import set_pycortex_config_file

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
group = sys.argv[4]
server_project = sys.argv[5]
memory_val = 30
nb_procs = 8
hour_proc = 10

# Load settings
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../../"))
settings_path = os.path.join(base_dir, project_dir, "settings.yml")
prf_settings_path = os.path.join(base_dir, project_dir, "prf-analysis.yml")
settings = load_settings([settings_path, prf_settings_path])
analysis_info = settings[0]

cluster_name  = analysis_info['cluster_name']
formats = analysis_info['formats']
extensions = analysis_info['extensions']
nCSF_task_names = analysis_info['nCSF_task_names'][0]
preproc_prep = analysis_info['preproc_prep']
filtering = analysis_info['filtering']
normalization = analysis_info['normalization']
avg_methods = analysis_info['avg_methods']

# # Set pycortex db and colormaps
# cortex_dir = "{}/{}/derivatives/pp_data/cortex".format(main_dir, project_dir)
# set_pycortex_config_file(cortex_dir)

# Define directories
pp_dir = "{}/{}/derivatives/pp_data".format(main_dir, project_dir)

# Define fns (filenames)
pp_fns = []
for avg_method in avg_methods:
    dct_avg_gii_fns = "{}/{}/fsnative/func/{}_{}_{}_{}/*_task-{}_*{}*.func.gii".format(
        pp_dir, subject, preproc_prep, filtering, normalization, avg_method, nCSF_task_names, avg_method)
    dct_avg_nii_fns = "{}/{}/170k/func/{}_{}_{}_{}/*_task-{}_*{}*.dtseries.nii".format(
        pp_dir, subject, preproc_prep, filtering, normalization, avg_method, nCSF_task_names, avg_method)

    # Accumulate the results
    pp_fns.extend(glob.glob(dct_avg_gii_fns))
    pp_fns.extend(glob.glob(dct_avg_nii_fns))

for fit_num, pp_fn in enumerate(pp_fns):
    if pp_fn.endswith('.nii'):
        ncsf_dir = "{}/{}/170k/ncsf".format(pp_dir, subject)
        os.makedirs(ncsf_dir, exist_ok=True)
        ncsf_jobs_dir = "{}/{}/170k/ncsf/jobs".format(pp_dir, subject)
        os.makedirs(ncsf_jobs_dir, exist_ok=True)
        ncsf_logs_dir = "{}/{}/170k/ncsf/log_outputs".format(pp_dir, subject)
        os.makedirs(ncsf_logs_dir, exist_ok=True)

    elif pp_fn.endswith('.gii'):
        ncsf_dir = "{}/{}/fsnative/ncsf".format(pp_dir, subject)
        os.makedirs(ncsf_dir, exist_ok=True)
        ncsf_jobs_dir = "{}/{}/fsnative/ncsf/jobs".format(pp_dir, subject)
        os.makedirs(ncsf_jobs_dir, exist_ok=True)
        ncsf_logs_dir = "{}/{}/fsnative/ncsf/log_outputs".format(pp_dir, subject)
        os.makedirs(ncsf_logs_dir, exist_ok=True)
    
    slurm_cmd = """\
#!/bin/bash
#SBATCH -p {cluster_name}
#SBATCH -A {server_project}
#SBATCH --nodes=1
#SBATCH --mem={memory_val}gb
#SBATCH --cpus-per-task={nb_procs}
#SBATCH --time={hour_proc}:00:00
#SBATCH -e {log_dir}/{subject}_ncsf_fit_%N_%j_%a.err
#SBATCH -o {log_dir}/{subject}_ncsf_fit_%N_%j_%a.out
#SBATCH -J {subject}_ncsf_fit
""".format(server_project=server_project, 
           cluster_name=cluster_name,
           nb_procs=nb_procs, 
           hour_proc=hour_proc, 
           subject=subject, 
           memory_val=memory_val, 
           log_dir=ncsf_logs_dir)

    # Define fit cmd
    fit_cmd = "python ncsf_fit.py {} {} {} {} {}".format(
        main_dir, project_dir, subject, pp_fn, nb_procs)
    
    # Create shs
    sh_fn = "{}/jobs/{}_ncsf_fit-{}.sh".format(ncsf_dir, subject, fit_num)

    of = open(sh_fn, 'w')
    of.write("{} \n{}".format(slurm_cmd, fit_cmd))
    of.close()

    # # Submit jobs
    print("Submitting {} to queue".format(sh_fn))
    os.system("sbatch {}".format(sh_fn))
    stop