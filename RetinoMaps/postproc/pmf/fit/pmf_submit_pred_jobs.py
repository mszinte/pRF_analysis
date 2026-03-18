"""
-----------------------------------------------------------------------------------------
pmf_submit_pred_jobs.py
-----------------------------------------------------------------------------------------
Goal of the script:
Create and submit jobscript to make a prf gaussian fit
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
>> cd ~/projects/pRF_analysis/analysis_code/postproc/prf/fit
2. run python command
python pmf_submit_pred_jobs.py [main directory] [project name] [subject]
                                [group] [server project]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/projects/pRF_analysis/RetinoMaps/postproc/pmf/fit
python pmf_submit_pred_jobs.py /scratch/mszinte/data RetinoMaps sub-01 327 b327
-----------------------------------------------------------------------------------------
Written by Martin Szinte (martin.szinte@gmail.com)
and Uriel Lascombes (uriel.lascombes@laposte.net)
adapted by Sina Kling (sina.kling@outlook.de)
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
from pycortex_utils import set_pycortex_config_file
from settings_utils import load_settings

# inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
group = sys.argv[4]
server_project = sys.argv[5]

memory_val = 30
hour_proc = 6
nb_procs = 8

# Load settings
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../../"))
settings_path = os.path.join(base_dir, project_dir, "pmf-settings.yml")
prf_settings_path = os.path.join(base_dir, project_dir, "prf-analysis.yml")
settings = load_settings([settings_path, prf_settings_path])
analysis_info = settings[0]

cluster_name  = analysis_info['cluster_name']
prf_task_names = analysis_info['prf_task_names']
pmf_task_names = analysis_info['pmf_task_names']   # 'SacLoc'
preproc_prep = analysis_info['preproc_prep']
filtering = analysis_info['filtering']
normalization = analysis_info['normalization']
avg_methods = analysis_info['avg_methods']

# Set pycortex db and colormaps
cortex_dir = "{}/{}/derivatives/pp_data/cortex".format(main_dir, project_dir)
set_pycortex_config_file(cortex_dir)

# Define directories
pp_dir = "{}/{}/derivatives/pp_data".format(main_dir, project_dir)

# define permission cmd
chmod_cmd = "chmod -Rf 771 {}/{}".format(main_dir, project_dir)
chgrp_cmd = "chgrp -Rf {} {}/{}".format(group, main_dir, project_dir)

# Define helper to find the matching pRF avg gauss fit file
def get_prf_avg_gauss_fit_fn(pp_dir, subject, pp_fn, preproc_prep, filtering, normalization):
    """
    Find the pRF avg gauss fit file matching the format and hemisphere of pp_fn.
    - For 170k (.nii): matches a single cifti file
    - For fsnative (.gii): extracts hemisphere (hemi-L or hemi-R) from pp_fn and matches accordingly
    """
    if pp_fn.endswith('.nii'):
        pattern = "{}/{}/170k/prf/fit/{}_task-{}_{}_{}_{}_avg_prf-gauss_fit*.nii".format(pp_dir, subject, subject, prf_task_names[0], preproc_prep, filtering, normalization)
        matches = glob.glob(pattern)

    elif pp_fn.endswith('.gii'):
        # Extract hemisphere token from the input filename
        fn_basename = os.path.basename(pp_fn)
        if 'hemi-L' in fn_basename:
            hemi = 'hemi-L'
        elif 'hemi-R' in fn_basename:
            hemi = 'hemi-R'
        else:
            raise ValueError("Could not determine hemisphere from filename: {}".format(fn_basename))

        pattern = "{}/{}/fsnative/prf/fit/{}_task-{}_{}_{}_{}_{}_avg_prf-gauss_fit*.gii".format(pp_dir, subject, subject, prf_task_names[0], hemi, preproc_prep, filtering, normalization)
        matches = glob.glob(pattern)

    if not matches:
        raise FileNotFoundError("No pRF avg gauss fit found for pattern: {}".format(pattern))

    return matches[0]  # expecting a single match per format/hemisphere


# Define fns (filenames): SacLoc task only, concat averaging method only
pp_fns = []
for avg_method in avg_methods:
    if avg_method != 'concat':
        continue  # only concat

    dct_avg_gii_fns = "{}/{}/fsnative/func/{}_{}_{}_{}/*_task-{}*{}*.func.gii".format(pp_dir, subject, preproc_prep, filtering, normalization, avg_method,pmf_task_names[0], avg_method)
    dct_avg_nii_fns = "{}/{}/170k/func/{}_{}_{}_{}/*_task-{}*{}*.dtseries.nii".format(pp_dir, subject, preproc_prep, filtering, normalization, avg_method,pmf_task_names[0], avg_method)

    pp_fns.extend(glob.glob(dct_avg_gii_fns))
    pp_fns.extend(glob.glob(dct_avg_nii_fns))

                
for fit_num, pp_fn in enumerate(pp_fns):

    if pp_fn.endswith('.nii'):
        pmf_dir = "{}/{}/170k/pmf".format(pp_dir, subject)
        os.makedirs(pmf_dir, exist_ok=True)
        pmf_jobs_dir = "{}/{}/170k/pmf/jobs".format(pp_dir, subject)
        os.makedirs(pmf_jobs_dir, exist_ok=True)
        prf_logs_dir = "{}/{}/170k/pmf/log_outputs".format(pp_dir, subject)
        os.makedirs(pmf_logs_dir, exist_ok=True)

    elif pp_fn.endswith('.gii'):
        pmf_dir = "{}/{}/fsnative/pmf".format(pp_dir, subject)
        os.makedirs(pmf_dir, exist_ok=True)
        pmf_jobs_dir = "{}/{}/fsnative/pmf/jobs".format(pp_dir, subject)
        os.makedirs(pmf_jobs_dir, exist_ok=True)
        pmf_logs_dir = "{}/{}/fsnative/pmf/log_outputs".format(pp_dir, subject)
        os.makedirs(pmf_logs_dir, exist_ok=True)

    # Find matching pRF avg gauss fit
    prf_gauss_fit_fn = get_prf_avg_gauss_fit_fn(pp_dir, subject, pp_fn, preproc_prep, filtering, normalization)

    slurm_cmd = """\
#!/bin/bash
#SBATCH -p {cluster_name}
#SBATCH -A {server_project}
#SBATCH --nodes=1
#SBATCH --mem={memory_val}gb
#SBATCH --cpus-per-task={nb_procs}
#SBATCH --time={hour_proc}:00:00
#SBATCH -e {log_dir}/{subject}_pmf_pred_%N_%j_%a.err
#SBATCH -o {log_dir}/{subject}_pmf_pred_%N_%j_%a.out
#SBATCH -J {subject}_pmfpred
""".format(server_project=server_project, cluster_name=cluster_name,
           nb_procs=nb_procs, hour_proc=hour_proc, 
           subject=subject, memory_val=memory_val, log_dir=pmf_logs_dir)

    # define fit cmd  
    fit_cmd = "python pmf_pred_from_prf.py {} {} {} {} {} {}".format(main_dir, project_dir, subject, pp_fn, prf_gauss_fit_fn, nb_procs)
    
    # create sh
    sh_fn = "{}/jobs/{}_pmf_pred-{}.sh".format(pmf_dir, subject, fit_num)
    of = open(sh_fn, 'w')
    of.write("{} \n{} \n{} \n{}".format(slurm_cmd, fit_cmd, chmod_cmd, chgrp_cmd))
    of.close()

    # submit jobs
    print("Submitting {} to queue".format(sh_fn))
    os.system("sbatch {}".format(sh_fn))