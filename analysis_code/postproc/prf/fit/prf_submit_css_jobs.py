"""
-----------------------------------------------------------------------------------------
prf_submit_css_jobs.py
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
>> cd ~/projects/pRF_analysis/analysis_code/postproc/prf/fit
2. run python command
python prf_submit_css_jobs.py [main directory] [project name] [subject] 
                              [group] [server project]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/projects/pRF_analysis/analysis_code/postproc/prf/fit
python prf_submit_css_jobs.py /scratch/mszinte/data RetinoMaps sub-01 327 b327
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
import json

# Personal imports
sys.path.append("{}/../../../utils".format(os.getcwd()))
from pycortex_utils import set_pycortex_config_file

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
group = sys.argv[4]
server_project = sys.argv[5]
memory_val = 30
nb_procs = 32

# Cluster settings
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../../"))
settings_path = os.path.join(base_dir, project_dir, "settings.json")

with open(settings_path) as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)
cluster_name  = analysis_info['cluster_name']
formats = analysis_info['formats']
extensions = analysis_info['extensions']
filter_rois = analysis_info['filter_rois']
if filter_rois: hour_proc = 8
else: hour_proc = 16
prf_task_names = analysis_info['prf_task_names']
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

# Define fns (filenames)
pp_fns = []
for avg_method in avg_methods:
    for prf_task_name in prf_task_names:
        dct_avg_gii_fns = "{}/{}/fsnative/func/{}_{}_{}_{}/*_task-{}*{}*.func.gii".format(
            pp_dir, subject, preproc_prep, filtering, normalization, avg_method, prf_task_name, avg_method)
        dct_avg_nii_fns = "{}/{}/170k/func/{}_{}_{}_{}/*_task-{}*{}*.dtseries.nii".format(
            pp_dir, subject, preproc_prep, filtering, normalization, avg_method, prf_task_name, avg_method)

        # Accumulate the results
        pp_fns.extend(glob.glob(dct_avg_gii_fns))
        pp_fns.extend(glob.glob(dct_avg_nii_fns))

for fit_num, pp_fn in enumerate(pp_fns):
    if pp_fn.endswith('.nii'):
        prf_dir = "{}/{}/170k/prf".format(pp_dir, subject)
        os.makedirs(prf_dir, exist_ok=True)
        prf_jobs_dir = "{}/{}/170k/prf/jobs".format(pp_dir, subject)
        os.makedirs(prf_jobs_dir, exist_ok=True)
        prf_logs_dir = "{}/{}/170k/prf/log_outputs".format(pp_dir, subject)
        os.makedirs(prf_logs_dir, exist_ok=True)

    elif pp_fn.endswith('.gii'):
        prf_dir = "{}/{}/fsnative/prf".format(pp_dir, subject)
        os.makedirs(prf_dir, exist_ok=True)
        prf_jobs_dir = "{}/{}/fsnative/prf/jobs".format(pp_dir, subject)
        os.makedirs(prf_jobs_dir, exist_ok=True)
        prf_logs_dir = "{}/{}/fsnative/prf/log_outputs".format(pp_dir, subject)
        os.makedirs(prf_logs_dir, exist_ok=True)
    
    slurm_cmd = """\
#!/bin/bash
#SBATCH -p {cluster_name}
#SBATCH -A {server_project}
#SBATCH --nodes=1
#SBATCH --mem={memory_val}gb
#SBATCH --cpus-per-task={nb_procs}
#SBATCH --time={hour_proc}:00:00
#SBATCH -e {log_dir}/{subject}_css_fit_%N_%j_%a.err
#SBATCH -o {log_dir}/{subject}_css_fit_%N_%j_%a.out
#SBATCH -J {subject}_css_fit
""".format(server_project=server_project, 
           cluster_name=cluster_name,
           nb_procs=nb_procs, 
           hour_proc=hour_proc, 
           subject=subject, 
           memory_val=memory_val, 
           log_dir=prf_logs_dir)

    # Define fit cmd
    fit_cmd = "python prf_cssfit.py {} {} {} {} {}".format(
        main_dir, project_dir, subject, pp_fn, nb_procs)
    
    # Create sh
    sh_fn = "{}/jobs/{}_prf_css_fit-{}.sh".format(prf_dir, subject, fit_num)

    of = open(sh_fn, 'w')
    of.write("{} \n{} \n{} \n{}".format(slurm_cmd, fit_cmd, 
                                        chmod_cmd, chgrp_cmd))
    of.close()

    # Submit jobs
    print("Submitting {} to queue".format(sh_fn))
    os.system("sbatch {}".format(sh_fn))