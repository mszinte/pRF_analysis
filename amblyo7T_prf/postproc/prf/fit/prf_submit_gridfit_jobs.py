"""
-----------------------------------------------------------------------------------------
prf_submit_gridfit_jobs.py
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
sys.argv[6]: OPTIONAL main analysis folder (e.g. prf_em_ctrl)
-----------------------------------------------------------------------------------------
Output(s):
.sh file to execute in server
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/pRF_analysis/amblyo7T_prf/postproc/prf/fit/
2. run python command
python prf_submit_gridfit_jobs.py [main directory] [project name] [subject] 
                                  [group] [server project] [analysis folder - optional]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/projects/pRF_analysis/amblyo7T_prf/postproc/prf/fit/
python prf_submit_gridfit_jobs.py /scratch/mszinte/data amblyo7T_prf sub-01 327 b327 pRFRightEye
python prf_submit_gridfit_jobs.py /scratch/mszinte/data amblyo7T_prf sub-01 327 b327 pRFLeftEye
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
import ipdb
deb = ipdb.set_trace 
# inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
group = sys.argv[4]
server_project = sys.argv[5]
if len(sys.argv) > 6: output_folder = sys.argv[6]
else: output_folder = "prf"

memory_val = 30
hour_proc = 6
nb_procs = 8

# cluster settings
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../../"))
settings_path = os.path.join(base_dir, project_dir, "settings.json")

with open(settings_path) as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)
cluster_name  = analysis_info['cluster_name']
prf_task_name = output_folder

# Define directories
pp_dir = "{}/{}/derivatives/pp_data".format(main_dir, project_dir)

# define permission cmd
chmod_cmd = "chmod -Rf 771 {}/{}".format(main_dir, project_dir)
chgrp_cmd = "chgrp -Rf {} {}/{}".format(group, main_dir, project_dir)

# Define fns (filenames)
dct_concat_gii_fns = "{}/{}/fsnative/func/fmriprep_dct_concat/*_task-{}_*concat*.func.gii".format(pp_dir, subject, prf_task_name)
dct_concat_nii_fns = "{}/{}/170k/func/fmriprep_dct_concat/*_task-{}_*concat*.dtseries.nii".format(pp_dir, subject, prf_task_name)
pp_fns = glob.glob(dct_concat_gii_fns) + glob.glob(dct_concat_nii_fns) 

for fit_num, pp_fn in enumerate(pp_fns):

    if pp_fn.endswith('.nii'):
        prf_dir = "{}/{}/170k/{}".format(pp_dir, subject, output_folder)
        os.makedirs(prf_dir, exist_ok=True)
        prf_jobs_dir = "{}/{}/170k/{}/jobs".format(pp_dir, subject, output_folder)
        os.makedirs(prf_jobs_dir, exist_ok=True)
        prf_logs_dir = "{}/{}/170k/{}/log_outputs".format(pp_dir, subject, output_folder)
        os.makedirs(prf_logs_dir, exist_ok=True)

    elif pp_fn.endswith('.gii'):
        prf_dir = "{}/{}/fsnative/{}".format(pp_dir, subject, output_folder)
        os.makedirs(prf_dir, exist_ok=True)
        prf_jobs_dir = "{}/{}/fsnative/{}/jobs".format(pp_dir, subject, output_folder)
        os.makedirs(prf_jobs_dir, exist_ok=True)
        prf_logs_dir = "{}/{}/fsnative/{}/log_outputs".format(pp_dir, subject, output_folder)
        os.makedirs(prf_logs_dir, exist_ok=True)
    
    slurm_cmd = """\
#!/bin/bash
#SBATCH -p {cluster_name}
#SBATCH -A {server_project}
#SBATCH --nodes=1
#SBATCH --mem={memory_val}gb
#SBATCH --cpus-per-task={nb_procs}
#SBATCH --time={hour_proc}:00:00
#SBATCH -e {log_dir}/{subject}_grid_gaussfit_%N_%j_%a.err
#SBATCH -o {log_dir}/{subject}_grid_gaussfit_%N_%j_%a.out
#SBATCH -J {subject}_grid_gaussfit
""".format(server_project=server_project, cluster_name=cluster_name,
           nb_procs=nb_procs, hour_proc=hour_proc, 
           subject=subject, memory_val=memory_val, log_dir=prf_logs_dir)

    # define fit cmd
    fit_cmd = "python prf_gridfit.py {} {} {} {} {} {}".format(
        main_dir, project_dir, subject, pp_fn, nb_procs, output_folder)
    
    # create sh
    sh_fn = "{}/jobs/{}_prf_gridfit-{}.sh".format(prf_dir, subject, fit_num)
    of = open(sh_fn, 'w')
    of.write("{} \n{} \n{} \n{}".format(slurm_cmd, fit_cmd, 
                                        chmod_cmd, chgrp_cmd))
    of.close()

    # submit jobs
    print("Submitting {} to queue".format(sh_fn))
    os.system("sbatch {}".format(sh_fn))