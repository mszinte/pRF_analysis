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
-----------------------------------------------------------------------------------------
Output(s):
.sh file to execute in server
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/RetinoMaps/analysis_code/postproc/prf/fit
2. run python command
python prf_submit_gridfit_jobs.py [main directory] [project name] [subject] 
                                  [group] [server project]
-----------------------------------------------------------------------------------------
Exemple:
python prf_submit_gridfit_jobs.py /scratch/mszinte/data RetinoMaps sub-21 327 b327
-----------------------------------------------------------------------------------------
Written by Martin Szinte (mail@martinszinte.net)
Edited by Uriel Lascombes (uriel.lascombes@laposte.net)
-----------------------------------------------------------------------------------------
"""

# stop warnings
import warnings
warnings.filterwarnings("ignore")

# general imports
import os
import json
import sys
import glob
import ipdb
deb = ipdb.set_trace

# inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
group = sys.argv[4]
server_project = sys.argv[5]
memory_val = 30
hour_proc = 2
nb_procs = 8

# cluster settings
with open('../../../settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)
cluster_name  = analysis_info['cluster_name']

# Define directories
pp_dir = "{}/{}/derivatives/pp_data".format(main_dir, project_dir)

# define permission cmd
chmod_cmd = "chmod -Rf 771 {main_dir}/{project_dir}".format(main_dir=main_dir, project_dir=project_dir)
chgrp_cmd = "chgrp -Rf {group} {main_dir}/{project_dir}".format(main_dir=main_dir, project_dir=project_dir, group=group)

# Define fns (filenames)
dct_avg_nii_fns = "{}/{}/170k/func/fmriprep_dct_avg/*_task-pRF_*avg*.dtseries.nii".format(pp_dir,subject)
dct_avg_gii_fns = "{}/{}/fsnative/func/fmriprep_dct_avg/*_task-pRF_*avg*.func.gii".format(pp_dir,subject)


pp_fns=  glob.glob(dct_avg_gii_fns) + glob.glob(dct_avg_nii_fns) 
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
#SBATCH -e {log_dir}/{subject}_grid_gaussfit_%N_%j_%a.err
#SBATCH -o {log_dir}/{subject}_grid_gaussfit_%N_%j_%a.out
#SBATCH -J {subject}_grid_gaussfit
""".format(server_project=server_project, cluster_name=cluster_name,
           nb_procs=nb_procs, hour_proc=hour_proc, 
           subject=subject, memory_val=memory_val, log_dir=prf_logs_dir)

    # define fit cmd
    fit_cmd = "python prf_gridfit.py {} {} {} {} {} ".format(
        main_dir, project_dir, subject, pp_fn, nb_procs )
    
    # create sh
    sh_fn = "{}/jobs/{}_prf_gridfit-{}.sh".format(prf_dir, subject, fit_num)
    of = open(sh_fn, 'w')
    of.write("{} \n{} \n{} \n{}".format(slurm_cmd, fit_cmd, 
                                        chmod_cmd, chgrp_cmd))
    of.close()

    # submit jobs
    print("Submitting {} to queue".format(sh_fn))
    os.system("sbatch {}".format(sh_fn))