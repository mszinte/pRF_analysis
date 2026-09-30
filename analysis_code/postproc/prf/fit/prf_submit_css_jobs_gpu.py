"""
-----------------------------------------------------------------------------------------
prf_submit_css_jobs_gpu.py
-----------------------------------------------------------------------------------------
Goal of the script:
Create and submit jobscript to make a css fit on GPU for pRF analysis (prfmodel)
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject name (e.g. sub-01)
sys.argv[4]: analysis name (e.g. prf)
sys.argv[5]: group (e.g. 327)
sys.argv[6]: server project (e.g. b327)
-----------------------------------------------------------------------------------------
Output(s):
.sh file to execute in server
-----------------------------------------------------------------------------------------
To run:
0. activate the prfmodel environment (it is passed on to the jobs)
>> conda activate prfmodel
1. cd to function
>> cd ~/projects/pRF_analysis/analysis_code/postproc/prf/fit
2. run python command
python prf_submit_css_jobs_gpu.py [main directory] [project name] [subject]
                                  [analysis name] [group] [server project]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/projects/pRF_analysis/analysis_code/postproc/prf/fit
python prf_submit_css_jobs_gpu.py /scratch/mszinte/data amsterdam24 sub-03 prf 327 b327
-----------------------------------------------------------------------------------------
Written by Sina Kling (sina.kling@outlook.de)
adapted from Uriel Lascombes (uriel.lascombes@laposte.net)
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
from settings_utils import load_settings
from pycortex_utils import set_pycortex_config_file
 
# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
analysis_name = sys.argv[4]
group = sys.argv[5]
server_project = sys.argv[6]
 
# GPU job settings
cluster_name = 'volta'      # GPU 
nb_gpus = 1                 # prfmodel uses a single GPU
nb_procs = 4                # CPUs for data loading / saving 
memory_val = 40             # RAM (GB): data, grid predictions and outputs on the host
hour_proc = 2               
fit_script = 'prf_cssfit_gpu.py'
 
# Load settings
base_dir = os.path.abspath(os.path.join(script_dir, "../../../../"))
general_settings_path = os.path.join(base_dir, project_dir, "settings.yml")
analysis_settings_path = os.path.join(base_dir, project_dir, f"{analysis_name}-analysis.yml")
settings = load_settings([general_settings_path, analysis_settings_path])
analysis_info = settings[0]
 
formats = analysis_info['formats']
extensions = analysis_info['extensions']
task_names = analysis_info['analysis_task_names']
preproc_prep = analysis_info['preproc_prep']
filtering = analysis_info['filtering']
normalization = analysis_info['normalization']
avg_methods = analysis_info['avg_methods'] # selects input type
output_folder = analysis_info["output_folder"]
dm_name = analysis_info["dm_name"]
 
# Set pycortex db and colormaps
cortex_dir = "{}/{}/derivatives/pp_data/cortex".format(main_dir, project_dir)
set_pycortex_config_file(cortex_dir)
 
# Define directories
pp_dir = "{}/{}/derivatives/pp_data".format(main_dir, project_dir)
 
# define permission cmd
chmod_cmd = "chmod -Rf 771 {}/{}".format(main_dir, project_dir)
chgrp_cmd = "chgrp -Rf {} {}/{}".format(group, main_dir, project_dir)

# Conda base of the submitting shell (e.g. ~/softwares/miniconda3), from $CONDA_EXE = <base>/bin/conda
conda_exe = os.environ.get('CONDA_EXE')
if conda_exe is None:
    raise EnvironmentError("CONDA_EXE not set: run this script from a shell where conda is initialised")
conda_dir = os.path.dirname(os.path.dirname(conda_exe))  
conda_env = 'prfmodel'
 
# Define environment cmd: activate the prfmodel env inside the job only
env_cmd = """\
source {conda_dir}/etc/profile.d/conda.sh
conda activate {conda_env}
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH
export KERAS_BACKEND=tensorflow
export TF_FORCE_GPU_ALLOW_GROWTH=true
cd {script_dir}
echo "python: $(which python)"
nvidia-smi""".format(conda_dir=conda_dir, conda_env=conda_env, script_dir=script_dir)
 
# Define fns (filenames)
pp_fns = []
for avg_method in avg_methods:
    for task_name in task_names:
        print(f"=========================")
        print(f"Running CSS GPU Fit with: data: {avg_method}, task: {task_name}, dm: {dm_name}, analysis: {analysis_name}")
        print(f"=========================\n")
        dct_avg_gii_fns = "{}/{}/fsnative/func/{}_{}_{}_{}/*_task-{}_*{}*.func.gii".format(
            pp_dir, subject, preproc_prep, filtering, normalization, avg_method, task_name, avg_method)
        dct_avg_nii_fns = "{}/{}/170k/func/{}_{}_{}_{}/*_task-{}_*{}*.dtseries.nii".format(
            pp_dir, subject, preproc_prep, filtering, normalization, avg_method, task_name, avg_method)
 
        # Accumulate the results
        pp_fns.extend(glob.glob(dct_avg_gii_fns))
        pp_fns.extend(glob.glob(dct_avg_nii_fns))
 
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
#SBATCH -e {log_dir}/{subject}_{avg_method}-{analysis_name}-css{dm_name}_gpu_fit_%N_%j_%a.err
#SBATCH -o {log_dir}/{subject}_{avg_method}-{analysis_name}-css{dm_name}_gpu_fit_%N_%j_%a.out
#SBATCH -J {subject}_{avg_method}-{analysis_name}-css{dm_name}_gpu_fit
""".format(server_project=server_project,
           cluster_name=cluster_name,
           nb_gpus=nb_gpus,
           nb_procs=nb_procs,
           hour_proc=hour_proc,
           subject=subject,
           memory_val=memory_val,
           log_dir=prf_logs_dir,
           avg_method=avg_method,
           analysis_name=analysis_name,
           dm_name=dm_name)
 
    # Define fit cmd (last argument = number of jobs, unused on GPU but kept for the fit script's inputs)
    fit_cmd = "python {} {} {} {} {} {} {}".format(
        fit_script, main_dir, project_dir, subject, pp_fn, analysis_name, nb_procs)
 
    # Create shs
    sh_fn = "{}/jobs/{}_{}-{}-css{}_gpu_fit-{}.sh".format(prf_dir, subject, avg_method, analysis_name, dm_name, fit_num)
 
    of = open(sh_fn, 'w')
    of.write("{} \n{} \n{} \n{} \n{}".format(slurm_cmd, env_cmd, fit_cmd,
                                             chmod_cmd, chgrp_cmd))
    of.close()
 
    # Submit jobs
    print("Submitting {} to queue".format(sh_fn))
    os.system("sbatch {}".format(sh_fn))
 