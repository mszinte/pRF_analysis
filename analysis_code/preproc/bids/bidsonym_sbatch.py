"""
-----------------------------------------------------------------------------------------
bidsonym_sbatch.py
-----------------------------------------------------------------------------------------
Goal of the script:
Deface anatomical images
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject name (e.g. sub-01)
sys.argv[4]: group (e.g. 327)
sys.argv[5]: server project (e.g. b327)
sys.argv[6]: Bidsoym singularity version (e.g. v0.0.4 or latest)
sys.argv[7]: defacing algorithm (e.g. pydeface, mri_deface, quickshear, mridefacer)
-----------------------------------------------------------------------------------------
Output(s):
sh file for running batch command
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/pRF_analysis/analysis_code/preproc/bids/
2. run python command
>> python bidsonym_sbatch.py [main directory] [project name] [subject num] [group] 
                             [server project] [bidsonym_version] [defacing_algorithm]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/projects/pRF_analysis/analysis_code/preproc/bids/
python bidsonym_sbatch.py /scratch/mszinte/data nCSF sub-01 327 b327 v0.0.4 pydeface
python bidsonym_sbatch.py /scratch/mszinte/data amblyo7T_prf sub-16 327 b327 v0.0.4 pydeface
python bidsonym_sbatch.py /scratch/ulascombes/data CRASH sub-01 327 b327 latest pydeface
-----------------------------------------------------------------------------------------
Note: 'nobrainer' brain extraction (used for QC only) was removed after v0.0.4,
'bet' is used for other versions (e.g. latest from ghcr.io/peerherholz/bidsonym)
-----------------------------------------------------------------------------------------
Written by Uriel Lascombes (uriel.lascombes@laposte.net)
edited by Martin Szinte (martin.szinte@gmail.com)
-----------------------------------------------------------------------------------------
"""

# Stop warnings
import warnings
warnings.filterwarnings("ignore")

# debug 
import ipdb 
deb = ipdb.set_trace

# general imports
import os
import sys

# Personal imports
sys.path.append("{}/../../../analysis_code/utils".format(os.getcwd()))
from settings_utils import load_settings

# inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
group = sys.argv[4]
server_project = sys.argv[5]
bydsonym_version = sys.argv[6]
defacing_algorithm = sys.argv[7]

subject_num = subject.split('-')[1]

# Define analysis parameters
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../"))
settings_path = os.path.join(base_dir, project_dir, "settings.yml")
settings = load_settings([settings_path])
analysis_info = settings[0]
cluster_name  = analysis_info['cluster_name']
nb_procs = 1
memory_val = 48
hour_proc = 1

# set folders
log_dir = "{}/{}/derivatives/log_outputs".format(main_dir, project_dir)
job_dir = "{}/{}/derivatives/jobs".format(main_dir, project_dir)
os.makedirs(log_dir, exist_ok=True)
os.makedirs(job_dir, exist_ok=True)
bids_dir = "{}/{}".format(main_dir, project_dir)
bidsonym_img = "{}/code/singularity/bidsonym-{}.simg".format(bids_dir, bydsonym_version)

slurm_cmd = """\
#!/bin/bash
#SBATCH -p {cluster_name}
#SBATCH -A {server_project}
#SBATCH --nodes=1
#SBATCH --mem={memory_val}gb
#SBATCH --cpus-per-task={nb_procs}
#SBATCH --time={hour_proc}:00:00
#SBATCH -e {log_dir}/{subject}_bidsonym_%N_%j_%a.err
#SBATCH -o {log_dir}/{subject}_bidsonym_%N_%j_%a.out
#SBATCH -J {subject}_bidsonym
""".format(server_project=server_project, cluster_name=cluster_name, subject=subject,
nb_procs=nb_procs, hour_proc=hour_proc, memory_val=memory_val, log_dir=log_dir)

# bidsonym_cmd = "singularity run --cleanenv -B {}:/bids_dataset {} /bids_dataset participant --participant-label {} --deid {} --brainextraction --bet_frac 0.5".format(bids_dir, bidsonym_img, subject_num, defacing_algorithm)
# brain extraction is only used for QC: nobrainer exists only up to v0.0.4, bet otherwise
if bydsonym_version == 'v0.0.4':
    brainextraction_opt = "--brainextraction nobrainer"
else:
    brainextraction_opt = "--brainextraction bet --bet_frac 0.5"

bidsonym_cmd = "singularity run --cleanenv -B {}:/bids_dataset {} /bids_dataset participant --participant_label {} --deid {} {} --skip_bids_validation".format(bids_dir, bidsonym_img, subject_num, defacing_algorithm, brainextraction_opt)

# create sh fn
sh_fn = "{}/{}_bidsonym.sh".format(job_dir, subject)
of = open(sh_fn, 'w')

of.write("{} \n{}".format(slurm_cmd, bidsonym_cmd))
of.close()

# Submit jobs
print("Submitting {} to queue".format(sh_fn))
os.system("sbatch {}".format(sh_fn))
