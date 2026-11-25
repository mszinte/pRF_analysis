"""
-----------------------------------------------------------------------------------------
siam_sbatch.py
-----------------------------------------------------------------------------------------
Goal of the script:
Run SIAM (Segment It All Model) head tissue segmentation on mesocentre using job mode
https://github.com/romainVala/SIAM
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject (e.g. sub-001)
sys.argv[4]: data group (e.g. 327)
sys.argv[5]: server project (e.g. b327)
-----------------------------------------------------------------------------------------
Output(s):
SIAM segmentations in derivatives/siam/{subject}
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/pRF_analysis/analysis_code/preproc/anatomical
2. run python command
python siam_sbatch.py [main] [project] [subject] [group] [server]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/projects/pRF_analysis/analysis_code/preproc/anatomical
python siam_sbatch.py /scratch/ulascombes/data CRASH sub-2100563718 327 b327
-----------------------------------------------------------------------------------------
Written by Uriel Lascombes (uriel.lascombes@laposte.net)
-----------------------------------------------------------------------------------------
"""
# General imports
import os
import sys
import glob
import ipdb
deb = ipdb.set_trace

# Personal imports
sys.path.append("{}/../../utils".format(os.getcwd()))
from settings_utils import load_settings

# inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
group = sys.argv[4]
server_project = sys.argv[5]

# File and folder
bids_dir = f"{main_dir}/{project_dir}"
siam_dir = f"{main_dir}/{project_dir}/derivatives/siam"
subject_dir = f"{siam_dir}/{subject}"
log_dir = f"{siam_dir}/log_outputs"
job_dir = f"{siam_dir}/jobs"
singularity_dir = f"{main_dir}/{project_dir}/code/singularity"
os.makedirs(subject_dir, exist_ok=True)
os.makedirs(log_dir, exist_ok=True)
os.makedirs(job_dir, exist_ok=True)

# Define cluster/server specific parameters
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../"))
settings_path = os.path.join(base_dir, project_dir, "settings.yml")
siam_settings_path = os.path.join(base_dir, project_dir, "siam-settings.yml")
settings = load_settings([settings_path, siam_settings_path])
analysis_info = settings[0]

cluster_name = analysis_info['cluster_name']
email_sbatch = analysis_info['email_sbatch']
memory_sbatch = analysis_info['siam_memory_sbatch']
procs_hours = analysis_info['siam_procs_hours']
nb_procs = analysis_info['siam_nb_procs']
siam_simg = analysis_info['siam_simg']
siam_model = analysis_info['siam_model']
siam_device = analysis_info['siam_device']
siam_gpu_cluster_name = analysis_info['siam_gpu_cluster_name']
siam_voxelsize = analysis_info['siam_voxelsize']

if siam_device == 'cuda':
    cluster_name_in = siam_gpu_cluster_name
    gpu_sbatch_in = "#SBATCH --gres=gpu:1\n"
    nv_in = '--nv'
    gpu_check_cmd = "nvidia-smi \n"
else:
    cluster_name_in = cluster_name
    gpu_sbatch_in = ''
    nv_in = ''
    gpu_check_cmd = ''

# Find T1w images of the subject (defaced by bidsonym)
t1w_fns = sorted(glob.glob(f"{bids_dir}/{subject}/ses-*/anat/{subject}_*_T1w.nii.gz") + \
                 glob.glob(f"{bids_dir}/{subject}/anat/{subject}_*T1w.nii.gz"))
if len(t1w_fns) == 0:
    sys.exit(f"No T1w image found for {subject} in {bids_dir}")

# define SLURM cmd
slurm_cmd = f"\
#!/bin/bash\n\
#SBATCH --mail-type=ALL\n\
#SBATCH -p {cluster_name_in}\n\
#SBATCH --mail-user={email_sbatch}\n\
#SBATCH -A {server_project}\n\
#SBATCH --nodes=1\n\
#SBATCH --mem={memory_sbatch}gb\n\
#SBATCH --cpus-per-task={nb_procs}\n\
{gpu_sbatch_in}\
#SBATCH --time={procs_hours}\n\
#SBATCH -e {log_dir}/{subject}_siam_%N_%j_%a.err\n\
#SBATCH -o {log_dir}/{subject}_siam_%N_%j_%a.out\n\
#SBATCH -J {subject}_siam\n\
#SBATCH --mail-type=BEGIN,END\n\n"

# copy T1w in derivatives (SIAM writes its outputs next to its input)
copy_cmd = ''
for t1w_fn in t1w_fns:
    copy_cmd += f"cp {t1w_fn} {subject_dir}/ \n"

# define singularity cmd (input folder: all T1w copied in derivatives/siam/{subject})
singularity_cmd = f"singularity exec --cleanenv {nv_in} -B {main_dir}:/work_dir {singularity_dir}/{siam_simg} siam-pred -i /work_dir/{project_dir}/derivatives/siam/{subject} -m {siam_model} -device {siam_device} -nbthread {nb_procs} -voxelsize {siam_voxelsize} \n"

# define permission cmd
chmod_cmd = f"chmod -Rf 771 {main_dir}/{project_dir} \n"
chgrp_cmd = f"chgrp -Rf {group} {main_dir}/{project_dir} \n"

# create sh file
sh_fn = f"{job_dir}/{subject}_siam.sh"
of = open(sh_fn, 'w')
of.write(f"{slurm_cmd}{gpu_check_cmd}{copy_cmd}{singularity_cmd}{chmod_cmd}{chgrp_cmd}")
of.close()

# Submit jobs
print(f"Submitting {sh_fn} to queue")
os.chdir(log_dir)
os.system(f"sbatch {sh_fn}")
