"""
-----------------------------------------------------------------------------------------
freesurfer_siam_pial.py
-----------------------------------------------------------------------------------------
Goal of the script:
FALLBACK, for subjects already processed by fMRIPrep without the SIAM brain mask
(the default is to give the SIAM mask to fMRIPrep, see fmriprep-settings.yml):
replace the freesurfer brainmask by the SIAM brain mask (instead of manual edition with 
pial_edits.sh) and run freesurfer pial steps with this new brainmask
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: freesurfer subject (e.g. sub-01 or sub-01_ses-01)
sys.argv[4]: group of shared data (e.g. 327)
sys.argv[5]: server project (e.g. b327)
-----------------------------------------------------------------------------------------
Output(s):
derivatives/siam/{sub}[/{ses}]/anat/{subject}_desc-brain_mask.nii.gz (SIAM brain mask, T1w space)
mri/brainmask_orig.mgz (copy of fmriprep brainmask), mri/brainmask_siam.mgz (SIAM mask in 
freesurfer space), mri/brainmask.mgz (new brainmask) and new freesurfer pial surfaces
-----------------------------------------------------------------------------------------
To run:
0. run after siam_sbatch.py and fmriprep anat-only, in a python environment with nibabel
   and scipy
1. cd to function
>> cd ~/projects/pRF_analysis/analysis_code/preproc/anatomical/
2. run python command
python freesurfer_siam_pial.py [main directory] [project name] [subject] [group] [server]
-----------------------------------------------------------------------------------------
Exemple:
python freesurfer_siam_pial.py /scratch/ulascombes/data CRASH sub-2100563718_ses-preseason 327 b327
-----------------------------------------------------------------------------------------
Written by Uriel Lascombes (uriel.lascombes@laposte.net)
-----------------------------------------------------------------------------------------
"""
# General imports
import os
import sys
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

# Define cluster/server specific parameters
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../"))
settings_path = os.path.join(base_dir, project_dir, "settings.yml")
siam_settings_path = os.path.join(base_dir, project_dir, "siam-settings.yml")
settings = load_settings([settings_path, siam_settings_path])
analysis_info = settings[0]
cluster_name = analysis_info['cluster_name']
email_sbatch = analysis_info['email_sbatch']
hour_proc = analysis_info['fs_pial_procs_hours']
nb_procs = analysis_info['fs_pial_nb_procs']
memory_val = analysis_info['fs_pial_memory_sbatch']

# define directory and freesurfer licence
script_dir = os.getcwd()
log_dir = f"{main_dir}/{project_dir}/derivatives/freesurfer_pial/log_outputs"
job_dir = f"{main_dir}/{project_dir}/derivatives/freesurfer_pial/jobs"
fs_dir = f"{main_dir}/{project_dir}/derivatives/fmriprep/freesurfer"
fs_mri_dir = f"{fs_dir}/{subject}/mri"
fs_licence = f"{main_dir}/{project_dir}/code/freesurfer/license.txt"
# freesurfer subject sub-XX or sub-XX_ses-YY (sessionwise) and corresponding SIAM brain mask
sub_name = subject.split('_')[0]
if '_ses-' in subject: 
    ses_name = subject.split('_')[1]
    siam_mask_fn = f"{main_dir}/{project_dir}/derivatives/siam/{sub_name}/{ses_name}/anat/{subject}_desc-brain_mask.nii.gz"
else: 
    siam_mask_fn = f"{main_dir}/{project_dir}/derivatives/siam/{sub_name}/anat/{subject}_desc-brain_mask.nii.gz"
os.makedirs(log_dir, exist_ok=True)
os.makedirs(job_dir, exist_ok=True)

# define SLURM cmd
slurm_cmd = f"""\
#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH -p {cluster_name}
#SBATCH --mail-user={email_sbatch}
#SBATCH -A {server_project}
#SBATCH --nodes=1
#SBATCH --mem={memory_val}gb
#SBATCH --cpus-per-task={nb_procs}
#SBATCH --time={hour_proc}
#SBATCH -e {log_dir}/{subject}_freesurfer-siam-pial_%N_%j_%a.err
#SBATCH -o {log_dir}/{subject}_freesurfer-siam-pial_%N_%j_%a.out
#SBATCH -J {subject}_freesurfer-siam-pial\n"""

# define permission cmd
chmod_cmd = f"chmod -Rf 771 {main_dir}/{project_dir}\n"
chgrp_cmd = f"chgrp -Rf {group} {main_dir}/{project_dir}\n"

# define SIAM brain mask cmd (T1w native space)
siam_mask_cmd = f"""\
cd {script_dir}
{sys.executable} siam_brainmask.py {main_dir} {project_dir} {sub_name}\n"""

# define freesurfer command
freesurfer_cmd = f"""\
export FREESURFER_HOME={main_dir}/{project_dir}/code/freesurfer
export SUBJECTS_DIR={fs_dir}
export FS_LICENSE={fs_licence}
source $FREESURFER_HOME/SetUpFreeSurfer.sh
if [ ! -f {fs_mri_dir}/brainmask_orig.mgz ]; then cp {fs_mri_dir}/brainmask.mgz {fs_mri_dir}/brainmask_orig.mgz; fi
mri_vol2vol --mov {siam_mask_fn} --targ {fs_mri_dir}/T1.mgz --regheader --nearest --o {fs_mri_dir}/brainmask_siam.mgz
mri_mask {fs_mri_dir}/T1.mgz {fs_mri_dir}/brainmask_siam.mgz {fs_mri_dir}/brainmask.mgz
recon-all -autorecon-pial -subjid {subject} -no-isrunning -openmp {nb_procs}\n"""

# create sh folder and file
sh_fn = f"{job_dir}/{subject}_freesurfer-siam-pial.sh"
of = open(sh_fn, 'w')
of.write(f"{slurm_cmd}\n{siam_mask_cmd}{freesurfer_cmd}{chmod_cmd}{chgrp_cmd}")
of.close()

# submit jobs
print(f"Submitting {sh_fn} to queue")
os.chdir(log_dir)
os.system(f"sbatch {sh_fn}")
