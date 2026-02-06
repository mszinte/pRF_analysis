"""
-----------------------------------------------------------------------------------------
flatten_sbatch.py
-----------------------------------------------------------------------------------------
Goal of the script:
Run mris_flatten on mesocentre using job mode
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: freesurfer subject name (e.g. sub-01 or sub-01_ses-01)
sys.argv[4]: group (e.g. 327)
sys.argv[5]: server project (e.g. b327)
-----------------------------------------------------------------------------------------
Output(s):
preprocessed files
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/pRF_analysis/analysis_code/preproc/anatomical/
2. run python command
python flatten_sbatch.py [main directory] [project name] [subject] [group] [server]
-----------------------------------------------------------------------------------------
Example:
python flatten_sbatch.py /scratch/mszinte/data RetinoMaps sub-01 327 b327
-----------------------------------------------------------------------------------------
Written by Martin Szinte (martin.szinte@gmail.com)
Edited by Uriel Lascombes (uriel.lascombes@laposte.net)
-----------------------------------------------------------------------------------------
"""

# imports modules
import sys
import os
import ipdb
deb = ipdb.set_trace

# Personal imports
sys.path.append("{}/../../utils".format(os.getcwd()))
from settings_utils import load_settings

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
sub_num = subject[-2:]
group = sys.argv[4]
server_project = sys.argv[5]
cluster_name  = 'skylake'
nb_procs = 8
memory_val = 48
hour_proc = 20
hemis = ['lh', 'rh']

# Load input
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../"))
settings_path = os.path.join(base_dir, project_dir, "settings.yml")
settings = load_settings([settings_path])
analysis_info = settings[0]
flattening_method = analysis_info['flattening_method']

# Define cluster/server specific parameters
log_dir = "{}/{}/derivatives/flatten/log_outputs".format(main_dir,project_dir)
fs_dir = "{}/{}/derivatives/fmriprep/freesurfer".format(main_dir, project_dir)
job_dir = "{}/{}/derivatives/flatten/jobs".format(main_dir,project_dir)
fs_licence = '{}/{}/code/freesurfer/license.txt'.format(main_dir, project_dir)
os.makedirs(log_dir, exist_ok=True)
os.makedirs(job_dir, exist_ok=True)

# Define slurm cmd
for hemi in hemis:
    slurm_cmd = f"\
#!/bin/bash\n\
#SBATCH -p skylake\n\
#SBATCH -A {server_project}\n\
#SBATCH --nodes=1\n\
#SBATCH --mem={memory_val}gb\n\
#SBATCH --cpus-per-task={nb_procs}\n\
#SBATCH --time={hour_proc}:00:00\n\
#SBATCH -e {log_dir}/{subject}_{hemi}_{flattening_method}_%N_%j_%a.err\n\
#SBATCH -o {log_dir}/{subject}_{hemi}_{flattening_method}_%N_%j_%a.out\n\
#SBATCH -J {subject}_{hemi}_{flattening_method}\n"

    freesurfer_cmd = f"\
export FREESURFER_HOME={main_dir}/{project_dir}/code/freesurfer\n\
export SUBJECTS_DIR={fs_dir}\n\
export FS_LICENSE={fs_licence}\n\
source $FREESURFER_HOME/SetUpFreeSurfer.sh\n\
cd {fs_dir}/{subject}/surf\n"
    
    chmod_cmd = f"chmod -Rf 771 {main_dir}/{project_dir}"
    chgrp_cmd = f"chgrp -Rf {group} {main_dir}/{project_dir}"

    # Define autoflatten cmd
    autoflatten_cmd = f"autoflatten {fs_dir}/{subject} --parallel --hemispheres {hemi} --overwrite"
    
    # Define mris_flatten cmd
    mris_flatten_cmd = f"mris_flatten {hemi}.full.patch.3d {hemi}.full.flat.patch.3d"

    # Create sh fn
    sh_fn = f"{job_dir}/{subject}_{hemi}_{flattening_method}.sh"
    
    # Write commands
    of = open(sh_fn, 'w')
    if flattening_method == 'autoflatten':
        of.write(f"{slurm_cmd}\n")
        of.write(f"{freesurfer_cmd}\n")
        of.write(f"{autoflatten_cmd}\n")
        of.write(f"{chmod_cmd}\n")
        of.write(f"{chgrp_cmd}")
        
    elif flattening_method == 'mrisflatten':
        of.write(f"{slurm_cmd}\n")
        of.write(f"{freesurfer_cmd}\n")
        of.write(f"{mris_flatten_cmd}\n")
        of.write(f"{chmod_cmd}\n")
        of.write(f"{chgrp_cmd}")
        
    of.close()

    # Submit jobs
    print("Submitting {} to queue".format(sh_fn))
    os.chdir(log_dir)
    os.system("sbatch {}".format(sh_fn))