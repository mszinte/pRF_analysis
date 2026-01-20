"""
-----------------------------------------------------------------------------------------
fmriprep_sbatch.py
-----------------------------------------------------------------------------------------
Goal of the script:
Run fMRIprep on mesocentre using job mode
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject (e.g. sub-001)
sys.argv[4]: data group (e.g. 327)
sys.argv[5]: server project (e.g. b327)
sys.argv[6]: anat only (1: yes; 0: no)
-----------------------------------------------------------------------------------------
Output(s):
preprocessed files
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/pRF_analysis/analysis_code/preproc/functional
2. run python command
python fmriprep_sbatch.py [main] [project] [subject] [group] [server] [anat_only] 
-----------------------------------------------------------------------------------------
Exemple:
cd ~/projects/pRF_analysis/analysis_code/preproc/functional
python fmriprep_sbatch.py /scratch/mszinte/data RetinoMaps sub-01 327 b327 1
-----------------------------------------------------------------------------------------
Written by Martin Szinte (martin.szinte@gmail.com)
Edited by Uriel Lascombes (uriel.lascombes@laposte.net)
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
sub_id = subject[4:]
group = sys.argv[4]
server_project = sys.argv[5]
anat_only = sys.argv[6]

# File and folder
log_dir = f"{main_dir}/{project_dir}/derivatives/fmriprep/log_outputs"
singularity_dir = "{}/{}/code/singularity".format(main_dir, project_dir)

# Define cluster/server specific parameters
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../"))
settings_path = os.path.join(base_dir, project_dir, "settings.yml")
fmriprep_settings_path = os.path.join(base_dir, project_dir, "fmriprep-settings.yml")
settings = load_settings([settings_path, fmriprep_settings_path])
analysis_info = settings[0]

cluster_name  = analysis_info['cluster_name']
fmriprep_simg = analysis_info['fmriprep_simg']
procs_hours = analysis_info['procs_hours']
procs_hours_anat = analysis_info['procs_hours_anat']
email_sbatch = analysis_info['email_sbatch']
memory_sbatch = analysis_info['memory_sbatch']
bold2anat_dof = analysis_info['bold2anat_dof']
output_spaces = analysis_info['output-spaces']
cifti_output = analysis_info['cifti-output']
if cifti_output=='': cifti_output_in = ''
else: cifti_output_in = f"cifti-output {cifti_output}"
mem_limit = analysis_info['mem_limit']
fs_no_resume = analysis_info['fs_no_resume']
if fs_no_resume: fs_no_resume_in = '--fs-no-resume'
else: fs_no_resume_in = ''
skip_bids_valid = analysis_info['skip_bids_valid']
if skip_bids_valid: skip_bids_valid_in = "--skip_bids_validation"
else: skip_bids_valid_in = ''

if anat_only:
    nb_procs = 1
    anat_only_in = '--anat-only'
    anat_only_end = '_anat_only'
    procs_hours_in = procs_hours
else: 
    nb_procs = 32
    anat_only_in = ''
    anat_only_end = ''
    procs_hours_in = procs_hours_anat

temp_dir = "{}/temp/{}/sub-{}".format(main_dir, project_dir, sub_num)
os.makedirs(temp_dir, exist_ok=True)

# define SLURM cmd
slurm_cmd = f"\
#!/bin/bash\n\
#SBATCH --mail-type=ALL\n\
#SBATCH -p {cluster_name}\n\
#SBATCH --mail-user={email_sbatch}\n\
#SBATCH -A {server_project}\n\
#SBATCH --nodes=1\n\
#SBATCH --mem={memory_sbatch}gb\n\
#SBATCH --cpus-per-task={nb_procs}\n\
#SBATCH --time={procs_hours_in}\n\
#SBATCH -e {log_dir}/{subject}_fmriprep{anat_only_end}_%N_%j_%a.err\n\
#SBATCH -o {log_dir}/{subject}_fmriprep{anat_only_end}_%N_%j_%a.out\n\
#SBATCH -J {subject}_fmriprep{anat_only_end}\n\
#SBATCH --mail-type=BEGIN,END\n\n"

#define singularity cmd
singularity_cmd = "singularity run --cleanenv {tf_bind} -B {main_dir}:/work_dir {simg} --fs-license-file /work_dir/{project_dir}/code/freesurfer/license.txt --fs-subjects-dir /work_dir/{project_dir}/derivatives/fmriprep/freesurfer/ /work_dir/{project_dir}/ /work_dir/{project_dir}/derivatives/fmriprep/fmriprep{aroma_end}/ participant --participant-label {sub_num} -w /work_dir/temp/{project_dir}/sub-{sub_num} --bold2anat-dof {dof} --output-spaces T1w fsnative {fsaverage} {hcp_cifti} --mem-mb {memory_val}000 --nthreads {nb_procs:.0f} {anat_only}{use_aroma}{use_fmapfree}{use_fs_no_resume}{use_skip_bids_val}".format(
        tf_bind=tf_bind, main_dir=main_dir, project_dir=project_dir,
        simg=singularity_dir, sub_num=sub_num, nb_procs=nb_procs, use_fs_no_resume=use_fs_no_resume,
        anat_only=anat_only, use_aroma=use_aroma, use_fmapfree=use_fmapfree,
        use_skip_bids_val=use_skip_bids_val, fsaverage = fsaverage,hcp_cifti=hcp_cifti, memory_val=memory_val,
        temp_dir=temp_dir,
        dof=dof, aroma_end=aroma_end)

# define permission cmd
chmod_cmd = f"chmod -Rf 771 {main_dir}/{project_dir} \n"
chgrp_cmd = f"chgrp -Rf {group} {main_dir}/{project_dir} \n"
chmod_cmd_temp = f"chmod -Rf 771 {main_dir}/temp \n"
chgrp_cmd_temp = f"chgrp -Rf {group} {main_dir}/temp \n"

# create sh folder and file
sh_fn = f"{main_dir}/{project_dir}/derivatives/fmriprep/jobs/sub-{sub_id}_fmriprep{anat_only_end}.sh"
os.makedirs(f"{main_dir}/{project_dir}/derivatives/fmriprep/jobs", exist_ok=True)
os.makedirs(f"{main_dir}/{project_dir}/derivatives/fmriprep/log_outputs", exist_ok=True)

of = open(sh_fn, 'w')
of.write(f"{slurm_cmd} {singularity_cmd} {chmod_cmd} {chgrp_cmd} {chmod_cmd_temp} {chgrp_cmd_temp}")
of.close()

# Submit jobs
print(f"Submitting {sh_fn} to queue")
os.chdir(log_dir)
os.system(f"sbatch {sh_fn}")