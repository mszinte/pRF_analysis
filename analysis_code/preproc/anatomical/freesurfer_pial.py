"""
-----------------------------------------------------------------------------------------
freesurfer_pial.py
-----------------------------------------------------------------------------------------
Goal of the script:
Run freesurfer with new brainmask manually edited
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject (e.g. sub-01)
sys.argv[4]: group of shared data (e.g. 327)
-----------------------------------------------------------------------------------------
Output(s):
new freesurfer segmentation files
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/RetinoMaps/analysis_code/preproc/anatomical/
2. run python command
python freesurfer_pial.py [main directory] [project name] [subject] [group] [proj_name]
-----------------------------------------------------------------------------------------
Exemple:
python freesurfer_pial.py /scratch/mszinte/data RetinoMaps sub-01 327 b327
-----------------------------------------------------------------------------------------
Written by Martin Szinte (mail@martinszinte.net)
-----------------------------------------------------------------------------------------
"""

# imports modules
import os
import sys
import json

# inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
group = sys.argv[4]
proj_name = sys.argv[5]

# Define cluster/server specific parameters
cluster_name = 'skylake'
hour_proc = 20
nb_procs = 8
memory_val = 48


# define directory and freesurfer licence
log_dir = "{}/{}/derivatives/freesurfer_pial/log_outputs".format(main_dir,project_dir)
job_dir = "{}/{}/derivatives/freesurfer_pial/jobs".format(main_dir,project_dir)
fs_dir = "{}/{}/derivatives/fmriprep/freesurfer/".format(main_dir, project_dir)
fs_licence = '{}/{}/code/freesurfer/license.txt'.format(main_dir, project_dir)
os.makedirs(log_dir, exist_ok=True)
os.makedirs(job_dir, exist_ok=True)

# define SLURM cmd
slurm_cmd = """\
#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH -p skylake
#SBATCH -A {proj_name}
#SBATCH --nodes=1
#SBATCH --mem={memory_val}gb
#SBATCH --cpus-per-task={nb_procs}
#SBATCH --time={hour_proc}:00:00
#SBATCH -e {log_dir}/{subject}_freesurfer-pial_%N_%j_%a.err
#SBATCH -o {log_dir}/{subject}_freesurfer-pial_%N_%j_%a.out
#SBATCH -J {subject}_freesurfer-pial\n\n""".format( nb_procs=nb_procs, subject=subject, memory_val = memory_val, 
                                                    hour_proc=hour_proc, log_dir=log_dir, proj_name=proj_name)

# define permission cmd
chmod_cmd = "chmod -Rf 771 {main_dir}/{project_dir}\n".format(main_dir=main_dir, project_dir=project_dir)
chgrp_cmd = "chgrp -Rf {group} {main_dir}/{project_dir}\n".format(main_dir=main_dir, project_dir=project_dir, group=group)

# define freesurfer command
freesurfer_cmd = """\
export FREESURFER_HOME={}/{}/code/freesurfer
export SUBJECTS_DIR={}\n\
export FS_LICENSE={}\n\
source $FREESURFER_HOME/SetUpFreeSurfer.sh
recon-all -autorecon-pial -subjid {} -no-isrunning\n""".format(main_dir, project_dir, fs_dir, fs_licence, subject)

# create sh folder and file
sh_dir = "{}/{}_freesurfer-pial.sh".format(job_dir, subject)

of = open(sh_dir, 'w')
of.write("{}{}{}{}".format(slurm_cmd, chmod_cmd, chgrp_cmd, freesurfer_cmd, chmod_cmd, chgrp_cmd))
of.close()

# submit jobs
print("Submitting {} to queue".format(sh_dir))
os.chdir(log_dir)
os.system("sbatch {}".format(sh_dir))