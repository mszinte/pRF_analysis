"""
-----------------------------------------------------------------------------------------
freesurfer_import_pycortex.py
-----------------------------------------------------------------------------------------
Goal of the script:
Load freesurfer and import subject in pycortex database
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject (e.g. sub-01)
sys.argv[4]: group (e.g. 327)
-----------------------------------------------------------------------------------------
Output(s):
None
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/[PROJECT]/analysis_code/preproc/functional/
2. run python command
python freesurfer_import_pycortex.py [main directory] [project name] [subject] [group]
-----------------------------------------------------------------------------------------
Executions:
cd ~/projects/pRF_analysis/analysis_code/preproc/functional/
python freesurfer_import_pycortex.py /scratch/mszinte/data MotConf sub-01 327
-----------------------------------------------------------------------------------------
Written by Martin Szinte (mail@martinszinte.net)
Edited by Uriel Lascombes (uriel.lascombes@laposte.net)
-----------------------------------------------------------------------------------------
"""
# stop warnings
import warnings
warnings.filterwarnings("ignore")

# Debug 
import ipdb
deb = ipdb.set_trace

# imports
import os
import sys

# functions import
sys.path.append("{}/../../utils".format(os.getcwd()))

# get inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
group = sys.argv[4]

# define directories
jobs_dir = "{}/{}/derivatives/pp_data/jobs".format(main_dir, project_dir)
fs_dir = "{}/{}/derivatives/fmriprep/freesurfer".format(main_dir, project_dir)
fs_license = "{}/{}/code/freesurfer/license.txt".format(main_dir, project_dir)
cortex_dir = "{}/{}/derivatives/pp_data/cortex".format(main_dir, project_dir)
os.makedirs(jobs_dir, exist_ok=True)

# define freesurfer command
freesurfer_cmd = """\
export FREESURFER_HOME={}/{}/code/freesurfer
export SUBJECTS_DIR={}\n\
export FS_LICENSE={}\n\
source $FREESURFER_HOME/SetUpFreeSurfer.sh""".format(main_dir, project_dir, fs_dir, fs_license)

#define pycortex cmd
py_cortex_cmd = "python pycortex_import.py {} {} {} {}".format(main_dir,project_dir,subject,group)

# create sh folder and file
sh_dir = "{}/{}_freesurfer_import_pycortex.sh".format(jobs_dir, subject)

# # Define permission cmd
chmod_cmd = "chmod -Rf 771 {main_dir}/{project_dir}".format(main_dir=main_dir, project_dir=project_dir)
chgrp_cmd = "chgrp -Rf {group} {main_dir}/{project_dir}".format(main_dir=main_dir, project_dir=project_dir, group=group)

of = open(sh_dir, 'w')
of.write("{}\n{}\n{}\n{}".format(freesurfer_cmd,py_cortex_cmd,chmod_cmd,chgrp_cmd))
of.close()

# Run freesurfer and pycortex
os.system("{}".format(sh_dir))