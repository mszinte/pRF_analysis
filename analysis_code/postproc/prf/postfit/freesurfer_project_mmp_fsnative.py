"""
-----------------------------------------------------------------------------------------
freesurfer_project_mmp_fsnative.py
-----------------------------------------------------------------------------------------
Goal of the script:
Load freesurfer and use it to project MMP from fsaverage to fsnative 
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: freesurfer subject name (e.g. sub-01)
sys.argv[4]: group of shared data (e.g. 327)
-----------------------------------------------------------------------------------------
Output(s):
None
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/pRF_analysis/analysis_code/postproc/prf/postfit/
2. run python command
python freesurfer_project_mmp_fsnative.py [main directory] [project name] [subject name] [group]
-----------------------------------------------------------------------------------------
Executions:
cd ~/projects/pRF_analysis/analysis_code/postproc/prf/postfit/
python freesurfer_project_mmp_fsnative.py /scratch/mszinte/data RetinoMaps sub-01 327
python freesurfer_project_mmp_fsnative.py /scratch/mszinte/data amblyo7T_prf sub-01_ses-01 327
-----------------------------------------------------------------------------------------
Written by Uriel Lascombes (uriel.lascombes@laposte.net)
-----------------------------------------------------------------------------------------
"""

# stop warnings
import warnings
warnings.filterwarnings("ignore")

# Debug 
import ipdb
deb = ipdb.set_trace

# General imports
import os
import sys
import subprocess
from pathlib import Path
import urllib.request

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
group = sys.argv[4]

# Directories
fs_dir = "{}/{}/derivatives/fmriprep/freesurfer".format(main_dir, project_dir)
fs_license = "{}/{}/code/freesurfer/license.txt".format(main_dir, project_dir)
fs_average_dir = '{}/fsaverage'.format(fs_dir)
lh_mmp_fn = Path('{}/label/lh.HCPMMP1.annot'.format(fs_average_dir))
rh_mmp_fn = Path('{}/label/rh.HCPMMP1.annot'.format(fs_average_dir))

# Check if annot files existes in fressurfer fsaverage folder 
if lh_mmp_fn.is_file() and rh_mmp_fn.is_file():
    print("Atlas found")
else: 
    print('\n\nAtlas files not found. Please download manually:')
    print('1. Go to: https://figshare.com/articles/dataset/HCP-MMP1_0_projected_on_fsaverage/3498446')
    print('2. Download lh.HCPMMP1.annot and rh.HCPMMP1.annot')
    print(f'3. Place them in: {lh_mmp_fn.parent}')
    print('4. Rerun this script\n')
    quit()
      
# Createcommand and Sh file 
jobs_dir = "{}/{}/derivatives/pp_data/jobs".format(main_dir, project_dir)
os.makedirs(jobs_dir, exist_ok=True)
sh_dir = "{}/{}_project_mmp_fsnative.sh".format(jobs_dir, subject)
hemis = ['lh', 'rh']

# define freesurfer command
freesurfer_cmd = """\
export FREESURFER_HOME={}/{}/code/freesurfer
export SUBJECTS_DIR={}\n\
export FS_LICENSE={}\n\
source $FREESURFER_HOME/SetUpFreeSurfer.sh\n\
""".format(main_dir, project_dir, fs_dir, fs_license)
of = open(sh_dir, 'w')
of.write("{}".format(freesurfer_cmd))
   
# Freesurfer cmd for both hemi
for hemi in hemis:
    # surf2surf command 
    surf2surf = """\
    mri_surf2surf \
        --srcsubject fsaverage \
        --trgsubject {subject} \
        --hemi {hemi} \
        --sval-annot $SUBJECTS_DIR/fsaverage/label/{hemi}.HCPMMP1.annot \
        --tval $SUBJECTS_DIR/{subject}/label/{hemi}.HCPMMP1.annot
    """.format(subject=subject, hemi=hemi)
    of.write("\n{}".format(surf2surf))

of.close()

# Run commands
subprocess.run(["bash", sh_dir], check=True)

# Try to set permissions (may fail without sudo)
try:
    print('\nChanging files permissions...')
    subprocess.run(["chmod", "-Rf", "771", f"{main_dir}/{project_dir}"])
    subprocess.run(["chgrp", "-Rf", group, f"{main_dir}/{project_dir}"])
except subprocess.CalledProcessError:
    print("Warning: Could not set permissions/group (may require admin rights)")