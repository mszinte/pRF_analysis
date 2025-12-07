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
sys.argv[5]: session name (optional, e.g. ses-01)
-----------------------------------------------------------------------------------------
Output(s):
None
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/pRF_analysis/analysis_code/preproc/functional/
2. run python command
python freesurfer_import_pycortex.py [main directory] [project name] [subject] [group] [session (optional)]
-----------------------------------------------------------------------------------------
Executions:
cd ~/projects/pRF_analysis/analysis_code/preproc/functional/
python freesurfer_import_pycortex.py /scratch/mszinte/data MotConf sub-01 327
python freesurfer_import_pycortex.py /scratch/mszinte/data amblyo7T_prf sub-01 327 ses-01
-----------------------------------------------------------------------------------------
Written by Martin Szinte (martin.szinte@gmail.com)
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
import urllib.request
import json

# functions import
sys.path.append("{}/../../utils".format(os.getcwd()))

# get inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
group = sys.argv[4]
session = sys.argv[5] if len(sys.argv) > 5 else None

# Handle session parameter
if session:
    subject_name = f"{subject}_{session}"
else:
    subject_name = subject

# define directories
jobs_dir = "{}/{}/derivatives/pp_data/jobs".format(main_dir, project_dir)
fs_dir = "{}/{}/derivatives/fmriprep/freesurfer".format(main_dir, project_dir)
fs_license = "{}/{}/code/freesurfer/license.txt".format(main_dir, project_dir)
cortex_dir = "{}/{}/derivatives/pp_data/cortex".format(main_dir, project_dir)
os.makedirs(jobs_dir, exist_ok=True)

# Check and setup pycortex directory structure
def setup_pycortex_dirs(cortex_dir):
    """Check for cortex/colormaps and cortex/db folders, create if missing and download colormaps"""
    colormaps_dir = os.path.join(cortex_dir, "colormaps")
    db_dir = os.path.join(cortex_dir, "db")
    
    # Create directories if they don't exist
    os.makedirs(colormaps_dir, exist_ok=True)
    os.makedirs(db_dir, exist_ok=True)
    
    # Check if colormaps directory is empty
    if not os.listdir(colormaps_dir):
        print("Downloading colormaps from GitHub...")
        # GitHub API URL to list files in the colormaps directory
        api_url = "https://api.github.com/repos/gallantlab/pycortex/contents/filestore/colormaps"
        
        try:
            with urllib.request.urlopen(api_url) as response:
                files = json.loads(response.read())
            
            # Download each colormap file
            for file_info in files:
                if file_info['type'] == 'file':
                    file_url = file_info['download_url']
                    file_name = file_info['name']
                    file_path = os.path.join(colormaps_dir, file_name)
                    
                    print(f"  Downloading {file_name}...")
                    urllib.request.urlretrieve(file_url, file_path)
            
            print("Colormaps downloaded successfully.")
        except Exception as e:
            print(f"Warning: Could not download colormaps: {e}")
    else:
        print("Colormaps directory already contains files.")

# Setup pycortex directories
setup_pycortex_dirs(cortex_dir)

# define freesurfer command
freesurfer_cmd = """\
export FREESURFER_HOME={}/{}/code/freesurfer
export SUBJECTS_DIR={}\n\
export FS_LICENSE={}\n\
source $FREESURFER_HOME/SetUpFreeSurfer.sh""".format(main_dir, project_dir, fs_dir, fs_license)

#define pycortex cmd
if session:
    py_cortex_cmd = "python pycortex_import.py {} {} {} {} {}".format(main_dir, project_dir, subject, group, session)
else:
    py_cortex_cmd = "python pycortex_import.py {} {} {} {}".format(main_dir, project_dir, subject, group)

# create sh folder and file
sh_dir = "{}/{}_freesurfer_import_pycortex.sh".format(jobs_dir, subject_name)

# # Define permission cmd
chmod_cmd = "chmod -Rf 771 {main_dir}/{project_dir}".format(main_dir=main_dir, project_dir=project_dir)
chgrp_cmd = "chgrp -Rf {group} {main_dir}/{project_dir}".format(main_dir=main_dir, project_dir=project_dir, group=group)

of = open(sh_dir, 'w')
of.write("{}\n{}\n{}\n{}".format(freesurfer_cmd,py_cortex_cmd,chmod_cmd,chgrp_cmd))
of.close()

# Run freesurfer and pycortex
os.system("{}".format(sh_dir))
