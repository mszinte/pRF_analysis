"""
-----------------------------------------------------------------------------------------
pycortex_import.py
-----------------------------------------------------------------------------------------
Goal of the script:
Import subject in pycortex database
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
>> cd ~/projects/[PROJECT]/analysis_code/preproc/functional/
2. run python command
python pycortex_import.py [main directory] [project name] [subject] [group] [session (optional)]
-----------------------------------------------------------------------------------------
Executions:
cd ~/projects/pRF_analysis/analysis_code/preproc/functional/
python pycortex_import.py /scratch/mszinte/data MotConf sub-01 327
python pycortex_import.py /scratch/mszinte/data amblyo7T_prf sub-01 327 ses-01
-----------------------------------------------------------------------------------------
Written by Martin Szinte (mail@martinszinte.net)
and Uriel Lascombes (uriel.lascombes@laposte.net)
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
import cortex
import importlib
import numpy as np

# functions import
sys.path.append("{}/../../utils".format(os.getcwd()))
from pycortex_utils import set_pycortex_config_file

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

# define directories and get fns
fmriprep_dir = "{}/{}/derivatives/fmriprep".format(main_dir, project_dir)
fs_dir = "{}/{}/derivatives/fmriprep/freesurfer".format(main_dir, project_dir)
cortex_dir = "{}/{}/derivatives/pp_data/cortex".format(main_dir, project_dir)
temp_dir = "{}/{}/derivatives/temp_data/{}_rand_ds/".format(main_dir, project_dir, subject_name)

# set pycortex db and colormaps
set_pycortex_config_file(cortex_dir)
importlib.reload(cortex)

# add participant to pycortex db
print('import subject in pycortex')
cortex.freesurfer.import_subj(subject_name, subject_name, fs_dir, 'smoothwm')

# add participant flat maps
print('import subject flatmaps')
try: cortex.freesurfer.import_flat(fs_subject=subject_name, cx_subject=subject_name, 
                                  freesurfer_subject_dir=fs_dir, patch='full', auto_overwrite=True)
except: pass

# create participant pycortex overlays
print('create subject pycortex overlays to check')
surfs = [cortex.polyutils.Surface(*d) for d in cortex.db.get_surf(subject_name, "fiducial")]
num_verts = surfs[0].pts.shape[0] + surfs[1].pts.shape[0]
rand_data = np.random.randn(num_verts)
vertex_data = cortex.Vertex(rand_data, subject_name)
ds = cortex.Dataset(rand=vertex_data)
cortex.webgl.make_static(outpath=temp_dir, data=ds)

# # Define permission cmd
# os.system("chmod -Rf 771 {main_dir}/{project_dir}".format(main_dir=main_dir, project_dir=project_dir))
# os.system("chgrp -Rf {group} {main_dir}/{project_dir}".format(main_dir=main_dir, project_dir=project_dir, group=group))
