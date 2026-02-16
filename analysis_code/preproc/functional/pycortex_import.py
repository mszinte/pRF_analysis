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
sys.argv[3]: freesurfer subject name (e.g. sub-01 or sub-01_ses-01)
sys.argv[4]: pycortex subject name (e.g. sub-01)
sys.argv[5]: group (e.g. 327)
-----------------------------------------------------------------------------------------
Output(s):
None
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/pRF_analysis/analysis_code/preproc/functional/
2. run python command
python pycortex_import.py [main directory] [project name] [fs_subject] 
                          [pycortex_subject] [group]
-----------------------------------------------------------------------------------------
Executions:
cd ~/projects/pRF_analysis/analysis_code/preproc/functional/
python pycortex_import.py /scratch/mszinte/data MotConf sub-01_ses-01 sub-01 327
python pycortex_import.py /scratch/mszinte/data RetinoMaps sub-01 sub-01 327
-----------------------------------------------------------------------------------------
Written by Martin Szinte (martin.szinte@gmail.com)
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
from pycortex_utils import set_pycortex_config_file, setup_pycortex_dirs
from settings_utils import load_settings

# get inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
fs_subject = sys.argv[3]
cx_subject = sys.argv[4]
group = sys.argv[5]

# Load input
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../"))
settings_path = os.path.join(base_dir, project_dir, "settings.yml")
settings = load_settings([settings_path])
analysis_info = settings[0]
flattening_method = analysis_info['flattening_method']
if flattening_method == 'autoflatten':
    patch_name = 'autoflatten'
elif flattening_method == 'mrisflatten':
    patch_name = 'full'

# define directories and get fns
fmriprep_dir = "{}/{}/derivatives/fmriprep".format(main_dir, project_dir)
fs_dir = "{}/{}/derivatives/fmriprep/freesurfer".format(main_dir, project_dir)
cortex_dir = "{}/{}/derivatives/pp_data/cortex".format(main_dir, project_dir)
temp_dir = "{}/{}/derivatives/temp_data/{}_rand_ds/".format(main_dir, project_dir, cx_subject)

# Setup pycortex directories
setup_pycortex_dirs(cortex_dir)

# Set pycortex db and colormaps
set_pycortex_config_file(cortex_dir)
importlib.reload(cortex)

# add participant to pycortex db
print('import subject in pycortex')
cortex.freesurfer.import_subj(fs_subject, cx_subject, fs_dir, 'smoothwm')
# cortex.freesurfer.import_subj(fs_subject, cx_subject, fs_dir, 'white')

# add participant flat maps
print('import subject flatmaps')
try: cortex.freesurfer.import_flat(fs_subject=fs_subject, cx_subject=cx_subject, 
                                  freesurfer_subject_dir=fs_dir, patch=patch_name, auto_overwrite=True)
except Exception as e:
    print(f"Warning: Flatmap import failed → {e}")

# create participant pycortex overlays
print('create subject pycortex overlays to check')
surfs = [cortex.polyutils.Surface(*d) for d in cortex.db.get_surf(cx_subject, "fiducial")]
num_verts = surfs[0].pts.shape[0] + surfs[1].pts.shape[0]
rand_data = np.random.randn(num_verts)
vertex_data = cortex.Vertex(rand_data, cx_subject)
ds = cortex.Dataset(rand=vertex_data)
cortex.webgl.make_static(outpath=temp_dir, data=ds)
