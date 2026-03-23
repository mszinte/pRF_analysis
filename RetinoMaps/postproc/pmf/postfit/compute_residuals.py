"""
-----------------------------------------------------------------------------------------
compute_residuals.py
-----------------------------------------------------------------------------------------
Goal of the script:
Compute residuals of voxels from pmf pred
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject name
sys.argv[4]: input file name (path to the bold data)
sys.argv[5]: prf pred file name (path to the prf fit data for same subject) 
sys.argv[6]: number of jobs
-----------------------------------------------------------------------------------------
Output(s):
fit tester numpy arrays
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/pRF_analysis/analysis_code/postproc/prf/fit
2. run python command
python pmf_pred_from_prf.py [main directory] [project name] [subject name] 
                       [inout file name] [prf fit file] [number of jobs]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/projects/pRF_analysis/RetinoMaps/postproc/pmf/fit
python compute_residuals.py /scratch/mszinte/data RetinoMaps sub-03 [file path] [file_path] 32  
-----------------------------------------------------------------------------------------
Written by Sina Kling (sina.kling@outlook.de)
-----------------------------------------------------------------------------------------
"""
# Stop warnings
import warnings
warnings.filterwarnings("ignore")
# Debug
import ipdb
deb = ipdb.set_trace
# General imports
import os
import sys
import datetime
import numpy as np
import nibabel as nb
# Personal imports
script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.abspath(os.path.join(script_dir, "../../../../analysis_code/utils")))
from surface_utils import make_surface_image, load_surface
from screen_utils import get_screen_settings
from settings_utils import load_settings

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
sub_num = subject[4:]
input_fn = sys.argv[4]
prf_pred_fn = sys.argv[5]
n_jobs = int(sys.argv[6])
n_batches = n_jobs
verbose = True
gauss_params_num = 8

# Load settings
base_dir = os.path.abspath(os.path.join(script_dir, "../../../../"))
settings_path = os.path.join(base_dir, project_dir, "pmf-settings.yml")
prf_settings_path = os.path.join(base_dir, project_dir, "prf-analysis.yml")
settings = load_settings([settings_path, prf_settings_path])
analysis_info = settings[0]

# Define directories and output filename
if input_fn.endswith('.nii'):
    residuals_dir = "{}/{}/derivatives/pp_data/{}/170k/pmf/pmf_residuals".format(
        main_dir, project_dir, subject)
elif input_fn.endswith('.gii'):
    residuals_dir = "{}/{}/derivatives/pp_data/{}/fsnative/pmf/pmf_residuals".format(
        main_dir, project_dir, subject)

os.makedirs(residuals_dir, exist_ok=True)

residual_fn = input_fn.split('/')[-1]
residual_fn = residual_fn.replace('bold', 'pmf-residuals')

# Load data
print('Loading pRF prediction data: {}'.format(prf_pred_fn))
prf_pred_img, prf_pred_data = load_surface(fn=prf_pred_fn)

print('Loading BOLD data: {}'.format(input_fn))
bold_img, bold_data = load_surface(fn=input_fn)

print('prf_pred_data shape: {}'.format(prf_pred_data.shape))
print('bold_data shape: {}'.format(bold_data.shape))

# NaN mask: keep only voxels clean in both datasets
nan_mask = (
    np.isnan(bold_data).any(axis=0) |
    np.isnan(prf_pred_data).any(axis=0)
)
clean_vox = np.where(~nan_mask)[0]

# Compute residuals on clean voxels only
residual_data = np.full_like(bold_data, np.nan)
residual_data[:, clean_vox] = bold_data[:, clean_vox] - prf_pred_data[:, clean_vox]

# Export residuals
residual_img = make_surface_image(data=residual_data, source_img=bold_img)
residual_fullfn = '{}/{}'.format(residuals_dir, residual_fn)

print('Saving residuals: {}'.format(residual_fullfn))
nb.save(residual_img, residual_fullfn)

print('Done: residuals saved to {}'.format(residual_fullfn))