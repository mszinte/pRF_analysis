"""
-----------------------------------------------------------------------------------------
prf_gridfit.py
-----------------------------------------------------------------------------------------
Goal of the script:
Prf fit computing gaussian grid fit
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject name
sys.argv[4]: input file name (path to the data to fit)
sys.argv[5]: number of jobs 
sys.argv[6]: OPTIONAL main analysis folder (e.g. prf_em_ctrl)
-----------------------------------------------------------------------------------------
Output(s):
fit tester numpy arrays
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/pRF_analysis/amblyo7T_prf/postproc/prf/fit/
2. run python command
python prf_gridfit.py [main directory] [project name] [subject name] 
[inout file name] [number of jobs] [analysis folder - optional]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/projects/pRF_analysis/amblyo7T_prf/postproc/prf/fit/
python prf_gridfit.py /scratch/mszinte/data amblyo7T_prf sub-01 
    /scratch/mszinte/data/amblyo7T_prf/derivatives/pp_data/sub-01/fsnative/func/fmriprep_dct_concat/sub-01_task-pRFRightEye_hemi-R_dct_concat_bold.func.gii
    8 pRFRightEye 
-----------------------------------------------------------------------------------------
Written by Martin Szinte (martin.szinte@gmail.com)
and Uriel Lascombes (uriel.lascombes@laposte.net)
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
import json
import datetime
import numpy as np

# MRI analysis imports
import nibabel as nb
from prfpy.stimulus import PRFStimulus2D
from prfpy.model import Iso2DGaussianModel 
from prfpy.fit import Iso2DGaussianFitter 

# Personal imports
sys.path.append("{}/../../../../analysis_code/utils".format(os.getcwd()))
from surface_utils import make_surface_image , load_surface
from screen_utils import get_screen_settings

# Get inputs
start_time = datetime.datetime.now()

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
sub_num = subject[4:]
input_fn = sys.argv[4]
n_jobs = int(sys.argv[5])
if len(sys.argv) > 6: output_folder = sys.argv[6]
else: output_folder = "prf"
n_batches = n_jobs
verbose = True
gauss_params_num = 8

# Analysis parameters
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../../"))
settings_path = os.path.join(base_dir, project_dir, 'settings.json')

with open(settings_path) as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)

# Load screen settings from subject dependend task-events.json
screen_size_cm, screen_distance_cm = get_screen_settings(main_dir,project_dir, sub_num, output_folder)

TR = analysis_info['TR']
vdm_width = analysis_info['vdm_size_pix'][0] 
vdm_height = analysis_info['vdm_size_pix'][1]
gauss_grid_nr = analysis_info['gauss_grid_nr']
max_ecc_size = analysis_info['max_ecc_size']
prf_task_name = output_folder

# Define directories
if input_fn.endswith('.nii'):
    prf_fit_dir = "{}/{}/derivatives/pp_data/{}/170k/{}/fit".format(
        main_dir, project_dir, subject, output_folder)
    os.makedirs(prf_fit_dir, exist_ok=True)

elif input_fn.endswith('.gii'):
    prf_fit_dir = "{}/{}/derivatives/pp_data/{}/fsnative/{}/fit".format(
        main_dir, project_dir, subject, output_folder)
    os.makedirs(prf_fit_dir, exist_ok=True)

fit_fn_gauss_gridfit = input_fn.split('/')[-1]
fit_fn_gauss_gridfit = fit_fn_gauss_gridfit.replace('bold', 'prf-fit_gauss_gridfit')

pred_fn_gauss_gridfit = input_fn.split('/')[-1]
pred_fn_gauss_gridfit = pred_fn_gauss_gridfit.replace('bold', 'prf-pred_gauss_gridfit')

# Get task specific visual design matrix
vdm_fn = '{}/{}/derivatives/vdm/vdm_{}_{}_{}.npy'.format(
    main_dir, project_dir, prf_task_name, vdm_width, vdm_height)
vdm = np.load(vdm_fn)

# defind model parameter grid range
sizes = max_ecc_size * np.linspace(0.1,1,gauss_grid_nr)**2
eccs = max_ecc_size * np.linspace(0.1,1,gauss_grid_nr)**2
polars = np.linspace(0, 2*np.pi, gauss_grid_nr)

print("\n===== PRF FIT PARAMETERS =====")
print(f"Screen Size (cm): {screen_size_cm}")
print(f"Screen Distance (cm): {screen_distance_cm}")
print(f"TR: {TR}")
print(f"Max eccentricity/size values: {max_ecc_size}")
print("==============================\n")

# load data
img, raw_data = load_surface(fn=input_fn)

# exlude nan voxel from the analysis 
valid_vertices = ~np.isnan(raw_data).any(axis=0)
valid_vertices_idx = np.where(valid_vertices)[0]
data = raw_data[:,valid_vertices]

# determine stimulus
stimulus = PRFStimulus2D(screen_size_cm=screen_size_cm[0],
                         screen_distance_cm=screen_distance_cm,
                         design_matrix=vdm, 
                         TR=TR)

print("\n===== PRF MODEL PARAMETERS =====")
print("Stimulus x min/max (deg):", np.nanmin(stimulus.x_coordinates ), np.nanmax(stimulus.x_coordinates ))
print("Stimulus y min/max (deg) :", np.nanmin(stimulus.y_coordinates), np.nanmax(stimulus.y_coordinates))
print("Eccentricity grid range:", np.min(eccs), np.max(eccs))
print("Size grid range:", np.min(sizes), np.max(sizes))
print("==============================\n")


# determine gaussian model
gauss_model = Iso2DGaussianModel(stimulus=stimulus)

# grid fit gauss model
gauss_fitter = Iso2DGaussianFitter(data=data.T, 
                                   model=gauss_model, 
                                   n_jobs=n_jobs)

gauss_fitter.grid_fit(ecc_grid=eccs, 
                      polar_grid=polars, 
                      size_grid=sizes, 
                      verbose=verbose, 
                      n_batches=n_batches)

# rearange result of Gauss model 
gauss_fit = gauss_fitter.gridsearch_params
gauss_fit_mat = np.zeros((raw_data.shape[1],gauss_params_num))
gauss_pred_mat = np.zeros_like(raw_data) 

for est,vert in enumerate(valid_vertices_idx):
    gauss_fit_mat[vert] = gauss_fit[est]
    gauss_pred_mat[:,vert] = gauss_model.return_prediction(mu_x=gauss_fit[est][0], 
                                                          mu_y=gauss_fit[est][1], 
                                                          size=gauss_fit[est][2], 
                                                          beta=gauss_fit[est][3], 
                                                          baseline=gauss_fit[est][4],
                                                          hrf_1=gauss_fit[est][5],
                                                          hrf_2=gauss_fit[est][6])

gauss_fit_mat = np.where(gauss_fit_mat == 0, np.nan, gauss_fit_mat)
gauss_pred_mat = np.where(gauss_pred_mat == 0, np.nan, gauss_pred_mat)

#export data from gauss model fit
maps_names = ['mu_x', 'mu_y', 'prf_size', 'prf_amplitude', 'bold_baseline', 
              'hrf_1','hrf_2', 'r_squared']
              
# export fit
img_gauss_gridfit_fit_mat = make_surface_image(data=gauss_fit_mat.T, source_img=img, maps_names=maps_names)
nb.save(img_gauss_gridfit_fit_mat,'{}/{}'.format(prf_fit_dir, fit_fn_gauss_gridfit)) 

# export pred
img_gauss_gridfit_pred_mat = make_surface_image(data=gauss_pred_mat, source_img=img)
nb.save(img_gauss_gridfit_pred_mat,'{}/{}'.format(prf_fit_dir, pred_fn_gauss_gridfit)) 

# Print duration
end_time = datetime.datetime.now()
print("\nStart time:\t{start_time}\nEnd time:\t{end_time}\nDuration:\t{dur}".format(
start_time=start_time, end_time=end_time, dur=end_time - start_time))