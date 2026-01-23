"""
-----------------------------------------------------------------------------------------
prf_gaussfit.py
-----------------------------------------------------------------------------------------
Goal of the script:
Prf fit using gaussian model
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject name
sys.argv[4]: input file name (path to the data to fit)
sys.argv[5]: number of jobs 
-----------------------------------------------------------------------------------------
Output(s):
fit tester numpy arrays
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/pRF_analysis/analysis_code/postproc/prf/fit
2. run python command
python prf_gaussfit.py [main directory] [project name] [subject name] 
                       [inout file name] [number of jobs]
-----------------------------------------------------------------------------------------
Exemple:
python prf_gaussfit.py /scratch/mszinte/data RetinoMaps sub-03 [file path] 32  
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
sys.path.append("{}/../../../utils".format(os.getcwd()))
from surface_utils import make_surface_image , load_surface
from screen_utils import get_screen_settings
from settings_utils import load_settings

# Get inputs
start_time = datetime.datetime.now()

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
sub_num = subject[4:]
input_fn = sys.argv[4]
n_jobs = int(sys.argv[5])
n_batches = n_jobs
verbose = True
gauss_params_num = 8

# Load settings
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../../"))
settings_path = os.path.join(base_dir, project_dir, "settings.yml")
prf_settings_path = os.path.join(base_dir, project_dir, "prf-analysis.yml")
settings = load_settings([settings_path, prf_settings_path])
analysis_info = settings[0]

TR = analysis_info['TR']
vdm_width = analysis_info['vdm_size_pix'][0] 
vdm_height = analysis_info['vdm_size_pix'][1]
gauss_grid_nr = analysis_info['gauss_grid_nr']
max_ecc_size = analysis_info['max_ecc_size']
rsq_iterative_th = analysis_info['rsq_iterative_th']
size_th = analysis_info['size_th']
prf_amp_th = analysis_info['prf_amp_th']

# Load screen settings from subject dependend task-events.json
prf_task_name = input_fn.split("task-")[1].split("_")[0] # from the file path
screen_size_cm, screen_distance_cm = get_screen_settings(main_dir,project_dir, sub_num, prf_task_name)

print("\n===== PRF FIT PARAMETERS =====")
print(f"Screen Size (cm): {screen_size_cm}")
print(f"Screen Distance (cm): {screen_distance_cm}")
print(f"TR: {TR}")
print(f"Max eccentricity/size values: {max_ecc_size}")
print("==============================\n")

# Define directories
if input_fn.endswith('.nii'):
    prf_fit_dir = "{}/{}/derivatives/pp_data/{}/170k/prf/fit".format(
        main_dir, project_dir, subject)
    os.makedirs(prf_fit_dir, exist_ok=True)

elif input_fn.endswith('.gii'):
    prf_fit_dir = "{}/{}/derivatives/pp_data/{}/fsnative/prf/fit".format(
        main_dir, project_dir, subject)
    os.makedirs(prf_fit_dir, exist_ok=True)

gauss_fit_fn  = input_fn.split('/')[-1]
gauss_fit_fn = gauss_fit_fn.replace('bold', 'prf-gauss_fit')

gauss_pred_fn = input_fn.split('/')[-1]
gauss_pred_fn = gauss_pred_fn.replace('bold', 'prf-gauss_pred')

# Get task specific visual design matrix
# TODO add option to have run specific vdm
# if vdm per run setting is true, check for different path
vdm_fn = '{}/{}/derivatives/vdm/task-{}_vdm.npy'.format(
    main_dir, project_dir, prf_task_name)
vdm = np.load(vdm_fn)

# define model parameter grid range
sizes = max_ecc_size * np.linspace(0.1, 1, gauss_grid_nr) ** 2
eccs = max_ecc_size * np.linspace(0.1, 1, gauss_grid_nr) ** 2
polars = np.linspace(0, 2*np.pi, gauss_grid_nr)

# load data
img, raw_data = load_surface(fn=input_fn)

# exlude nan voxel from the analysis 
valid_vertices = ~np.isnan(raw_data).any(axis=0)
valid_vertices_idx = np.where(valid_vertices)[0]
data = raw_data[:,valid_vertices]

# determine stimulus
stimulus = PRFStimulus2D(screen_size_cm=screen_size_cm[1],
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

# Iterative fit gauss model
gauss_bounds = [(-max_ecc_size, max_ecc_size), # x
                (-max_ecc_size, max_ecc_size), # y
                (size_th[0], size_th[1]), # prf size
                (prf_amp_th[0], prf_amp_th[1]), # prf amplitude
                (-2, 2), # bold baseline
                (0, 10), # hrf1
                (0, 0) # hrf2
                ]

print("\n===== PRF FIT BOUNDS =====")
print("Gauss bounds:")
print("  x range:", gauss_bounds[0])
print("  y range:", gauss_bounds[1])
print("  prf size range:", gauss_bounds[2])
print("  prf amplitude range:", gauss_bounds[3])
print("  bold baseline range:", gauss_bounds[4])
print("  hrf1 range:", gauss_bounds[5])
print("  hrf2 range:", gauss_bounds[6])

gauss_fitter.iterative_fit(rsq_threshold=rsq_iterative_th, 
                           bounds=gauss_bounds,
                           verbose=verbose)
gauss_fit = gauss_fitter.iterative_search_params

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
img_gauss_fit_mat = make_surface_image(data=gauss_fit_mat.T, source_img=img, maps_names=maps_names)
nb.save(img_gauss_fit_mat,'{}/{}'.format(prf_fit_dir, gauss_fit_fn)) 

# export pred
img_gauss_pred_mat = make_surface_image(data=gauss_pred_mat, source_img=img)
nb.save(img_gauss_pred_mat,'{}/{}'.format(prf_fit_dir, gauss_pred_fn)) 

# Print duration
end_time = datetime.datetime.now()
print("\nStart time:\t{start_time}\nEnd time:\t{end_time}\nDuration:\t{dur}".format(
start_time=start_time, end_time=end_time, dur=end_time - start_time))