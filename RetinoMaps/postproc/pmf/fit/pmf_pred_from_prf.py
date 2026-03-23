"""
-----------------------------------------------------------------------------------------
pmf_pred_from_prf.py
-----------------------------------------------------------------------------------------
Goal of the script:
Prf fit using gaussian model
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject name
sys.argv[4]: input file name (path to the bold data)
sys.argv[5]: prf fit file name (path to the prf fit data for same subject) 
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
python pmf_pred_from_prf.py /scratch/mszinte/data RetinoMaps sub-03 [file path] [file_path] 32  
-----------------------------------------------------------------------------------------
Written by Martin Szinte (martin.szinte@gmail.com)
and Uriel Lascombes (uriel.lascombes@laposte.net)
edited by Sina Kling (sina.kling@outlook.de)
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

# MRI analysis imports
import nibabel as nb
from prfpy.stimulus import PRFStimulus2D
from prfpy.model import Iso2DGaussianModel 
from prfpy.fit import Iso2DGaussianFitter 


# Personal imports
script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.abspath(os.path.join(script_dir, "../../../../analysis_code/utils")))
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
prf_fit_fn = sys.argv[5]
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

TR = analysis_info['TR']
gauss_grid_nr = analysis_info['gauss_grid_nr']
max_ecc_size = analysis_info['max_ecc_size']
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
    prf_fit_dir = "{}/{}/derivatives/pp_data/{}/170k/pmf/fit".format(
        main_dir, project_dir, subject)
    os.makedirs(prf_fit_dir, exist_ok=True)

elif input_fn.endswith('.gii'):
    prf_fit_dir = "{}/{}/derivatives/pp_data/{}/fsnative/pmf/fit".format(
        main_dir, project_dir, subject)
    os.makedirs(prf_fit_dir, exist_ok=True)

gauss_pred_fn = input_fn.split('/')[-1]
gauss_pred_fn = gauss_pred_fn.replace('bold', 'pmf-gauss_pred')


# Find vdm: check subject-specific directory first, then general vdm directory
vdm_base_dir = '{}/{}/derivatives/vdm'.format(main_dir, project_dir)
vdm_fn_subject = '{}/sub-{}/sub-{}_task-{}_vdm.npy'.format(vdm_base_dir, sub_num, sub_num, prf_task_name)
vdm_fn_general = '{}/task-{}_vdm.npy'.format(vdm_base_dir, prf_task_name)

if os.path.isfile(vdm_fn_subject):
    vdm_fn = vdm_fn_subject
elif os.path.isfile(vdm_fn_general):
    vdm_fn = vdm_fn_general
else:
    raise FileNotFoundError(
        f"No VDM found for task '{prf_task_name}'.\n"
        f"  Checked: {vdm_fn_subject}\n"
        f"  Checked: {vdm_fn_general}"
    )

print(f"Loading VDM from: {vdm_fn}")
vdm = np.load(vdm_fn)
vdm_norm = vdm * (1.0 / (vdm.max() + 1e-8))     # normalize to [0,1] first
vdm = vdm_norm * 100                            # then scale up


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


# load prf fit 
prf_fit_img, prf_fit_data = load_surface(fn=prf_fit_fn)
# extract params from fit data corresponding to: mu_x, mu_Y, prf_size, amplidue, baseline, hrf_1, hrf_2, rsq
prev_gridsearch_params = np.column_stack([prf_fit_data[0, :], prf_fit_data[1, :], prf_fit_data[2, :], prf_fit_data[3, :], prf_fit_data[4, :], prf_fit_data[5, :], prf_fit_data[6, :], prf_fit_data[7, :]])

# Build a param mask
params_valid = (
    np.isfinite(prf_fit_data[0, :]) &                          # mu_x
    np.isfinite(prf_fit_data[1, :]) &                          # mu_y
    np.isfinite(prf_fit_data[2, :]) &                          # size
    (prf_fit_data[3, :] >= prf_amp_th[0]) &                    # amplitude lower bound
    (prf_fit_data[3, :] <= prf_amp_th[1]) &                    # amplitude upper bound
    np.isfinite(prf_fit_data[4, :]) &                          # baseline
    np.isfinite(prf_fit_data[5, :]) &                          # hrf_1
    (prf_fit_data[2, :] > 0)        &                          # size > 0
    (prf_fit_data[7, :] > 0.1)                                 # r²
)


# Intersect with existing valid_vertices_idx
valid_vertices_idx_prf = np.array([
    v for v in valid_vertices_idx if params_valid[v]
])

# determine stimulus
stimulus = PRFStimulus2D(screen_size_cm=screen_size_cm[1],
                         screen_distance_cm=screen_distance_cm,
                         design_matrix=vdm, 
                         TR=TR)

print("\n===== PRF MODEL PARAMETERS =====")
print("Stimulus x min/max (deg):", np.nanmin(stimulus.x_coordinates ), np.nanmax(stimulus.x_coordinates))
print("Stimulus y min/max (deg) :", np.nanmin(stimulus.y_coordinates), np.nanmax(stimulus.y_coordinates))
print("Eccentricity grid range:", np.min(eccs), np.max(eccs))
print("Size grid range:", np.min(sizes), np.max(sizes))
print("==============================\n")

# determine gaussian model
gauss_model = Iso2DGaussianModel(stimulus=stimulus) # make model with SacLoc vdm 

gauss_fit_mat  = np.full((raw_data.shape[1], gauss_params_num), np.nan, dtype=np.float32)
gauss_pred_mat = np.full(raw_data.shape, np.nan, dtype=np.float32)

eps = 1e-6

for vert in valid_vertices_idx_prf:
    params = prf_fit_data[:, vert]
    bold_tc = raw_data[:, vert]  # (TRs,)

    # one voxel at a time
    fitter = Iso2DGaussianFitter(
        data=bold_tc[np.newaxis, :],   # (1, TRs)
        model=gauss_model,
        n_jobs=1
    )

    # inject this voxel's pRF params as the grid result
    fitter.gridsearch_params = np.array([[
        params[0],  # mu_x
        params[1],  # mu_y
        params[2],  # size
        params[3],  # beta     — starting point
        params[4],  # baseline — starting point
        params[5],  # hrf_1
        0.0,        # hrf_2
        params[7],  # rsq
    ]])

    # freeze everything except beta and baseline
    gauss_bounds = [
        (params[0] - eps, params[0] + eps),  # mu_x     frozen
        (params[1] - eps, params[1] + eps),  # mu_y     frozen
        (params[2] - eps, params[2] + eps),  # size     frozen
        (prf_amp_th[0], prf_amp_th[1]),       # beta     free
        (-2, 2),                              # baseline free
        (params[5] - eps, params[5] + eps),  # hrf_1    frozen
        (0, 0),                               # hrf_2    frozen
    ]

    fitter.iterative_fit(
        rsq_threshold=0.1,   
        bounds=gauss_bounds,
        verbose=False
    )

    fit = fitter.iterative_search_params[0]  # (8,)
    gauss_fit_mat[vert] = fit
    gauss_pred_mat[:, vert] = gauss_model.return_prediction(
        mu_x=fit[0], mu_y=fit[1], size=fit[2],
        beta=fit[3], baseline=fit[4]
    ).squeeze()

for i, vert in enumerate(valid_vertices_idx_prf):
    if i % 1000 == 0:
        print(f"  {i} / {len(valid_vertices_idx_prf)} voxels done...")


# prediction step (no fitter needed)
# gauss_pred_mat = np.full(raw_data.shape, np.nan, dtype=np.float32)

# for vert in valid_vertices_idx_prf:
#     params = prf_fit_data[:, vert]  # shape (8,) — index by voxel directly
#     gauss_pred_mat[:, vert] = gauss_model.return_prediction(
#         mu_x=params[0],
#         mu_y=params[1],
#         size=params[2],
#         beta=params[3],
#         baseline=params[4]
#     ).squeeze()


# #only NaN-out voxels that were never predicted
# predicted_mask = np.zeros(raw_data.shape[1], dtype=bool)
# predicted_mask[valid_vertices_idx_prf] = True
# gauss_pred_mat[:, ~predicted_mask] = np.nan

              
# # export pred
# img_gauss_pred_mat = make_surface_image(data=gauss_pred_mat, source_img=img)
# nb.save(img_gauss_pred_mat,'{}/{}'.format(prf_fit_dir, gauss_pred_fn)) 
# print(f"Saved: {prf_fit_dir}/{gauss_pred_fn}")

# # Print duration
# end_time = datetime.datetime.now()
# print("\nStart time:\t{start_time}\nEnd time:\t{end_time}\nDuration:\t{dur}".format(
# start_time=start_time, end_time=end_time, dur=end_time - start_time))