# Skeleton python script for prfmodel css fit 
"""
-----------------------------------------------------------------------------------------
prf_cssfit_gpu.py
-----------------------------------------------------------------------------------------
Goal of the script:
Prf fit computing css fit with gpu usage using prfmodel toolbox: 
grid (polar model) → least squares (= β/baseline) 
→ iterative Gaussian fit (bounded SGD) →  CSS SGD
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject name
sys.argv[4]: input file name (path to the data to fit)
sys.argv[5]: analysis task name (ex. prf)
sys.argv[6]: number of jobs 
-----------------------------------------------------------------------------------------
Output(s):
fit 
prediction 
-----------------------------------------------------------------------------------------
To run:
conda activate prfmodel
1. cd to function
>> cd ~/projects/pRF_analysis/analysis_code/postproc/prf/fit
2. run python command
python prf_cssfit_gpu.py [main directory] [project name] [subject name] 
                     [input file name] [analysis task name] [number of jobs]
-----------------------------------------------------------------------------------------
Exemple:
for direct gpu test on mesocentre: 
srun -p volta --gres=gpu:1 --cpus-per-task=4 --time=02:00:00 --pty bash -i
python prf_cssfit_gpu.py /scratch/mszinte/data amsterdam24 sub-03 /scratch/mszinte/data/amsterdam24/derivatives/pp_data/sub-03/fsnative/func/fmriprep_dct_z-score_avg/sub-03_task-pRF_hemi-L_fmriprep_dct_z-score_avg_bold.func.gii prf 32  
python prf_cssfit_gpu.py /scratch/mszinte/data amsterdam24 sub-03 /scratch/mszinte/data/amsterdam24/derivatives/pp_data/sub-03/170k/func/fmriprep_dct_z-score_avg/sub-03_task-pRF_fmriprep_dct_z-score_avg_bold.dtseries.nii prf 32

-----------------------------------------------------------------------------------------
Written by Sina Kling (sina.kling@outlook.de)
adapted from prf_css_fit.py by Uriel Lascombes 
-----------------------------------------------------------------------------------------
"""
# Stop warnings
import warnings
warnings.filterwarnings("ignore")

# Debug
import ipdb
deb = ipdb.set_trace

import os
from importlib.util import find_spec
import pandas as pd
import numpy as np
import datetime
import sys
from keras import ops
import nibabel as nb
import keras

#os.environ["KERAS_BACKEND"] = "torch" #only for local usage
from prfmodel.stimuli import PRFStimulus
from prfmodel.impulse import DerivativeTwoGammaImpulse
from prfmodel.models.prf import Gaussian2DPRFModel
from prfmodel.fitters import GridFitter
from prfmodel.fitters import LeastSquaresFitter
from prfmodel.models.prf import init_css_from_gaussian
from prfmodel.fitters import SGDFitter
from prfmodel.fitters.adapter import Adapter, ParameterTransform
from prfmodel.utils import batched
from prfmodel.models.compression import CompressiveEncoder
from prfmodel.models.prf import PRFStimulusEncoder
import tensorflow as tf
from prfmodel.fitters.losses import CorrelationLoss

tf.config.list_physical_devices('GPU')
for gpu in tf.config.list_physical_devices("GPU"):
    tf.config.experimental.set_memory_growth(gpu, True)

# Personal imports
sys.path.append("{}/../../../utils".format(os.getcwd()))
from maths_utils import r2_score_surf
from settings_utils import load_settings
from screen_utils import get_screen_settings
from pycortex_utils import set_pycortex_config_file
from surface_utils import load_surface ,make_surface_image
from prfmodel.models.prf.canonical import CanonicalPRFModel
from prfmodel_utils import PolarGaussian2DPRFTuning, polar_to_cartesian, bounded_transform, clip_to_bounds

# Get inputs
start_time = datetime.datetime.now()

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
sub_num = subject[4:]
input_fn = sys.argv[4]
analysis_name = sys.argv[5]
n_jobs = int(sys.argv[6])

n_batches = n_jobs
verbose = True
css_params_num = 9

# Load settings
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../../"))
general_settings_path = os.path.join(base_dir, project_dir, "settings.yml")
analysis_settings_path = os.path.join(base_dir, project_dir, f"{analysis_name}-analysis.yml")
settings = load_settings([general_settings_path, analysis_settings_path])
analysis_info = settings[0]

tr = analysis_info['TR']
gauss_grid_nr = analysis_info['gauss_grid_nr']
max_ecc_size = analysis_info['max_ecc_size']
css_grid_nr = analysis_info['css_grid_nr']
css_exponent_bounds = analysis_info['n_th']

#NEW SETTINGS 
#sgd_num_steps = #1000 
#batch_size = = #100?

# Load screen settings from subject dependent task-events.json
task_name = input_fn.split("task-")[1].split("_")[0]  # from the file path
screen_size_cm, screen_distance_cm = get_screen_settings(main_dir, project_dir, sub_num, task_name)

# Override with fake screen settings if defined in analysis yml
if 'fake_screen_size_cm' in analysis_info:
    print("[INFO] Fake screen size found in settings — overriding real screen size.")
    screen_size_cm = analysis_info['fake_screen_size_cm']

print(f"Fitting data: {input_fn}\n")
print("\n===== PRF FIT PARAMETERS =====")
print(f"fit : {input_fn}")
print(f"Screen Size (cm): {screen_size_cm} {'[FAKE]' if 'fake_screen_size_cm' in analysis_info else '[REAL]'}")
print(f"Screen Distance (cm): {screen_distance_cm}")
print(f"TR: {tr}")
print(f"Max eccentricity/size values: {max_ecc_size}") 
print(f"Inputing x/y = [-{max_ecc_size * 0.5}, {max_ecc_size * 0.5}]")
print("==============================\n")
# Set pycortex db and colormaps
cortex_dir = "{}/{}/derivatives/pp_data/cortex".format(main_dir, project_dir)
set_pycortex_config_file(cortex_dir)

# Get task specific (visual) design matrix
# Find dm: check subject-specific directory first, then general vdm directory
dm_name = analysis_info['dm_name']
dm_name_tofind = dm_name[1:] if dm_name else "vdm"
dm_base_dir = '{}/{}/derivatives/vdm'.format(main_dir, project_dir)
dm_fn_subject = '{}/sub-{}/sub-{}_task-{}_{}.npy'.format(dm_base_dir, sub_num, sub_num, task_name, dm_name_tofind)
dm_fn_general = '{}/task-{}_{}.npy'.format(dm_base_dir, task_name, dm_name_tofind)


if os.path.isfile(dm_fn_subject):
    dm_fn = dm_fn_subject
elif os.path.isfile(dm_fn_general):
    dm_fn = dm_fn_general
else:
    raise FileNotFoundError(
        f"No DM found for task '{task_name}'.\n"
        f"  Checked: {dm_fn_subject}\n"
        f"  Checked: {dm_fn_general}"
    )

print(f"Loading DM from: {dm_fn}")
dm = np.load(dm_fn)

# Define directories and files names (fn)
output_folder = analysis_info["output_folder"]
if input_fn.endswith('.nii'):
    prf_fit_dir = "{}/{}/derivatives/pp_data/{}/170k/{}/fit".format(
        main_dir, project_dir, subject, output_folder)
    os.makedirs(prf_fit_dir, exist_ok=True)

elif input_fn.endswith('.gii'):
    prf_fit_dir = "{}/{}/derivatives/pp_data/{}/fsnative/{}/fit".format(
        main_dir, project_dir, subject,  output_folder)
    os.makedirs(prf_fit_dir, exist_ok=True)

css_fit_fn = input_fn.split('/')[-1]
css_fit_fn = css_fit_fn.replace('bold', f'{analysis_name}-css{dm_name}_fit') 

css_pred_fn = input_fn.split('/')[-1] 
css_pred_fn = css_pred_fn.replace('bold', f'{analysis_name}-css{dm_name}_pred')

#------------LOAD DATA----------------
img, data = load_surface(fn=input_fn) #shape: (n times, n vertices)
# Exclude vertices with all-NaN timeseries to avoid errors during fitting
valid_vertices = ~np.isnan(data).any(axis=0)
valid_vertices_idx = np.where(valid_vertices)[0]
# Filter data to only include valid vertices
data_clean = data[:, valid_vertices]
n_excluded = data.shape[1] - data_clean.shape[1]
if n_excluded > 0:
    print(f"Excluded {n_excluded} vertices with all-NaN values")
print(f"Fitting {valid_vertices.sum()} valid vertices")


#----------GAUSSIAN GRID -----------------
rows, columns = dm.shape[0], dm.shape[1]                                               # design size (50 x 50)
design = np.moveaxis(dm, -1, 0)                                                        # time first, then rows, then columns

screen_deg = 2 * np.degrees(np.arctan(screen_size_cm[1] / (2 * screen_distance_cm)))   # screen height in dva
pixel_size = screen_deg / rows                                                         # each design pixel covers screen_deg / rows dva

x = (np.arange(columns) - (columns - 1) / 2) * pixel_size                              # pixel centres, left -> right
y = -(np.arange(rows) - (rows - 1) / 2) * pixel_size                                   # pixel centres, top -> bottom (row 0 = +y, as in prfpy)
xv, yv = np.meshgrid(x, y)
# y comes first because design axis 1 is the row (height) axis
grid_coordinates = np.stack((yv, xv), axis=-1)                                         # shape (rows, columns, 2)

stimulus = PRFStimulus(design = design, grid = grid_coordinates)

# create 2D Gaussian prf model 
impulse_model = DerivativeTwoGammaImpulse(resolution=tr)
# Only weight_deriv is set; the two-gamma parameters use the model's default Glover HRF values
impulse_default_params = pd.DataFrame({
    "weight_deriv": [-0.5],
})

# Define pRF model with custom impulse response submodel
prf_model = Gaussian2DPRFModel(
    impulse_model=impulse_model,
)
# Same Gaussian + HRF, but centre given as (ecc, polar): used only for the grid search
polar_model = CanonicalPRFModel(prf_model=PolarGaussian2DPRFTuning(), impulse_model=impulse_model)

# Define model parameter grid range (as in prf_cssfit.py)
sizes  = max_ecc_size * np.linspace(0.1, 1, gauss_grid_nr)**2
eccs   = max_ecc_size * np.linspace(0.1, 1, gauss_grid_nr)**2
polars = np.linspace(0, 2 * np.pi, gauss_grid_nr)

grid_fitter = GridFitter(model=polar_model, stimulus=stimulus, compile_step=True)
_, grid_params = grid_fitter.fit(data=data_clean.T, batch_size=20, parameter_values={
    "ecc": eccs, "polar": polars, "sigma": sizes,
    "weight_deriv": [-0.5], "baseline": [0.0], "amplitude": [1.0]})

grid_params = polar_to_cartesian(grid_params)                                         # (ecc, polar) -> (mu_x, mu_y) for least squares / SGD


# -----------LEAST SQUARES FIT ------------
# Amplitude and baseline for the grid result 
print("Running least squares fit....")
ls_fitter = LeastSquaresFitter(model=prf_model, stimulus=stimulus)
_, gaussian_params = ls_fitter.fit(data=data_clean.T, parameters=grid_params,        #override the gaussian parameters
                                   slope_name="amplitude", intercept_name="baseline", batch_size=10000)


# ---------- ITERATIVE GAUSSIAN FIT --------------
gauss_sgd = SGDFitter(model=prf_model, stimulus=stimulus, compile_step=True,
                      adapter=Adapter([ParameterTransform(["sigma"], ops.log, ops.exp)]),
                      optimizer=keras.optimizers.Adam(learning_rate=0.01))
gauss_history, gaussian_params = gauss_sgd.fit(data=data_clean.T, init_parameters=gaussian_params,
                                               fixed_parameters=["weight_deriv"], num_steps=1000)


# -----------CSS FIT WITH SGD ----------------
css_adapter = Adapter([
    ParameterTransform(["sigma", "n"], ops.log, ops.exp)                             #adapter that log-transforms sigma and n during parameter optimization, forcing them to stay positive.
])

# 1. CSS model (compression here)
css_model = Gaussian2DPRFModel(
    encoding_model=CompressiveEncoder(encoding_model=PRFStimulusEncoder()),
    impulse_model=impulse_model,
)

# --- CSS grid over n (prfpy: exponent grid, size scaled by sqrt(n)) ---
corr_loss = CorrelationLoss(reduction="none")                                       # returns -r per vertex

n_grid = np.linspace(css_exponent_bounds[0], css_exponent_bounds[1], css_grid_nr)
best_r, best_n = np.full(len(data_clean.T), -np.inf), np.full(len(data_clean.T), n_grid[-1])
for n_val in n_grid:
    params_n = gaussian_params.assign(gain=1.0, n=n_val, sigma=gaussian_params["sigma"] * np.sqrt(n_val))
    pred_n = batched(css_model)(stimulus, params_n, batch_size=20000)
    r = -ops.convert_to_numpy(corr_loss(data_clean.T, pred_n))  # Pearson r per vertex
    is_better = r > best_r
    best_r[is_better], best_n[is_better] = r[is_better], n_val

init_params = gaussian_params.assign(gain=1.0, n=best_n, sigma=gaussian_params["sigma"] * np.sqrt(best_n))
_, init_params = LeastSquaresFitter(model=css_model, stimulus=stimulus).fit(
    data=data_clean.T, parameters=init_params, slope_name="amplitude", intercept_name="baseline", batch_size=10000)

# --- Iterative CSS fit ---
n_lo, n_hi = css_exponent_bounds                                                   # n_th = [0.01, 1]
css_adapter = Adapter([ParameterTransform(["sigma"], ops.log, ops.exp),
                       bounded_transform(["n"], n_lo, n_hi)])
init_params = clip_to_bounds(init_params, {"n": (n_lo, n_hi)})                     # the n grid ends exactly at 1.0

css_sgd = SGDFitter(model=css_model, stimulus=stimulus, adapter=css_adapter, compile_step=True,
                    optimizer=keras.optimizers.Adam(learning_rate=0.01))
css_history, css_params = css_sgd.fit(data=data_clean.T, init_parameters=init_params,
                                      fixed_parameters=["weight_deriv", "gain"], num_steps=1000)

print("All fitted!!")

# ------------ PREDICTION FROM CSS AND SAVING ------------------
# 3. Prediction = same CSS model + SGD parameters -> compressed response
css_pred = batched(css_model)(stimulus, css_params, batch_size=200)                # the CSS model, NOT the Gaussian prf_model # np.ndarray, (n_valid_vertices, n_frames)


# Rearrange CSS results into the prf_cssfit.py layout (columns selected by name)
p = css_params
css_fit = np.column_stack([
    p["mu_x"], p["mu_y"], p["sigma"],
    p["amplitude"] * p["gain"],                                                    # effective scale; equals amplitude if gain is fixed at 1
    p["baseline"], p["n"],
    p["weight_deriv"],                                                             # hrf_1 slot: NOT on prfpy's scale (see below)
    np.zeros(len(p)),                                                              # hrf_2 slot: fixed at 0, as in your prfpy bounds
])

css_fit_mat = np.full((data.shape[1], css_params_num), np.nan, dtype=float)        # (n_vertices, 9)
css_pred_mat = np.full_like(data, np.nan, dtype=float)                             # (188, n_vertices)

css_fit_mat[valid_vertices, :8] = css_fit
css_pred_mat[:, valid_vertices] = css_pred.T
css_fit_mat[:, 8] = r2_score_surf(bold_signal=data, model_prediction=css_pred_mat)

if 'loo-avg' in input_fn:
    loo_bold_fn = input_fn.replace('loo-avg-', 'loo-')
    loo_img, loo_bold = load_surface(fn=loo_bold_fn)
    loo_r2 = r2_score_surf(bold_signal=loo_bold, model_prediction=css_pred_mat)
    css_fit_mat = np.column_stack((css_fit_mat, loo_r2))
    maps_names = ['mu_x', 'mu_y', 'prf_size', 'prf_amplitude', 'bold_baseline',
                  'n', 'hrf_1', 'hrf_2', 'r_squared', 'loo_r_squared']
else:
    maps_names = ['mu_x', 'mu_y', 'prf_size', 'prf_amplitude', 'bold_baseline',
                  'n', 'hrf_1', 'hrf_2', 'r_squared']

assert css_fit_mat.shape == (data.shape[1], len(maps_names))
assert css_pred_mat.shape == data.shape

img_css_fit_mat = make_surface_image(data=css_fit_mat.T, source_img=img, maps_names=maps_names)
nb.save(img_css_fit_mat, '{}/{}'.format(prf_fit_dir, css_fit_fn))
print(f"Saved fit to: {prf_fit_dir}/{css_fit_fn}")
img_css_pred_mat = make_surface_image(data=css_pred_mat, source_img=img)           # no .T
nb.save(img_css_pred_mat, '{}/{}'.format(prf_fit_dir, css_pred_fn))
print(f"Saved prediction to: {prf_fit_dir}/{css_pred_fn}")

# Print duration
end_time = datetime.datetime.now()
print("\nStart time:\t{start_time}\nEnd time:\t{end_time}\nDuration:\t{dur}".format(start_time=start_time, 
                                                                                    end_time=end_time, 
                                                                                    dur=end_time - start_time))