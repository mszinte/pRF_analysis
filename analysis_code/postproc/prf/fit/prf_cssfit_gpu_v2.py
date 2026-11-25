"""
-----------------------------------------------------------------------------------------
prf_cssfit_gpu_v2.py
-----------------------------------------------------------------------------------------
Goal of the script:
pRF fit (CSS model) on GPU with prfmodel >= 0.2.0, following the procedure of the
prfmodel workshop (Coimbra, 2026-11-05) and of prf_cssfit_gpu.py:
grid (polar Gaussian) -> least squares (amplitude/baseline) -> Gaussian SGD
-> CSS least squares -> CSS SGD (n bounded to [n_th])

Changes compared with prf_cssfit_gpu.py (see "[V2]" comments):
1. [V2-SCREEN] screen geometry unchanged (design assumed to span the screen height, or the
   fake_screen_size_cm of the analysis yml, as in RetinoMaps); new check that warns when the
   stimulated rows of the design match the aspect ratio of a non-square screen (design
   spanning the WIDTH, e.g. amblyo7T: set fake_screen_size_cm to [width, width]).
2. [V2-GRID] Gaussian grid: same polar grid as prfpy (quadratic eccentricity spacing: more
   points near the fovea). Default: same ranges as before (max_ecc_size). If grid_mu_extent is
   set in the analysis yml: ranges of prfmodel.fitters.grid_values_2d_prf (eccentricity up to
   grid_mu_extent x the half-extent of the stimulus, sigma log-spaced from the stimulus
   resolution to its extent). Polar angles without the duplicate 2*pi in both cases.
3. [V2-N] CSS start at n = 0.5 for every vertex (as prfmodel.init_css_from_gaussian), sigma
   scaled by sqrt(n) so that the starting effective size sigma/sqrt(n) equals the Gaussian
   size. The previous grid over n included the upper bound n = 1: in a simulation with the
   real design matrix, 96 % of the vertices started at n = 1 stayed exactly at 1 after SGD
   (sigmoid-bounded parameter started at the edge), including vertices with a true n of 0.5.
4. [V2-MODEL] CSS model = prfmodel.models.prf.Gaussian2DCSSPRFModel.
5. [V2-GAUSS] the Gaussian fit is also saved (same layout, n = 1).
6. Removed: unused imports, a first css_adapter that was overwritten, unused n_batches.
-----------------------------------------------------------------------------------------
Input(s): same as prf_cssfit_gpu.py
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject name
sys.argv[4]: input file name (path to the data to fit)
sys.argv[5]: analysis task name (ex. prf)
sys.argv[6]: number of jobs (kept for compatibility with the submission script)
-----------------------------------------------------------------------------------------
Output(s): same files as prf_cssfit_gpu.py (CSS fit and prediction) + Gaussian fit
-----------------------------------------------------------------------------------------
Written by Sina Kling (sina.kling@outlook.de), adapted from prf_css_fit.py by Uriel Lascombes
V2 proposal: adapted to prfmodel 0.2.0 (not yet run on the cluster: test on one subject first)
-----------------------------------------------------------------------------------------
"""
# Stop warnings
import warnings
warnings.filterwarnings("ignore")

import os
os.environ["KERAS_BACKEND"] = "torch"                                               # must be set before keras/prfmodel are imported
import sys
import datetime
from importlib.util import find_spec

PRFMODEL_URL = "https://github.com/popylar-org/prfmodel"
missing_packages = [package for package in ['prfmodel', 'keras', 'torch'] if find_spec(package) is None]
if missing_packages:
    raise ImportError("Missing python package(s) for the GPU pRF fit: {}\nInstall prfmodel (>= 0.2.0), see {}\n"
                      "(current python: {})".format(', '.join(missing_packages), PRFMODEL_URL, sys.executable))

import numpy as np
import pandas as pd
import nibabel as nb
import keras
import torch
from keras import ops
torch._dynamo.config.automatic_dynamic_shapes = False                               # see prf_cssfit_gpu.py (torch 2.15 recompilation crash)
from prfmodel.stimuli import PRFStimulus
from prfmodel.impulse import DerivativeTwoGammaImpulse
from prfmodel.models.prf import Gaussian2DPRFModel, Gaussian2DCSSPRFModel, init_css_from_gaussian   # [V2-MODEL]
from prfmodel.models.prf.canonical import CanonicalPRFModel
from prfmodel.fitters import GridFitter, LeastSquaresFitter, SGDFitter, grid_values_2d_prf
from prfmodel.fitters.adapter import Adapter, ParameterTransform
from prfmodel.utils import batched

# Personal imports
sys.path.append("{}/../../../utils".format(os.getcwd()))
from maths_utils import r2_score_surf
from settings_utils import load_settings
from screen_utils import get_screen_settings
from pycortex_utils import set_pycortex_config_file
from surface_utils import load_surface, make_surface_image
from prfmodel_utils import PolarGaussian2DPRFTuning, polar_to_cartesian, bounded_transform, clip_to_bounds, host_safe

start_time = datetime.datetime.now()

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
sub_num = subject[4:]
input_fn = sys.argv[4]
analysis_name = sys.argv[5]
n_jobs = int(sys.argv[6])
css_params_num = 9

# Settings
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../../"))
analysis_info = load_settings([os.path.join(base_dir, project_dir, "settings.yml"),
                               os.path.join(base_dir, project_dir, f"{analysis_name}-analysis.yml")])[0]
tr = analysis_info['TR']
gauss_grid_nr = analysis_info['gauss_grid_nr']
max_ecc_size = analysis_info['max_ecc_size']
grid_mu_extent = analysis_info.get('grid_mu_extent', None)                         # [V2-GRID] None = previous grid (max_ecc_size)
n_lo, n_hi = analysis_info['n_th']                                                 # CSS exponent bounds, e.g. [0.01, 1]
n_init = analysis_info.get('css_n_init', 0.5)                                      # [V2-N] starting n for every vertex
sgd_steps = analysis_info.get('sgd_num_steps', 1000)
sgd_lr = analysis_info.get('sgd_learning_rate', 0.01)

# Screen settings from the task-<task>_events.json files (optional override in the analysis yml)
task_name = input_fn.split("task-")[1].split("_")[0]
screen_size_cm, screen_distance_cm = get_screen_settings(main_dir, project_dir, sub_num, task_name)
real_screen_size_cm = list(screen_size_cm)                                         # kept for the geometry check
if 'fake_screen_size_cm' in analysis_info:
    print("[INFO] Fake screen size found in settings — overriding real screen size.")
    screen_size_cm = analysis_info['fake_screen_size_cm']

set_pycortex_config_file("{}/{}/derivatives/pp_data/cortex".format(main_dir, project_dir))

# Design matrix (subject-specific first, then general)
dm_name = analysis_info['dm_name']
dm_name_tofind = dm_name[1:] if dm_name else "vdm"
dm_base_dir = '{}/{}/derivatives/vdm'.format(main_dir, project_dir)
dm_fn_subject = '{}/sub-{}/sub-{}_task-{}_{}.npy'.format(dm_base_dir, sub_num, sub_num, task_name, dm_name_tofind)
dm_fn_general = '{}/task-{}_{}.npy'.format(dm_base_dir, task_name, dm_name_tofind)
dm_fn = dm_fn_subject if os.path.isfile(dm_fn_subject) else dm_fn_general
if not os.path.isfile(dm_fn):
    raise FileNotFoundError(f"No DM found for task '{task_name}'.\n  Checked: {dm_fn_subject}\n  Checked: {dm_fn_general}")
print(f"Loading DM from: {dm_fn}")
dm = np.load(dm_fn)                                                                # (rows, columns, time)

# Output files (same names as prf_cssfit_gpu.py + Gaussian fit)
output_folder = analysis_info["output_folder"]
fmt = '170k' if input_fn.endswith('.nii') else 'fsnative'
prf_fit_dir = "{}/{}/derivatives/pp_data/{}/{}/{}/fit".format(main_dir, project_dir, subject, fmt, output_folder)
os.makedirs(prf_fit_dir, exist_ok=True)
base_fn = input_fn.split('/')[-1]
css_fit_fn = base_fn.replace('bold', f'{analysis_name}-css{dm_name}_fit')
css_pred_fn = base_fn.replace('bold', f'{analysis_name}-css{dm_name}_pred')
gauss_fit_fn = base_fn.replace('bold', f'{analysis_name}-gauss{dm_name}_fit')     # [V2-GAUSS]

# ------------ LOAD DATA ------------
img, data = load_surface(fn=input_fn)                                              # (n times, n vertices)
valid_vertices = ~np.isnan(data).any(axis=0)
data_clean = data[:, valid_vertices]
print(f"Fitting {valid_vertices.sum()} valid vertices ({data.shape[1] - valid_vertices.sum()} excluded: NaN)")

# ------------ STIMULUS ------------
rows, columns = dm.shape[0], dm.shape[1]
design = np.moveaxis(dm, -1, 0)                                                    # (time, rows, columns)

# Pixel size as before: the square design is assumed to span the screen HEIGHT (screen_size_cm[1]), or the
# fake_screen_size_cm of the analysis yml (used by RetinoMaps, and by amblyo7T: [67.9, 67.9], see below).
screen_deg = 2 * np.degrees(np.arctan(screen_size_cm[1] / (2 * screen_distance_cm)))
pixel_size = screen_deg / rows
print(f"Design spans {screen_deg:.2f} deg ({'fake' if 'fake_screen_size_cm' in analysis_info else 'screen height'}); "
      f"pixel = {pixel_size:.3f} deg")

# [V2-SCREEN] Check. If the stimulus video was padded to a square (black bands above and below), the design spans the
# screen WIDTH and only (height / width) of its rows are ever stimulated (amblyo7T: 1024 x 768 screen, 76 % of the rows).
# The height convention then underestimates every position and size by height / width.
on = dm.max(axis=-1) > 0
frac_rows, frac_cols = on.any(axis=1).mean(), on.any(axis=0).mean()
real_aspect = real_screen_size_cm[1] / real_screen_size_cm[0]
print(f"Stimulated part of the design: {100 * frac_rows:.0f} % of the rows, {100 * frac_cols:.0f} % of the columns "
      f"(screen height / width = {real_aspect:.2f})")
if ('fake_screen_size_cm' not in analysis_info and real_aspect < 0.95 and frac_rows < 0.9 * frac_cols
        and abs(frac_rows / frac_cols - real_aspect) < 0.1):
    print("[WARNING] The stimulated rows match the screen aspect ratio: the design probably spans the screen WIDTH. "
          f"Set fake_screen_size_cm: [{real_screen_size_cm[0]}, {real_screen_size_cm[0]}] in the analysis yml.")

x = (np.arange(columns) - (columns - 1) / 2) * pixel_size                          # pixel centres, left -> right
y = -(np.arange(rows) - (rows - 1) / 2) * pixel_size                               # pixel centres, top -> bottom (row 0 = +y)
xv, yv = np.meshgrid(x, y)
grid_coordinates = np.stack((yv, xv), axis=-1)                                     # prfmodel: grid[..., 0] = y, grid[..., 1] = x
stimulus = PRFStimulus(design=design, grid=grid_coordinates, frame_duration=tr)

impulse_model = DerivativeTwoGammaImpulse(resolution=tr)                           # default Glover two-gamma + derivative
prf_model = Gaussian2DPRFModel(impulse_model=impulse_model)
polar_model = CanonicalPRFModel(prf_model=PolarGaussian2DPRFTuning(), impulse_model=impulse_model)   # grid only

# ------------ GAUSSIAN GRID ------------
# [V2-GRID] polar grid with quadratic eccentricity spacing (prfpy: more grid points near the fovea)
polars = np.linspace(0, 2 * np.pi, gauss_grid_nr, endpoint=False)                  # 0 and 2*pi are the same angle
if grid_mu_extent is None:
    # previous ranges (prf_cssfit_gpu.py / prf_cssfit.py)
    eccs = max_ecc_size * np.linspace(0.1, 1, gauss_grid_nr) ** 2
    sizes = max_ecc_size * np.linspace(0.1, 1, gauss_grid_nr) ** 2
else:
    # ranges as prfmodel.fitters.grid_values_2d_prf: eccentricity up to grid_mu_extent x the half-extent of the
    # stimulus grid, sigma log-spaced from the grid spacing to the stimulus extent (grid_values_2d_prf itself)
    eccs = grid_mu_extent * 0.5 * max(np.ptp(x), np.ptp(y)) * np.linspace(0.1, 1, gauss_grid_nr) ** 2
    sizes = grid_values_2d_prf(stimulus, num_sigma=gauss_grid_nr)["sigma"]
print(f"Grid: eccentricity {eccs.min():.2f}-{eccs.max():.2f} deg (quadratic), sigma {sizes.min():.2f}-{sizes.max():.2f} deg (log), "
      f"{gauss_grid_nr} polar angles")
grid_fitter = GridFitter(model=polar_model, stimulus=stimulus, compile_step=True)
_, grid_params = grid_fitter.fit(data=data_clean.T, batch_size=216, parameter_values={
    "ecc": eccs, "polar": polars, "sigma": sizes, "weight_deriv": [-0.5], "baseline": [0.0], "amplitude": [1.0]})
grid_params = polar_to_cartesian(grid_params)

# ------------ GAUSSIAN LEAST SQUARES + SGD ------------
ls_fitter = LeastSquaresFitter(model=prf_model, stimulus=stimulus)
_, gaussian_params = ls_fitter.fit(data=data_clean.T, parameters=grid_params, slope_name="amplitude",
                                   intercept_name="baseline", batch_size=None)
gauss_sgd = SGDFitter(model=prf_model, stimulus=stimulus, compile_step=True,
                      adapter=Adapter([ParameterTransform(["sigma"], host_safe(ops.log), host_safe(ops.exp))]),
                      optimizer=keras.optimizers.Adam(learning_rate=sgd_lr))
gauss_history, gaussian_params = gauss_sgd.fit(data=data_clean.T, init_parameters=gaussian_params,
                                               fixed_parameters=["weight_deriv"], num_steps=sgd_steps)

# ------------ CSS: START, LEAST SQUARES, SGD ------------
css_model = Gaussian2DCSSPRFModel(impulse_model=impulse_model)                     # [V2-MODEL]
# [V2-N] start at n = n_init for every vertex; sigma * sqrt(n) keeps the Gaussian effective size sigma / sqrt(n)
css_init = init_css_from_gaussian(gaussian_params, n=n_init)                       # as the prfmodel examples (adds n)
css_init["sigma"] = css_init["sigma"] * np.sqrt(n_init)                            # prfpy convention: keep sigma / sqrt(n)
css_init = clip_to_bounds(css_init, {"n": (n_lo, n_hi)})
_, css_init = LeastSquaresFitter(model=css_model, stimulus=stimulus).fit(
    data=data_clean.T, parameters=css_init, slope_name="amplitude", intercept_name="baseline", batch_size=None)

# n bounded to the open interval (n_lo, n_hi) with a scaled sigmoid (prfmodel ParameterConstraint takes only one
# bound at a time: with upper=1 alone, n could become <= 0, which the CSS model does not allow)
css_adapter = Adapter([ParameterTransform(["sigma"], host_safe(ops.log), host_safe(ops.exp)),
                       bounded_transform(["n"], n_lo, n_hi)])
css_sgd = SGDFitter(model=css_model, stimulus=stimulus, adapter=css_adapter, compile_step=True,
                    optimizer=keras.optimizers.Adam(learning_rate=sgd_lr))
css_history, css_params = css_sgd.fit(data=data_clean.T, init_parameters=css_init,
                                      fixed_parameters=["weight_deriv"], num_steps=sgd_steps)
print("All fitted!!")

# Diagnostics: convergence and parameters at the bounds
print("Gaussian SGD loss:", [round(float(v), 5) for v in gauss_history.history["loss"][::100]])
print("CSS SGD loss:     ", [round(float(v), 5) for v in css_history.history["loss"][::100]])
eps = 0.01 * (n_hi - n_lo)
print(f"n: median {css_params['n'].median():.2f}; {100 * (css_params['n'] >= n_hi - eps).mean():.1f} % at the upper bound, "
      f"{100 * (css_params['n'] <= n_lo + eps).mean():.1f} % at the lower bound")
print(f"Gaussian sigma at the largest grid size: {100 * (gaussian_params['sigma'] >= sizes.max()).mean():.1f} %")

# ------------ PREDICTION AND SAVING ------------
css_pred = batched(css_model)(stimulus, css_params, batch_size=len(css_params))
gauss_pred = batched(prf_model)(stimulus, gaussian_params, batch_size=len(gaussian_params))


def to_layout(p, n):
    """prf_cssfit.py layout: mu_x, mu_y, sigma, amplitude, baseline, n, hrf_1 (= weight_deriv, not prfpy's scale), hrf_2 (0)."""
    return np.column_stack([p["mu_x"], p["mu_y"], p["sigma"], p["amplitude"], p["baseline"], n,
                            p["weight_deriv"], np.zeros(len(p))])


maps_names = ['mu_x', 'mu_y', 'prf_size', 'prf_amplitude', 'bold_baseline', 'n', 'hrf_1', 'hrf_2', 'r_squared']
for fit_fn, pred_fn, params, pred, n_col in [(css_fit_fn, css_pred_fn, css_params, css_pred, css_params["n"]),
                                             (gauss_fit_fn, None, gaussian_params, gauss_pred, np.ones(len(gaussian_params)))]:
    fit_mat = np.full((data.shape[1], css_params_num), np.nan, dtype=float)
    pred_mat = np.full_like(data, np.nan, dtype=float)
    fit_mat[valid_vertices, :8] = to_layout(params, n_col)
    pred_mat[:, valid_vertices] = ops.convert_to_numpy(pred).T
    fit_mat[:, 8] = r2_score_surf(bold_signal=data, model_prediction=pred_mat)
    names = list(maps_names)
    if 'loo-avg' in input_fn:
        _, loo_bold = load_surface(fn=input_fn.replace('loo-avg-', 'loo-'))
        fit_mat = np.column_stack((fit_mat, r2_score_surf(bold_signal=loo_bold, model_prediction=pred_mat)))
        names.append('loo_r_squared')
    nb.save(make_surface_image(data=fit_mat.T, source_img=img, maps_names=names), '{}/{}'.format(prf_fit_dir, fit_fn))
    print(f"Saved fit to: {prf_fit_dir}/{fit_fn}")
    if pred_fn is not None:
        nb.save(make_surface_image(data=pred_mat, source_img=img), '{}/{}'.format(prf_fit_dir, pred_fn))
        print(f"Saved prediction to: {prf_fit_dir}/{pred_fn}")

if torch.cuda.is_available():
    print(f"Peak GPU memory: {torch.cuda.max_memory_allocated() / 1e9:.1f} GB "
          f"/ {torch.cuda.get_device_properties(0).total_memory / 1e9:.1f} GB")
end_time = datetime.datetime.now()
print("\nStart time:\t{}\nEnd time:\t{}\nDuration:\t{}".format(start_time, end_time, end_time - start_time))
