"""
-----------------------------------------------------------------------------------------
ncsf_fit_gpu.py
-----------------------------------------------------------------------------------------
Goal of the script:
nCSF fit on GPU with prfmodel >= 0.2.0 (replaces ncsf_fit.py, based on prfpy_csenf),
following practical 7 of the prfmodel workshop (Coimbra, 2026-11-05) and the inputs /
outputs of analysis_code/postproc/prf/fit/prf_cssfit_gpu_v2.py:
grid (cs_peak, sf_peak, [width_l], width_r, slope_crf) -> least squares (amplitude/baseline)
-> SGD (Adam) with every CSF parameter kept inside (lower, upper) by bounded_transform
(prfmodel_utils.py; same as box_transform of the prfmodel nCSF practical: prfmodel 0.2.0's
ParameterConstraint takes only one bound per parameter)

Changes compared with ncsf_fit.py (prfpy_csenf):
1. [PRFMODEL] model = prfmodel.models.csf.CSFModel (aLP CSF + Naka-Rushton CRF, Roelofzen et al. 2025).
   Parameter names: CSp -> cs_peak, SFp -> sf_peak, crf_exp -> slope_crf, amp_1 -> amplitude,
   bold_baseline -> baseline. prfmodel fixes width_l (0.68): set ncsf_fit_width_l: True in
   nCSF-analysis.yml to fit it as before (CSFTuningFreeWidthL, nCSF/utils/ncsf_utils.py).
2. [HRF] HRF = prfmodel DerivativeTwoGammaImpulse with weight_deriv fixed (ncsf_weight_deriv),
   instead of prfpy hrf=[1, 1, 0]. The hrf_1 map now stores weight_deriv, hrf_2 = 0.
3. [IO] same inputs and folders as prf_cssfit_gpu_v2.py: analysis name argument (loads
   {analysis_name}-analysis.yml), task name read from the input file name (the yml key
   'nCSF_task_name' no longer exists), output folder = output_folder of the yml,
   non-fitted vertices = NaN (were 0), r2 computed with r2_score_surf.
4. [DERIV] auc, normalize_auc and SFmax computed with nCSF/utils/ncsf_utils.py (no prfpy_csenf).
   Map layout unchanged (maps_names_ncsf), so the postfit scripts read the same columns.
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject name
sys.argv[4]: input file name (path to the data to fit)
sys.argv[5]: analysis name (ex. nCSF, loads nCSF-analysis.yml)
sys.argv[6]: number of jobs (kept for compatibility with the submission script)
-----------------------------------------------------------------------------------------
Output(s):
{output_folder}/fit/*_ncsf_fit.{func.gii,dtseries.nii}: maps_names_ncsf (+ ncsf_loo_rsq)
{output_folder}/fit/*_ncsf_pred.{func.gii,dtseries.nii}: predicted time series
-----------------------------------------------------------------------------------------
To run:
0. use a python environment where prfmodel is installed (https://github.com/popylar-org/prfmodel)
1. cd to function
>> cd ~/projects/pRF_analysis/nCSF/postproc/nCSF/fit
2. run python command
python ncsf_fit_gpu.py [main directory] [project name] [subject name]
                       [input file name] [analysis name] [number of jobs]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/projects/pRF_analysis/nCSF/postproc/nCSF/fit
python ncsf_fit_gpu.py /scratch/mszinte/data nCSF sub-01 /scratch/mszinte/data/nCSF/derivatives/pp_data/sub-01/fsnative/func/fmriprep_dct_z-score_avg/sub-01_task-nCSF_hemi-R_fmriprep_dct_z-score_avg_bold.func.gii nCSF 4
-----------------------------------------------------------------------------------------
Written by Martin Szinte (martin.szinte@gmail.com)
and Uriel Lascombes (uriel.lascombes@laposte.net)
prfmodel version adapted from prf_cssfit_gpu_v2.py (Sina Kling)
-----------------------------------------------------------------------------------------
"""
# Stop warnings
import warnings
warnings.filterwarnings("ignore")

import os
os.environ.setdefault("KERAS_BACKEND", "torch")                                    # must be set before keras/prfmodel are imported
import re
import sys
import glob
import datetime
from importlib.util import find_spec

PRFMODEL_URL = "https://github.com/popylar-org/prfmodel"
missing_packages = [package for package in ['prfmodel', 'keras', 'torch'] if find_spec(package) is None]
if missing_packages:
    raise ImportError("Missing python package(s) for the GPU nCSF fit: {}\nInstall prfmodel (>= 0.2.0), see {}\n"
                      "(current python: {})".format(', '.join(missing_packages), PRFMODEL_URL, sys.executable))

import numpy as np
import nibabel as nb
import keras
import torch
from keras import ops
torch._dynamo.config.automatic_dynamic_shapes = False                               # see prf_cssfit_gpu.py (torch recompilation crash)
from prfmodel.stimuli import CSFStimulus
from prfmodel.impulse import DerivativeTwoGammaImpulse
from prfmodel.models.csf import CSFModel
from prfmodel.models.csf.canonical import CanonicalCSFModel
from prfmodel.fitters import GridFitter, LeastSquaresFitter, SGDFitter
from prfmodel.fitters.adapter import Adapter
from prfmodel.utils import batched

# Personal imports
script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.abspath(os.path.join(script_dir, "../../../../analysis_code/utils")))
sys.path.append(os.path.abspath(os.path.join(script_dir, "../../../utils")))
from maths_utils import r2_score_surf
from settings_utils import load_settings
from surface_utils import load_surface, make_surface_image
from prfmodel_utils import bounded_transform, clip_to_bounds
from ncsf_utils import (ncsf_stimulus_levels, ncsf_sequence_from_events, make_csf_tuning_free_width_l,
                        ncsf_aulcsf, ncsf_sfmax)

start_time = datetime.datetime.now()

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
input_fn = sys.argv[4]
analysis_name = sys.argv[5]
n_jobs = int(sys.argv[6])

# Settings: prf-analysis.yml (if present) only provides the shared ROI definitions,
# {analysis_name}-analysis.yml is loaded last and overrides it
base_dir = os.path.abspath(os.path.join(script_dir, "../../../../"))
settings_fns = [os.path.join(base_dir, project_dir, "settings.yml"),
                os.path.join(base_dir, project_dir, "prf-analysis.yml"),
                os.path.join(base_dir, project_dir, f"{analysis_name}-analysis.yml")]
analysis_info = load_settings([fn for fn in settings_fns if os.path.isfile(fn)])[0]

tr = analysis_info['TR']
output_folder = analysis_info['output_folder']
nCSF_session = analysis_info['nCSF_session']
maps_names_ncsf = list(analysis_info['maps_names_ncsf'])
run_ncsf_fit_on_rois = analysis_info['run_ncsf_fit_on_rois']
fit_width_l = analysis_info.get('ncsf_fit_width_l', False)
width_l_fixed = analysis_info.get('ncsf_width_l', 0.68)
weight_deriv = analysis_info.get('ncsf_weight_deriv', -0.5)
grid_settings = analysis_info['ncsf_grid']                                         # {param: [min, max, n]}, geomspace
bounds = {k: tuple(v) for k, v in analysis_info['ncsf_bounds'].items()}            # {param: [lower, upper]} for SGD
sgd_steps = analysis_info.get('sgd_num_steps', 1000)
sgd_lr = analysis_info.get('sgd_learning_rate', 0.01)
grid_batch_size = analysis_info.get('ncsf_grid_batch_size', 100)
max_sfmax = analysis_info.get('ncsf_max_sfmax', 50)
sf_levels, con_levels = ncsf_stimulus_levels(analysis_info['sf_minFreq'], analysis_info['sf_maxFreq'],
                                             analysis_info['sf_filtNum'], analysis_info['minCont'],
                                             analysis_info['maxCont'], analysis_info['contNum'])
csf_params = ["cs_peak", "sf_peak"] + (["width_l"] if fit_width_l else []) + ["width_r", "slope_crf"]
if not fit_width_l:
    grid_settings = {k: v for k, v in grid_settings.items() if k != 'width_l'}
    bounds = {k: v for k, v in bounds.items() if k != 'width_l'}

# Output files (same folders as prf_cssfit_gpu_v2.py; names kept for the postfit scripts)
fmt = '170k' if input_fn.endswith('.nii') else 'fsnative'
ncsf_fit_dir = "{}/{}/derivatives/pp_data/{}/{}/{}/fit".format(main_dir, project_dir, subject, fmt, output_folder)
os.makedirs(ncsf_fit_dir, exist_ok=True)
base_fn = input_fn.split('/')[-1]
ncsf_fit_fn = base_fn.replace('bold', 'ncsf_fit')
ncsf_pred_fn = base_fn.replace('bold', 'ncsf_pred')

# ------------ STIMULUS ------------
# Task name from the input file name (as prf_cssfit_gpu_v2.py), events of the first run of the session
task_name = input_fn.split("task-")[1].split("_")[0]
events_dir = '{}/{}/{}/{}/func'.format(main_dir, project_dir, subject, nCSF_session)
events_fns = sorted(glob.glob('{}/{}_{}_task-{}_*events.tsv'.format(events_dir, subject, nCSF_session, task_name)))
if not events_fns:
    raise FileNotFoundError(f"No events file for task '{task_name}' in {events_dir}")
sf_seq, con_seq = ncsf_sequence_from_events(events_fns[0], sf_levels, con_levels)
for other_fn in events_fns[1:]:                                                    # the fit assumes the same sequence in every run
    sf_other, con_other = ncsf_sequence_from_events(other_fn, sf_levels, con_levels)
    if sf_other.shape != sf_seq.shape or not (np.array_equal(sf_other, sf_seq) and np.array_equal(con_other, con_seq)):
        print(f"[WARNING] stimulus sequence of {os.path.basename(other_fn)} differs from {os.path.basename(events_fns[0])}")
print(f"Stimulus from: {events_fns[0]} ({len(sf_seq)} TRs, {np.mean(sf_seq == 0) * 100:.0f} % blank)")

# ------------ LOAD DATA ------------
img, raw_data = load_surface(fn=input_fn)                                          # (n times, n vertices)
if raw_data.shape[0] != len(sf_seq):
    raise ValueError(f"Data have {raw_data.shape[0]} TRs but the events file {len(sf_seq)} rows")
valid_vertices = ~np.isnan(raw_data).any(axis=0) & (np.nanstd(raw_data, axis=0) > 0)

if run_ncsf_fit_on_rois:
    from pycortex_utils import set_pycortex_config_file, get_rois
    set_pycortex_config_file("{}/{}/derivatives/pp_data/cortex".format(main_dir, project_dir))
    if fmt == '170k':
        pycortex_subject, hemi = analysis_info['pycortex_subject_template'], None
    else:
        pycortex_subject, hemi = subject, re.search(r'hemi-[LR]', input_fn).group()
    rois_method_format = analysis_info['rois_methods'][fmt][0]
    rois = analysis_info[rois_method_format]
    rois = list(rois.keys()) if isinstance(rois, dict) else rois
    print('Running fit only on: {}'.format(rois))
    roi_verts_dict = get_rois(pycortex_subject, surf_format=fmt, rois_type=rois_method_format,
                              mask=False, rois=rois, hemis=hemi)
    in_rois = np.zeros(raw_data.shape[1], dtype=bool)
    in_rois[np.unique(np.concatenate(list(roi_verts_dict.values()))).astype(int)] = True
    valid_vertices &= in_rois

data_clean = raw_data[:, valid_vertices].T                                         # prfmodel: (n vertices, n times)
print(f"Fitting {valid_vertices.sum()} vertices ({raw_data.shape[1] - valid_vertices.sum()} excluded: NaN, flat or outside ROIs)")

# ------------ MODEL ------------
stimulus = CSFStimulus(sf=sf_seq, contrast=con_seq, frame_duration=tr)             # contrast in % Michelson
impulse_model = DerivativeTwoGammaImpulse(resolution=tr)
if fit_width_l:
    ncsf_model = CanonicalCSFModel(csf_model=make_csf_tuning_free_width_l(), impulse_model=impulse_model)
else:
    ncsf_model = CSFModel(width_l=width_l_fixed, impulse_model=impulse_model)
print(f"Model parameters: {ncsf_model.parameter_names}")

# ------------ GRID ------------
grid_values = {name: np.geomspace(lo, hi, int(n)) for name, (lo, hi, n) in grid_settings.items()}
for name, (lower, upper) in bounds.items():                                        # start far from the bounds: the sigmoid is flat
    if not lower < np.min(grid_values[name]) <= np.max(grid_values[name]) < upper:  # there and a parameter starting on a bound stays
        raise ValueError(f"Grid of {name} ({np.min(grid_values[name])}-{np.max(grid_values[name])}) "
                         f"must lie strictly inside its bounds ({lower}, {upper})")
grid_values.update({"weight_deriv": [weight_deriv], "baseline": [0.0], "amplitude": [1.0]})
print("Grid: " + ", ".join(f"{k} {np.min(v):.2f}-{np.max(v):.2f} ({len(v)})" for k, v in grid_values.items()
                           if k in csf_params) + f" = {int(np.prod([len(grid_values[k]) for k in csf_params]))} predictions")
grid_fitter = GridFitter(model=ncsf_model, stimulus=stimulus, compile_step=True)
_, grid_params = grid_fitter.fit(data=data_clean, parameter_values=grid_values, batch_size=grid_batch_size)

# ------------ LEAST SQUARES (amplitude, baseline) ------------
ls_fitter = LeastSquaresFitter(model=ncsf_model, stimulus=stimulus)
_, ls_params = ls_fitter.fit(data=data_clean, parameters=grid_params, slope_name="amplitude",
                             intercept_name="baseline", batch_size=None)

# ------------ SGD ------------
# Adam + Adapter as in the prfmodel workshop. Every CSF parameter is kept strictly inside (lower, upper): Adam works on
# logit((x - lower) / (upper - x)), mapped back with lower + (upper - lower) * sigmoid before every prediction.
# Start values clipped strictly inside the bounds (a value on a bound has an infinite logit)
ls_params = clip_to_bounds(ls_params, bounds)
adapter = Adapter([bounded_transform([name], lo, hi) for name, (lo, hi) in bounds.items()])
sgd_fitter = SGDFitter(model=ncsf_model, stimulus=stimulus, adapter=adapter, compile_step=True,
                       optimizer=keras.optimizers.Adam(learning_rate=sgd_lr))
sgd_history, ncsf_params = sgd_fitter.fit(data=data_clean, init_parameters=ls_params,
                                          fixed_parameters=["weight_deriv"], num_steps=sgd_steps)
print("All fitted!!")

# Diagnostics: convergence, parameters at the bounds, negative amplitudes
print("SGD loss:", [round(float(v), 5) for v in sgd_history.history["loss"][::max(1, sgd_steps // 10)]])
for name, (lo, hi) in bounds.items():                                             # at a bound = not constrained by the data
    eps = 0.01 * (hi - lo)
    print(f"{name}: median {ncsf_params[name].median():.2f}; {100 * (ncsf_params[name] >= hi - eps).mean():.1f} % at "
          f"the upper bound, {100 * (ncsf_params[name] <= lo + eps).mean():.1f} % at the lower bound")
print(f"amplitude <= 0: {100 * (ncsf_params['amplitude'] <= 0).mean():.1f} % of the vertices "
      "(response decreasing with contrast, exclude them in postfit)")

# ------------ DERIVATIVES ------------
width_l = ncsf_params["width_l"].to_numpy() if fit_width_l else np.full(len(ncsf_params), width_l_fixed)
sf_min, sf_max = sf_levels[sf_levels > 0].min(), sf_levels.max()
auc = ncsf_aulcsf(ncsf_params["cs_peak"], ncsf_params["sf_peak"], width_l, ncsf_params["width_r"], sf_min, sf_max)
normalize_auc = ncsf_aulcsf(ncsf_params["cs_peak"], ncsf_params["sf_peak"], width_l, ncsf_params["width_r"],
                            sf_min, sf_max, normalize=True, cs_max=100 / analysis_info['minCont'])
sfmax = ncsf_sfmax(ncsf_params["cs_peak"], ncsf_params["sf_peak"], ncsf_params["width_r"], max_sfmax=max_sfmax)

# ------------ PREDICTION AND SAVING ------------
ncsf_pred = batched(ncsf_model)(stimulus, ncsf_params, batch_size=len(ncsf_params))
pred_mat = np.full_like(raw_data, np.nan, dtype=float)
pred_mat[:, valid_vertices] = ops.convert_to_numpy(ncsf_pred).T

# same layout as ncsf_fit.py (maps_names_ncsf): width_r, SFp, CSp, width_l, crf_exp, amp_1, bold_baseline,
# hrf_1 (= weight_deriv), hrf_2 (0), ncsf_rsq, auc, normalize_auc, SFmax
layout = {'width_r': ncsf_params["width_r"], 'SFp': ncsf_params["sf_peak"], 'CSp': ncsf_params["cs_peak"],
          'width_l': width_l, 'crf_exp': ncsf_params["slope_crf"], 'amp_1': ncsf_params["amplitude"],
          'bold_baseline': ncsf_params["baseline"], 'hrf_1': ncsf_params["weight_deriv"],
          'hrf_2': np.zeros(len(ncsf_params)), 'auc': auc, 'normalize_auc': normalize_auc, 'SFmax': sfmax}
fit_mat = np.full((raw_data.shape[1], len(maps_names_ncsf)), np.nan, dtype=float)
for idx, name in enumerate(maps_names_ncsf):
    if name == 'ncsf_rsq':
        fit_mat[:, idx] = r2_score_surf(bold_signal=raw_data, model_prediction=pred_mat)
    else:
        fit_mat[valid_vertices, idx] = np.asarray(layout[name], dtype=float)

if 'loo-avg' in input_fn:
    _, loo_bold = load_surface(fn=input_fn.replace('loo-avg-', 'loo-'))
    fit_mat = np.column_stack((fit_mat, r2_score_surf(bold_signal=loo_bold, model_prediction=pred_mat)))
    maps_names_ncsf = maps_names_ncsf + ['ncsf_loo_rsq']

nb.save(make_surface_image(data=fit_mat.T, source_img=img, maps_names=maps_names_ncsf),
        '{}/{}'.format(ncsf_fit_dir, ncsf_fit_fn))
print(f"Saved fit to: {ncsf_fit_dir}/{ncsf_fit_fn}")
nb.save(make_surface_image(data=pred_mat, source_img=img), '{}/{}'.format(ncsf_fit_dir, ncsf_pred_fn))
print(f"Saved prediction to: {ncsf_fit_dir}/{ncsf_pred_fn}")

if torch.cuda.is_available():
    print(f"Peak GPU memory: {torch.cuda.max_memory_allocated() / 1e9:.1f} GB "
          f"/ {torch.cuda.get_device_properties(0).total_memory / 1e9:.1f} GB")
end_time = datetime.datetime.now()
print("\nStart time:\t{}\nEnd time:\t{}\nDuration:\t{}".format(start_time, end_time, end_time - start_time))
