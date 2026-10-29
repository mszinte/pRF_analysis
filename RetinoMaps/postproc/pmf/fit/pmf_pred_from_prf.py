"""
-----------------------------------------------------------------------------------------
pmf_pred_from_prf.py
-----------------------------------------------------------------------------------------
Goal of the script:
Predict the eye-movement task timecourse (e.g. SacLoc) from the visual pRF of each vertex.
The pRF parameters (x, y, size, n) come from the retinotopy task CSS fit and stay fixed;
only amplitude and baseline are estimated on the task data (prfmodel LeastSquaresFitter).
Model (CSS + SPM HRF) and grid are the same as in prf_cssfit_gpu.py.
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject name
sys.argv[4]: input file name (path to the task bold data, e.g. SacLoc concat)
-----------------------------------------------------------------------------------------
Output(s):
{subject}_task-{task}_..._pmf-prfpred_pred : predicted timecourses (NaN where no valid pRF)
{subject}_task-{task}_..._pmf-prfpred_fit  : amplitude, baseline, r_squared
-----------------------------------------------------------------------------------------
To run:
cd ~/projects/pRF_analysis/RetinoMaps/postproc/pmf/fit
python pmf_pred_from_prf.py [main directory] [project name] [subject name] [input file name]
-----------------------------------------------------------------------------------------
Exemple:
python pmf_pred_from_prf.py /scratch/mszinte/data RetinoMaps sub-03 \
    /scratch/mszinte/data/RetinoMaps/derivatives/pp_data/sub-03/fsnative/func/fmriprep_dct_z-score_concat/sub-03_task-SacLoc_hemi-L_fmriprep_dct_z-score_concat_bold.func.gii
-----------------------------------------------------------------------------------------
Written by Sina Kling (sina.kling@outlook.de)
-----------------------------------------------------------------------------------------
"""

# Stop warnings
import warnings
warnings.filterwarnings("ignore")

# General imports
import os
import sys
import datetime
import numpy as np
import pandas as pd
import nibabel as nb

# prfmodel imports
from prfmodel.stimuli import PRFStimulus
from prfmodel.impulse import DerivativeTwoGammaImpulse
from prfmodel.models.prf import Gaussian2DPRFModel, PRFStimulusEncoder
from prfmodel.models.compression import CompressiveEncoder
from prfmodel.fitters import LeastSquaresFitter
from prfmodel.utils import batched

# Personal imports
script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.abspath(os.path.join(script_dir, "../../../../analysis_code/utils")))
from surface_utils import make_surface_image, load_surface
from settings_utils import load_settings

start_time = datetime.datetime.now()

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
sub_num = subject[4:]
input_fn = sys.argv[4]

PRF_RSQ_TH = 0.1     # only vertices with a usable visual pRF are predicted


# Load settings (pmf and prf settings separately: they define different avg_methods)
base_dir = os.path.abspath(os.path.join(script_dir, "../../../../"))
general_settings_path = os.path.join(base_dir, project_dir, "settings.yml")
pmf_info = load_settings([general_settings_path, os.path.join(base_dir, project_dir, "pmf-analysis.yml")])[0]
prf_info = load_settings([general_settings_path, os.path.join(base_dir, project_dir, "prf-analysis.yml")])[0]

TR = pmf_info['TR']
dm_range = pmf_info['dm_range']
task_name = input_fn.split("task-")[1].split("_")[0]

# [NORMALIZATION] read from the input file name (e.g. ..._fmriprep_dct_psc_concat_bold.func.gii -> 'psc')
input_fn_parts = os.path.basename(input_fn).split('_')
normalization_found = [n for n in pmf_info['prf_amp_th'] if n in input_fn_parts]
if len(normalization_found) != 1:
    raise ValueError('Cannot identify normalization {} in {}'.format(list(pmf_info['prf_amp_th']), input_fn))
normalization = normalization_found[0]
prf_amp_th = pmf_info['prf_amp_th'][normalization]

# ── pRF derivatives of the retinotopy task (same file make_corr_tsv.py reads) ─────────
file_format = "170k" if input_fn.endswith(".nii") else "fsnative"
hemi_tag = "" if file_format == "170k" else ("_hemi-L" if "hemi-L" in input_fn else "_hemi-R")
ext = "dtseries.nii" if file_format == "170k" else "func.gii"
prf_avg = prf_info['avg_methods'][0]
prf_maps_names = prf_info['maps_names_css_loo'] if "loo" in prf_avg else prf_info['maps_names_css']
prf_deriv_fn = (f"{main_dir}/{project_dir}/derivatives/pp_data/{subject}/{file_format}/prf/prf_derivatives/"
                f"{subject}_task-pRF{hemi_tag}_{prf_info['preproc_prep']}_{prf_info['filtering']}_"
                f"{prf_info['normalization']}_{prf_avg}_prf-css_deriv.{ext}")

# Output directory and names
prf_fit_dir = f"{main_dir}/{project_dir}/derivatives/pp_data/{subject}/{file_format}/pmf/fit"
os.makedirs(prf_fit_dir, exist_ok=True)
suffix = "pmf-prfpred"
pred_fn = input_fn.split('/')[-1].replace('bold', f'{suffix}_pred')
fit_fn = input_fn.split('/')[-1].replace('bold', f'{suffix}_fit')

# Retinal design matrix of the task (row 0 = +y, spans ±dm_range)
vdm_fn = f"{main_dir}/{project_dir}/derivatives/vdm/{subject}/{subject}_task-{task_name}_rdm.npy"
vdm = np.load(vdm_fn)

print("\n===== PREDICTION FROM pRF =====")
print(f"Task data:     {input_fn}")
print(f"pRF params:    {prf_deriv_fn}")
print(f"Design matrix: {vdm_fn}  shape={vdm.shape}, ±{dm_range} dva")
print("===============================\n")

# Load task data (n_TRs, n_vertices) and pRF derivatives (n_maps, n_vertices)
img, raw_data = load_surface(fn=input_fn)
_, prf_deriv = load_surface(fn=prf_deriv_fn)
prf = {name: prf_deriv[i] for i, name in enumerate(prf_maps_names)}

# Vertices with finite data and a usable pRF
valid = (np.isfinite(raw_data).all(axis=0)
         & np.isfinite(prf['prf_x']) & np.isfinite(prf['prf_y'])
         & (prf['prf_size'] > 0) & (prf['prf_n'] > 0) & (prf['prf_rsq'] > PRF_RSQ_TH))
v = np.where(valid)[0]
print(f"{len(v)} of {raw_data.shape[1]} vertices predicted (pRF R2 > {PRF_RSQ_TH})")

# ── Stimulus grid (same as prf_cssfit_gpu.py) ────────────────────────────────────────
rows, cols = vdm.shape[0], vdm.shape[1]
pixel_size = 2 * dm_range / rows
x = (np.arange(cols) - (cols - 1) / 2) * pixel_size
y = -(np.arange(rows) - (rows - 1) / 2) * pixel_size       # row 0 = +y
xv, yv = np.meshgrid(x, y)
stimulus = PRFStimulus(design=np.moveaxis(vdm, -1, 0), grid=np.stack((yv, xv), axis=-1))

# ── CSS model with the pRF-fit HRF ────────────────────────────────────────────────────
impulse_model = DerivativeTwoGammaImpulse(resolution=TR, default_parameters="spm_hrf")
css_model = Gaussian2DPRFModel(encoding_model=CompressiveEncoder(encoding_model=PRFStimulusEncoder()),
                               impulse_model=impulse_model)

# Fixed pRF parameters; amplitude and baseline are the only ones estimated
params = pd.DataFrame({
    "mu_x": prf['prf_x'][v],
    "mu_y": prf['prf_y'][v],
    "sigma": prf['prf_size'][v],
    "n": prf['prf_n'][v],
    "gain": 1.0,
    "amplitude": 1.0,
    "baseline": 0.0,
    "weight_deriv": -1.0,
})

ls_fitter = LeastSquaresFitter(model=css_model, stimulus=stimulus)
_, ls_params = ls_fitter.fit(data=raw_data[:, v].T, parameters=params,
                             slope_name="amplitude", intercept_name="baseline", batch_size=10000)
pred = np.asarray(batched(css_model)(stimulus, ls_params, batch_size=1000))        # (n_vertices, n_TRs)

# R2 of the prediction
data_v = raw_data[:, v].T
rsq = 1 - ((data_v - pred) ** 2).sum(1) / ((data_v - data_v.mean(1, keepdims=True)) ** 2).sum(1)
print(f"median R2 = {np.median(rsq):.3f}, median amplitude = {np.median(ls_params['amplitude']):.3f}, "
      f"negative amplitudes: {(ls_params['amplitude'] < 0).mean():.1%}")

# ── Save ──────────────────────────────────────────────────────────────────────────────
pred_mat = np.full(raw_data.shape, np.nan, dtype=np.float32)
pred_mat[:, v] = pred.T
nb.save(make_surface_image(data=pred_mat, source_img=img), f"{prf_fit_dir}/{pred_fn}")
print(f"Saved: {prf_fit_dir}/{pred_fn}")

fit_mat = np.full((3, raw_data.shape[1]), np.nan, dtype=np.float32)
fit_mat[:, v] = np.vstack([ls_params['amplitude'], ls_params['baseline'], rsq])
nb.save(make_surface_image(data=fit_mat, source_img=img, maps_names=['prf_amplitude', 'bold_baseline', 'r_squared']),
        f"{prf_fit_dir}/{fit_fn}")
print(f"Saved: {prf_fit_dir}/{fit_fn}")

# print duration
end_time = datetime.datetime.now()
print("\nStart time:\t{}\nEnd time:\t{}\nDuration:\t{}".format(start_time, end_time, end_time - start_time))
