"""
-----------------------------------------------------------------------------------------
prf_cssfit.py
-----------------------------------------------------------------------------------------
Goal of the script:
Prf fit computing css fit
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject name
sys.argv[4]: input file name (path to the data to fit)
sys.argv[5]: number of jobs 
sys.argv[6]: OPTIONAL main analysis folder (e.g. prf_em_ctrl)
sys.argv[7]: OPTIONAL session number for freesurfer (e.g. ses-01)
sys.argv[8]: OPTIONAL filter_rois (0 or 1, default=1) - whether to filter NaN vertices
-----------------------------------------------------------------------------------------
Output(s):
fit tester numpy arrays
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/pRF_analysis/amblyo7T_prf/postproc/prf/fit/
2. run python command
python prf_cssfit.py [main directory] [project name] [subject name] 
[inout file name] [number of jobs] [analysis folder - optional] [session - optional] [filter_rois - optional]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/projects/pRF_analysis/amblyo7T_prf/postproc/prf/fit/
python prf_cssfit.py /scratch/mszinte/data amblyo7T_prf sub-01  /scratch/mszinte/data/amblyo7T_prf/derivatives/pp_data/sub-01/fsnative/func/fmriprep_dct_concat/sub-01_task-pRFRightEye_hemi-R_dct_concat_bold.func.gii 32 pRFRightEye ses-01
python prf_cssfit.py /scratch/mszinte/data amblyo7T_prf sub-01  /scratch/mszinte/data/amblyo7T_prf/derivatives/pp_data/sub-01/fsnative/func/fmriprep_dct_concat/sub-01_task-pRFRightEye_hemi-R_dct_concat_bold.func.gii 32 pRFRightEye ses-01 0
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
from prfpy.model import Iso2DGaussianModel, CSS_Iso2DGaussianModel
from prfpy.fit import Iso2DGaussianFitter, CSS_Iso2DGaussianFitter

# Personal imports
sys.path.append("{}/../../../../analysis_code/utils".format(os.getcwd()))
from surface_utils import make_surface_image, load_surface
from pycortex_utils import data_from_rois, set_pycortex_config_file
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
if len(sys.argv) > 7: session = sys.argv[7]
else: session = None
if len(sys.argv) > 8: filter_rois = bool(int(sys.argv[8]))
else: filter_rois = True

# Handle session parameter for pycortex subject name
if session:
    pycortex_subject = f"{subject}_{session}"
else:
    pycortex_subject = subject


n_batches = n_jobs
verbose = True
css_params_num = 9

# Analysis parameters
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../../"))
settings_path = os.path.join(base_dir, project_dir, 'settings.json')

with open(settings_path) as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)

# Load screen settings from subject dependend task-events.json
screen_size_cm, screen_distance_cm = get_screen_settings(main_dir, project_dir, sub_num, output_folder)

TR = analysis_info['TR']
vdm_width = analysis_info['vdm_size_pix'][0] 
vdm_height = analysis_info['vdm_size_pix'][1]
gauss_grid_nr = analysis_info['gauss_grid_nr']
max_ecc_size = analysis_info['max_ecc_size']
n_th = analysis_info['n_th']
rsq_iterative_th = analysis_info['rsq_iterative_th']
css_grid_nr = analysis_info['css_grid_nr']
ecc_th = analysis_info['ecc_th']
size_th = analysis_info['size_th']
prf_task_name = output_folder

# Set pycortex db and colormaps
cortex_dir = "{}/{}/derivatives/pp_data/cortex".format(main_dir, project_dir)
set_pycortex_config_file(cortex_dir)

# Define directories
if input_fn.endswith('.nii'):
    prf_fit_dir = "{}/{}/derivatives/pp_data/{}/170k/{}/fit".format(
        main_dir, project_dir, subject, output_folder)
    os.makedirs(prf_fit_dir, exist_ok=True)
    rois = analysis_info['mmp_rois']

elif input_fn.endswith('.gii'):
    prf_fit_dir = "{}/{}/derivatives/pp_data/{}/fsnative/{}/fit".format(
        main_dir, project_dir, subject, output_folder)
    os.makedirs(prf_fit_dir, exist_ok=True)
    rois = analysis_info['rois']

fit_fn_css = input_fn.split('/')[-1]
fit_fn_css = fit_fn_css.replace('bold', 'prf-fit_css')

pred_fn_css = input_fn.split('/')[-1]
pred_fn_css = pred_fn_css.replace('bold', 'prf-pred_css')

# Get task specific visual design matrix
vdm_fn = '{}/{}/derivatives/vdm/vdm_{}_{}_{}.npy'.format(
    main_dir, project_dir, prf_task_name, vdm_width, vdm_height)
vdm = np.load(vdm_fn)

# Define model parameter grid range
sizes = max_ecc_size * np.linspace(0.1, 1, gauss_grid_nr)**2
eccs = max_ecc_size * np.linspace(0.1, 1, gauss_grid_nr)**2
polars = np.linspace(0, 2 * np.pi, gauss_grid_nr)
exponent_css_grid = np.linspace(n_th[0], n_th[1], css_grid_nr)

print("\n===== PRF FIT PARAMETERS =====")
print(f"Screen Size (cm): {screen_size_cm}")
print(f"Screen Distance (cm): {screen_distance_cm}")
print(f"TR: {TR}")
print(f"Max eccentricity/size values: {max_ecc_size}")
print(f"CSS exponent grid: {exponent_css_grid}")
print("==============================\n")


# Load data
img, data, data_roi, roi_idx = data_from_rois(fn=input_fn, 
                                              subject=pycortex_subject,  # Changed from subject
                                              rois=rois,
                                              filter_rois=filter_rois)

print('roi extraction done')

# Exclude vertices with all-NaN timeseries to avoid errors during fitting
valid_vertices = ~np.isnan(data_roi).any(axis=0)
valid_vertices_idx = np.where(valid_vertices)[0]
original_roi_idx = roi_idx.copy()
roi_idx = original_roi_idx[valid_vertices_idx]
data_roi = data_roi[:, valid_vertices]
n_excluded = len(valid_vertices_idx) - valid_vertices.sum()
if n_excluded > 0:
    print(f"Excluded {n_excluded} vertices with all-NaN values")
print(f"Fitting {valid_vertices.sum()} valid vertices")

# Determine stimulus
stimulus = PRFStimulus2D(screen_size_cm=screen_size_cm[0],
                         screen_distance_cm=screen_distance_cm,
                         design_matrix=vdm, 
                         TR=TR)

print("\n===== PRF MODEL PARAMETERS =====")
print("Stimulus x min/max (deg):", np.nanmin(stimulus.x_coordinates), np.nanmax(stimulus.x_coordinates))
print("Stimulus y min/max (deg):", np.nanmin(stimulus.y_coordinates), np.nanmax(stimulus.y_coordinates))
print("Eccentricity grid range:", np.min(eccs), np.max(eccs))
print("Size grid range:", np.min(sizes), np.max(sizes))
print("==============================\n")

# Gauss fit
# ---------
gauss_bounds = [(-max_ecc_size, max_ecc_size),  # x
                (-max_ecc_size, max_ecc_size),  # y
                (size_th[0], size_th[1]),  # prf size
                (-6, 6),  # prf amplitude
                (-1, 1),# bold baseline
                (0, 10),  # hrf1
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


# Determine gaussian model
gauss_model = Iso2DGaussianModel(stimulus=stimulus)

# Grid fit gauss model
gauss_fitter = Iso2DGaussianFitter(data=data_roi.T, 
                                   model=gauss_model, 
                                   n_jobs=n_jobs)

gauss_fitter.grid_fit(ecc_grid=eccs, 
                      polar_grid=polars, 
                      size_grid=sizes, 
                      verbose=verbose, 
                      n_batches=n_batches)

# Iterative fit gauss model
gauss_fitter.iterative_fit(rsq_threshold=rsq_iterative_th, 
                           bounds=gauss_bounds,
                           verbose=verbose)
gauss_fit = gauss_fitter.iterative_search_params

# CSS fit
# -------
css_bounds = [(-max_ecc_size, max_ecc_size),  # x
              (-max_ecc_size, max_ecc_size),  # y
              (size_th[0], size_th[1]),  # prf size
              (-6, 6),  # prf amplitude
              (-1, 1),  # bold baseline 
              (n_th[0], n_th[1]),  # n
              (0, 10),  # hrf1
              (0, 0) # hrf2
              ] 

print("\nCSS bounds:")
print("  x range:", css_bounds[0])
print("  y range:", css_bounds[1])
print("  prf size range:", css_bounds[2])
print("  prf amplitude range:", css_bounds[3])
print("  bold baseline range:", css_bounds[4])
print("  n exponent range:", css_bounds[5])
print("  hrf1 range:", css_bounds[6])
print("  hrf2 range:", css_bounds[7])
print("==========================\n")

# Define CSS model
css_model = CSS_Iso2DGaussianModel(stimulus=stimulus)

# Grid fit CSS model
css_fitter = CSS_Iso2DGaussianFitter(data=data_roi.T, 
                                     model=css_model, 
                                     n_jobs=n_jobs,
                                     use_previous_gaussian_fitter_hrf=False,
                                     previous_gaussian_fitter=gauss_fitter)

# Grid fit css model
css_fitter.grid_fit(exponent_grid=exponent_css_grid,
                    verbose=verbose,
                    n_batches=n_batches)

# Iterative fit css model
css_fitter.iterative_fit(rsq_threshold=rsq_iterative_th, 
                         verbose=verbose, 
                         bounds=css_bounds,
                         xtol=1e-4, 
                         ftol=1e-4)

css_fit = css_fitter.iterative_search_params

# Rearrange result of CSS model
css_fit_mat = np.full((data.shape[1], css_params_num), np.nan, dtype=float)
css_pred_mat = np.full_like(data, np.nan, dtype=float) 

for est, vert in enumerate(roi_idx):
    css_fit_mat[vert] = css_fit[est]
    css_pred_mat[:,vert] = css_model.return_prediction(mu_x=css_fit[est][0],
                                                       mu_y=css_fit[est][1], 
                                                       size=css_fit[est][2], 
                                                       beta=css_fit[est][3], 
                                                       baseline=css_fit[est][4],
                                                       n=css_fit[est][5],
                                                       hrf_1=css_fit[est][6],
                                                       hrf_2=css_fit[est][7])

# Export data from CSS model fit
maps_names = ['mu_x', 'mu_y', 'prf_size', 'prf_amplitude', 'bold_baseline',
              'n', 'hrf_1', 'hrf_2', 'r_squared']

# Export fit
img_css_fit_mat = make_surface_image(data=css_fit_mat.T, source_img=img, maps_names=maps_names)
nb.save(img_css_fit_mat, '{}/{}'.format(prf_fit_dir, fit_fn_css)) 

# Export pred
img_css_pred_mat = make_surface_image(data=css_pred_mat, source_img=img)
nb.save(img_css_pred_mat, '{}/{}'.format(prf_fit_dir, pred_fn_css)) 

# Print duration
end_time = datetime.datetime.now()
print("\nStart time:\t{start_time}\nEnd time:\t{end_time}\nDuration:\t{dur}".format(
    start_time=start_time, end_time=end_time, dur=end_time - start_time))