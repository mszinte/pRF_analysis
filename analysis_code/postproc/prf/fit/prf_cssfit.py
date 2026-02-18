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
-----------------------------------------------------------------------------------------
Output(s):
fit tester numpy arrays
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/pRF_analysis/analysis_code/postproc/prf/fit
2. run python command
python prf_cssfit.py [main directory] [project name] [subject name] 
                     [input file name] [number of jobs]
-----------------------------------------------------------------------------------------
Exemple:
python prf_cssfit.py /scratch/mszinte/data RetinoMaps sub-03 [file path] 32  
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
import datetime
import numpy as np

# MRI analysis imports
import nibabel as nb
from prfpy.stimulus import PRFStimulus2D
from prfpy.model import Iso2DGaussianModel, CSS_Iso2DGaussianModel
# from prfpy.fit import Iso2DGaussianFitter, CSS_Iso2DGaussianFitter

# Personal imports
sys.path.append("{}/../../../utils".format(os.getcwd()))
from maths_utils import r2_score_surf
from settings_utils import load_settings
from screen_utils import get_screen_settings
from pycortex_utils import set_pycortex_config_file
from surface_utils import load_surface ,make_surface_image
from prfpy_utils import Iso2DGaussianFitter, CSS_Iso2DGaussianFitter

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
css_params_num = 9

# Load settings
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../../"))
settings_path = os.path.join(base_dir, project_dir, "settings.yml")
prf_settings_path = os.path.join(base_dir, project_dir, "prf-analysis.yml")
settings = load_settings([settings_path, prf_settings_path])
analysis_info = settings[0]

TR = analysis_info['TR']
gauss_grid_nr = analysis_info['gauss_grid_nr']
max_ecc_size = analysis_info['max_ecc_size']
n_th = analysis_info['n_th']
rsq_iterative_th = analysis_info['rsq_iterative_th']
css_grid_nr = analysis_info['css_grid_nr']
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

# Set pycortex db and colormaps
cortex_dir = "{}/{}/derivatives/pp_data/cortex".format(main_dir, project_dir)
set_pycortex_config_file(cortex_dir)

# Define directories and files names (fn)
if input_fn.endswith('.nii'):
    prf_fit_dir = "{}/{}/derivatives/pp_data/{}/170k/prf/fit".format(
        main_dir, project_dir, subject)
    os.makedirs(prf_fit_dir, exist_ok=True)

elif input_fn.endswith('.gii'):
    prf_fit_dir = "{}/{}/derivatives/pp_data/{}/fsnative/prf/fit".format(
        main_dir, project_dir, subject)
    os.makedirs(prf_fit_dir, exist_ok=True)

css_fit_fn = input_fn.split('/')[-1]
css_fit_fn = css_fit_fn.replace('bold', 'prf-css_fit')

css_pred_fn = input_fn.split('/')[-1] 
css_pred_fn = css_pred_fn.replace('bold', 'prf-css_pred')

# Get task specific visual design matrix
vdm_fn = '{}/{}/derivatives/vdm/task-{}_vdm.npy'.format(
    main_dir, project_dir, prf_task_name)
vdm = np.load(vdm_fn)

# Define model parameter grid range
sizes = max_ecc_size * np.linspace(0.1, 1, gauss_grid_nr)**2
eccs = max_ecc_size * np.linspace(0.1, 1, gauss_grid_nr)**2
polars = np.linspace(0, 2 * np.pi, gauss_grid_nr)
exponent_css_grid = np.linspace(n_th[0], n_th[1], css_grid_nr)

# Load data
img, data = load_surface(fn=input_fn)

# Exclude vertices with all-NaN timeseries to avoid errors during fitting
valid_vertices = ~np.isnan(data).any(axis=0)
valid_vertices_idx = np.where(valid_vertices)[0]

# Filter data to only include valid vertices
data = data[:, valid_vertices]

n_excluded = len(valid_vertices_idx) - valid_vertices.sum()
if n_excluded > 0:
    print(f"Excluded {n_excluded} vertices with all-NaN values")
print(f"Fitting {valid_vertices.sum()} valid vertices")

# Determine visual design
stimulus = PRFStimulus2D(screen_size_cm=screen_size_cm[1],
                         screen_distance_cm=screen_distance_cm,
                         design_matrix=vdm, 
                         TR=TR)

print("\n===== PRF MODEL PARAMETERS =====")
print("Stimulus x min/max (deg):", np.nanmin(stimulus.x_coordinates ), np.nanmax(stimulus.x_coordinates ))
print("Stimulus y min/max (deg):", np.nanmin(stimulus.y_coordinates), np.nanmax(stimulus.y_coordinates))
print("Eccentricity grid range:", np.min(eccs), np.max(eccs))
print("CSS exponent grid:", exponent_css_grid)
print("Size grid range:", np.min(sizes), np.max(sizes))
print("==============================\n")

# Gauss fit
# ---------
gauss_bounds = [(-max_ecc_size, max_ecc_size),  # x
                (-max_ecc_size, max_ecc_size),  # y
                (size_th[0], size_th[1]),  # prf size
                (prf_amp_th[0], prf_amp_th[1]), # prf amplitude
                (-2, 2),# bold baseline
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

# CSS fit
# -------
css_bounds = [(-max_ecc_size, max_ecc_size),  # x
              (-max_ecc_size, max_ecc_size),  # y
              (size_th[0], size_th[1]),  # prf size
              (prf_amp_th[0], prf_amp_th[1]), # prf amplitude
              (-2, 2),  # bold baseline 
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

# Define gauss model
gauss_model = Iso2DGaussianModel(stimulus=stimulus)

# Grid fit gauss model
gauss_fitter = Iso2DGaussianFitter(data=data.T, 
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

# Define CSS model
css_model = CSS_Iso2DGaussianModel(stimulus=stimulus)

# Grid fit CSS model
css_fitter = CSS_Iso2DGaussianFitter(data=data.T, 
                                     model=css_model, 
                                     n_jobs=n_jobs,
                                     use_previous_gaussian_fitter_hrf=False,
                                     previous_gaussian_fitter=gauss_fitter)
# Run css grid fit
css_fitter.grid_fit(exponent_grid=exponent_css_grid,
                    verbose=verbose,
                    n_batches=n_batches)

# Run iterative fit
css_fitter.iterative_fit(rsq_threshold=rsq_iterative_th, 
                         verbose=verbose, 
                         bounds=css_bounds,
                         xtol=1e-4, 
                         ftol=1e-4)

css_fit = css_fitter.iterative_search_params

# Rearrange result of CSS model
css_fit_mat = np.full((data.shape[1], css_params_num), np.nan, dtype=float)
css_pred_mat = np.full_like(data, np.nan, dtype=float) 

for est, vert in enumerate(valid_vertices_idx):
    css_fit_mat[vert] = css_fit[est]
    css_pred_mat[:,vert] = css_model.return_prediction(mu_x=css_fit[est][0],
                                                      mu_y=css_fit[est][1], 
                                                      size=css_fit[est][2], 
                                                      beta=css_fit[est][3], 
                                                      baseline=css_fit[est][4],
                                                      n=css_fit[est][5],
                                                      hrf_1=css_fit[est][6],
                                                      hrf_2=css_fit[est][7])
if 'loo-avg' in input_fn:
    # Compute loo r2
    loo_bold_fn = input_fn.replace('loo-avg-', 'loo-')
    loo_img, loo_bold = load_surface(fn=loo_bold_fn)
    loo_r2 = r2_score_surf(bold_signal=loo_bold, model_prediction=css_pred_mat)
    
    # Add loo r2 css_fit_mat
    css_fit_mat = np.column_stack((css_fit_mat, loo_r2))

    # Export data
    maps_names = ['mu_x', 'mu_y', 'prf_size', 'prf_amplitude', 'bold_baseline',
                  'n', 'hrf_1','hrf_2', 'r_squared', 'loo_r_squared']
    
else:
    # Export data from CSS model fit
    maps_names = ['mu_x', 'mu_y', 'prf_size', 'prf_amplitude', 'bold_baseline',
                  'n', 'hrf_1', 'hrf_2', 'r_squared']
    
# Export fit
img_css_fit_mat = make_surface_image(data=css_fit_mat.T, source_img=img, maps_names=maps_names)
nb.save(img_css_fit_mat,'{}/{}'.format(prf_fit_dir, css_fit_fn))

# Export prediction
img_css_pred_mat = make_surface_image(data=css_pred_mat, source_img=img)
nb.save(img_css_pred_mat,'{}/{}'.format(prf_fit_dir, css_pred_fn))

# Print duration
end_time = datetime.datetime.now()
print("\nStart time:\t{start_time}\nEnd time:\t{end_time}\nDuration:\t{dur}".format(start_time=start_time, 
                                                                                    end_time=end_time, 
                                                                                    dur=end_time - start_time))