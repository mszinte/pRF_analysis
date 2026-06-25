"""
-----------------------------------------------------------------------------------------
compute_residuals.py
-----------------------------------------------------------------------------------------
Goal of the script:
Compute residuals between BOLD data and PMF gauss predictions.
Finds all necessary files automatically from settings — no separate job submission needed.
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject name (e.g. sub-01)
sys.argv[4]: analysis name (e.g. pmf)
sys.argv[5]: group (e.g. 327)
-----------------------------------------------------------------------------------------
Output(s):
Residual files saved to: .../func/fmriprep_dct_z-score_concat-residuals/
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/pRF_analysis/RetinoMaps/postproc/pmf/postfit
2. run python command
>> python compute_residuals.py [main directory] [project name] [subject] [analysis name] [group]
-----------------------------------------------------------------------------------------
Example:
python compute_residuals.py /scratch/mszinte/data RetinoMaps sub-01 pmf 327
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
import glob
import numpy as np
import nibabel as nb

# Personal imports
script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.abspath(os.path.join(script_dir, "../../../../analysis_code/utils")))
from surface_utils import make_surface_image, load_surface
from pycortex_utils import set_pycortex_config_file
from settings_utils import load_settings

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
analysis_name = sys.argv[4]
group = sys.argv[5]

# Load settings
base_dir = os.path.abspath(os.path.join(script_dir, "../../../../"))
general_settings_path = os.path.join(base_dir, project_dir, "settings.yml")
analysis_settings_path = os.path.join(base_dir, project_dir, f"{analysis_name}-analysis.yml")
settings = load_settings([general_settings_path, analysis_settings_path])
analysis_info = settings[0]

formats = analysis_info['formats']
extensions = analysis_info['extensions']
task_names = analysis_info['analysis_task_names']
preproc_prep = analysis_info['preproc_prep']
filtering = analysis_info['filtering']
normalization = analysis_info['normalization']
avg_methods = analysis_info['avg_methods']
output_folder = analysis_info['output_folder']
dm_name = analysis_info['dm_name']

# Set pycortex config
cortex_dir = "{}/{}/derivatives/pp_data/cortex".format(main_dir, project_dir)
set_pycortex_config_file(cortex_dir)

pp_dir = "{}/{}/derivatives/pp_data".format(main_dir, project_dir)

# Only process concat avg_method
for avg_method in avg_methods:
    if avg_method != 'concat':
        continue

    for format_, extension in zip(formats, extensions):

        func_dir = f'{pp_dir}/{subject}/{format_}/func/{preproc_prep}_{filtering}_{normalization}_{avg_method}'
        fit_dir  = f'{pp_dir}/{subject}/{format_}/{output_folder}/fit'

        for task_name in task_names:
            print(f'\n{"="*72}')
            print(f'Processing: {avg_method} - {format_} - {task_name}')

            # Find BOLD files 
            if format_ == 'fsnative':
                bold_fns = glob.glob(f'{func_dir}/*task-{task_name}_hemi-L*_{avg_method}_bold.{extension}')
                bold_fns += glob.glob(f'{func_dir}/*task-{task_name}_hemi-R*_{avg_method}_bold.{extension}')
            elif format_ == '170k':
                bold_fns = glob.glob(f'{func_dir}/*task-{task_name}_*_{avg_method}_bold*.{extension}')

            if not bold_fns:
                print(f'No BOLD files found in: {func_dir}')
                continue

            # Find pred files 
            pred_fns = glob.glob(f'{fit_dir}/*task-{task_name}_*_{avg_method}*{analysis_name}-css{dm_name}_pred.{extension}')

            if not pred_fns:
                print(f'No pred files found in: {fit_dir}')
                continue

            for bold_fn in bold_fns:
                print(f'\nBOLD: {bold_fn}')

                # Match pred to bold by hemisphere (fsnative) or directly (170k)
                if format_ == 'fsnative':
                    hemi = 'hemi-L' if 'hemi-L' in os.path.basename(bold_fn) else 'hemi-R'
                    matched_preds = [fn for fn in pred_fns if hemi in os.path.basename(fn)]
                elif format_ == '170k':
                    hemi = None
                    matched_preds = pred_fns

                if not matched_preds:
                    print(f'No matching pred file found for: {bold_fn}')
                    continue

                pred_fn = matched_preds[0]
                print(f'Pred: {pred_fn}')

                # Build output path and filename
                residuals_dir = f'{pp_dir}/{subject}/{format_}/func/{preproc_prep}_{filtering}_{normalization}_concat-residuals'
                os.makedirs(residuals_dir, exist_ok=True)

                if format_ == 'fsnative':
                    residual_fn = f'{subject}_task-{task_name}_{hemi}_{preproc_prep}_{filtering}_{normalization}_concat-residuals_bold.{extension}'
                elif format_ == '170k':
                    residual_fn = f'{subject}_task-{task_name}_{preproc_prep}_{filtering}_{normalization}_concat-residuals_bold.{extension}'

                residual_fullfn = os.path.join(residuals_dir, residual_fn)
                print(f'Output: {residual_fullfn}')

                # Load data
                print('Loading pred data...')
                pred_img, pred_data = load_surface(fn=pred_fn)
                print('Loading BOLD data...')
                bold_img, bold_data = load_surface(fn=bold_fn)

                print(f'  BOLD shape: {bold_data.shape}')
                print(f'  Pred shape: {pred_data.shape}')

                # Sanity check
                if bold_data.shape != pred_data.shape:
                    raise ValueError(
                        f'Shape mismatch — BOLD {bold_data.shape} vs pred {pred_data.shape}.'
                    )

                # NaN mask
                nan_mask = (
                    np.isnan(bold_data).any(axis=0) |
                    np.isnan(pred_data).any(axis=0)
                )
                clean_vox = np.where(~nan_mask)[0]
                print(f'  Clean voxels: {len(clean_vox)} / {bold_data.shape[1]}')

                # Compute residuals
                residual_data = np.full_like(bold_data, np.nan)
                residual_data[:, clean_vox] = bold_data[:, clean_vox] - pred_data[:, clean_vox]

                # Save
                residual_img = make_surface_image(data=residual_data, source_img=bold_img)
                nb.save(residual_img, residual_fullfn)
                print(f'Saved: {residual_fullfn}')

print('\nDone — all residuals computed.')

# Permissions
print('Changing files permissions in {}/{}'.format(main_dir, project_dir))
os.system("chmod -Rf 771 {}/{}".format(main_dir, project_dir))
os.system("chgrp -Rf {} {}/{}".format(group, main_dir, project_dir))