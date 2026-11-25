"""
-----------------------------------------------------------------------------------------
make_od_img.py
-----------------------------------------------------------------------------------------
Goal of the script:
Compute the ocular dominance (OD) between eyes from the pRF CSS fits of each eye
and save it as a surface image.
Each eye was scanned in separate runs and has its own pRF prediction.
1. Common regressor: the pRF predictions of both eyes are z-scored on the whole
   concatenation (same scale for both eyes) and averaged.
2. One intercept per run: the run mean is removed from the regressor and from the data.
3. GLM of each eye (one regressor, all runs): beta = sum(x*y) / sum(x^2) and its
   standard error.
Then:
- od_le_re = (beta_le - beta_re) / (beta_le + beta_re) (only where both betas > 0)
- od_re_le = (beta_re - beta_le) / (beta_le + beta_re)
- z_le_re  = (beta_le - beta_re) / sqrt(se_le^2 + se_re^2), z_re_le = -z_le_re
No threshold is applied (thresholds are applied in the plotting scripts).
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject name (e.g. sub-01)
sys.argv[4]: analysis task name (ex. prf)
sys.argv[5]: group (e.g. 327)
-----------------------------------------------------------------------------------------
Output(s):
OD surface image (maps_names_css_od) in the prf_derivatives folder
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/pRF_analysis/amblyo7T_prf/postproc/prf/postfit/
2. run python command
>> python make_od_img.py [main directory] [project name]
                         [subject num] [analysis name] [group]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/projects/pRF_analysis/amblyo7T_prf/postproc/prf/postfit/
python make_od_img.py /scratch/mszinte/data amblyo7T_prf sub-01 prf 327
python make_od_img.py ~/disks/meso_shared amblyo7T_prf sub-20 prf 327
-----------------------------------------------------------------------------------------
Written by Uriel Lascombes (uriel.lascombes@laposte.net)
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
import numpy as np
import nibabel as nb

# Personal imports
sys.path.append("{}/../../../../analysis_code/utils".format(os.getcwd()))
from surface_utils import make_surface_image, load_surface
from settings_utils import load_settings

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
analysis_name = sys.argv[4]
group = sys.argv[5]

# Load settings
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../../"))
general_settings_path = os.path.join(base_dir, project_dir, "settings.yml")
analysis_settings_path = os.path.join(base_dir, project_dir, f"{analysis_name}-analysis.yml")
settings = load_settings([general_settings_path, analysis_settings_path])
analysis_info = settings[0]

formats = analysis_info['formats']
extensions = analysis_info['extensions']
preproc_prep = analysis_info['preproc_prep']
filtering = analysis_info['filtering']
normalizations = analysis_info['normalization']
avg_methods = analysis_info['avg_methods']
output_folder = analysis_info["output_folder"]
dm_name = analysis_info["dm_name"]
TRs = analysis_info['TRs']
od_task_names = analysis_info['od_task_names']
maps_names = analysis_info['maps_names_css_od']
eyes = ['LeftEye', 'RightEye']

# Define folders
pp_dir = "{}/{}/derivatives/pp_data".format(main_dir, project_dir)

for normalization in normalizations:
    for avg_method in avg_methods:
        for format_, extension in zip(formats, extensions):

            # define/create folders
            prf_fit_dir = '{}/{}/{}/{}/fit'.format(
                pp_dir, subject, format_, output_folder)
            prf_deriv_dir = "{}/{}/{}/{}/prf_derivatives".format(
                pp_dir, subject, format_, output_folder)
            func_dir = "{}/{}/{}/func/{}_{}_{}_{}".format(
                pp_dir, subject, format_, preproc_prep, filtering, normalization, avg_method)

            # hemispheres (fsnative) or whole brain (170k)
            if format_ == 'fsnative': hemis = ['_hemi-L', '_hemi-R']
            elif format_ == '170k': hemis = ['']

            for task_name in od_task_names:
                for hemi in hemis:
                    print(f'{normalization} - {avg_method} - {format_} - {task_name}{hemi}')

                    # Define pRF predictions and data files of both eyes
                    pred_fns, bold_fns = {}, {}
                    for eye in eyes:
                        fn_spec = f'{subject}_task-{task_name}{eye}{hemi}_{preproc_prep}_{filtering}_{normalization}_{avg_method}'
                        pred_fns[eye] = f'{prf_fit_dir}/{fn_spec}_{analysis_name}-css{dm_name}_pred.{extension}'
                        bold_fns[eye] = f'{func_dir}/{fn_spec}_bold.{extension}'

                    missing_fns = [fn for fn in list(pred_fns.values()) + list(bold_fns.values()) if not os.path.isfile(fn)]
                    if missing_fns:
                        print(f'[SKIP] missing files: {missing_fns}')
                        continue

                    # Load pRF predictions and data of both eyes
                    pred_data, bold_data = {}, {}
                    for eye in eyes:
                        print(f'Loading: {pred_fns[eye]}')
                        pred_img, pred_data[eye] = load_surface(pred_fns[eye])
                        print(f'Loading: {bold_fns[eye]}')
                        _, bold_data[eye] = load_surface(bold_fns[eye])
                        pred_data[eye] = pred_data[eye].astype(float)
                        bold_data[eye] = bold_data[eye].astype(float)

                    # Number of runs in the concatenation
                    n_trs, n_vertices = bold_data['LeftEye'].shape
                    if n_trs % TRs:
                        print(f'[SKIP] {n_trs} TRs is not a multiple of {TRs} TRs per run')
                        continue
                    n_runs = n_trs // TRs

                    # Common regressor: mean of both eyes predictions, each z-scored on the whole concatenation
                    pred_z = []
                    for eye in eyes:
                        with np.errstate(divide='ignore', invalid='ignore'):
                            pred_z.append(np.nan_to_num((pred_data[eye] - pred_data[eye].mean(axis=0)) / pred_data[eye].std(axis=0)))
                    common_reg = (pred_z[0] + pred_z[1]) / 2

                    # One intercept per run: run mean removed from the regressor
                    common_reg = common_reg.reshape(n_runs, TRs, n_vertices)
                    common_reg = (common_reg - common_reg.mean(axis=1, keepdims=True)).reshape(n_trs, n_vertices)
                    sum_xx = (common_reg ** 2).sum(axis=0)

                    # GLM of each eye: beta and standard error
                    dof = n_trs - n_runs - 1
                    betas, ses = {}, {}
                    for eye in eyes:
                        # One intercept per run: run mean removed from the data
                        bold = bold_data[eye].reshape(n_runs, TRs, n_vertices)
                        bold = (bold - bold.mean(axis=1, keepdims=True)).reshape(n_trs, n_vertices)
                        sum_xy = (common_reg * bold).sum(axis=0)
                        sum_yy = (bold ** 2).sum(axis=0)
                        with np.errstate(divide='ignore', invalid='ignore'):
                            betas[eye] = sum_xy / sum_xx
                            ses[eye] = np.sqrt(np.clip(sum_yy - betas[eye] * sum_xy, 0, None) / dof / sum_xx)

                    # Ocular dominance index (both betas > 0) and z
                    beta_le, beta_re = betas['LeftEye'], betas['RightEye']
                    se_le, se_re = ses['LeftEye'], ses['RightEye']
                    with np.errstate(divide='ignore', invalid='ignore'):
                        od_le_re = (beta_le - beta_re) / (beta_le + beta_re)
                        z_le_re = (beta_le - beta_re) / np.sqrt(se_le ** 2 + se_re ** 2)
                    od_le_re[~((beta_le > 0) & (beta_re > 0))] = np.nan
                    od_re_le = -od_le_re
                    z_re_le = -z_le_re

                    # Export results (order of maps_names_css_od)
                    od_maps = {'od_beta_le': beta_le, 'od_beta_re': beta_re,
                               'od_se_le': se_le, 'od_se_re': se_re,
                               'od_le_re': od_le_re, 'od_re_le': od_re_le,
                               'od_z_le_re': z_le_re, 'od_z_re_le': z_re_le}
                    od_mat = np.vstack([od_maps[maps_name] for maps_name in maps_names])

                    od_fn = f'{prf_deriv_dir}/{subject}_task-{task_name}{hemi}_{preproc_prep}_{filtering}_{normalization}_{avg_method}_{analysis_name}-css{dm_name}_od.{extension}'
                    od_img = make_surface_image(data=od_mat, source_img=pred_img, maps_names=maps_names)
                    os.makedirs(prf_deriv_dir, exist_ok=True)
                    nb.save(od_img, od_fn)
                    print(f'Saving OD: {od_fn}')

# # Define permission cmd
# print('Changing files permissions in {}/{}'.format(main_dir, project_dir))
# os.system("chmod -Rf 771 {}/{}".format(main_dir, project_dir))
# os.system("chgrp -Rf {} {}/{}".format(group, main_dir, project_dir))
