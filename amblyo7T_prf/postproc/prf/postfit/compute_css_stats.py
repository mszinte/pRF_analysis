"""
-----------------------------------------------------------------------------------------
compute_css_stats.py
-----------------------------------------------------------------------------------------
Goal of the script:
Compute the linear regression between the CSS pRF predictions and the bold signal
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject name (e.g. sub-01)
sys.argv[4]: group (e.g. 327)
sys.argv[5]: OPTIONAL main analysis folder (e.g. prf_em_ctrl)
-----------------------------------------------------------------------------------------
Output(s):
results of linear regression 
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/pRF_analysis/amblyo7T_prf/postproc/prf/postfit/
2. run python command
>> python compute_css_stats.py [main directory] [project name] 
                               [subject num] [group] [analysis folder - optional]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/projects/pRF_analysis/amblyo7T_prf/postproc/prf/postfit/
python compute_css_stats.py /scratch/mszinte/data amblyo7T_prf sub-01 327 pRFRightEye
-----------------------------------------------------------------------------------------
Written by Uriel Lascombes (uriel.lascombes@laposte.net)
Edited by Martin Szinte (martin.szinte@gmail.com)
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
import re
import sys
import glob
import json
import numpy as np
import nibabel as nb
from scipy import stats

# Personal imports
sys.path.append("{}/../../../../analysis_code/utils".format(os.getcwd()))
from pycortex_utils import set_pycortex_config_file
from surface_utils import make_surface_image, load_surface
from maths_utils import linear_regression_surf, multipletests_surface, median_subject_template

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
group = sys.argv[4]
if len(sys.argv) > 5: output_folder = sys.argv[5]
else: output_folder = "prf"

# load settings
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../../"))
settings_path = os.path.join(base_dir, project_dir, "settings.json")

with open(settings_path) as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)
fdr_alpha = analysis_info['fdr_alpha']
formats = analysis_info['formats']
extensions = analysis_info['extensions']
maps_names = analysis_info['maps_names_css_stats']
prf_task_name = analysis_info['prf_task_name']
subjects = analysis_info['subjects']

# Set pycortex db and colormaps
cortex_dir = "{}/{}/derivatives/pp_data/cortex".format(main_dir, project_dir)
set_pycortex_config_file(cortex_dir)

# Index
slope_idx, intercept_idx, rvalue_idx, pvalue_idx, stderr_idx, \
    trs_idx, corr_pvalue_5pt_idx, corr_pvalue_1pt_idx = 0, 1, 2, 3, 4, 5, 6, 7

# sub-170k exception
if subject != 'sub-170k':
    pp_dir = "{}/{}/derivatives/pp_data".format(main_dir, project_dir)
    
    for format_, extension in zip(formats, extensions): 
        print(format_)
        
        # Find pRF fit files 
        prf_fit_dir = '{}/{}/{}/{}/fit'.format(pp_dir, subject, format_, output_folder)
        prf_bold_dir = '{}/{}/{}/func/fmriprep_dct_concat'.format(pp_dir, subject, format_)
        
        prf_pred_fns_list = glob.glob('{}/*task-{}*_prf-pred_css.{}'.format(
            prf_fit_dir, prf_task_name, extension))

        for prf_pred_fn in prf_pred_fns_list:
            # Find the corresponding bold signal
            if format_ == 'fsnative': 
                hemi = re.search(r'hemi-(\w)', prf_pred_fn).group(1)
                prf_bold_fn = '{}/{}_task-{}_hemi-{}_dct_concat_bold.{}'.format(
                    prf_bold_dir, subject, output_folder, hemi, extension)
            elif format_ == '170k':
                prf_bold_fn = '{}/{}_task-{}_dct_concat_bold.{}'.format(
                    prf_bold_dir, subject, output_folder, extension)

            # load data  
            pred_img, pred_data = load_surface(prf_pred_fn)
            bold_img, bold_data = load_surface(prf_bold_fn)
            
            # Compute linear regression 
            results = linear_regression_surf(bold_signal=bold_data, 
                                             model_prediction=pred_data, 
                                             correction='fdr_tsbh', 
                                             alpha=fdr_alpha)
            
            # Compute two-sided corrected p-values
            t_statistic = results[slope_idx, :] / results[stderr_idx, :]
            degrees_of_freedom = np.nanmax(results[trs_idx, :]) - 2
            p_values = 2 * (1 - stats.t.cdf(abs(t_statistic), df=degrees_of_freedom)) 
            corrected_p_values = multipletests_surface(pvals=p_values, 
                                                       correction='fdr_tsbh', 
                                                       alpha=fdr_alpha)
            results[pvalue_idx, :] = p_values
            results[corr_pvalue_5pt_idx, :] = corrected_p_values[0,:]
            results[corr_pvalue_1pt_idx, :] = corrected_p_values[1,:]
            
            # Save results
            prf_deriv_dir = "{}/{}/{}/{}/prf_derivatives".format(
                pp_dir, subject, format_, output_folder)
            os.makedirs(prf_deriv_dir, exist_ok=True)
            
            stat_prf_fn = prf_pred_fn.split('/')[-1].replace('pred_css', 'stats')
            stat_prf_img = make_surface_image(data=results, 
                                              source_img=bold_img, 
                                              maps_names=maps_names)
            print('Saving: {}/{}'.format(prf_deriv_dir, stat_prf_fn))
            nb.save(stat_prf_img, '{}/{}'.format(prf_deriv_dir, stat_prf_fn))

# Sub-170k median across subjects
elif subject == 'sub-170k':
    print('sub-170k, computing median prf stats across subjects...')
    
    # find all the subject prf derivatives
    subjects_stats = []
    for subj in subjects: 
        stats_fn = "{}/{}/derivatives/pp_data/{}/170k/{}/prf_derivatives/{}_task-{}_fmriprep_dct_concat_prf-stats.dtseries.nii".format(
            main_dir, project_dir, subj, output_folder, subj, prf_task_name)
        if os.path.exists(stats_fn):
            subjects_stats.append(stats_fn)

    if not subjects_stats:
        print("No subject stats files found for sub-170k median computation")
        sys.exit(0)

    # Computing median across subjects
    img, data_stat_median = median_subject_template(fns=subjects_stats)
    
    # Compute two sided corrected p-values
    t_statistic = data_stat_median[slope_idx, :] / data_stat_median[stderr_idx, :]
    degrees_of_freedom = np.nanmax(data_stat_median[trs_idx, :]) - 2
    p_values = 2 * (1 - stats.t.cdf(abs(t_statistic), df=degrees_of_freedom)) 
    corrected_p_values = multipletests_surface(pvals=p_values, 
                                               correction='fdr_tsbh', 
                                               alpha=fdr_alpha)
    data_stat_median[pvalue_idx, :] = p_values
    data_stat_median[corr_pvalue_5pt_idx, :] = corrected_p_values[0,:]
    data_stat_median[corr_pvalue_1pt_idx, :] = corrected_p_values[1,:]
        
    # Export results
    sub_170k_stats_dir = "{}/{}/derivatives/pp_data/sub-170k/170k/{}/prf_derivatives/".format(
        main_dir, project_dir, output_folder)
    os.makedirs(sub_170k_stats_dir, exist_ok=True)
    
    sub_170k_stat_fn = "{}/sub-170k_task-{}_fmriprep_dct_concat_prf-stats.dtseries.nii".format(
        sub_170k_stats_dir, prf_task_name)
    print("save: {}".format(sub_170k_stat_fn))
    sub_170k_stat_img = make_surface_image(data=data_stat_median, 
                                           source_img=img, 
                                           maps_names=maps_names)
    nb.save(sub_170k_stat_img, sub_170k_stat_fn)

# # Define permission cmd
# print('Changing files permissions in {}/{}'.format(main_dir, project_dir))
# os.system("chmod -Rf 771 {}/{}".format(main_dir, project_dir))
# os.system("chgrp -Rf {} {}/{}".format(group, main_dir, project_dir))
