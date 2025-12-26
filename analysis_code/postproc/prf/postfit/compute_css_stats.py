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
-----------------------------------------------------------------------------------------
Output(s):
results of linear regression 
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/pRF_analysis/analysis_code/postproc/prf/postfit/
2. run python command
>> python compute_css_stats.py [main directory] [project name] 
                               [subject num] [group]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/projects/pRF_analysis/analysis_code/postproc/prf/postfit/
python compute_css_stats.py /scratch/mszinte/data MotConf sub-01 327
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
sys.path.append("{}/../../../utils".format(os.getcwd()))
from pycortex_utils import set_pycortex_config_file
from surface_utils import make_surface_image , load_surface
from maths_utils import linear_regression_surf, multipletests_surface, median_subject_template

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
group = sys.argv[4]

# load settings
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../../"))
settings_path = os.path.join(base_dir, project_dir, "settings.json")

with open(settings_path) as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)
fdr_alpha = analysis_info['fdr_alpha']
formats = analysis_info['formats']
extensions = analysis_info['extensions']
TRs = analysis_info['TRs']
maps_names = analysis_info['maps_names_css_stats']
prf_task_names = analysis_info['prf_task_names']
subjects = analysis_info['subjects']
preproc_prep = analysis_info['preproc_prep']
filtering = analysis_info['filtering']
normalization = analysis_info['normalization']
avg_methods = analysis_info['avg_methods']

# Set pycortex db and colormaps
cortex_dir = "{}/{}/derivatives/pp_data/cortex".format(main_dir, project_dir)
set_pycortex_config_file(cortex_dir)

# Index
slope_idx, intercept_idx, rvalue_idx, pvalue_idx, stderr_idx, \
    trs_idx, corr_pvalue_5pt_idx, corr_pvalue_1pt_idx = 0, 1, 2, 3, 4, 5, 6, 7

# Define preprocessing folder
pp_dir = "{}/{}/derivatives/pp_data".format(main_dir, project_dir)

# sub-170k exeption
if subject != 'sub-170k':
    for avg_method in avg_methods:
        for format_, extension in zip(formats, extensions): 
            
            # define/create folders
            prf_fit_dir = '{}/{}/{}/prf/fit'.format(
                    pp_dir, subject, format_)
            prf_func_dir = '{}/{}/{}/func/{}_{}_{}_{}'.format(
                    pp_dir, subject, format_, 
                    preproc_prep, filtering, normalization, avg_method)
            prf_deriv_dir = "{}/{}/{}/prf/prf_derivatives".format(
                pp_dir, subject, format_)
            os.makedirs(prf_deriv_dir, exist_ok=True)

            for prf_task_name in prf_task_names:
                print(f'{avg_method} - {format_} - {prf_task_name}')
                
                # Find pRF func/pred files
                prf_pred_fns = glob.glob('{}/*task-{}*_{}*_prf-css_pred.{}'.format(
                    prf_fit_dir, prf_task_name, avg_method, extension))
    
                for prf_pred_fn in prf_pred_fns :
                    if 'loo' in prf_pred_fn:
                        loo_number = re.search(r'loo-avg-(\d+)', prf_pred_fn).group(1)
                        if format_ == 'fsnative': 
                            hemi = re.search(r'hemi-(\w)', prf_pred_fn).group(1)
                            prf_bold_fn = glob.glob('{}/*task-{}_hemi-{}*_loo-{}_bold.{}'.format(
                                prf_func_dir, prf_task_name, hemi, loo_number, extension))[0]
                        elif format_ == '170k':
                            prf_bold_fn = glob.glob('{}/*task-{}*_loo-{}_bold.{}'.format(
                                prf_func_dir, prf_task_name, loo_number, extension))[0]
                    else:
                        if format_ == 'fsnative': 
                            hemi = re.search(r'hemi-(\w)', prf_pred_fn).group(1)
                            prf_bold_fn = glob.glob('{}/*task-{}_hemi-{}*_{}_bold.{}'.format(
                                prf_func_dir, prf_task_name, hemi, avg_method, extension))[0]
                        elif format_ == '170k':
                            prf_bold_fn = glob.glob('{}/*task-{}*_{}_bold.{}'.format(
                                prf_func_dir, prf_task_name, avg_method, extension))[0]
                        
                    # load data  
                    print(f'Loading pred: {prf_pred_fn}') 
                    bold_img, bold_data = load_surface(prf_bold_fn)
                    print(f'Loading bold: {prf_bold_fn}')
                    pred_img, pred_data = load_surface(prf_pred_fn)
                    
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
                    prf_stats_fn = prf_pred_fn.split('/')[-1].replace('prf-css_pred', 'prf-css_stats')
                    prf_stats_img = make_surface_image(data=results, 
                                                       source_img=pred_img, 
                                                       maps_names=maps_names)
                    
                    print('Saving: {}/{}'.format(prf_deriv_dir, prf_stats_fn))
                    nb.save(prf_stats_img, '{}/{}'.format(prf_deriv_dir, prf_stats_fn))
            
                # Compute median across leave-one-out fit
                if 'loo-avg' in avg_method:
                    print('Computing median across LOO')

                    # Get LOO files (excluding any with "median" in the name)
                    loo_prf_stats_fns = glob.glob(f"{prf_deriv_dir}/*task-{prf_task_name}*loo-avg-*_prf-css_stats.{extension}")

                    # Group files by hemisphere/format
                    loo_prf_stats_fsnative_hemi_L_fns = [fn for fn in loo_prf_stats_fns if "hemi-L" in fn]
                    loo_prf_stats_fsnative_hemi_R_fns = [fn for fn in loo_prf_stats_fns if "hemi-R" in fn]
                    loo_prf_stats_170k_fns = [fn for fn in loo_prf_stats_fns if "hemi-L" not in fn and "hemi-R" not in fn]
                    
                    # Process each group
                    for group_files, hemi in [(loo_prf_stats_fsnative_hemi_L_fns, "_hemi-L"),
                                              (loo_prf_stats_fsnative_hemi_R_fns, "_hemi-R"),
                                              (loo_prf_stats_170k_fns, "")]:
                        
                        if len(group_files)>0:
    
                            # Load first file to initialize median array and define fn
                            stats_img, stats_data = load_surface(group_files[0])
                            loo_prf_stats = np.zeros_like(stats_data)
                            loo_prf_stats_fn =  '{}/{}_task-{}{}_{}_{}_{}_loo-avg_prf-css_stats.{}'.format(
                                prf_deriv_dir, subject, prf_task_name, hemi, 
                                preproc_prep, filtering, normalization, extension)
                        
                            # Compute median across LOO runs
                            for n_run, loo_stats_fn in enumerate(group_files):
                                print(f'Loadding loo stats: {loo_stats_fn}')
                                _, loo_prf_stats_run_data = load_surface(loo_stats_fn)
                                if n_run == 0: loo_prf_stats = np.copy(loo_prf_stats_run_data)
                                else: loo_prf_stats = np.nanmedian(np.array([loo_prf_stats, loo_prf_stats_run_data]), axis=0)
                        
                            # Recalculate p-values for median data
                            t_statistic = loo_prf_stats[slope_idx, :] / loo_prf_stats[stderr_idx, :]
                            degrees_of_freedom = np.nanmax(loo_prf_stats[trs_idx, :]) - 2
                            p_values = 2 * (1 - stats.t.cdf(np.abs(t_statistic), df=degrees_of_freedom))
                            corrected_p_values = multipletests_surface(p_values, correction="fdr_tsbh", alpha=fdr_alpha)
                            
                            # Update median data with recalculated p-values
                            loo_prf_stats[pvalue_idx, :] = p_values
                            loo_prf_stats[corr_pvalue_5pt_idx, :] = corrected_p_values[0, :]
                            loo_prf_stats[corr_pvalue_1pt_idx, :] = corrected_p_values[1, :]
    
                            # Save median results
                            loo_prf_stats_img = make_surface_image(loo_prf_stats, stats_img, maps_names)
                            nb.save(loo_prf_stats_img, loo_prf_stats_fn)
                            print(f"Saving median: {loo_prf_stats_fn}")

# # Sub-170k median
# elif subject == 'sub-170k':
#     print('sub-170, computing median prf stats across subject...')
    
#     # find all the subject prf derivatives
#     subjects_stats = []
#     for subject in subjects: 
#         subjects_stats += ["{}/{}/derivatives/pp_data/{}/170k/prf/prf_derivatives/{}_task-{}_fmriprep_dct_avg_prf-stats_loo-median.dtseries.nii".format(
#                 main_dir, project_dir, subject, subject, prf_task_name)]

#     # Computing median across subject
#     img, data_stat_median = median_subject_template(fns=subjects_stats)
    
#     # Compute two sided corrected p-values
#     t_statistic = data_stat_median[slope_idx, :] / data_stat_median[stderr_idx, :]
#     degrees_of_freedom = np.nanmax(data_stat_median[trs_idx, :]) - 2
#     p_values = 2 * (1 - stats.t.cdf(abs(t_statistic), df=degrees_of_freedom)) 
#     corrected_p_values = multipletests_surface(pvals=p_values, 
#                                                correction='fdr_tsbh', 
#                                                alpha=fdr_alpha)
#     data_stat_median[pvalue_idx, :] = p_values
#     data_stat_median[corr_pvalue_5pt_idx, :] = corrected_p_values[0,:]
#     data_stat_median[corr_pvalue_1pt_idx, :] = corrected_p_values[1,:]
        
#     # Export results
#     sub_170k_stats_dir = "{}/{}/derivatives/pp_data/sub-170k/170k/prf/prf_derivatives/".format(
#             main_dir, project_dir)
#     os.makedirs(sub_170k_stats_dir, exist_ok=True)
    
#     sub_170k_stat_fn = "{}/sub-170k_task-{}_fmriprep_dct_avg_prf-stats_loo-median.dtseries.nii".format(sub_170k_stats_dir, prf_task_name)
#     print("save: {}".format(sub_170k_stat_fn))
#     sub_170k_stat_img = make_surface_image(data=data_stat_median, 
#                                            source_img=img, 
#                                            maps_names=maps_names)
#     nb.save(sub_170k_stat_img, sub_170k_stat_fn)

# Define permission cmd
print('Changing files permissions in {}/{}'.format(main_dir, project_dir))
os.system("chmod -Rf 771 {}/{}".format(main_dir, project_dir))
os.system("chgrp -Rf {} {}/{}".format(group, main_dir, project_dir))