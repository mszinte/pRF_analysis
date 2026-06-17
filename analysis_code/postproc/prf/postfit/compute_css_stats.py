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
sys.argv[4]: analysis name (e.g prf)
sys.argv[5]: group (e.g. 327)
-----------------------------------------------------------------------------------------
Output(s):
results of linear regression 
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/pRF_analysis/analysis_code/postproc/prf/postfit/
2. run python command
>> python compute_css_stats.py [main directory] [project name] [analysis name]
                               [subject num] [group]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/projects/pRF_analysis/analysis_code/postproc/prf/postfit/
python compute_css_stats.py /scratch/mszinte/data RetinoMaps sub-01 prf 327
python compute_css_stats.py /scratch/mszinte/data RetinoMaps template_avg prf 327
python compute_css_stats.py /scratch/mszinte/data amblyo7T_prf sub-04 prf 327
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
import numpy as np
import nibabel as nb
from scipy import stats

# Personal imports
sys.path.append("{}/../../../utils".format(os.getcwd()))
from pycortex_utils import set_pycortex_config_file
from surface_utils import make_surface_image , load_surface
from maths_utils import linear_regression_surf, multipletests_surface, median_subject_template
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

fdr_alpha = analysis_info['fdr_alpha']
formats = analysis_info['formats']
extensions = analysis_info['extensions']
TRs = analysis_info['TRs']
maps_names = analysis_info['maps_names_css_stats']
task_names = analysis_info['analysis_task_names']
subjects = analysis_info['subjects']
preproc_prep = analysis_info['preproc_prep']
filtering = analysis_info['filtering']
normalization = analysis_info['normalization']
avg_methods = analysis_info['avg_methods']
averaging_templates = analysis_info['averaging_templates']
output_folder = analysis_info["output_folder"]
dm_name = analysis_info["dm_name"]

# Set pycortex db and colormaps
cortex_dir = f"{main_dir}/{project_dir}/derivatives/pp_data/cortex"
set_pycortex_config_file(cortex_dir)

# Index
slope_idx, intercept_idx, rvalue_idx, pvalue_idx, stderr_idx, \
    trs_idx, corr_pvalue_5pt_idx, corr_pvalue_1pt_idx = 0, 1, 2, 3, 4, 5, 6, 7

# Define preprocessing folder
pp_dir = f"{main_dir}/{project_dir}/derivatives/pp_data"


# template_avg exeption
if subject != 'template_avg':
    for avg_method in avg_methods:
        print(avg_method)
        for format_, extension in zip(formats, extensions): 
            
            # define/create folders
            prf_fit_dir = f'{pp_dir}/{subject}/{format_}/{output_folder}/fit'
            prf_func_dir = f'{pp_dir}/{subject}/{format_}/func/{preproc_prep}_{filtering}_{normalization}_{avg_method}'
            prf_deriv_dir = f'{pp_dir}/{subject}/{format_}/{output_folder}/prf_derivatives'
            os.makedirs(prf_deriv_dir, exist_ok=True)

            for task_name in task_names:
                print(f'{avg_method} - {format_} - {task_name}')
                # Find pRF func/pred files
                prf_pred_fns = glob.glob(f'{prf_fit_dir}/*task-{task_name}_*_{avg_method}*{analysis_name}-css{dm_name}_pred.{extension}')
                
                # NOT WORKING 
                # for prf_pred_fn in prf_pred_fns:
                #     if 'loo' in prf_pred_fn:
                #         loo_number = re.search(r'loo-avg-(\d+)', prf_pred_fn).group(1)
                #         if format_ == 'fsnative':
                #             hemi = re.search(r'hemi-(\w)', prf_pred_fn).group(1)
                #             bold_matches = glob.glob(f'{prf_func_dir}/*task-{task_name}_hemi-{hemi}*_loo-avg-{loo_number}.{extension}')
                #             if not bold_matches:
                #                 bold_matches = glob.glob(f'{prf_func_dir}/*task-{task_name}_hemi-{hemi}*_loo-{loo_number}.{extension}')
                #         elif format_ == '170k':
                #             bold_matches = glob.glob(f'{prf_func_dir}/*task-{task_name}_*_loo-avg-{loo_number}.{extension}')
                #             if not bold_matches:
                #                 bold_matches = glob.glob(f'{prf_func_dir}/*task-{task_name}_*_loo-{loo_number}.{extension}')
                #     else:                                          # <-- same level as 'if loo'
                #         if format_ == 'fsnative':
                #             hemi = re.search(r'hemi-(\w)', prf_pred_fn).group(1)
                #             bold_matches = glob.glob(f'{prf_func_dir}/*task-{task_name}_hemi-{hemi}*_{avg_method}.{extension}')
                #             if not bold_matches:
                #                 bold_matches = glob.glob(f'{prf_func_dir}/*task-{task_name}_hemi-{hemi}*_{avg_method}_bold.{extension}')
                #         elif format_ == '170k':
                #             bold_matches = glob.glob(f'{prf_func_dir}/*task-{task_name}_*_{avg_method}.{extension}')
                #             if not bold_matches:
                #                 bold_matches = glob.glob(f'{prf_func_dir}/*task-{task_name}_*_{avg_method}_bold.{extension}')

                #     prf_bold_fn = bold_matches[0]  

                
                for prf_pred_fn in prf_pred_fns :
                    if 'loo' in prf_pred_fn:
                        loo_number = re.search(r'loo-avg-(\d+)', prf_pred_fn).group(1)
                        if format_ == 'fsnative': 
                            hemi = re.search(r'hemi-(\w)', prf_pred_fn).group(1)
                            prf_bold_fn = glob.glob('{}/*task-{}_hemi-{}*_loo-{}_bold*.{}'.format(
                                prf_func_dir, task_name, hemi, loo_number, extension))[0]
                        elif format_ == '170k':
                            prf_bold_fn = glob.glob('{}/*task-{}_*_loo-{}_bold.{}'.format(
                                prf_func_dir, task_name, loo_number, extension))[0]
                    else:
                        if format_ == 'fsnative': 
                            hemi = re.search(r'hemi-(\w)', prf_pred_fn).group(1)
                            prf_bold_fn = glob.glob('{}/*task-{}_hemi-{}*_{}_bold.{}'.format(
                                prf_func_dir, task_name, hemi, avg_method, extension))[0]
                        elif format_ == '170k':
                            prf_bold_fn = glob.glob('{}/*task-{}_*_{}_bold*.{}'.format(
                                prf_func_dir, task_name, avg_method, extension))[0]

                    # load data
                    print(f'Loading bold: {prf_bold_fn}')
                    bold_img, bold_data = load_surface(prf_bold_fn)
                    print(f'Loading pred: {prf_pred_fn}')
                    pred_img, pred_data = load_surface(prf_pred_fn)
                    
                    # Compute linear regression 
                    results = linear_regression_surf(bold_signal=bold_data, 
                                                     model_prediction=pred_data, 
                                                     alternative='two-sided',
                                                     correction='fdr_tsbh', 
                                                     alpha=fdr_alpha)
                    
                    # Save results
                    prf_stats_fn = prf_pred_fn.split('/')[-1].replace(f'{analysis_name}-css{dm_name}_pred', f'{analysis_name}-css{dm_name}_stats')
                    prf_stats_img = make_surface_image(data=results, 
                                                       source_img=pred_img, 
                                                       maps_names=maps_names)

                    print(f'Saving: {prf_deriv_dir}/{prf_stats_fn}')
                    nb.save(prf_stats_img, f'{prf_deriv_dir}/{prf_stats_fn}')
                
                # Compute median across leave-one-out fit
                if 'loo-avg' in avg_method:
                    print('Computing median across LOO')

                    # Get LOO files (excluding any with "median" in the name) # sub-13_task-pRF_hemi-R_fmriprep_dct_z-score_loo-avg-1_prf-css_deriv.func.gii
                    loo_prf_stats_fns = glob.glob(f"{prf_deriv_dir}/*task-{task_name}_*loo-avg-*_{analysis_name}-css{dm_name}_stats.{extension}")
                    
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
                            loo_prf_stats_fn =  f'{prf_deriv_dir}/{subject}_task-{task_name}{hemi}_{preproc_prep}_{filtering}_{normalization}_loo-avg_{analysis_name}-css{dm_name}_stats.{extension}'
                            
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

# template_avg median
elif subject == 'template_avg':
    for averaging_template_name, averaging_template_format in averaging_templates.items(): 
        print('{}, Median corr across subject...'. format(averaging_template_name))
    
        for task_name in task_names:
            
            for avg_method in avg_methods:
                
                # find all the subject prf stats
                prf_stats_fns = []
                for subject in subjects: 
                    prf_deriv_dir = f"{main_dir}/{project_dir}/derivatives/pp_data/{subject}/{averaging_template_format}/{output_folder}/prf_derivatives"
                    prf_stats_fns += [f"{prf_deriv_dir}/{subject}_task-{task_name}_{preproc_prep}_{filtering}_{normalization}_{avg_method}_{analysis_name}-css{dm_name}_stats.dtseries.nii"]
        
                # Computing  across subject
                img, data_stat_median = median_subject_template(fns=prf_stats_fns)
                
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
                template_stats_dir = f"{main_dir}/{project_dir}/derivatives/pp_data/{averaging_template_name}/{averaging_template_format}/{output_folder}/prf_derivatives"
                os.makedirs(template_stats_dir, exist_ok=True)
                
                template_stat_fn = f"{template_stats_dir}/{averaging_template_name}_task-{task_name}_{preproc_prep}_{filtering}_{normalization}_{avg_method}_{analysis_name}-css{dm_name}_stats.dtseries.nii"
                print("saving: {}".format(template_stat_fn))
                template_stat_img = make_surface_image(data=data_stat_median, 
                                                       source_img=img, 
                                                       maps_names=maps_names)
                nb.save(template_stat_img, template_stat_fn)

# Define permission cmd
print('Changing files permissions in {}/{}'.format(main_dir, project_dir))
os.system("chmod -Rf 771 {}/{}".format(main_dir, project_dir))
os.system("chgrp -Rf {} {}/{}".format(group, main_dir, project_dir))