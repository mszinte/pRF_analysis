"""
-----------------------------------------------------------------------------------------
compute_ncsf_stats.py
-----------------------------------------------------------------------------------------
Goal of the script:
Compute the linear regression between the nCSF  predictions and the bold signal
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
>> cd ~/projects/pRF_analysis/nCSF/postproc/nCSF/postfit/
2. run python command
>> python compute_ncsf_stats.py [main directory] [project name] 
                               [subject num] [group]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/projects/pRF_analysis/nCSF/postproc/nCSF/postfit/
python compute_ncsf_stats.py /scratch/mszinte/data nCSF sub-01 327
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
sys.path.append("{}/../../../../analysis_code/utils".format(os.getcwd()))
from pycortex_utils import set_pycortex_config_file
from surface_utils import make_surface_image , load_surface
from maths_utils import linear_regression_surf, multipletests_surface, median_subject_template
from settings_utils import load_settings

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
group = sys.argv[4]

# Load settings
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../../"))
settings_path = os.path.join(base_dir, project_dir, "settings.yml")
nCSF_settings_path = os.path.join(base_dir, project_dir, "nCSF-analysis.yml")
prf_settings_path = os.path.join(base_dir, project_dir, "prf-analysis.yml")
settings = load_settings([settings_path, nCSF_settings_path, prf_settings_path])
analysis_info = settings[0]

fdr_alpha = analysis_info['fdr_alpha']
formats = analysis_info['formats']
extensions = analysis_info['extensions']
TRs = analysis_info['TRs']
maps_names = analysis_info['maps_names_ncsf_stats']
ncsf_task_name = analysis_info['nCSF_task_name']
subjects = analysis_info['subjects']
preproc_prep = analysis_info['preproc_prep']
filtering = analysis_info['filtering']
normalization = analysis_info['normalization']
avg_methods = analysis_info['avg_methods']
averaging_templates = analysis_info['averaging_templates']

# Set pycortex db and colormaps
cortex_dir = "{}/{}/derivatives/pp_data/cortex".format(main_dir, project_dir)
set_pycortex_config_file(cortex_dir)

# Index
slope_idx, intercept_idx, rvalue_idx, pvalue_idx, stderr_idx, \
    trs_idx, corr_pvalue_5pt_idx, corr_pvalue_1pt_idx = 0, 1, 2, 3, 4, 5, 6, 7

# Define preprocessing folder
pp_dir = "{}/{}/derivatives/pp_data".format(main_dir, project_dir)


# template_avg exeption
if subject != 'template_avg':
    for avg_method in avg_methods:
        for format_, extension in zip(formats, extensions): 
            
            # define/create folders
            ncsf_fit_dir = '{}/{}/{}/ncsf/fit'.format(
                    pp_dir, subject, format_)
            ncsf_func_dir = '{}/{}/{}/func/{}_{}_{}_{}'.format(
                    pp_dir, subject, format_, 
                    preproc_prep, filtering, normalization, avg_method)
            ncsf_deriv_dir = "{}/{}/{}/ncsf/ncsf_derivatives".format(
                pp_dir, subject, format_)
            os.makedirs(ncsf_deriv_dir, exist_ok=True)
            
            # Find ncsf func/pred files
            ncsf_pred_fns = glob.glob('{}/*task-{}_*_{}*_ncsf_pred.{}'.format(
                ncsf_fit_dir, ncsf_task_name, avg_method, extension))
            for ncsf_pred_fn in ncsf_pred_fns :
                if 'loo' in ncsf_pred_fn:
                    loo_number = re.search(r'loo-avg-(\d+)', ncsf_pred_fn).group(1)
                    if format_ == 'fsnative': 
                        hemi = re.search(r'hemi-(\w)', ncsf_pred_fn).group(1)
                        ncsf_bold_fn = glob.glob('{}/*task-{}_hemi-{}*_loo-{}_bold*.{}'.format(
                            ncsf_func_dir, ncsf_task_name, hemi, loo_number, extension))[0]
                    elif format_ == '170k':
                        ncsf_bold_fn = glob.glob('{}/*task-{}_*_loo-{}_bold.{}'.format(
                            ncsf_func_dir, ncsf_task_name, loo_number, extension))[0]
                else:
                    if format_ == 'fsnative': 
                        hemi = re.search(r'hemi-(\w)', ncsf_pred_fn).group(1)
                        ncsf_bold_fn = glob.glob('{}/*task-{}_hemi-{}*_{}_bold.{}'.format(
                            ncsf_func_dir, ncsf_task_name, hemi, avg_method, extension))[0]
                    elif format_ == '170k':
                        ncsf_bold_fn = glob.glob('{}/*task-{}_*_{}_bold*.{}'.format(
                            ncsf_func_dir, ncsf_task_name, avg_method, extension))[0]

                # load data  
                print(f'Loading pred: {ncsf_pred_fn}') 
                bold_img, bold_data = load_surface(ncsf_bold_fn)
                print(f'Loading bold: {ncsf_bold_fn}')
                pred_img, pred_data = load_surface(ncsf_pred_fn)
                
                # Compute linear regression 
                results = linear_regression_surf(bold_signal=bold_data, 
                                                 model_prediction=pred_data, 
                                                 alternative='two-sided',
                                                 correction='fdr_tsbh', 
                                                 alpha=fdr_alpha)
                
                # Save results
                ncsf_stats_fn = ncsf_pred_fn.split('/')[-1].replace('ncsf_pred', 'ncsf_stats')
                ncsf_stats_img = make_surface_image(data=results, 
                                                   source_img=pred_img, 
                                                   maps_names=maps_names)

                print('Saving: {}/{}'.format(ncsf_deriv_dir, ncsf_stats_fn))
                nb.save(ncsf_stats_img, '{}/{}'.format(ncsf_deriv_dir, ncsf_stats_fn))
            
                # Compute median across leave-one-out fit
                if 'loo-avg' in avg_method:
                    print('Computing median across LOO')

                    # Get LOO files (excluding any with "median" in the name)
                    loo_ncsf_stats_fns = glob.glob(f"{ncsf_deriv_dir}/*task-{ncsf_task_name}_*loo-avg-*_ncsf_stats.{extension}")

                    # Group files by hemisphere/format
                    loo_ncsf_stats_fsnative_hemi_L_fns = [fn for fn in loo_ncsf_stats_fns if "hemi-L" in fn]
                    loo_ncsf_stats_fsnative_hemi_R_fns = [fn for fn in loo_ncsf_stats_fns if "hemi-R" in fn]
                    loo_ncsf_stats_170k_fns = [fn for fn in loo_ncsf_stats_fns if "hemi-L" not in fn and "hemi-R" not in fn]
                    
                    # Process each group
                    for group_files, hemi in [(loo_ncsf_stats_fsnative_hemi_L_fns, "_hemi-L"),
                                              (loo_ncsf_stats_fsnative_hemi_R_fns, "_hemi-R"),
                                              (loo_ncsf_stats_170k_fns, "")]:
                        
                        if len(group_files)>0:
    
                            # Load first file to initialize median array and define fn
                            stats_img, stats_data = load_surface(group_files[0])
                            loo_ncsf_stats = np.zeros_like(stats_data)
                            loo_ncsf_stats_fn =  '{}/{}_task-{}{}_{}_{}_{}_loo-avg_ncsf_stats.{}'.format(
                                ncsf_deriv_dir, subject, ncsf_task_name, hemi, 
                                preproc_prep, filtering, normalization, extension)
                        
                            # Compute median across LOO runs
                            for n_run, loo_stats_fn in enumerate(group_files):
                                print(f'Loadding loo stats: {loo_stats_fn}')
                                _, loo_ncsf_stats_run_data = load_surface(loo_stats_fn)
                                if n_run == 0: loo_ncsf_stats = np.copy(loo_ncsf_stats_run_data)
                                else: loo_ncsf_stats = np.nanmedian(np.array([loo_ncsf_stats, loo_ncsf_stats_run_data]), axis=0)
                        
                            # Recalculate p-values for median data
                            t_statistic = loo_ncsf_stats[slope_idx, :] / loo_ncsf_stats[stderr_idx, :]
                            degrees_of_freedom = np.nanmax(loo_ncsf_stats[trs_idx, :]) - 2
                            p_values = 2 * (1 - stats.t.cdf(np.abs(t_statistic), df=degrees_of_freedom))
                            corrected_p_values = multipletests_surface(p_values, correction="fdr_tsbh", alpha=fdr_alpha)
                            
                            # Update median data with recalculated p-values
                            loo_ncsf_stats[pvalue_idx, :] = p_values
                            loo_ncsf_stats[corr_pvalue_5pt_idx, :] = corrected_p_values[0, :]
                            loo_ncsf_stats[corr_pvalue_1pt_idx, :] = corrected_p_values[1, :]
    
                            # Save median results
                            loo_ncsf_stats_img = make_surface_image(loo_ncsf_stats, stats_img, maps_names)
                            nb.save(loo_ncsf_stats_img, loo_ncsf_stats_fn)
                            print(f"Saving median: {loo_ncsf_stats_fn}")

# template_avg median
elif subject == 'template_avg':
    for averaging_template_name, averaging_template_format in averaging_templates.items(): 
        print('{}, Median corr across subject...'. format(averaging_template_name))
                
        for avg_method in avg_methods:
            
            # find all the subject ncsf stats
            ncsf_stats_fns = []
            for subject in subjects: 
                ncsf_deriv_dir = "{}/{}/derivatives/pp_data/{}/{}/ncsf/ncsf_derivatives".format(
                    main_dir, project_dir, subject, averaging_template_format)
                ncsf_stats_fns += ["{}/{}_task-{}_{}_{}_{}_{}_ncsf_stats.dtseries.nii".format(
                    ncsf_deriv_dir, subject, ncsf_task_name,
                    preproc_prep, filtering, normalization, avg_method)]
    
            # Computing  across subject
            img, data_stat_median = median_subject_template(fns=ncsf_stats_fns)
            
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
            template_stats_dir = "{}/{}/derivatives/pp_data/{}/{}/ncsf/ncsf_derivatives".format(
                    main_dir, project_dir, averaging_template_name, averaging_template_format)
            os.makedirs(template_stats_dir, exist_ok=True)
            
            template_stat_fn = "{}/{}_task-{}_{}_{}_{}_{}_ncsf_stats.dtseries.nii".format(
                template_stats_dir, averaging_template_name, ncsf_task_name, 
                preproc_prep, filtering, normalization, avg_method)
            print("saving: {}".format(template_stat_fn))
            template_stat_img = make_surface_image(data=data_stat_median, 
                                                   source_img=img, 
                                                   maps_names=maps_names)
            nb.save(template_stat_img, template_stat_fn)

# Define permission cmd
print('Changing files permissions in {}/{}'.format(main_dir, project_dir))
os.system("chmod -Rf 771 {}/{}".format(main_dir, project_dir))
os.system("chgrp -Rf {} {}/{}".format(group, main_dir, project_dir))