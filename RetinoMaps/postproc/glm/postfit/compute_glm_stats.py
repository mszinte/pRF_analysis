"""
-----------------------------------------------------------------------------------------
compute_glm_stats.py
-----------------------------------------------------------------------------------------
Goal of the script:
Compute the linear regression between the GLM prediction and the bold signal
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
>> cd ~/projects/pRF_analysis/RetinoMaps/postproc/glm/postfit
2. run python command
>> python compute_glm_stats.py [main directory] [project name] [subject num] [group]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/projects/pRF_analysis/RetinoMaps/postproc/glm/postfit
python compute_glm_stats.py /scratch/mszinte/data RetinoMaps sub-01 327
python compute_glm_stats.py /scratch/mszinte/data RetinoMaps template_avg 327
-----------------------------------------------------------------------------------------
Written by Uriel Lascombes (uriel.lascombes@laposte.net)
Edited by Martin Szinte (martin.szinte@gmail.com) 
-----------------------------------------------------------------------------------------
"""
# Stop warnings
import warnings
warnings.filterwarnings("ignore")

# debug 
import ipdb
deb = ipdb.set_trace

# General imports
import os
import re
import sys
import glob
import numpy as np
import nibabel as nb
# [CORRECTION #20] pandas necessaire pour lire la matrice de design (signe du beta de la tache)
# AVANT : (ligne ajoutee)
import pandas as pd
from scipy import stats

# personal imports
sys.path.append("{}/../../../../analysis_code/utils".format(os.getcwd()))
from settings_utils import load_settings
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
settings_path = os.path.join(base_dir, project_dir, "settings.yml")
prf_settings_path = os.path.join(base_dir, project_dir, "prf-analysis.yml")
glm_settings_path = os.path.join(base_dir, project_dir, "glm-analysis.yml")
settings = load_settings([settings_path, prf_settings_path, glm_settings_path])
analysis_info = settings[0]
    
TRs = analysis_info['TRs']
tasks = analysis_info['task_glm']
formats = analysis_info['formats']
subjects = analysis_info['subjects']
fdr_alpha = analysis_info['fdr_alpha']
extensions = analysis_info['extensions']
maps_names = analysis_info['maps_names_css_stats']
averaging_templates = analysis_info['averaging_templates']
pycortex_subject_template = analysis_info['pycortex_subject_template']
# [NORMALIZATION] preprocessing labels loaded (were hard-coded in the file names)
# AVANT : (lignes ajoutees)
preproc_prep = analysis_info['preproc_prep']
filtering = analysis_info['filtering']
normalizations = analysis_info['normalization']

# Set pycortex db and colormaps
cortex_dir = "{}/{}/derivatives/pp_data/cortex".format(main_dir, project_dir)
set_pycortex_config_file(cortex_dir)

# Index
slope_idx, intercept_idx, rvalue_idx, pvalue_idx, stderr_idx, \
    trs_idx, corr_pvalue_5pt_idx, corr_pvalue_1pt_idx = 0, 1, 2, 3, 4, 5, 6, 7

# template_avg exeption
if subject != 'template_avg':
    # [NORMALIZATION] loop over normalizations (inside the if/elif: the template branch reassigns `subject`)
    for normalization in normalizations:
        for format_, extension in zip(formats, extensions): 
            print(format_)
            for task in tasks:
                # Find glm fit files 
                glm_fit_dir = '{}/{}/derivatives/pp_data/{}/{}/glm/fit'.format(
                    main_dir, project_dir, subject, format_)
                # [NORMALIZATION] hard-coded 'fmriprep_dct_z-score' replaced by the settings labels
                # AVANT : glm_bold_dir = '{}/{}/derivatives/pp_data/{}/{}/func/fmriprep_dct_z-score_loo-avg'.format(
                glm_bold_dir = '{}/{}/derivatives/pp_data/{}/{}/func/{preproc_prep}_{filtering}_{normalization}_loo-avg'.format(
                    main_dir, project_dir, subject, format_, preproc_prep=preproc_prep, filtering=filtering, normalization=normalization)
                # [NORMALIZATION] glob restricted to the current normalization
                # AVANT : glm_pred_loo_fns_list = glob.glob('{}/*task-{}*loo-avg*_glm-pred.{}'.format(
                # AVANT :     glm_fit_dir, task, extension))
                glm_pred_loo_fns_list = glob.glob('{}/*task-{}*_{}_{}_{}_loo-avg*_glm-pred.{}'.format(
                    glm_fit_dir, task, preproc_prep, filtering, normalization, extension))

                for glm_pred_loo_fn in glm_pred_loo_fns_list : 
                    # Find the correponding bold signal to the loo prediction
                    loo_number = re.search(r'loo-avg-(\d+)', glm_pred_loo_fn).group(1)
                    if format_ == 'fsnative': 
                        hemi = re.search(r'hemi-(\w)', glm_pred_loo_fn).group(1)
                        glm_bold_fn = glob.glob('{}/*task-{}*_hemi-{}*_loo-{}_bold.{}'.format(
                            glm_bold_dir, task, hemi, loo_number, extension))[0]
                    elif format_ == '170k':
                        glm_bold_fn = glob.glob('{}/*task-{}*_loo-{}_bold.{}'.format(
                            glm_bold_dir, task, loo_number, extension))[0]

                    # load data  
                    print(f'Loading pred: {glm_pred_loo_fn}') 
                    pred_img, pred_data = load_surface(glm_pred_loo_fn)
                    print(f'Loading bold: {glm_bold_fn}')
                    bold_img, bold_data = load_surface(glm_bold_fn)
                
                    # [CORRECTION #20] on multiplie la prediction par le signe du beta de la tache : un vertex desactive (beta < 0) donne alors une pente
                    # [CORRECTION #20] negative et n'est plus compte comme 'actif' avec le test unilateral ci-dessous
                    # AVANT : (bloc ajoute)
                    # Sign of the task beta (prediction = beta_task * task regressor + constant)
                    dm_fn = '{}/{}/derivatives/pp_data/{}/{}/glm/designe_matrix/{}_task-{}_designe_matrix.tsv'.format(
                        main_dir, project_dir, subject, format_, subject, task)
                    # [CORRECTION #40] verification : la matrice de design doit contenir exactement un regresseur de tache + 'constant',
                    # [CORRECTION #40] sinon iloc[:, 0] prendrait un regresseur arbitraire pour le signe du beta
                    # AVANT : task_reg = pd.read_table(dm_fn).drop(columns=['constant']).iloc[:, 0].values
                    task_dm = pd.read_table(dm_fn)
                    if task_dm.shape[1] != 2 or 'constant' not in task_dm.columns:
                        raise ValueError('{} : expected 1 task regressor + constant, found columns {}'.format(
                            dm_fn, list(task_dm.columns)))
                    task_reg = task_dm.drop(columns=['constant']).iloc[:, 0].values
                    beta_sign = np.sign(np.dot(task_reg - task_reg.mean(), pred_data - np.nanmean(pred_data, axis=0)))
                    pred_data = pred_data * beta_sign
                
                    # Compute linear regression 
                    results = linear_regression_surf(bold_signal=bold_data, 
                                                     model_prediction=pred_data, 
                                                     # [CORRECTION #20] test unilateral : pente > 0 (vertex active par la tache et prediction repliquee sur le run laisse de cote)
                                                     # AVANT : alternative='two-sided',
                                                     alternative='greater',
                                                     correction='fdr_tsbh', 
                                                     alpha=fdr_alpha)
    
                    # Save results 
                    glm_deriv_dir = '{}/{}/derivatives/pp_data/{}/{}/glm/glm_derivatives'.format(
                        main_dir, project_dir, subject, format_)
                    os.makedirs(glm_deriv_dir, exist_ok=True)
                    stat_glm_loo_fn = glm_pred_loo_fn.split('/')[-1].replace('pred', 'stats')
                    stat_glm_loo_img = make_surface_image(data=results, 
                                                          source_img=bold_img, 
                                                          maps_names=maps_names)
                    print('Saving: {}/{}'.format(glm_deriv_dir, stat_glm_loo_fn))
                    nb.save(stat_glm_loo_img, '{}/{}'.format(glm_deriv_dir, stat_glm_loo_fn))
          
        # Compute median across LOO
        # Get files
        glm_stats_loo_fns_list = []
        for format_, extension in zip(formats, extensions):
            # [NORMALIZATION] glob restricted to the current normalization
            # AVANT : list_ = glob.glob("{}/{}/derivatives/pp_data/{}/{}/glm/glm_derivatives/*loo-*_glm-stats.{}".format(
            # AVANT :     main_dir, project_dir, subject, format_, extension))
            list_ = glob.glob("{}/{}/derivatives/pp_data/{}/{}/glm/glm_derivatives/*_{}_{}_{}_loo-*_glm-stats.{}".format(
                main_dir, project_dir, subject, format_, preproc_prep, filtering, normalization, extension))
            # [CORRECTION #19] on ne garde que les plis loo-avg-N : le fichier median loo-avg deja ecrit etait repris a chaque relance
            # AVANT : list_ = [item for item in list_ if "loo-median" not in item]
            list_ = [item for item in list_ if re.search(r'loo-avg-\d+_glm-stats', item)]
            glm_stats_loo_fns_list.extend(list_)
            
        # split filtered files  depending of their nature
        stats_fsnative_hemi_L, stats_fsnative_hemi_R, stats_170k = [], [], []
        for subtype in glm_stats_loo_fns_list:
            if "hemi-L" in subtype: stats_fsnative_hemi_L.append(subtype)
            elif "hemi-R" in subtype: stats_fsnative_hemi_R.append(subtype)
            else : stats_170k.append(subtype)
        loo_stats_fns_list = [stats_fsnative_hemi_L, stats_fsnative_hemi_R, stats_170k]
        hemi_data_median = {'hemi-L': [], 'hemi-R': [], '170k': []}
    
        # Compute median across leave-one-out fit
        for loo_stats_fns in loo_stats_fns_list:
            for task in tasks:
                # defind output files names 
                loo_stats_fns_task = [file for file in loo_stats_fns if task in file]
            
                # Check if task files exist
                if not loo_stats_fns_task:
                    print('No files for {}'.format(task))
                    continue
            
                if loo_stats_fns_task[0].find('hemi-L') != -1: hemi = 'hemi-L'
                elif loo_stats_fns_task[0].find('hemi-R') != -1: hemi = 'hemi-R'
                else: hemi = None
        
                # Median 
                stats_img, stats_data = load_surface(fn=loo_stats_fns[0])
                loo_stats_data_median = np.zeros(stats_data.shape)
            
                for n_run, loo_stats_fn in enumerate(loo_stats_fns_task):
                    loo_stats_median_fn = loo_stats_fn.split('/')[-1]
                    loo_stats_median_fn = re.sub(r'loo-avg-\d+', 'loo-avg', loo_stats_median_fn)
  
                    # load data 
                    print('adding {} to averaging'.format(loo_stats_fn))
                    loo_stats_img, loo_stats_data = load_surface(fn=loo_stats_fn)
                
                    # Median
                    if n_run == 0: loo_stats_data_median = np.copy(loo_stats_data)
                    else: loo_stats_data_median = np.nanmedian(np.array([loo_stats_data_median, loo_stats_data]), axis=0)
                
                # [CORRECTION #20] commentaire mis a jour : p-values unilaterales (t.sf), plus bilaterales
                # AVANT : # Compute two sided corrected p-values
                # Compute one sided (slope > 0) corrected p-values
                t_statistic = loo_stats_data_median[slope_idx, :] / loo_stats_data_median[stderr_idx, :]
                degrees_of_freedom = np.nanmax(loo_stats_data_median[trs_idx, :]) - 2
                # [CORRECTION #20] p unilateral (pente > 0), coherent avec alternative='greater'
                # AVANT : p_values = 2 * (1 - stats.t.cdf(abs(t_statistic), df=degrees_of_freedom))
                p_values = stats.t.sf(t_statistic, df=degrees_of_freedom)
                corrected_p_values = multipletests_surface(pvals=p_values, 
                                                           correction='fdr_tsbh', 
                                                           alpha=fdr_alpha)
                loo_stats_data_median[pvalue_idx, :] = p_values
                loo_stats_data_median[corr_pvalue_5pt_idx, :] = corrected_p_values[0,:]
                loo_stats_data_median[corr_pvalue_1pt_idx, :] = corrected_p_values[1,:]
                
                if hemi:
                    median_fn = '{}/{}/derivatives/pp_data/{}/fsnative/glm/glm_derivatives/{}'.format(
                        main_dir, project_dir, subject, loo_stats_median_fn)
                    hemi_data_median[hemi] = loo_stats_data_median
                else:
                    median_fn = '{}/{}/derivatives/pp_data/{}/170k/glm/glm_derivatives/{}'.format(
                        main_dir, project_dir, subject, loo_stats_median_fn)
                    hemi_data_median['170k'] = loo_stats_data_median
            
                # export averaged data in surface format 
                loo_stats_img = make_surface_image(data=loo_stats_data_median, 
                                                   source_img=loo_stats_img, 
                                                   maps_names=maps_names)
                print('Saving {}'.format(median_fn))
                nb.save(loo_stats_img, median_fn)
            
# template_avg median
elif subject == 'template_avg':
    # [NORMALIZATION] loop over normalizations (inside the if/elif: the template branch reassigns `subject`)
    for normalization in normalizations:
        for averaging_template_name, averaging_template_format in averaging_templates.items(): 
            print('{}, Median corr across subject...'. format(averaging_template_name))
    
            for task in tasks:
                # find all the subject prf derivatives
                subjects_stats_task = []
                for subject in subjects:
                    # [NORMALIZATION] hard-coded 'fmriprep_dct_z-score' replaced by the settings labels
                    # AVANT : subjects_stats_task += ["{}/{}/derivatives/pp_data/{}/170k/glm/glm_derivatives/{}_task-{}_fmriprep_dct_z-score_loo-avg_glm-stats.dtseries.nii".format(
                    subjects_stats_task += ["{}/{}/derivatives/pp_data/{}/170k/glm/glm_derivatives/{}_task-{}_{preproc_prep}_{filtering}_{normalization}_loo-avg_glm-stats.dtseries.nii".format(
                            main_dir, project_dir, subject, subject, task, preproc_prep=preproc_prep, filtering=filtering, normalization=normalization)]
        
                # Computing median across subject
                img, data_stat_median = median_subject_template(fns=subjects_stats_task)
            
                # [CORRECTION #20] commentaire mis a jour : p-values unilaterales (t.sf), plus bilaterales
                # AVANT : # Compute two sided corrected p-values
                # Compute one sided (slope > 0) corrected p-values
                t_statistic = data_stat_median[slope_idx, :] / data_stat_median[stderr_idx, :]
                degrees_of_freedom = np.nanmax(data_stat_median[trs_idx, :]) - 2
                # [CORRECTION #20] p unilateral (pente > 0), coherent avec alternative='greater'
                # AVANT : p_values = 2 * (1 - stats.t.cdf(abs(t_statistic), df=degrees_of_freedom))
                p_values = stats.t.sf(t_statistic, df=degrees_of_freedom)
                corrected_p_values = multipletests_surface(pvals=p_values, 
                                                           correction='fdr_tsbh', 
                                                           alpha=fdr_alpha)
                data_stat_median[pvalue_idx, :] = p_values
                data_stat_median[corr_pvalue_5pt_idx, :] = corrected_p_values[0,:]
                data_stat_median[corr_pvalue_1pt_idx, :] = corrected_p_values[1,:]
                
                # Export results
                template_stats_dir = "{}/{}/derivatives/pp_data/{}/170k/glm/glm_derivatives/".format(
                        main_dir, project_dir, pycortex_subject_template)
                os.makedirs(template_stats_dir, exist_ok=True)
    
                # [NORMALIZATION] hard-coded 'fmriprep_dct_z-score' replaced by the settings labels
                # AVANT : template_stat_fn = "{}/{}_task-{}_fmriprep_dct_z-score_loo-avg_glm-stats.dtseries.nii".format(template_stats_dir, pycortex_subject_template, task)
                template_stat_fn = "{}/{}_task-{}_{preproc_prep}_{filtering}_{normalization}_loo-avg_glm-stats.dtseries.nii".format(template_stats_dir, pycortex_subject_template, task, preproc_prep=preproc_prep, filtering=filtering, normalization=normalization)
                print("save: {}".format(template_stat_fn))
                template_stat_img = make_surface_image(data=data_stat_median, 
                                                       source_img=img, 
                                                       maps_names=maps_names)
                nb.save(template_stat_img, template_stat_fn)

# # Define permission cmd
# print('Changing files permissions in {}/{}'.format(main_dir, project_dir))
# os.system("chmod -Rf 771 {main_dir}/{project_dir}".format(main_dir=main_dir, project_dir=project_dir))
# os.system("chgrp -Rf {group} {main_dir}/{project_dir}".format(main_dir=main_dir, project_dir=project_dir, group=group))