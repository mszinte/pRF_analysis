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
>> cd ~/projects/pRF_analysis/RetinoMaps/glm/postfit
2. run python command
>> python compute_glm_stats.py [main directory] [project name] [subject num] [group]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/projects/pRF_analysis/RetinoMaps/glm/postfit
python compute_glm_stats.py /scratch/mszinte/data RetinoMaps sub-01 327
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
import json
import numpy as np
import nibabel as nb
from scipy import stats

# personal imports
sys.path.append("{}/../../../analysis_code/utils".format(os.getcwd()))
from pycortex_utils import set_pycortex_config_file
from surface_utils import make_surface_image , load_surface
from maths_utils import linear_regression_surf, multipletests_surface

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
group = sys.argv[4]

# load settings
with open('../../settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)
fdr_alpha = analysis_info['fdr_alpha']
formats = analysis_info['formats']
extensions = analysis_info['extensions']
maps_names = analysis_info['maps_names_css_stats']
tasks = analysis_info['task_glm']
TRs = analysis_info['TRs']

# Set pycortex db and colormaps
cortex_dir = "{}/{}/derivatives/pp_data/cortex".format(main_dir, project_dir)
set_pycortex_config_file(cortex_dir)

# Index
slope_idx, intercept_idx, rvalue_idx, pvalue_idx, stderr_idx, \
    trs_idx, corr_pvalue_5pt_idx, corr_pvalue_1pt_idx = 0, 1, 2, 3, 4, 5, 6, 7

for format_, extension in zip(formats, extensions): 
    print(format_)
    for task in tasks:
        # Find glm fit files 
        glm_fit_dir = '{}/{}/derivatives/pp_data/{}/{}/glm/glm_fit'.format(
            main_dir, project_dir, subject, format_)
        glm_bold_dir = '{}/{}/derivatives/pp_data/{}/{}/func/fmriprep_dct_loo_avg'.format(
            main_dir, project_dir, subject, format_)
        glm_pred_loo_fns_list = glob.glob('{}/*task-{}*loo-*_glm-pred.{}'.format(
            glm_fit_dir, task, extension))
        
        for glm_pred_loo_fn in glm_pred_loo_fns_list : 
            # Find the correponding bold signal to the loo prediction
            loo_number = re.search(r'loo-(\d+)', glm_pred_loo_fn).group(1)
            if format_ == 'fsnative': 
                rois = analysis_info['rois']
                hemi = re.search(r'hemi-(\w)', glm_pred_loo_fn).group(1)
                glm_bold_fn = '{}/{}_task-{}_hemi-{}_fmriprep_dct_loo-{}_bold.{}'.format(
                    glm_bold_dir, subject, task, hemi, loo_number, extension)
            elif format_ == '170k':
                rois = analysis_info['mmp_rois']
                glm_bold_fn = '{}/{}_task-{}_fmriprep_dct_loo-{}_bold.{}'.format(
                    glm_bold_dir, subject, task, loo_number, extension)
            
            # load data  
            pred_img, pred_data = load_surface(glm_pred_loo_fn)
            bold_img, bold_data = load_surface(glm_bold_fn)
            
            # Compute linear regression 
            results = linear_regression_surf(
                bold_signal=bold_data, model_prediction=pred_data, correction='fdr_tsbh', alpha=fdr_alpha)

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
    list_ = glob.glob("{}/{}/derivatives/pp_data/{}/{}/glm/glm_derivatives/*loo-*_glm-stats.{}".format(
        main_dir, project_dir, subject, format_, extension))
    list_ = [item for item in list_ if "loo-median" not in item]
    glm_stats_loo_fns_list.extend(list_)
        
# split filtered files  depending of their nature
stats_fsnative_hemi_L, stats_fsnative_hemi_R, stats_170k = [], [], []
for subtype in glm_stats_loo_fns_list:
    if "hemi-L" in subtype: stats_fsnative_hemi_L.append(subtype)
    elif "hemi-R" in subtype: stats_fsnative_hemi_R.append(subtype)
    else : stats_170k.append(subtype)
loo_stats_fns_list = [stats_fsnative_hemi_L, stats_fsnative_hemi_R, stats_170k]
hemi_data_median = {'hemi-L': [], 'hemi-R': [], '170k': []}

# Median
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
            loo_stats_median_fn = re.sub(r'avg_loo-\d+_glm-stats', 'avg_glm-stats_loo-median', loo_stats_median_fn)
            
            # load data 
            print('adding {} to averaging'.format(loo_stats_fn))
            loo_stats_img, loo_stats_data = load_surface(fn=loo_stats_fn)
            
            # Median
            if n_run == 0: loo_stats_data_median = np.copy(loo_stats_data)
            else: loo_stats_data_median = np.nanmedian(np.array([loo_stats_data_median, loo_stats_data]), axis=0)
            
        # Compute two sided corrected p-values
        t_statistic = loo_stats_data_median[slope_idx, :] / loo_stats_data_median[stderr_idx, :]
        degrees_of_freedom = np.nanmax(loo_stats_data_median[trs_idx, :]) - 2
        p_values = 2 * (1 - stats.t.cdf(abs(t_statistic), df=degrees_of_freedom)) 
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
        nb.save(loo_stats_img, median_fn)

# Define permission cmd
print('Changing files permissions in {}/{}'.format(main_dir, project_dir))
os.system("chmod -Rf 771 {main_dir}/{project_dir}".format(main_dir=main_dir, project_dir=project_dir))
os.system("chgrp -Rf {group} {main_dir}/{project_dir}".format(main_dir=main_dir, project_dir=project_dir, group=group))