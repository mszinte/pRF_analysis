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
>> cd ~/projects/[PROJECT]/analysis_code/postproc/prf/postfit/
2. run python command
>> python compute_css_stats.py [main directory] [project name] 
                               [subject num] [group] [analysis folder - optional]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/projects/pRF_analysis/analysis_code/postproc/prf/postfit/

python compute_css_stats.py /scratch/mszinte/data MotConf sub-01 327
python compute_css_stats.py /scratch/mszinte/data MotConf sub-170k 327

python compute_css_stats.py /scratch/mszinte/data RetinoMaps sub-01 327
python compute_css_stats.py /scratch/mszinte/data RetinoMaps sub-170k 327

python compute_css_stats.py /scratch/mszinte/data amblyo_prf sub-01 327
python compute_css_stats.py /scratch/mszinte/data amblyo_prf sub-01 327 prf_em_ctrl
python compute_css_stats.py /scratch/mszinte/data amblyo_prf sub-170k 327
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
TRs = analysis_info['TRs']
maps_names = analysis_info['maps_names_css_stats']
prf_task_name = analysis_info['prf_task_name']
subjects = analysis_info['subjects']

# Set pycortex db and colormaps
cortex_dir = "{}/{}/derivatives/pp_data/cortex".format(main_dir, project_dir)
set_pycortex_config_file(cortex_dir)

# Index
slope_idx, intercept_idx, rvalue_idx, pvalue_idx, stderr_idx, \
    trs_idx, corr_pvalue_5pt_idx, corr_pvalue_1pt_idx = 0, 1, 2, 3, 4, 5, 6, 7

# sub-170k exeption
if subject != 'sub-170k':
    pp_dir = "{}/{}/derivatives/pp_data".format(main_dir, project_dir)
    for format_, extension in zip(formats, extensions): 
        print(format_)
        
        # Find pRF fit files 
        prf_fit_dir = '{}/{}/{}/{}/fit'.format(
            pp_dir, subject, format_, output_folder)
        prf_bold_dir = '{}/{}/{}/func/fmriprep_dct_loo_avg'.format(
             pp_dir, subject, format_)
        prf_pred_loo_fns_list = glob.glob('{}/*task-{}*loo-*_prf-pred_css.{}'.format(
            prf_fit_dir, prf_task_name, extension))
        
        for prf_pred_loo_fn in prf_pred_loo_fns_list : 
            # Find the correponding bold signal to the loo prediction
            loo_number = re.search(r'loo-(\d+)', prf_pred_loo_fn).group(1)
            if format_ == 'fsnative': 
                hemi = re.search(r'hemi-(\w)', prf_pred_loo_fn).group(1)
                prf_bold_fn = '{}/{}_task-{}_hemi-{}_fmriprep_dct_loo-{}_bold.{}'.format(
                    prf_bold_dir, subject, prf_task_name, hemi, loo_number, extension)
            elif format_ == '170k':
                prf_bold_fn = '{}/{}_task-{}_fmriprep_dct_loo-{}_bold.{}'.format(
                    prf_bold_dir, subject, prf_task_name, loo_number, extension)
            
            # load data  
            pred_img, pred_data = load_surface(prf_pred_loo_fn)
            bold_img, bold_data = load_surface(prf_bold_fn)
            
            # Compute linear regression 
            results = linear_regression_surf(bold_signal=bold_data, 
                                             model_prediction=pred_data, 
                                             correction='fdr_tsbh', 
                                             alpha=fdr_alpha)
            
            # Save results
            prf_deriv_dir = "{}/{}/{}/{}/prf_derivatives".format(
                pp_dir, subject, format_, output_folder)
            stat_prf_loo_fn = prf_pred_loo_fn.split('/')[-1].replace('pred_css', 'stats')
            stat_prf_loo_img = make_surface_image(data=results, 
                                                  source_img=bold_img, 
                                                  maps_names=maps_names)
            print('Saving: {}/{}'.format(prf_deriv_dir, stat_prf_loo_fn))
            nb.save(stat_prf_loo_img, '{}/{}'.format(prf_deriv_dir, stat_prf_loo_fn))
            
    # Compute median across LOO
    print('Compute median across LOO')
    
    # Get files
    prf_stats_loo_fns_list = []
    for format_, extension in zip(formats, extensions):
        list_ = glob.glob("{}/{}/{}/{}/prf_derivatives/*loo-*_prf-stats.{}".format(
            pp_dir, subject, format_, output_folder, extension))
        list_ = [item for item in list_ if "loo-median" not in item]
        prf_stats_loo_fns_list.extend(list_)
                
    # Split files depending of their nature
    stats_fsnative_hemi_L, stats_fsnative_hemi_R, stats_170k = [], [], []
    for subtype in prf_stats_loo_fns_list:
        if "hemi-L" in subtype: stats_fsnative_hemi_L.append(subtype)
        elif "hemi-R" in subtype: stats_fsnative_hemi_R.append(subtype)
        else : stats_170k.append(subtype)
    loo_stats_fns_list = [stats_fsnative_hemi_L, stats_fsnative_hemi_R, stats_170k]
    hemi_data_median = {'hemi-L': [], 'hemi-R': [], '170k': []}
    
    # Computing median
    for loo_stats_fns in loo_stats_fns_list:
        if not loo_stats_fns:
            continue
        if loo_stats_fns[0].find('hemi-L') != -1:  hemi = 'hemi-L'
        elif loo_stats_fns[0].find('hemi-R') != -1: hemi = 'hemi-R'
        else: hemi = None
    
        # Computing median
        stats_img, stats_data = load_surface(fn=loo_stats_fns[0])
        loo_stats_data_median = np.zeros(stats_data.shape)
        
        for n_run, loo_stats_fn in enumerate(loo_stats_fns):
            loo_stats_median_fn = loo_stats_fn.split('/')[-1]
            loo_stats_median_fn = re.sub(r'avg_loo-\d+_prf-stats', 'avg_prf-stats_loo-median', loo_stats_median_fn)
    
            # Load data 
            print('adding {} to computing median'.format(loo_stats_fn))
            loo_stats_img, loo_stats_data = load_surface(fn=loo_stats_fn)
    
            # median
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
            median_fn = '{}/{}/fsnative/{}/prf_derivatives/{}'.format(
                pp_dir, subject, output_folder, loo_stats_median_fn)
            hemi_data_median[hemi] = loo_stats_data_median
        else:
            median_fn = '{}/{}/170k/{}/prf_derivatives/{}'.format(
                pp_dir, subject, output_folder, loo_stats_median_fn)
            hemi_data_median['170k'] = loo_stats_data_median
    
        # Saving data in surface format
        loo_stats_img = make_surface_image(data=loo_stats_data_median, 
                                           source_img=loo_stats_img, 
                                           maps_names=maps_names)
        print('Saving median: {}'.format(median_fn))
        nb.save(loo_stats_img, median_fn)
        
# Sub-170k median
elif subject == 'sub-170k':
    print('sub-170, computing median prf stats across subject...')
    
    # find all the subject prf derivatives
    subjects_stats = []
    for subject in subjects: 
        subjects_stats += ["{}/{}/derivatives/pp_data/{}/170k/{}/prf_derivatives/{}_task-{}_fmriprep_dct_avg_prf-stats_loo-median.dtseries.nii".format(
                main_dir, project_dir, subject, output_folder, subject, prf_task_name)]

    # Computing median across subject
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
    
    sub_170k_stat_fn = "{}/sub-170k_task-{}_fmriprep_dct_avg_prf-stats_loo-median.dtseries.nii".format(sub_170k_stats_dir, prf_task_name)
    print("save: {}".format(sub_170k_stat_fn))
    sub_170k_stat_img = make_surface_image(data=data_stat_median, 
                                           source_img=img, 
                                           maps_names=maps_names)
    nb.save(sub_170k_stat_img, sub_170k_stat_fn)

# Define permission cmd
print('Changing files permissions in {}/{}'.format(main_dir, project_dir))
os.system("chmod -Rf 771 {}/{}".format(main_dir, project_dir))
os.system("chgrp -Rf {} {}/{}".format(group, main_dir, project_dir))