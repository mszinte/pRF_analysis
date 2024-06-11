"""
-----------------------------------------------------------------------------------------
170k_stats_averaging.py
-----------------------------------------------------------------------------------------
Goal of the script:
Average all the subject of the studie on the 170k space.
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: group (e.g. 327)
sys.argv[5]: model (e.g. css, gauss)
-----------------------------------------------------------------------------------------
Output(s):
sh file for running batch command
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/RetinoMaps/analysis_code/postproc/stats/
2. run python command
>> python 170gaussgridfit_averaging.py [main directory] [project name] [group] [model]
    [server project]
-----------------------------------------------------------------------------------------
Exemple:
python 170k_stats_averaging.py /scratch/mszinte/data RetinoMaps 327 
-----------------------------------------------------------------------------------------
Written by Martin Szinte (martin.szinte@gmail.com)
Edited by Uriel Lascombes (uriel.lascombes@laposte.net)
-----------------------------------------------------------------------------------------
"""
# stop warnings
import warnings
warnings.filterwarnings("ignore")

# general imports
import os
import sys
import json
import ipdb
import glob
import numpy as np
import pandas as pd
import nibabel as nb
from scipy import stats
deb = ipdb.set_trace

# Personal imports
sys.path.append("{}/../../utils".format(os.getcwd()))
from pycortex_utils import get_rois
from surface_utils import load_surface , make_surface_image
from maths_utils import multipletests_surface

# inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
group = sys.argv[3]

with open('../../settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)
subjects = analysis_info['subjects']
tasks = analysis_info['task_glm'] + ['pRF']
fdr_alpha = analysis_info['fdr_alpha'][1]
TRs = analysis_info['TRs']

fdr_p_map_idx = 6


slope_idx, intercept_idx, rvalue_idx, pvalue_idx, stderr_idx  = 0,1,2,3,4

for n_task, task in enumerate(tasks) :
    print(task) 
    if task == 'pRF':
        analysis = 'prf'
        avg_170k_stats_dir = '{}/{}/derivatives/pp_data/sub-170k/170k/prf/stats'.format(main_dir, project_dir)
        os.makedirs(avg_170k_stats_dir, exist_ok=True)  
    else:
        analysis = 'glm'
        avg_170k_stats_dir = '{}/{}/derivatives/pp_data/sub-170k/170k/glm/stats'.format(main_dir, project_dir)
        os.makedirs(avg_170k_stats_dir, exist_ok=True)  
    for n_subject, subject in enumerate(subjects) :
        print('adding {} to averaging'.format(subject)) 
        
        if task == 'pRF':
            # define directory and file name
            stats_dir = '{}/{}/derivatives/pp_data/{}/170k/prf/stats'.format(main_dir, project_dir, subject)
            stats_fn = '{}_task-{}_fmriprep_dct_loo-avg_{}-stats.dtseries.nii'.format(subject, task,analysis)
        
        else:
            # define directory and file name
            stats_dir = '{}/{}/derivatives/pp_data/{}/170k/glm/stats'.format(main_dir, project_dir, subject)
            stats_fn = '{}_task-{}_fmriprep_dct_loo-avg_{}-stats.dtseries.nii'.format(subject, task,analysis)


        # Load data
        img, data = load_surface(fn='{}/{}'.format(stats_dir, stats_fn))
    
        # Average without considering nan 
        if n_subject == 0:
            data_avg = np.copy(data)
        else:
            data_avg = np.nanmean(np.array([data_avg, data]), axis=0)
            
    # Compute p-values en base om t-satistic and fdr-corrected p-values for averaged loo runs 
    t_statistic = data_avg[slope_idx, :] / data_avg[stderr_idx, :]
    
    # compute two sided p-values
    degrees_of_freedom = TRs - 2 
    p_values = 2 * (1 - stats.t.cdf(abs(t_statistic), df=degrees_of_freedom)) 
    
    corrected_p_values = multipletests_surface(pvals=p_values, correction='fdr_tsbh', alpha=fdr_alpha)
    slope_idx, intercept_idx, rvalue_idx, pvalue_idx, stderr_idx 
    
    
    loo_stats_data_avg = np.vstack((data_avg[slope_idx,:], 
                                    data_avg[intercept_idx,:], 
                                    data_avg[rvalue_idx,:], 
                                    p_values, 
                                    data_avg[stderr_idx,:], 
                                    corrected_p_values))
            
    #  export results 
    avg_170k_stats_fn = 'sub-170k_task-{}_fmriprep_dct_loo-avg_{}-stats.dtseries.nii'.format(task, analysis)


    print('saving {}/{}'.format(avg_170k_stats_dir, avg_170k_stats_fn))

    avg_img = make_surface_image(data=data_avg, source_img=img)
    nb.save(avg_img,'{}/{}'.format(avg_170k_stats_dir, avg_170k_stats_fn))



#  make the final map for sub 170k 
# find all the stats files 
glm_stats_fns = glob.glob("{}/{}/derivatives/pp_data/sub-170k/170k/glm/stats/*loo-avg*stats.dtseries.nii".format(main_dir, project_dir))
prf_stats_fns = glob.glob("{}/{}/derivatives/pp_data/sub-170k/170k/prf/stats/*loo-avg*stats.dtseries.nii".format(main_dir, project_dir))
    
stats_fns = glm_stats_fns + prf_stats_fns


img, data = load_surface(fn=stats_fns[0]) 
final_map = np.zeros((8,data[fdr_p_map_idx].shape[0]))
for stats_file in stats_fns:

    

        


    # make a final map with both the tasks 
    if 'PurLoc' in stats_file:
        task_idx = 1

    elif 'SacLoc' in stats_file:
        task_idx = 2
  
    elif 'pRF' in stats_file:
        task_idx = 4
    

        
    # load data 
    stats_img_task, stats_data_task = load_surface(fn=stats_file)

    fdr_p_map = stats_data_task[fdr_p_map_idx, :]
    
    for vert, fdr_value in enumerate(fdr_p_map):
        if fdr_value < fdr_alpha:
            final_map[task_idx,vert] += task_idx
                    
final_map[0,:] = np.sum(final_map, axis=0)                
#  Make specifique maps 
for vert, final_value in enumerate(final_map[0,:]):

    if final_value == 3 :
        final_map[3,vert] = 3
           
    elif final_value == 5 :
        final_map[5,vert] = 5
        
    elif final_value == 6 :
        final_map[6,vert] = 6
        
    elif final_value == 7 :
        final_map[7,vert] = 7
        

# Export finals map
final_stats_dir = '{}/{}/derivatives/pp_data/sub-170k/170k/final_stats/results/'.format(main_dir, project_dir)
os.makedirs(final_stats_dir, exist_ok=True)
final_stats_fn = 'sub-170k_final-stats.dtseries.nii'


maps_names = ['all','pursuit','saccade', 'pursuit_and_saccade', 'vision', 'vision_and_pursuit', 'vision_and_saccade', 'vision_and_pursuit_and_saccade']
final_img = make_surface_image(data=final_map, source_img=img, maps_names=maps_names)
nb.save(final_img, '{}/{}'.format(final_stats_dir, final_stats_fn))


#Export data in TSV
code_name = ['non_responding', 'pursuit', 'saccade', 'pursuit_and_saccade', 'vision', 'vision_and_pursuit', 'vision_and_saccade', 'vision_and_pursuit_and_saccade']                
concat_rois_list = [analysis_info['mmp_rois'], analysis_info['rois']]
for n_list, rois_list in enumerate(concat_rois_list):
    rois = rois_list
    if 'LO' in rois_list:
        atlas_name = 'mmp_group'
        tsv_suffix = 'derivatives_group'
    else:
        atlas_name = 'mmp'
        tsv_suffix = 'derivatives'
        
    roi_verts_dict = get_rois('sub-170k', 
                              return_concat_hemis=True, 
                              rois=rois, 
                              mask=True, 
                              atlas_name=atlas_name, 
                              surf_size='170k')

    prf_tsv_fn = '{}/{}/derivatives/pp_data/sub-170k/170k/prf/tsv/sub-170k_css-prf_{}.tsv'.format(main_dir, project_dir, tsv_suffix)
    tsv_df = pd.read_table(prf_tsv_fn)
    
    if 'stats_final' not in tsv_df.columns:
        tsv_df['stats_final'] = np.nan


    for roi in roi_verts_dict.keys():
        data_roi = final_map[0, roi_verts_dict[roi]]
        tsv_df.loc[tsv_df['rois'] == roi, 'stats_final'] = data_roi
        
    # replace the value by the code name in tsv_df

    # tsv_df['stats_final'] = tsv_df['stats_final'].apply(lambda x: code_name[int(x)] if pd.notnull(x) else x)
    tsv_df['stats_final'] = tsv_df['stats_final'].apply(lambda x: code_name[int(x)] if isinstance(x, float) and pd.notnull(x) else x)
    tsv_df.to_csv(prf_tsv_fn, sep="\t", na_rep='NaN', index=False)
   

 


# Define permission cmd
os.system("chmod -Rf 771 {main_dir}/{project_dir}".format(main_dir=main_dir, project_dir=project_dir))
os.system("chgrp -Rf {group} {main_dir}/{project_dir}".format(main_dir=main_dir, project_dir=project_dir, group=group))