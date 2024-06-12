"""
-----------------------------------------------------------------------------------------
stats_final.py
-----------------------------------------------------------------------------------------
Goal of the script:
Make final img with significant vertex 
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: group of shared data (e.g. 327)
-----------------------------------------------------------------------------------------
Output(s):
img with significant vertex 
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/RetinoMaps/analysis_code/postproc/stats/
2. run python command
python stats_final.py [main directory] [project name] [subject] [group]
-----------------------------------------------------------------------------------------
Exemple:
python stats_final.py /scratch/mszinte/data RetinoMaps sub-01 327
-----------------------------------------------------------------------------------------
Written by Martin Szinte (mail@martinszinte.net)
Edited by Uriel Lascombes (uriel.lascombes@laposte.net)
-----------------------------------------------------------------------------------------
"""
# stop warnings
import warnings
warnings.filterwarnings("ignore")

# Debug
import ipdb
deb = ipdb.set_trace

# General imports
import os
import sys
import json
import glob
import nibabel as nb
import numpy as np
import pandas as pd

# Personal import
sys.path.append("{}/../../utils".format(os.getcwd()))
from surface_utils import make_surface_image , load_surface
from pycortex_utils import set_pycortex_config_file, get_rois

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
group = sys.argv[4]


fdr_p_map_idx = 6

# load settings
with open('../../settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)
formats = analysis_info['formats']
extensions = analysis_info['extensions']
tasks = analysis_info['task_glm'] + ['pRF']
fdr_alpha = analysis_info['fdr_alpha'][1]

# find all the stats files 
glm_stats_fns = []
prf_stats_fns = []
for format_, extension in zip(formats, extensions):
    list_glm = glob.glob("{}/{}/derivatives/pp_data/{}/{}/glm/stats/*loo-avg*stats.{}".format(main_dir, project_dir, subject, format_, extension))
    list_prf = glob.glob("{}/{}/derivatives/pp_data/{}/{}/prf/prf_derivatives/*loo-avg*stats.{}".format(main_dir, project_dir, subject, format_, extension))
    
    glm_stats_fns.extend(list_glm)
    prf_stats_fns.extend(list_prf)   


stats_fns = glm_stats_fns + prf_stats_fns

# split filtered files  depending of their nature
stats_fsnative_hemi_L, stats_fsnative_hemi_R, stats_170k = [], [], []
for subtype in stats_fns:
    if "hemi-L" in subtype:
        stats_fsnative_hemi_L.append(subtype)
    elif "hemi-R" in subtype:
        stats_fsnative_hemi_R.append(subtype)
    elif "170k" in subtype:
        stats_170k.append(subtype)
        
stats_files_list = [stats_fsnative_hemi_L, 
                    stats_fsnative_hemi_R, 
                    stats_170k]

 

for stats_files in stats_files_list:
    img, data = load_surface(fn=stats_files[0])
    final_map = np.zeros((8,data[fdr_p_map_idx].shape[0]))
    

    if stats_files[0].find('hemi-L') != -1: hemi = 'hemi-L'
    elif stats_files[0].find('hemi-R') != -1: hemi = 'hemi-R'
    else: hemi = None
    
    for task in tasks:
        # defind output files names 
        stats_files_tasks = [file for file in stats_files if task in file]
        
        # make a final map with both the tasks 
        if task == 'PurLoc':
            task_idx = 1

        elif task == 'SacLoc':
            task_idx = 2
  
        elif task == 'pRF':
            task_idx = 4
        
        for stats_file in stats_files_tasks:
            
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
            
    code_name = ['non_responding', 'pursuit', 'saccade', 'pursuit_and_saccade', 'vision', 'vision_and_pursuit', 'vision_and_saccade', 'vision_and_pursuit_and_saccade']                
    # Export finals map
    if hemi:
        final_stats_dir = '{}/{}/derivatives/pp_data/{}/fsnative/final_stats/results/'.format(main_dir, project_dir, subject)
        os.makedirs(final_stats_dir, exist_ok=True)
        final_stats_fn = '{}_{}_final-stats.func.gii' .format(subject, hemi)
        
        # #Export data in TSV
        # # get a dicst with the surface vertices contained in each ROI
        # rois = analysis_info['rois']
        # roi_verts_L, roi_verts_R  = get_rois(subject, return_concat_hemis=False, rois=rois, mask=True, atlas_name=None, surf_size=None)
        
        # if hemi == 'hemi-L':
        #     roi_verts = roi_verts_L
        # elif hemi == 'hemi-R':
        #     roi_verts = roi_verts_R

        # prf_tsv_fn = '{}/{}/derivatives/pp_data/{}/fsnative/prf/tsv/{}_css-prf_derivatives.tsv'.format(main_dir, project_dir, subject, subject)
        # tsv_df = pd.read_table(prf_tsv_fn)
        
        
        # if 'stats_final' not in tsv_df.columns:
        #     tsv_df['stats_final'] = np.nan

    
        # for roi in roi_verts.keys():
        #     data_roi = final_map[0, roi_verts[roi]]
        #     tsv_df.loc[(tsv_df['rois'] == roi) & (tsv_df['hemi'] == hemi)  , 'stats_final'] = data_roi
            
          
        # # tsv_df.loc[tsv_df['hemi'] == hemi  , 'stats_final'].apply(lambda x: code_name[int(x)] if pd.notnull(x) else x)
        # tsv_df['stats_final'] = tsv_df['stats_final'].apply(lambda x: code_name[int(x)] if isinstance(x, float) and pd.notnull(x) else x)


        # tsv_df.to_csv(prf_tsv_fn, sep="\t", na_rep='NaN', index=False)
        
    else:
        final_stats_dir = '{}/{}/derivatives/pp_data/{}/170k/final_stats/results/'.format(main_dir, project_dir, subject)
        os.makedirs(final_stats_dir, exist_ok=True)
        final_stats_fn = '{}_final-stats.dtseries.nii' .format(subject)
        
        # #Export data in TSV
        # concat_rois_list = [analysis_info['mmp_rois'], analysis_info['rois']]
        # for n_list, rois_list in enumerate(concat_rois_list):
        #     rois = rois_list
        #     if 'LO' in rois_list:
        #         atlas_name = 'mmp_group'
        #         tsv_suffix = 'derivatives_group'
        #     else:
        #         atlas_name = 'mmp'
        #         tsv_suffix = 'derivatives'
                
        #     roi_verts_dict = get_rois('sub-170k', 
        #                               return_concat_hemis=True, 
        #                               rois=rois, 
        #                               mask=True, 
        #                               atlas_name=atlas_name, 
        #                               surf_size='170k')
        
        #     prf_tsv_fn = '{}/{}/derivatives/pp_data/{}/170k/prf/tsv/{}_css-prf_{}.tsv'.format(main_dir, project_dir, subject, subject, tsv_suffix)
        #     tsv_df = pd.read_table(prf_tsv_fn)
            
        #     if 'stats_final' not in tsv_df.columns:
        #         tsv_df['stats_final'] = np.nan

    
        #     for roi in roi_verts_dict.keys():
        #         data_roi = final_map[0, roi_verts_dict[roi]]
        #         tsv_df.loc[tsv_df['rois'] == roi, 'stats_final'] = data_roi
                
            
            
        #     # tsv_df['stats_final'] = tsv_df['stats_final'].apply(lambda x: code_name[int(x)] if pd.notnull(x) else x)
        #     tsv_df['stats_final'] = tsv_df['stats_final'].apply(lambda x: code_name[int(x)] if isinstance(x, float) and pd.notnull(x) else x)
        #     tsv_df.to_csv(prf_tsv_fn, sep="\t", na_rep='NaN', index=False)
       
    maps_names = ['all','pursuit','saccade', 'pursuit_and_saccade', 'vision', 'vision_and_pursuite', 'vision_and_saccade', 'vision_and_pursuite_and_saccade']
    final_img = make_surface_image(data=final_map, source_img=img, maps_names=maps_names)
    nb.save(final_img, '{}/{}'.format(final_stats_dir, final_stats_fn))
        
# Define permission cmd
os.system("chmod -Rf 771 {}/{}".format(main_dir, project_dir))
os.system("chgrp -Rf {} {}/{}".format(group, main_dir, project_dir)) 
    
    
