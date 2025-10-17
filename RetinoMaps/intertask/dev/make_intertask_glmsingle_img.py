"""
-----------------------------------------------------------------------------------------
make_intertask_glmsingle_img.py
-----------------------------------------------------------------------------------------
Goal of the script:
Make final img with significant vertex for pRF, SacLoc and PurLoc
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
>> cd ~/projects/pRF_analysis/RetinoMaps/intertask/dev
2. run python command
python stats_final.py [main directory] [project name] [subject] [group]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/projects/pRF_analysis/RetinoMaps/intertask/dev
python make_intertask_glmsingle_img.py /scratch/mszinte/data RetinoMaps sub-01 327
python make_intertask_glmsingle_img.py /scratch/mszinte/data RetinoMaps sub-170k 327
-----------------------------------------------------------------------------------------
Written by Uriel Lascombes (uriel.lascombes@laposte.net)
Edited by Martin Szinte (mail@martinszinte.net) 
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
import numpy as np
import nibabel as nb

# Personal import
sys.path.append("{}/../../../analysis_code/utils".format(os.getcwd()))
from surface_utils import make_surface_image , load_surface

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
group = sys.argv[4]

# Index
slope_idx, intercept_idx, rvalue_idx, pvalue_idx, stderr_idx, \
    trs_idx, corr_pvalue_5pt_idx, corr_pvalue_1pt_idx = 0, 1, 2, 3, 4, 5, 6, 7

# load settings
with open('../../settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)
if subject == 'sub-170k': 
    formats = ['170k']
    extensions = ['dtseries.nii']
else: 
    formats = analysis_info['formats']
    extensions = analysis_info['extensions']
group_tasks = analysis_info['task_intertask']
fdr_alpha = analysis_info['stats_th']
glm_code_names = analysis_info['glm_code_names']
maps_names_inter_task = analysis_info['maps_names_intertask']

#Set treshold
if fdr_alpha == 0.05: fdr_p_map_idx = corr_pvalue_5pt_idx
elif fdr_alpha == 0.01: fdr_p_map_idx = corr_pvalue_1pt_idx

# find all the stats files 
glm_stats_fns = []
prf_stats_fns = []
for format_, extension in zip(formats, extensions):
    list_glm = glob.glob("{}/{}/derivatives/pp_data/{}/{}/glm/glm_single_derivatives/*stats_loo-median.{}".format(
        main_dir, project_dir, subject, format_, extension))
    list_prf = glob.glob("{}/{}/derivatives/pp_data/{}/{}/prf/prf_derivatives/*stats_loo-median.{}".format(
        main_dir, project_dir, subject, format_, extension))
    
    glm_stats_fns.extend(list_glm)
    prf_stats_fns.extend(list_prf)   

stats_fns = glm_stats_fns + prf_stats_fns

# split filtered files  depending of their nature
if subject != 'sub-170k':
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
else: stats_files_list = [stats_fns]

# loop for different group of tasks (sac/pur, sacVE/purVE)
for tasks in group_tasks: 
    for stats_files in stats_files_list:
        #  Load file to get size
        img, data = load_surface(fn=stats_files[0])
        final_map = np.zeros((8,data[fdr_p_map_idx].shape[0]))
        
        if stats_files[0].find('hemi-L') != -1: hemi = 'hemi-L'
        elif stats_files[0].find('hemi-R') != -1: hemi = 'hemi-R'
        else: hemi = None
        
        for task in tasks:
            print(task)
            # defind output files names 
            stats_files_tasks = [file for file in stats_files if task in file]
            
            # make a final map with all  tasks 
            if 'Pur' in task: task_idx = glm_code_names['pursuit']
            elif 'Sac' in task: task_idx = glm_code_names['saccade']
            elif 'pRF' in task: task_idx = glm_code_names['vision']

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
            if final_value == glm_code_names['pursuit_and_saccade'] : final_map[3, vert] = glm_code_names['pursuit_and_saccade']
            elif final_value == glm_code_names['vision_and_pursuit'] : final_map[5, vert] = glm_code_names['vision_and_pursuit']
            elif final_value == glm_code_names['vision_and_saccade'] : final_map[6, vert] = glm_code_names['vision_and_saccade']
            elif final_value == glm_code_names['vision_and_pursuit_and_saccade'] : final_map[7, vert] = glm_code_names['vision_and_pursuit_and_saccade']
                    
        # Export finals map
        if 'SacVELoc' in tasks: suffix = 'SacVE_PurVE'
        else : suffix = 'Sac_Pur'
        if hemi:
            inter_task_dir = '{}/{}/derivatives/pp_data/{}/fsnative/intertask_glmsingle'.format(main_dir, project_dir, subject)
            os.makedirs(inter_task_dir, exist_ok=True)
            inter_task_fn = '{}_{}_intertask_{}.func.gii' .format(subject, hemi, suffix)    
        else:
            inter_task_dir = '{}/{}/derivatives/pp_data/{}/170k/intertask_glmsingle'.format(main_dir, project_dir, subject)
            os.makedirs(inter_task_dir, exist_ok=True)
            inter_task_fn = '{}_intertask_{}.dtseries.nii' .format(subject, suffix)
            
        # Save img 
        print('Save {}/{}'.format(inter_task_dir, inter_task_fn))
        final_img = make_surface_image(data=final_map, source_img=img, maps_names=maps_names_inter_task)
        nb.save(final_img, '{}/{}'.format(inter_task_dir, inter_task_fn))
            
# # Define permission cmd
# print('Changing files permissions in {}/{}'.format(main_dir, project_dir))
# os.system("chmod -Rf 771 {}/{}".format(main_dir, project_dir))
# os.system("chgrp -Rf {} {}/{}".format(group, main_dir, project_dir)) 
    
    