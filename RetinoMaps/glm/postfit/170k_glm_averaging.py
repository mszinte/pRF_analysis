"""
-----------------------------------------------------------------------------------------
170k_glm_averaging.py
-----------------------------------------------------------------------------------------
Goal of the script:
Average all all the subject of the studie on the 170k space.
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
>> cd ~/projects/RetinoMaps/analysis_code/postproc/glm/
2. run python command
>> python 170gaussgridfit_averaging.py [main directory] [project name] [group] [model]
    [server project]
-----------------------------------------------------------------------------------------
Exemple:
python 170k_glm_averaging.py /scratch/mszinte/data RetinoMaps 327 
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
import numpy as np
import nibabel as nb
deb = ipdb.set_trace

# Personal imports
sys.path.append("{}/../../utils".format(os.getcwd()))
from surface_utils import load_surface , make_surface_image

# inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
group = sys.argv[3]

with open('../../settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)
subjects = analysis_info['subjects']
tasks = analysis_info['task_glm']
glm_alpha = analysis_info['glm_alpha']

fdr_p_map_idx = 3

avg_170k_glm_dir = '{}/{}/derivatives/pp_data/sub-170k/170k/glm/glm_derivatives'.format(main_dir, project_dir)
os.makedirs(avg_170k_glm_dir, exist_ok=True)


for n_task, task in enumerate(tasks) :
    print(task)       
    for n_subject, subject in enumerate(subjects) :
        print('adding {} to averaging'.format(subject)) 
        
        # define directory and file name
        glm_dir = '{}/{}/derivatives/pp_data/{}/170k/glm/glm_derivatives'.format(main_dir, project_dir, subject)
        glm_fn = '{}_task-{}_space-fsLR_den-170k_dct_glm-fit_avg.dtseries.nii'.format(subject, task)

        # Load data
        img, data = load_surface(fn='{}/{}'.format(glm_dir, glm_fn))
    
        # Average without considering nan 
        if n_subject == 0:
            data_avg = np.copy(data)
        else:
            data_avg = np.nanmean(np.array([data_avg, data]), axis=0)
            
            if task == 'PurLoc':
                task_idx = 1

    
 
    #  export results 
    avg_170k_glm_fn = 'sub-170k_task-{}_space-fsLR_den-170k_dct_glm-fit_avg.dtseries.nii'.format(task)
    maps_names = ['z_map','z_p_map','fdr','fdr_p_map', 'r_map']
    
    print('saving {}/{}'.format(avg_170k_glm_dir, avg_170k_glm_fn))
    
    avg_img = make_surface_image(data=data_avg, source_img=img, maps_names=maps_names)
    nb.save(avg_img,'{}/{}'.format(avg_170k_glm_dir, avg_170k_glm_fn))
    
    #  make final map for sub-170k
    if n_task == 0:
        final_map = np.zeros((1,data[fdr_p_map_idx].shape[0])) 
        
    if task == 'PurLoc':
        task_idx = 1

    elif task == 'SacLoc':
        task_idx = 2
    
    for vert, fdr_value in enumerate(data_avg[fdr_p_map_idx,:]):
        if fdr_value < glm_alpha:
            final_map[:,vert] += task_idx
            
final_fn = "sub-170k_task-eyes-mvt_space-fsLR_den-170k_dct_glm-significant_map.dtseries.nii"
    
maps_names_2 = ['significative_map']
final_img = make_surface_image(data=final_map, source_img=img, maps_names=maps_names_2)
nb.save(final_img, '{}/{}'.format(avg_170k_glm_dir, final_fn))

# # Define permission cmd
# os.system("chmod -Rf 771 {main_dir}/{project_dir}".format(main_dir=main_dir, project_dir=project_dir))
# os.system("chgrp -Rf {group} {main_dir}/{project_dir}".format(main_dir=main_dir, project_dir=project_dir, group=group))