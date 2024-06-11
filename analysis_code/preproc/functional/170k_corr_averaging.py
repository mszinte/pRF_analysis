"""
-----------------------------------------------------------------------------------------
170k_corr_averaging.py
-----------------------------------------------------------------------------------------
Goal of the script:
Average correlations of the project tasks on the 170k format.
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: group (e.g. 327)
-----------------------------------------------------------------------------------------
Output(s):
sh file for running batch command
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/[PROJECT]/analysis_code/preproc/functional/
2. run python command
>> python 170k_corr_averaging.py [main directory] [project name] [group]
-----------------------------------------------------------------------------------------
Exemple:
python 170k_corr_averaging.py /scratch/mszinte/data RetinoMaps 327 
-----------------------------------------------------------------------------------------
Written by Uriel Lascombes (uriel.lascombes@laposte.net)
Edited by Martin Szinte (martin.szinte@gmail.com)
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
tasks = analysis_info['task_names']

avg_170k_corr_dir = '{}/{}/derivatives/pp_data/sub-170k/170k/corr/fmriprep_dct_corr'.format(
    main_dir, project_dir)
os.makedirs(avg_170k_corr_dir, exist_ok=True)
 
for task in tasks :
    print(task)       
    
    for n_subject, subject in enumerate(subjects) :
        print('adding {}...'.format(subject))
        
        corr_dir = '{}/{}/derivatives/pp_data/{}/170k/corr/fmriprep_dct_corr'.format(
            main_dir, project_dir, subject)
        corr_fn = '{}_task-{}_fmriprep_dct_corr_bold.dtseries.nii'.format(subject, task)    
        img, data = load_surface(fn='{}/{}'.format(corr_dir, corr_fn))
    
        # Average without nan
        if n_subject == 0:
            data_avg = np.copy(data)
        else:
            data_avg = np.nanmean(np.array([data_avg, data]), axis=0)
    
    # export results 
    avg_170k_corr_fn = 'sub-170k_task-{}_fmriprep_dct_corr_bold.dtseries.nii'.format(task)
    maps_names = ['runs_correlations']
    
    print('saving {}/{}'.format(avg_170k_corr_dir, avg_170k_corr_fn))
    avg_img = make_surface_image(data=data_avg, source_img=img, maps_names=maps_names)
    nb.save(avg_img,'{}/{}'.format(avg_170k_corr_dir, avg_170k_corr_fn))

# Define permission cmd
os.system("chmod -Rf 771 {main_dir}/{project_dir}".format(
    main_dir=main_dir, project_dir=project_dir))
os.system("chgrp -Rf {group} {main_dir}/{project_dir}".format(
    main_dir=main_dir, project_dir=project_dir, group=group))