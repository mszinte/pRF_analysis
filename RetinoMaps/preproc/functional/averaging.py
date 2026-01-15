"""
-----------------------------------------------------------------------------------------
averaging.py
-----------------------------------------------------------------------------------------
Goal of the script:
High-pass filter, z-score, pick anat files
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject name
sys.argv[4]: group of shared data (e.g. 327)
-----------------------------------------------------------------------------------------
Output(s):
# Preprocessed and averaged timeseries files
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/pRF_analysis/RetinoMaps/preproc/functional/
2. run python command
python averaging.py [main directory] [project name] [subject] [group]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/projects/pRF_analysis/RetinoMaps/preproc/functional/
python averaging.py /scratch/mszinte/data RetinoMaps sub-02 327
-----------------------------------------------------------------------------------------
Written by Martin Szinte (martin.szinte@gmail.com)
and Uriel Lascombes (uriel.lascombes@laposte.net)
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
import sys
import glob
import json
import shutil
import datetime
import numpy as np
import nibabel as nb
import itertools as it
from nilearn import signal
from nilearn.glm.first_level.design_matrix import _cosine_drift

# Personal imports
sys.path.append("{}/../../../analysis_code/utils".format(os.getcwd()))
from surface_utils import load_surface , make_surface_image
from pycortex_utils import set_pycortex_config_file
from settings_utils import load_settings

# Time
start_time = datetime.datetime.now()

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
group = sys.argv[4]

# Set pycortex db and colormaps
cortex_dir = "{}/{}/derivatives/pp_data/cortex".format(main_dir, project_dir)
set_pycortex_config_file(cortex_dir)

# Load settings
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../"))
settings_path = os.path.join(base_dir, project_dir, "settings.yml")
prf_settings_path = os.path.join(base_dir, project_dir, "prf-analysis.yml")
analysis_info, _, _, _ = load_settings([settings_path, prf_settings_path])
tasks = analysis_info['task_names']
sessions = analysis_info['sessions']
formats = analysis_info['formats']
extensions = analysis_info['extensions']
preproc_prep = analysis_info['preproc_prep']
filtering = analysis_info['filtering']
normalization = analysis_info['normalization']
avg_methods = analysis_info['avg_methods']

# Make extension folders
for format_, extension in zip(formats, extensions):
    for avg_method in avg_methods:
        os.makedirs("{}/{}/derivatives/pp_data/{}/{}/func/{}_{}_{}_{}".format(
            main_dir, project_dir, subject, format_, preproc_prep, 
            filtering, normalization, avg_method), exist_ok=True)


# Find all the preprocessed files 
preproc_fns = []

for format_, extension in zip(formats, extensions):
    list_ = glob.glob("{}/{}/derivatives/pp_data/{}/{}/func/{}_{}_{}/*_*.{}".format(
        main_dir, project_dir, subject, format_, preproc_prep, filtering, normalization, extension))
    preproc_fns.extend(list_)

# Split filtered files  depending of their nature
preproc_fsnative_hemi_L, preproc_fsnative_hemi_R, preproc_170k = [], [], []
for subtype in preproc_fns:
    if "hemi-L" in subtype:
        preproc_fsnative_hemi_L.append(subtype)
    elif "hemi-R" in subtype:
        preproc_fsnative_hemi_R.append(subtype)
    elif "170k" in subtype:
        preproc_170k.append(subtype)
        
preproc_files_list = [preproc_fsnative_hemi_L,
                      preproc_fsnative_hemi_R,
                      preproc_170k]

# Loop over files
for preproc_files in preproc_files_list:
    
    for task in tasks:
        # Defind output files names 
        preproc_files_task = [file for file in preproc_files if 'task-{}'.format(task) in file]

        if not preproc_files_task:
            print('No files for {}'.format(task))
            continue

        if preproc_files_task[0].find('hemi-L') != -1: hemi = 'hemi-L'
        elif preproc_files_task[0].find('hemi-R') != -1: hemi = 'hemi-R'
        else: hemi = None

        for avg_method in avg_methods:

            if avg_method == 'avg':
                
                # Averaging computation
                preproc_img, preproc_data = load_surface(fn=preproc_files_task[0])
                data_avg = np.zeros(preproc_data.shape)
                for preproc_file in preproc_files_task:
                    preproc_img, preproc_data = load_surface(fn=preproc_file)
                    data_avg += preproc_data/len(preproc_files_task)
            
                # Export averaged data
                if hemi:
                    avg_fn = "{}/{}/derivatives/pp_data/{}/fsnative/func/{}_{}_{}_avg/{}_task-{}_{}_{}_{}_{}_avg_bold.func.gii".format(
                        main_dir, project_dir, subject, preproc_prep, filtering, normalization, 
                        subject, task, hemi, preproc_prep, filtering, normalization)
                else:
                    avg_fn = "{}/{}/derivatives/pp_data/{}/170k/func/{}_{}_{}_avg/{}_task-{}_{}_{}_{}_avg_bold.dtseries.nii".format(
                        main_dir, project_dir, subject, preproc_prep, filtering, normalization, 
                        subject, task, preproc_prep, filtering, normalization)
        
                print('avg save: {}'.format(avg_fn))
                avg_img = make_surface_image(data=data_avg, source_img=preproc_img)
                nb.save(avg_img, avg_fn)    
        
            elif avg_method == 'loo-avg':
                
                # Leave-one-out averaging computation
                if len(preproc_files_task):
                    combi = []
                    combi = list(it.combinations(preproc_files_task, len(preproc_files_task)-1))
        
                for loo_num, avg_runs in enumerate(combi):
                    # Load data and make the loo_avg object
                    preproc_img, preproc_data = load_surface(fn=preproc_files_task[0])
                    data_loo_avg = np.zeros(preproc_data.shape)
                
                    # Compute leave-one-out averaging
                    for avg_run in avg_runs:
                        print('loo-avg-{} add: {}'.format(loo_num+1, avg_run))
                        preproc_img, preproc_data = load_surface(fn=avg_run)
                        data_loo_avg += preproc_data/len(avg_runs)
                        
                    # Export leave one out file 
                    if hemi:
                        loo_avg_fn = "{}/{}/derivatives/pp_data/{}/fsnative/func/{}_{}_{}_loo-avg/{}_task-{}_{}_{}_{}_{}_loo-avg-{}_bold.func.gii".format(
                            main_dir, project_dir, subject, preproc_prep, filtering, normalization, 
                            subject, task, hemi, preproc_prep, filtering, normalization, loo_num+1)
                        loo_fn = "{}/{}/derivatives/pp_data/{}/fsnative/func/{}_{}_{}_loo-avg/{}_task-{}_{}_{}_{}_{}_loo-{}_bold.func.gii".format(
                            main_dir, project_dir, subject, preproc_prep, filtering, normalization, 
                            subject, task, hemi, preproc_prep, filtering, normalization, loo_num+1)
                    else:
                        loo_avg_fn = "{}/{}/derivatives/pp_data/{}/170k/func/{}_{}_{}_loo-avg/{}_task-{}_{}_{}_{}_loo-avg-{}_bold.dtseries.nii".format(
                            main_dir, project_dir, subject, preproc_prep, filtering, normalization, 
                            subject, task, preproc_prep, filtering, normalization, loo_num+1)
                        loo_fn = "{}/{}/derivatives/pp_data/{}/170k/func/{}_{}_{}_loo-avg/{}_task-{}_{}_{}_{}_loo-{}_bold.dtseries.nii".format(
                            main_dir, project_dir, subject, preproc_prep, filtering, normalization, 
                            subject, task, preproc_prep, filtering, normalization, loo_num+1)

                    print('loo_avg save: {}'.format(loo_avg_fn))
                    loo_avg_img = make_surface_image(data = data_loo_avg, source_img=preproc_img)
                    nb.save(loo_avg_img, loo_avg_fn)
                
                    for loo in preproc_files_task:
                        if loo not in avg_runs:
                            print('loo_avg left: {}'.format(loo))
                            print("loo save: {}".format(loo_fn))
                            shutil.copyfile(loo, loo_fn)
                    

# Define permission cmd
print('Changing files permissions in {}/{}'.format(main_dir, project_dir))
os.system("chmod -Rf 771 {}/{}".format(main_dir, project_dir))
os.system("chgrp -Rf {} {}/{}".format(group, main_dir, project_dir))

# Time
end_time = datetime.datetime.now()
print("\nStart time:\t{start_time}\nEnd time:\t{end_time}\nDuration:\t{dur}".format(
        start_time=start_time,
        end_time=end_time,
        dur=end_time - start_time))