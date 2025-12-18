"""
-----------------------------------------------------------------------------------------
preproc_end.py
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
>> cd ~/projects/pRF_analysis/analysis_code/preproc/functional/
2. run python command
python preproc_end.py [main directory] [project name] [subject] [group]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/projects/pRF_analysis/analysis_code/preproc/functional/
python preproc_end.py /scratch/mszinte/data MotConf sub-01 327
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
sys.path.append("{}/../../utils".format(os.getcwd()))
from surface_utils import load_surface , make_surface_image
from pycortex_utils import set_pycortex_config_file

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
settings_path = os.path.join(base_dir, project_dir, "settings.json")

with open(settings_path) as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)
TR = analysis_info['TR']
tasks = analysis_info['task_names']
high_pass_threshold = analysis_info['high_pass_threshold'] 
sessions = analysis_info['sessions']
anat_session = analysis_info['anat_session'][0]
formats = analysis_info['formats']
extensions = analysis_info['extensions']

preproc_prep = analysis_info['preproc_prep']
normalization = analysis_info['normalization']
filtering = analysis_info['filtering']
# gauss_avg = analysis_info['gauss_avg']
# css_avg = analysis_info['css_avg']

# Make extension folders
for format_, extension in zip(formats, extensions):
    os.makedirs("{}/{}/derivatives/pp_data/{}/{}/func/{}_{}".format(
        main_dir, project_dir, subject, format_, preproc_prep, filtering), exist_ok=True)
    os.makedirs("{}/{}/derivatives/pp_data/{}/{}/corr/{}_{}_corr".format(
        main_dir, project_dir, subject, format_, preproc_prep, filtering), exist_ok=True)
    # os.makedirs("{}/{}/derivatives/pp_data/{}/{}/func/{}_{}_{}".format(
    #         main_dir, project_dir, subject, format_, preproc_prep, filtering, gauss_avg), exist_ok=True)
    # os.makedirs("{}/{}/derivatives/pp_data/{}/{}/func/{}_{}_{}".format(
    #         main_dir, project_dir, subject, format_, preproc_prep, filtering, css_avg), exist_ok=True)
    
# High pass filtering
for format_, extension in zip(formats, extensions):
    print('high pass filtering : {}'.format(format_))
    for session in sessions:
        # Find outputs from fMRIprep
        fmriprep_func_fns = glob.glob("{}/{}/derivatives/fmriprep/fmriprep/{}/{}/func/*{}*.{}".format(
            main_dir, project_dir, subject, session, format_, extension)) 

        if not fmriprep_func_fns:
            print('No files for {}'.format(session))
            continue
            
        for func_fn in fmriprep_func_fns :

            # Make output filtered filenames
            filtered_data_fn_end = func_fn.split('/')[-1].replace('bold', '{}_{}_bold'.format(
                filtering, normalization))

            # Load data
            surf_img, surf_data = load_surface(fn=func_fn)
           
            # High pass filtering 
            nb_tr = surf_data.shape[0]
            ft = np.linspace(0.5 * TR, (nb_tr + 0.5) * TR, nb_tr, endpoint=False)
            high_pass_set = _cosine_drift(high_pass_threshold, ft)
            surf_data = signal.clean(surf_data, 
                                      detrend=False,
                                      standardize=False, 
                                      confounds=high_pass_set)
           
            # Compute the Z-score 
            surf_data = (surf_data - np.mean(surf_data, axis=0)) / np.std(surf_data, axis=0)
            
            # Make an image with the preproceced data
            filtered_img = make_surface_image(data=surf_data, source_img=surf_img)

            # Save surface
            filtered_fn = "{}/{}/derivatives/pp_data/{}/{}/func/{}_{}_{}/{}".format(
                main_dir, project_dir, subject, format_, preproc_prep, filtering, normalization, filtered_data_fn_end)

            nb.save(filtered_img, filtered_fn)


# # Averaging

# # Find all the filtered files 
# preproc_fns = []
# for format_, extension in zip(formats, extensions):
#     list_ = glob.glob("{}/{}/derivatives/pp_data/{}/{}/func/{}_{}/*_*.{}".format(
#             main_dir, project_dir, subject, format_, preproc_prep, filtering, extension))
#     preproc_fns.extend(list_)

# # Split filtered files  depending of their nature
# preproc_fsnative_hemi_L, preproc_fsnative_hemi_R, preproc_170k = [], [], []
# for subtype in preproc_fns:
#     if "hemi-L" in subtype:
#         preproc_fsnative_hemi_L.append(subtype)
#     elif "hemi-R" in subtype:
#         preproc_fsnative_hemi_R.append(subtype)
#     elif "170k" in subtype:
#         preproc_170k.append(subtype)
        
# preproc_files_list = [preproc_fsnative_hemi_L,
#                       preproc_fsnative_hemi_R,
#                       preproc_170k]

# # Loop over files
# for preproc_files in preproc_files_list:
#     for task in tasks:
#         # Defind output files names 
#         preproc_files_task = [file for file in preproc_files if 'task-{}'.format(task) in file]

#         if not preproc_files_task:
#             print('No files for {}'.format(task))
#             continue

#         if preproc_files_task[0].find('hemi-L') != -1: hemi = 'hemi-L'
#         elif preproc_files_task[0].find('hemi-R') != -1: hemi = 'hemi-R'
#         else: hemi = None

#         # Averaging computation
#         preproc_img, preproc_data = load_surface(fn=preproc_files_task[0])
#         data_avg = np.zeros(preproc_data.shape)
#         for preproc_file in preproc_files_task:
#             preproc_img, preproc_data = load_surface(fn=preproc_file)
#             data_avg += preproc_data/len(preproc_files_task)
    
#         # Export averaged data
#         if hemi:
#             avg_fn = "{}/{}/derivatives/pp_data/{}/fsnative/func/{}_{}_{}/{}_task-{}_{}_{}_{}_{}_bold.func.gii".format(
#                 main_dir, project_dir, subject, preproc_prep, filtering, gauss_avg, subject, task, hemi, preproc_prep, filtering, gauss_avg)
#         else:
#             avg_fn = "{}/{}/derivatives/pp_data/{}/170k/func/{}_{}_{}/{}_task-{}_{}_{}_{}_bold.dtseries.nii".format(
#                 main_dir, project_dir, subject, preproc_prep, filtering, gauss_avg, subject, task, preproc_prep, filtering, gauss_avg)

#         print('avg save: {}'.format(avg_fn))
#         avg_img = make_surface_image(data=data_avg, source_img=preproc_img)
#         nb.save(avg_img, avg_fn)

#         # Leave-one-out averaging
#         if len(preproc_files_task):
#             combi = []
#             combi = list(it.combinations(preproc_files_task, len(preproc_files_task)-1))

#         for loo_num, avg_runs in enumerate(combi):
#             # Load data and make the loo_avg object
#             preproc_img, preproc_data = load_surface(fn=preproc_files_task[0])
#             data_loo_avg = np.zeros(preproc_data.shape)
        
#             # Compute leave-one-out averaging
#             for avg_run in avg_runs:
#                 print('loo_avg-{} add: {}'.format(loo_num+1, avg_run))
#                 preproc_img, preproc_data = load_surface(fn=avg_run)
#                 data_loo_avg += preproc_data/len(avg_runs)
                
#             # Export leave one out file 
#             preproc_prep = analysis_info['preproc_prep']
#             filtering = analysis_info['filtering']
#             gauss_avg = analysis_info['gauss_avg']
#             css_avg = analysis_info['css_avg']
    
            
#             if hemi:
#                 loo_avg_fn = "{}/{}/derivatives/pp_data/{}/fsnative/func/{}_{}_{}/{}_task-{}_{}_{}_{}_avg_loo-{}_bold.func.gii".format(
#                     main_dir, project_dir, subject, preproc_prep, filtering, css_avg, subject, task, hemi, preproc_prep, filtering, loo_num+1)
#                 loo_fn = "{}/{}/derivatives/pp_data/{}/fsnative/func/{}_{}_{}/{}_task-{}_{}_{}_{}_loo-{}_bold.func.gii".format(
#                     main_dir, project_dir, subject, preproc_prep, filtering, css_avg, subject, task, hemi, preproc_prep, filtering, loo_num+1)
#             else:
#                 loo_avg_fn = "{}/{}/derivatives/pp_data/{}/170k/func/{}_{}_{}/{}_task-{}_{}_{}_avg_loo-{}_bold.dtseries.nii".format(
#                     main_dir, project_dir, subject, preproc_prep, filtering, css_avg, subject, task, preproc_prep, filtering, loo_num+1)
#                 loo_fn = "{}/{}/derivatives/pp_data/{}/170k/func/{}_{}_{}/{}_task-{}_{}_{}_loo-{}_bold.dtseries.nii".format(
#                     main_dir, project_dir, subject, preproc_prep, filtering, css_avg, subject, task, preproc_prep, filtering, loo_num+1)
#             print('loo_avg save: {}'.format(loo_avg_fn))
#             loo_avg_img = make_surface_image(data = data_loo_avg, source_img=preproc_img)
#             nb.save(loo_avg_img, loo_avg_fn)
        
#             for loo in preproc_files_task:
#                 if loo not in avg_runs:
#                     print('loo_avg left: {}'.format(loo))
#                     print("loo save: {}".format(loo_fn))
#                     shutil.copyfile(loo, loo_fn)
                    
# Anatomy
for format_, pycortex_subject in zip(formats, [subject, 'sub-170k']):
    # define folders
    pycortex_flat_dir = '{}/{}/derivatives/pp_data/cortex/db/{}/surfaces'.format(
        main_dir, project_dir, pycortex_subject)
    dest_dir_anat = "{}/{}/derivatives/pp_data/{}/{}/anat".format(
        main_dir, project_dir, subject, format_)
    os.makedirs(dest_dir_anat, exist_ok=True)

    # Import surface anat data
    print("Copying anatomy {}".format(format_))
    for hemi in ['lh', 'rh']:
        if hemi == 'lh':
            anatomical_structure_primary = 'CortexLeft'
            save_hemi = 'hemi-L'
        elif hemi == 'rh':
            anatomical_structure_primary = 'CortexRight'
            save_hemi = 'hemi-R'

        for surf in ['pia', 'inflated', 'wm', 'flat']:
            if surf == 'pia':
                save_surf = 'pial'
                geometric_type = 'Anatomical'
                anatomical_structure_secondary = 'Pial'
            elif surf == 'wm':
                save_surf = 'smoothwm'
                geometric_type = 'Anatomical'
                anatomical_structure_secondary = 'GrayWhite'
            elif surf == 'inflated':
                save_surf = 'inflated'
                geometric_type = 'Inflated'
                anatomical_structure_secondary = None
            elif surf == 'flat':
                save_surf = 'flat'
                geometric_type = 'Flat'
                anatomical_structure_secondary = None
                
            img = nb.load('{}/{}_{}.gii'.format(pycortex_flat_dir, surf, hemi))
            img.darrays[0].meta['AnatomicalStructurePrimary'] = anatomical_structure_primary
            if anatomical_structure_secondary is not None:
                img.darrays[0].meta['AnatomicalStructureSecondary'] = anatomical_structure_secondary
            img.darrays[0].meta['GeometricType'] = geometric_type
            img.darrays[1].datatype = 'NIFTI_TYPE_FLOAT32'
            nb.save(img, '{}/{}_{}_{}.surf.gii'.format(
                dest_dir_anat, subject, save_hemi, save_surf), mode="compat")
            
# Time
end_time = datetime.datetime.now()
print("\nStart time:\t{start_time}\nEnd time:\t{end_time}\nDuration:\t{dur}".format(
        start_time=start_time,
        end_time=end_time,
        dur=end_time - start_time))

# Define permission cmd
print('Changing files permissions in {}/{}'.format(main_dir, project_dir))
os.system("chmod -Rf 771 {}/{}".format(main_dir, project_dir))
os.system("chgrp -Rf {} {}/{}".format(group, main_dir, project_dir))