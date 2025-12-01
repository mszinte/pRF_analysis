"""
-----------------------------------------------------------------------------------------
preproc_end.py
-----------------------------------------------------------------------------------------
Goal of the script:
High-pass filter, z-score, average, loo average and pick anat files
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
>> cd ~/projects/pRF_analysis/amblyo7T_prf/preproc/functional/
2. run python command
python preproc_end.py [main directory] [project name] [subject] [group]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/projects/pRF_analysis/amblyo7T_prf/preproc/functional/
python preproc_end.py /scratch/mszinte/data amblyo7T_prf sub-01 327
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
import ipdb
deb = ipdb.set_trace

# Personal imports
sys.path.append("{}/../../../analysis_code/utils".format(os.getcwd()))
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

# Make extension folders
for format_, extension in zip(formats, extensions):
    os.makedirs("{}/{}/derivatives/pp_data/{}/{}/func/fmriprep_dct".format(
        main_dir, project_dir, subject, format_), exist_ok=True)
    os.makedirs("{}/{}/derivatives/pp_data/{}/{}/func/fmriprep_dct_concat".format(
        main_dir, project_dir, subject, format_), exist_ok=True)
    
# High pass filtering
for format_, extension in zip(formats, extensions):
    print('High pass filtering + z-score : {}'.format(format_))
    for session in sessions:
        # Find outputs from fMRIprep
        fmriprep_func_fns = glob.glob("{}/{}/derivatives/fmriprep/fmriprep/{}/{}/func/*{}*.{}".format(
            main_dir, project_dir, subject, session, format_, extension)) 

        if not fmriprep_func_fns:
            print('No files for {}'.format(session))
            continue
            
        for func_fn in fmriprep_func_fns :

            # Make output filtered filenames
            filtered_data_fn_end = func_fn.split('/')[-1].replace('bold', 'dct_bold')

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
            filtered_fn = "{}/{}/derivatives/pp_data/{}/{}/func/fmriprep_dct/{}".format(
                main_dir, project_dir, subject, format_, filtered_data_fn_end)

            nb.save(filtered_img, filtered_fn)


# Find all the filtered files 
preproc_fns = []
for format_, extension in zip(formats, extensions):
    list_ = glob.glob("{}/{}/derivatives/pp_data/{}/{}/func/fmriprep_dct/*_*.{}".format(
            main_dir, project_dir, subject, format_, extension))
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

# Concatenate files by eye condition
print("Concatenating runs by eye condition")
for preproc_files in preproc_files_list:
    
    if preproc_files[0].find('hemi-L') != -1: 
        hemi = 'hemi-L'
    elif preproc_files[0].find('hemi-R') != -1: 
        hemi = 'hemi-R'
    else: 
        hemi = None
    
    # Separate by eye condition
    left_eye_files = sorted([f for f in preproc_files if 'LeftEye' in f])
    right_eye_files = sorted([f for f in preproc_files if 'RightEye' in f])
    
    for eye_files, eye_label in [(left_eye_files, 'LeftEye'), (right_eye_files, 'RightEye')]:
        
        if not eye_files:
            continue
            
        print(f"  Concatenating {len(eye_files)} files for {eye_label}")
        
        # Load all files and concatenate along time axis
        all_data = []
        concat_img = None
        
        for fn in eye_files:
            print(f"    Loading: {fn}")
            img, data = load_surface(fn=fn)
            all_data.append(data)
            if concat_img is None:
                concat_img = img
        
        # Concatenate along time (axis 0)
        concat_data = np.concatenate(all_data, axis=0)
        
        # Create output directory and filename based on hemi
        if hemi:
            output_dir = "{}/{}/derivatives/pp_data/{}/fsnative/func/fmriprep_dct_concat".format(
                main_dir, project_dir, subject)
            os.makedirs(output_dir, exist_ok=True)
            concat_fn = "{}/{}_task-pRF{}_{}_dct_concat_bold.func.gii".format(q
                output_dir, subject, eye_label, hemi)
        else:
            output_dir = "{}/{}/derivatives/pp_data/{}/170k/func/fmriprep_dct_concat".format(
                main_dir, project_dir, subject)
            os.makedirs(output_dir, exist_ok=True)
            concat_fn = "{}/{}_task-pRF{}_dct_concat_bold.dtseries.nii".format(
                output_dir, subject, eye_label)
        
        print(f"    Saved: {concat_fn}")
        concat_img_final = make_surface_image(data=concat_data, source_img=concat_img)
        nb.save(concat_img_final, concat_fn)

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