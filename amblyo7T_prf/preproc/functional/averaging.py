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
python preproc_end.py [main directory] [project name] [subject] [group] [session (optional)]
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
print("\nConcatenating runs by eye condition")
for preproc_files in preproc_files_list:
    
    if not preproc_files:
        continue
    
    if preproc_files[0].find('hemi-L') != -1: 
        hemi = 'hemi-L'
        mask_key = 'fsnative_hemi-L'
    elif preproc_files[0].find('hemi-R') != -1: 
        hemi = 'hemi-R'
        mask_key = 'fsnative_hemi-R'
    else: 
        hemi = None
        mask_key = '170k'
    
    # Load the consistent mask
    if mask_key in saved_masks:
        mask_img, mask_data = load_surface(fn=saved_masks[mask_key])
        consistent_mask = mask_data[0, :] > 0.5  # Convert back to boolean
    else:
        print(f"Warning: No mask found for {mask_key}")
        continue
    
    # Separate by eye condition
    left_eye_files = sorted([f for f in preproc_files if 'LeftEye' in f])
    right_eye_files = sorted([f for f in preproc_files if 'RightEye' in f])
    
    for eye_files, eye_label in [(left_eye_files, 'LeftEye'), (right_eye_files, 'RightEye')]:
        
        if not eye_files:
            continue
            
        print(f"  Concatenating {len(eye_files)} files for {eye_label} ({mask_key})")
        
        # Load all files and concatenate along time axis
        all_data = []
        concat_img = None
        
        for fn in eye_files:
            print(f"    Loading: {fn}")
            img, data = load_surface(fn=fn)
            # Apply consistent mask
            data[:, ~consistent_mask] = np.nan
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
            concat_fn = "{}/{}_task-pRF{}_{}_dct_concat_bold.func.gii".format(
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
for format_, pycortex_subject in zip(formats, [pycortex_subject_name, 'sub-170k']):
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
# os.system("chmod -Rf 771 {}/{}".format(main_dir, project_dir))
# os.system("chgrp -Rf {} {}/{}".format(group, main_dir, project_dir))










