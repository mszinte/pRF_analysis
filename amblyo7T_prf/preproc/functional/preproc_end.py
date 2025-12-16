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
sys.argv[5]: session name (optional, e.g. ses-01)
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
python preproc_end.py /scratch/mszinte/data amblyo7T_prf sub-01 327 ses-01
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
session = sys.argv[5] if len(sys.argv) > 5 else None

# Handle session parameter for pycortex subject name
if session:
    pycortex_subject_name = f"{subject}_{session}"
else:
    pycortex_subject_name = subject

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

preproc_prep = 'fmriprep'
filtering = 'dct'
run_grouping = 'concat'

# Make extension folders
for format_, extension in zip(formats, extensions):
    os.makedirs("{}/{}/derivatives/pp_data/{}/{}/func/{}_{}".format(
        main_dir, project_dir, subject, format_, preproc_prep, filtering), exist_ok=True)
    os.makedirs("{}/{}/derivatives/pp_data/{}/{}/func/{}_{}_{}".format(
        main_dir, project_dir, subject, format_, preproc_prep, filtering, run_grouping), exist_ok=True)
    
# High pass filtering
# Dictionary to collect masks per format and hemisphere
masks_dict = {'fsnative_hemi-L': [], 'fsnative_hemi-R': [], '170k': []}

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
            filtered_data_fn_end = func_fn.split('/')[-1].replace('bold', '{}_bold'.format(filtering))

            # Load data
            surf_img, surf_data = load_surface(fn=func_fn)
            
            # Identify valid vertices (exclude non-covered areas)
            vertex_means = np.mean(surf_data, axis=0)
            valid_mask = vertex_means > 10000  # Threshold between negatives and valid range
            
            print('Valid vertices: {}/{} ({:.1f}%)'.format(
                np.sum(valid_mask), len(valid_mask), 100*np.sum(valid_mask)/len(valid_mask)))
            
            # Store mask for later intersection
            if 'hemi-L' in func_fn:
                masks_dict['fsnative_hemi-L'].append(valid_mask)
            elif 'hemi-R' in func_fn:
                masks_dict['fsnative_hemi-R'].append(valid_mask)
            elif '170k' in func_fn:
                masks_dict['170k'].append(valid_mask)
           
            # High pass filtering 
            nb_tr = surf_data.shape[0]
            ft = np.linspace(0.5 * TR, (nb_tr + 0.5) * TR, nb_tr, endpoint=False)
            high_pass_set = _cosine_drift(high_pass_threshold, ft)
            surf_data = signal.clean(surf_data, 
                                      detrend=False,
                                      standardize=False, 
                                      confounds=high_pass_set)
           
            # Compute the Z-score only on valid vertices
            surf_data[:, valid_mask] = (surf_data[:, valid_mask] - np.mean(surf_data[:, valid_mask], axis=0)) / np.std(surf_data[:, valid_mask], axis=0)
            
            # Set invalid vertices to NaN
            surf_data[:, ~valid_mask] = np.nan
            
            # Make an image with the preproceced data
            filtered_img = make_surface_image(data=surf_data, source_img=surf_img)

            # Save surface
            filtered_fn = "{}/{}/derivatives/pp_data/{}/{}/func/{}_{}/{}".format(
                main_dir, project_dir, subject, format_, preproc_prep, filtering, filtered_data_fn_end)

            nb.save(filtered_img, filtered_fn)

# Create and save consistent masks
print("\nCreating consistent masks across runs")

saved_masks = {}
for key, mask_list in masks_dict.items():
    if not mask_list:
        continue
    
    # Create intersection mask
    consistent_mask = np.logical_and.reduce(mask_list)
    print(f'{key}: {np.sum(consistent_mask)}/{len(consistent_mask)} ({100*np.sum(consistent_mask)/len(consistent_mask):.1f}%)')
    
    # Create mask data (1 for valid, 0 for invalid)
    mask_data = consistent_mask.astype(np.float32)
    
    # Save mask as GIFTI or CIFTI
    if 'fsnative' in key:
        hemi = key.split('_')[1]
        mask_dir = "{}/{}/derivatives/pp_data/{}/fsnative/func/mask".format(main_dir, project_dir, subject)
        os.makedirs(mask_dir, exist_ok=True)
        # Load a reference file to get the structure
        ref_files = glob.glob("{}/{}/derivatives/pp_data/{}/fsnative/func/{}_{}/*{}*.func.gii".format(
            main_dir, project_dir, subject, preproc_prep, filtering, hemi))
        if ref_files:
            ref_img = nb.load(ref_files[0])
            mask_img = make_surface_image(data=mask_data[np.newaxis, :], source_img=ref_img)
            mask_fn = "{}/{}_{}_{}_mask.func.gii".format(mask_dir, subject, hemi, 'space-fsnative')
            nb.save(mask_img, mask_fn)
            saved_masks[key] = mask_fn
            print(f'  Saved: {mask_fn}')
    elif '170k' in key:
        mask_dir = "{}/{}/derivatives/pp_data/{}/170k/func/mask".format(main_dir, project_dir, subject)
        os.makedirs(mask_dir, exist_ok=True)
        ref_files = glob.glob("{}/{}/derivatives/pp_data/{}/170k/func/{}_{}/*170k*.dtseries.nii".format(
            main_dir, project_dir, subject, preproc_prep, filtering))
        if ref_files:
            ref_img = nb.load(ref_files[0])
            mask_img = make_surface_image(data=mask_data[np.newaxis, :], source_img=ref_img)
            mask_fn = "{}/{}_space-fsLR_den-170k_mask.dtseries.nii".format(mask_dir, subject)
            nb.save(mask_img, mask_fn)
            saved_masks[key] = mask_fn
            print(f'  Saved: {mask_fn}')


# Find all the filtered files 
preproc_fns = []
for format_, extension in zip(formats, extensions):
    list_ = glob.glob("{}/{}/derivatives/pp_data/{}/{}/func/{}_{}/*_*.{}".format(
            main_dir, project_dir, subject, format_, preproc_prep, filtering, extension))
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
            output_dir = "{}/{}/derivatives/pp_data/{}/fsnative/func/{}_{}_{}".format(
                main_dir, project_dir, subject, preproc_prep, filtering, run_grouping)
            os.makedirs(output_dir, exist_ok=True)
            concat_fn = "{}/{}_task-pRF{}_{}_{}_{}_{}_bold.func.gii".format(
                output_dir, subject, eye_label, hemi, preproc_prep, filtering, run_grouping)
        else:
            output_dir = "{}/{}/derivatives/pp_data/{}/170k/func/fmriprep_dct_concat".format(
                main_dir, project_dir, subject)
            os.makedirs(output_dir, exist_ok=True)
            concat_fn = "{}/{}_task-pRF{}_{}_{}_{}_bold.dtseries.nii".format(
                output_dir, subject, eye_label, preproc_prep, filtering, run_grouping)
        
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
os.system("chmod -Rf 771 {}/{}".format(main_dir, project_dir))
os.system("chgrp -Rf {} {}/{}".format(group, main_dir, project_dir))