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
python preproc_end.py /scratch/mszinte/data RetinoMaps sub-01 327
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
import yaml
import shutil
import datetime
import numpy as np
import nibabel as nb
import itertools as it
from nilearn import signal
from nilearn.glm.first_level import make_first_level_design_matrix

# Personal imports
sys.path.append("{}/../../utils".format(os.getcwd()))
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
settings = load_settings([settings_path])
analysis_info = settings[0]

TR = analysis_info['TR']
tasks = analysis_info['task_names']
high_pass_threshold = analysis_info['high_pass_threshold'] 
sessions = analysis_info['sessions']
formats = analysis_info['formats']
extensions = analysis_info['extensions']
preproc_prep = analysis_info['preproc_prep']
normalization = analysis_info['normalization']
filtering = analysis_info['filtering']
partial_scan = analysis_info['partial_scan']
pycortex_subject_template = analysis_info['pycortex_subject_template']

# Make extension folders
for format_, extension in zip(formats, extensions):
    os.makedirs("{}/{}/derivatives/pp_data/{}/{}/func/{}_{}_{}".format(
        main_dir, project_dir, subject, format_, preproc_prep, filtering, normalization), exist_ok=True)

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
            filtered_data_fn_end = func_fn.split('/')[-1].replace('bold', '{}_{}_{}_bold'.format(
                preproc_prep, filtering, normalization))

            # Load data
            surf_img, surf_data = load_surface(fn=func_fn)
            
            # Identify valid vertices (exclude non-covered areas)
            if partial_scan:
                vertex_means = np.mean(surf_data, axis=0)
                valid_mask = vertex_means > 10000  # Threshold between negatives and valid range
                print('Valid vertices: {}/{} ({:.1f}%)'.format(
                np.sum(valid_mask), len(valid_mask), 100*np.sum(valid_mask)/len(valid_mask)))
            else:
                valid_mask = np.ones(surf_data.shape[1], dtype=bool)
                            
            # Store mask for later intersection
            if 'hemi-L' in func_fn:
                masks_dict['fsnative_hemi-L'].append(valid_mask)
            elif 'hemi-R' in func_fn:
                masks_dict['fsnative_hemi-R'].append(valid_mask)
            elif '170k' in func_fn:
                masks_dict['170k'].append(valid_mask)


            # High pass filtering 
            # Create design matrix with cosine drift
            nb_tr = surf_data.shape[0]
            design_matrix = make_first_level_design_matrix(frame_times=np.arange(nb_tr) * TR,
                                                           drift_model='cosine',
                                                           high_pass=high_pass_threshold)
            cosine_drift = design_matrix.values[:, :-1]

            surf_data = signal.clean(surf_data, 
                                     detrend=False,
                                     standardize=False, 
                                     confounds=cosine_drift)
           
            # Compute the Z-score only on valid vertices
            surf_data[:, valid_mask] = (surf_data[:, valid_mask] - np.mean(surf_data[:, valid_mask], axis=0)) / np.std(surf_data[:, valid_mask], axis=0)
            
            # Set invalid vertices to NaN
            surf_data[:, ~valid_mask] = np.nan
            
            # Make an image with the preproceced data
            filtered_img = make_surface_image(data=surf_data, source_img=surf_img)

            # Save surface
            filtered_fn = "{}/{}/derivatives/pp_data/{}/{}/func/{}_{}_{}/{}".format(
                main_dir, project_dir, subject, format_, preproc_prep, filtering, normalization, filtered_data_fn_end)

            nb.save(filtered_img, filtered_fn)

# Anatomy
for format_, pycortex_subject in zip(formats, [subject, pycortex_subject_template]):
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