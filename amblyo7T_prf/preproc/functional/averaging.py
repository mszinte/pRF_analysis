"""
-----------------------------------------------------------------------------------------
averaging.py
-----------------------------------------------------------------------------------------
Goal of the script:
Averaging function project dependant
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
python averaging.py [main directory] [project name] [subject] [group]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/projects/pRF_analysis/amblyo7T_prf/preproc/functional/
python averaging.py /scratch/mszinte/data amblyo7T_prf sub-01 327
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
from surface_utils import load_surface, make_surface_image

# Time
start_time = datetime.datetime.now()

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
group = sys.argv[4]

# Load settings
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../"))
settings_path = os.path.join(base_dir, project_dir, "settings.json")

with open(settings_path) as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)

tasks = analysis_info['task_names']
formats = analysis_info['formats']
extensions = analysis_info['extensions']
preproc_prep = analysis_info['preproc_prep']
filtering = analysis_info['filtering']
normalization = analysis_info['normalization']
avg_methods = analysis_info['avg_methods']
prf_task_names = analysis_info['prf_task_names']

# Make extension folders
for format_, extension in zip(formats, extensions):
    for avg_method in avg_methods:
        os.makedirs("{}/{}/derivatives/pp_data/{}/{}/func/{}_{}_{}_{}".format(
            main_dir, project_dir, subject, format_, preproc_prep,
            filtering, normalization, avg_method), exist_ok=True)

# Find all preprocessed files
preproc_fns = []
for format_, extension in zip(formats, extensions):
    list_ = glob.glob("{}/{}/derivatives/pp_data/{}/{}/func/{}_{}_{}/*_*.{}".format(
            main_dir, project_dir, subject, format_, preproc_prep, filtering, normalization, extension))
    preproc_fns.extend(list_)

# Split filtered files depending on their nature
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
for prf_task_name in prf_task_names:
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
        
        if prf_task_name.endswith('RightEye'): eye_label = 'RightEye'
        elif prf_task_name.endswith('LeftEye'): eye_label = 'LeftEye'
        
        if 'BarsBarsRingsWedges' in prf_task_name:
            bold_fns = sorted(
                [f for f in preproc_files if eye_label in f and ('Bars' in f or 'Rings' in f or 'Wedges' in f)],
                key=lambda x: (
                    0 if 'Bars' in x else 1 if 'Rings' in x else 2,
                    0 if 'run-01' in x else 1 if 'run-02' in x else 2))
        elif 'BarsBars' in prf_task_name:
            bold_fns = sorted(
                [f for f in preproc_files if eye_label in f and 'Bars' in f],
                key=lambda x: (0 if 'run-01' in x else 1 if 'run-02' in x else 2))
        elif 'RingsWedges' in prf_task_name:
            bold_fns = sorted(
                [f for f in preproc_files if eye_label in f and ('Rings' in f or 'Wedges' in f)],
                key=lambda x: (0 if 'Rings' in x else 1))
    
        print(f"  Concatenating {len(bold_fns)} files for {prf_task_name} ({mask_key})")
    
        # Load all files and concatenate along time axis
        all_data = []

        for fn in bold_fns:
            print(f"    Loading: {fn}")
            img, data = load_surface(fn=fn)
            all_data.append(data)
            concat_img = img

        # Concatenate along time (axis 0)
        concat_data = np.concatenate(all_data, axis=0)

        # Check for vertices with all NaN values in the concatenated data
        nan_vertices = np.all(np.isnan(concat_data), axis=0)

        # If any vertex is all NaN, set the entire time series for that vertex to NaN
        if np.any(nan_vertices):
            concat_data[:, nan_vertices] = np.nan

        # Create output directory and filename based on hemi
        if hemi:
            output_dir = "{}/{}/derivatives/pp_data/{}/fsnative/func/{}_{}_{}_concat".format(
                main_dir, project_dir, subject, preproc_prep, filtering, normalization)
            os.makedirs(output_dir, exist_ok=True)
            concat_fn = "{}/{}_task-{}_{}_{}_{}_{}_concat_bold.func.gii".format(
                output_dir, subject, prf_task_name, hemi, preproc_prep, filtering, normalization)
        else:
            output_dir = "{}/{}/derivatives/pp_data/{}/170k/func/{}_{}_{}_concat".format(
                main_dir, project_dir, subject, preproc_prep, filtering, normalization)
            os.makedirs(output_dir, exist_ok=True)
            concat_fn = "{}/{}_task-{}_{}_{}_{}_concat_bold.dtseries.nii".format(
                output_dir, subject, prf_task_name, preproc_prep, filtering, normalization)

        print(f"    Saved: {concat_fn}")
        concat_img_final = make_surface_image(data=concat_data, source_img=concat_img)
        nb.save(concat_img_final, concat_fn)

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
