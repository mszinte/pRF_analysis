"""
-----------------------------------------------------------------------------------------
eyetrack_preproc.py
-----------------------------------------------------------------------------------------
Goal of the script:
Preprocess BIDS formatted eyetracking data 
- blinks are removed by excluding samples at which the pupil was lost entirely, 
excluding data before and after each occurrence. 
- removing slow signal drift by linear detrending
- median-centering of the gaze position time series (X and Y) 
(assumes that the median gaze position corresponds to central fixation)
- smoothing using a 50-ms running average
- downsampling
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject name
sys.argv[4]: group of shared data (e.g. 327)
-----------------------------------------------------------------------------------------
Output(s):
Cleaned timeseries data per run 
-----------------------------------------------------------------------------------------
To run:
cd ~/projects/pRF_analysis/RetinoMaps/eyetracking/preproc/
python eyetrack_preproc.py /scratch/mszinte/data RetinoMaps sub-02 327
-----------------------------------------------------------------------------------------
Written by Sina Kling
Edited by Uriel Lascombes (uriel.lascombes@laposte.net)
-----------------------------------------------------------------------------------------
"""
# Debug
import ipdb
deb = ipdb.set_trace

import os
import sys
import numpy as np
import pandas as pd

sys.path.append("{}/../../../analysis_code/utils".format(os.getcwd()))
from eyetrack_utils import *
from eyetrack_utils import extract_data
from settings_utils import load_settings

# inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
group = sys.argv[4]

# load settings
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../"))
eye_tracking_settings_path = os.path.join(base_dir, project_dir, "eyetracking-analysis.yml")
settings_path = os.path.join(base_dir, project_dir, "settings.yml")
settings = load_settings([eye_tracking_settings_path, settings_path])
analysis_info = settings[0]  

tasks = analysis_info['eye-tracking_task_names']

# Prepare save directory 
file_dir_save = f"{main_dir}/{project_dir}/derivatives/pp_data/{subject}/eyetracking/timeseries"
os.makedirs(file_dir_save, exist_ok=True)
design_dir_save = f"{main_dir}/{project_dir}/derivatives/pp_data/{subject}/eyetracking/exp_design"

for task in tasks :
    # Load tasks settings
    base_dir = os.path.abspath(os.path.join(os.getcwd(), "../"))
    task_settings_path = os.path.join(base_dir, "eye-tracking_{}.yml".format(task))
    task_settings = load_settings([task_settings_path, eye_tracking_settings_path])[0] 
    
    if subject == 'sub-01':
        if task == 'pRF': ses = 'ses-02'
        else: ses = 'ses-01'
    else: ses = task_settings['session']
    
    num_run = task_settings['num_run']
    eye = task_settings['eye']    
    
    # Load data
    df_event_runs = extract_data(main_dir, project_dir, subject, task, ses, num_run, eye, file_type="physioevents")
    df_data_runs = extract_data(main_dir, project_dir, subject, task, ses, num_run, eye, file_type="physio")
    
    # Preprocessing for each run
    for run_idx, (df_event, df_data) in enumerate(zip(df_event_runs, df_data_runs)):
    
        eye_data_run, time_start_eye, time_end_eye = extract_eye_data_and_triggers(df_event, df_data, task_settings['first_trial_pattern'], task_settings['last_trial_pattern'])
       
        # Apply preprocessing steps based on settings
        # --------- remove blinks ------------------
        eye_data_run = remove_blinks(eye_data_run, task_settings['blinks_remove'], task_settings['eyetrack_sampling'])
        # ------ convert to dva and center ---------
        eye_data_run = convert_to_dva(eye_data_run, task_settings['center'], task_settings['ppd'])
        # ------------ interpolate -----------------
        eye_data_run_x = interpolate_nans(eye_data_run[:,1])
        eye_data_run_y = interpolate_nans(eye_data_run[:,2])
        eye_data_run_p = interpolate_nans(eye_data_run[:,3])
        # ------- normalise pupil data -------------
        eye_data_run_p = normalize_data(eye_data_run_p)
    
        # ------------- detrending ----------------
        if task == "rest":
            eye_data_run_x = detrending_resting_state(eye_data_run_x)
            eye_data_run_y = detrending_resting_state(eye_data_run_y)
        else: 
            if task_settings.get('drift_corr'):
                eye_data_run_x = detrending(eye_data_run_x, subject, ses, run_idx, task_settings["fixation_column"], task, design_dir_save)
                eye_data_run_y = detrending(eye_data_run_y, subject, ses, run_idx, task_settings["fixation_column"], task, design_dir_save)
    
        
        eye_data_run = np.stack((eye_data_run[:,0],
                                 eye_data_run_x, 
                                 eye_data_run_y, 
                                 eye_data_run_p), axis=1)
        # ------------ downsampling ----------------
        if task_settings.get('downsampling'):
            eye_data_run   = downsample_data(eye_data_run, task_settings["eyetrack_sampling"], task_settings["desired_sampling_rate"])
    
        eye_data_run = pd.DataFrame(eye_data_run, columns=['timestamp', 'x', 'y', 'pupil_size'])
    
        # ------------ smoothing ------------------
        if task_settings.get('smoothing'):
            eye_data_run = apply_smoothing(eye_data_run, task_settings['smoothing'], task_settings)
        
        # Save the preprocessed data as tsv.gz
        tsv_file_path = f'{file_dir_save}/{subject}_task-{task}_run_0{run_idx+1}_eyedata.tsv.gz'        
        eye_data_run.to_csv(tsv_file_path, sep='\t', index=False, compression='gzip')
    
        print(f"Run {run_idx+1} preprocessed and saved at {tsv_file_path}")

# Define permission cmd
print('Changing files permissions in {}/{}'.format(main_dir, project_dir))
os.system("chmod -Rf 771 {}/{}".format(main_dir, project_dir))
os.system("chgrp -Rf {} {}/{}".format(group, main_dir, project_dir))