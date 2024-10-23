"""
-----------------------------------------------------------------------------------------
eyetrack_preproc.py
-----------------------------------------------------------------------------------------
Goal of the script:
- Preprocess BIDS formatted eyetracking data 
blinks are removed by excluding samples at which the pupil was lost entirely, excising data 100 ms before and 150 ms after each occurrence. 
removing slow signal drift by linear detrending and median-centering of the gaze position time series (X and Y) (assumes that the median gaze position corresponds to central fixation)
smoothing using a 50-ms running average. 
Any eyetracking run containing fewer than 1/3 valid samples after pre-processing is excluded from further analysis
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: subject 
sys.argv[2]: task
sys.argv[3]: session 
sys.argv[4]: eye
-----------------------------------------------------------------------------------------
Output(s):
Cleaned timeseries data per run 
Tsv trial trigger timestamps 
-----------------------------------------------------------------------------------------
To run:
cd /projects/prf_analysis/RetinoMaps/eyetracking/dev
python eyetrack_preproc.py sub-01 PurLoc ses-01 eye1
-----------------------------------------------------------------------------------------
"""
import pandas as pd
import json
import numpy as np
import re
import matplotlib.pyplot as plt
import glob 
import os
from sklearn.preprocessing import MinMaxScaler
import sys
from statistics import median
from scipy.signal import detrend
# path of utils folder  
sys.path.insert(0, "/Users/sinakling/projects/pRF_analysis/analysis_code/utils") #TODO make path general
from eyetrack_utils import *

# --------------------- Load settings and inputs -------------------------------------

def load_settings(settings_file):
    with open(settings_file) as f:
        settings = json.load(f)
    return settings

def load_inputs():
    return sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]

def ensure_save_dir(base_dir, subject):
    save_dir = f"{base_dir}/{subject}/eyetracking/timeseries"
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    return save_dir

def load_events(main_dir, subject, ses, task): 
    data_events = load_event_files(main_dir, subject, ses, task)
    return data_events 

# --------------------- Data Extraction ------------------------------------------------

def extract_event_and_physio_data(main_dir, subject, task, ses, num_run, eye):
    df_event_runs = extract_data(main_dir, subject, task, ses, num_run, eye, file_type="physioevents")
    df_data_runs = extract_data(main_dir, subject, task, ses, num_run, eye, file_type="physio")
    return df_event_runs, df_data_runs

# --------------------- Preprocessing Methods -----------------------------------------

def remove_blinks(data, method, sampling_rate):
    if method == 'pupil_off':
        return blinkrm_pupil_off(data, sampling_rate)
    else:
        print("No blink removal method specified")
        return data

def drift_correction(data, method, fixation_periods):
    if method == "linear":
        return detrend(data, type='linear')
    elif method == 'median':
        fixation_median = median([item for row in fixation_periods for item in row])
        return np.array([elem - fixation_median for elem in data])
    else:
        print("No drift correction method specified")
        return data

def interpolate_nans(data):
    return interpol_nans(data)

def normalize_data(data):
    print('- normalizing pupil data')
    scaler = MinMaxScaler(feature_range=(-1, 1))
    return scaler.fit_transform(data.reshape(-1, 1)).flatten()

def downsample_data(data, original_rate, target_rate):
    return downsample_to_targetrate(data, original_rate, target_rate)

def apply_smoothing(data, method, settings):
    if method == 'moving_avg':
        return moving_average_smoothing(data, 1000, 50)
    elif method == 'gaussian':
        sigma = settings.get("sigma")
        return gaussian_smoothing(data, 'x_coordinate', sigma)
    else:
        print(f"Unknown smoothing method: {method}")
        return data

# --------------------- Data Saving ----------------------------------------------------

def save_preprocessed_data(data, file_path):
    data.to_csv(file_path, sep='\t', index=False, compression='gzip')

# --------------------- Main Preprocessing Pipeline ------------------------------------

def main_preprocessing_pipeline():
    # Load inputs and settings
    subject, task, ses, eye = load_inputs()
    settings = load_settings('/Users/sinakling/projects/pRF_analysis/RetinoMaps/eyetracking/dev/PurLoc_SacLoc/behavior_settings.json')  #TODO make path general
    
    # Prepare save directory
    main_dir = settings.get('main_dir_mac')
    file_dir_save = ensure_save_dir(f'{main_dir}/derivatives/pp_data', subject)
    
    # Load data
    df_event_runs, df_data_runs = extract_event_and_physio_data(main_dir, subject, task, ses, settings['num_run'], eye)
    print(df_event_runs[0].head()) 
    print(df_data_runs[0].head())  
    
    # Preprocessing loop for each run
    for run_idx, (df_event, df_data) in enumerate(zip(df_event_runs, df_data_runs)):
        # from here per run 
        eye_data_run, time_start_eye, time_end_eye = extract_eye_data_and_triggers(df_event, df_data,settings['first_trial_pattern'], settings['last_trial_pattern']) 
        
        # Apply preprocessing steps based on settings
        eye_data_run = remove_blinks(eye_data_run, settings['blinks_remove'], settings['eyetrack_sampling'])
        eye_data_run = convert_to_dva(eye_data_run, settings['center'], settings['ppd'])
        eye_data_run_x = interpolate_nans(eye_data_run[:,1])
        eye_data_run_y = interpolate_nans(eye_data_run[:,2])
        eye_data_run_p = interpolate_nans(eye_data_run[:,3])
        eye_data_run_p = normalize_data(eye_data_run_p)

        if settings.get('drift_corr'):
<<<<<<< HEAD
            eye_data_run_x = detrending(eye_data_run_x)
            eye_data_run_y = detrending(eye_data_run_y)
=======
            eye_data_run_x = linear_detrending(eye_data_run_x)
            eye_data_run_y = linear_detrending(eye_data_run_y)
>>>>>>> 832e18d9c505f2eda02517588e5917df87f6b980
            
        
        eye_data_run = np.stack((eye_data_run[:,0],
                                 eye_data_run_x, 
                                 eye_data_run_y, 
                                 eye_data_run_p), axis=1)

        if settings.get('downsampling'):
            eye_data_run   = downsample_data(eye_data_run, 1000, 100)

        eye_data_run = pd.DataFrame(eye_data_run, columns=['timestamp', 'x', 'y', 'pupil_size'])

        
        if settings.get('smoothing'):
            eye_data_run = apply_smoothing(eye_data_run, settings['smoothing'], settings)
        
        # Save the preprocessed data
        tsv_file_path = f'{file_dir_save}/{subject}_task-{task}_run_0{run_idx+1}_eyedata.tsv.gz'
        save_preprocessed_data(eye_data_run, tsv_file_path)

        
        print(f"Run {run_idx+1} preprocessed and saved at {tsv_file_path}")

if __name__ == "__main__":
    main_preprocessing_pipeline()
