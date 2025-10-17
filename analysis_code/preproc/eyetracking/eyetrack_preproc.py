"""
-----------------------------------------------------------------------------------------
eyetrack_preproc.py
-----------------------------------------------------------------------------------------
Goal of the script:
- Preprocess BIDS formatted eyetracking data 
blinks are removed by excluding samples at which the pupil was lost entirely, excising 
data 100 ms before and 150 ms after each occurrence. 
removing slow signal drift by linear detrending and median-centering of the gaze position 
time series (X and Y) (assumes that the median gaze position corresponds to central 
fixation)
smoothing using a 50-ms running average. 
Any eyetracking run containing fewer than 1/3 valid samples after pre-processing is 
excluded from further analysis
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject name
sys.argv[4]: group of shared data (e.g. 327)
-----------------------------------------------------------------------------------------
Output(s):
Cleaned timeseries data per run 
Tsv trial trigger timestamps 
-----------------------------------------------------------------------------------------
To run:
cd ~/projects/pRF_analysis/analysis_code/preproc/eyetracking/
python eyetrack_preproc.py /scratch/mszinte/data RetinoMaps sub-01 327
-----------------------------------------------------------------------------------------
Written by Sina Kling
Edited by Uriel Lascombes (uriel.lascombes@laposte.net)
-----------------------------------------------------------------------------------------
"""
# Stop warnings
import warnings
warnings.filterwarnings("ignore")

# Debug
import ipdb
deb = ipdb.set_trace

# Imports
import os
import json
import sys
import numpy as np
import pandas as pd
from statistics import median
from scipy.signal import detrend
from sklearn.preprocessing import MinMaxScaler

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
group = sys.argv[4]

# Personal imports
sys.path.append("{}/../../utils".format(os.getcwd()))
from eyetrack_utils import *

# Load general settings
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../"))
settings_path = os.path.join(base_dir, project_dir, "settings.json")

with open(settings_path) as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)
tasks = analysis_info['task_intertask'][0]
prf_task_name = analysis_info['prf_task_name']

# Execption for subject 1 with no data for eye tracking
if subject == 'sub-01':
    if prf_task_name in tasks:
        tasks.remove(prf_task_name)

# --------------------- Data Extraction ------------------------------------------------
def extract_event_and_physio_data(main_dir, project_dir, subject, task, ses, num_run, eye):
    df_event_runs = extract_data(main_dir, project_dir, subject, task, ses, num_run, eye, file_type="physioevents")
    df_data_runs = extract_data(main_dir, project_dir, subject, task, ses, num_run, eye, file_type="physio")
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

# --------------------- Main Preprocessing Pipeline ------------------------------------

for task in tasks :
    print('Processing {} ...'.format(task))
    # Load task settings
    base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../"))
    settings_path = os.path.join(base_dir, project_dir, '{}_settings.json'.format(task))
    with open(settings_path) as f:
        settings = json.load(f)       
    if subject == 'sub-01':
        ses = 'ses-01'
    else: ses = settings['session']
    eye = settings['eye']
    fixation_column = settings['fixation_column']
    
    # Defind directories
    eye_tracking_dir = '{}/{}/derivatives/pp_data/{}/eyetracking'.format(main_dir, project_dir, subject)
    timeseries_dir = "{}/timeseries".format(eye_tracking_dir)
    os.makedirs(timeseries_dir, exist_ok=True)
    
    # Load data
    df_event_runs, df_data_runs = extract_event_and_physio_data(main_dir,  project_dir, subject, task, ses, settings['num_run'], eye)
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
            
            design_dir_save = '{}/{}/derivatives/pp_data/{}/eyetracking/exp_design'.format(main_dir, project_dir, subject)
            
            eye_data_run_x = detrending(eye_data_run_x, 
                                        subject, 
                                        ses, 
                                        run_idx, 
                                        fixation_column, 
                                        task, 
                                        design_dir_save)
            
            eye_data_run_y = detrending(eye_data_run_y, 
                                        subject, 
                                        ses, 
                                        run_idx, 
                                        fixation_column, 
                                        task, 
                                        design_dir_save)
                                        
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
        tsv_file_path = f'{timeseries_dir}/{subject}_task-{task}_run_0{run_idx+1}_eyedata.tsv.gz'
        eye_data_run.to_csv(tsv_file_path, sep='\t', index=False, compression='gzip')

        print(f"Run {run_idx+1} preprocessed and saved at {tsv_file_path}")

# # Define permission cmd
# print('Changing files permissions in {}/{}'.format(main_dir, project_dir))
# os.system("chmod -Rf 771 {}/{}".format(main_dir, project_dir))
# os.system("chgrp -Rf {} {}/{}".format(group, main_dir, project_dir))
