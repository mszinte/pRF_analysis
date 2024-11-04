"""
-----------------------------------------------------------------------------------------
extract_triggers.py
-----------------------------------------------------------------------------------------
Goal of the script:
- extract timestamps of experiment for saccade analysis 
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: subject 
sys.argv[2]: task
-----------------------------------------------------------------------------------------
Output(s):
Hdf5 file per run with all timestamps
tsv file with events and timestamps
-----------------------------------------------------------------------------------------
To run:
cd /projects/pRF_analysis/RetinoMaps/eyetracking/dev/PurLoc_SacLoc
python extract_triggers.py sub-01 PurLoc ses-01 eye1
-----------------------------------------------------------------------------------------
"""
#%%
import pandas as pd
import json
import numpy as np
import re
import matplotlib.pyplot as plt
import glob 
import os
from sklearn.preprocessing import MinMaxScaler
import sys
import math 
import h5py
# path of utils folder  
script_dir = os.path.dirname(os.path.abspath(__file__))  # Directory of the current script
utils_path = os.path.join(script_dir, "../../../analysis_code/utils")
sys.path.insert(0, utils_path)
from eyetrack_utils import *

# --------------------- Load settings and inputs -------------------------------------

def load_settings(settings_file):
    with open(settings_file) as f:
        settings = json.load(f)
    return settings

def load_events(main_dir, subject, ses, task): 
    data_events = load_event_files(main_dir, subject, ses, task)
    return data_events 

def load_inputs():
    return sys.argv[1], sys.argv[2]

def ensure_save_dir(base_dir, subject):
    save_dir = f"{base_dir}/{subject}/eyetracking/stats"
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    return save_dir

subject, task = load_inputs()
settings = load_settings('{task}_behavior_settings.json')
ses = settings['session']
eye = settings['eye']
# Define file list
main_dir = settings.get('main_dir_mac')

# Load main experiment settings 
num_run = settings['num_run']
num_seq = settings['num_seq']
seq_trs = settings['seq_trs']
eye_mov_seq = settings['eye_mov_seq']
trials_seq = settings['trials_seq']
rads = settings['rads']
pursuits_tr = np.arange(0,seq_trs,2)
saccades_tr = np.arange(1,seq_trs,2)
eyetracking_sampling = settings.['eyetrack_sampling']
screen_size = settings['screen_size']
ppd = settings['ppd']

file_dir_save = ensure_save_dir(f'{main_dir}/derivatives/pp_data', subject)


# ------------- Trigger extraction -------------------------
# Extract data from physio and physioevents as dataframes 
df_event_runs = extract_data(main_dir, subject, task, ses, num_run, eye, file_type = "physioevents")
df_data_runs = extract_data(main_dir, subject, task, ses, num_run, eye, file_type = "physio")

# Extract triggers
eye_data_runs_list = []

time_start_eye = np.zeros((1, num_run))
time_end_eye = np.zeros((1, num_run))
time_start_seq = np.zeros((num_seq, num_run))
time_end_seq = np.zeros((num_seq, num_run))
time_start_trial = np.zeros((seq_trs, num_seq, num_run))
time_end_trial = np.zeros((seq_trs, num_seq, num_run))

# Lists to collect record lines
record_lines = []
ongoing_trials = {}

# Regex patterns for matching
record_start_pattern = r'RECORD_START'
record_stop_pattern = r'RECORD_STOP'
seq_start_pattern = r'sequence\s(\d+)\sstarted'
seq_stop_pattern = r'sequence\s(\d+)\sstopped'
trial_onset_pattern = r'trial\s(\d+)\sonset'
trial_offset_pattern = r'trial\s(\d+)\soffset'
seq_9_stop_pattern = settings.get("last_trial_pattern")
seq_1_start_pattern = settings.get("first_trial_pattern")

# Loop through the 'messages' column to extract the patterns
for run_idx, df in enumerate(df_event_runs):
    for index, row in df.iterrows():
        message = row['message']
        
        if pd.isna(message):
            continue  # Skip if NaN
        
        # Check for RECORD_START and RECORD_STOP
        if re.search(record_start_pattern, message):
            record_lines.append(row['onset'])
        elif re.search(record_stop_pattern, message):
            record_lines.append(row['onset'])
        
        # Check for sequence 1 started
        if re.search(seq_1_start_pattern, message):
            time_start_eye[0, run_idx] = row['onset']  
            

        # Check for sequence start
        seq_start_match = re.search(seq_start_pattern, message)
        if seq_start_match:
            seq_num = int(seq_start_match.group(1))
            last_seq_num = seq_num  # Save sequence number for trial matching
            time_start_seq[seq_num - 1, run_idx] = row['onset']  

        # Check for sequence stop
        seq_stop_match = re.search(seq_stop_pattern, message)
        if seq_stop_match:
            seq_num = int(seq_stop_match.group(1))
            time_end_seq[seq_num - 1, run_idx] = row['onset']  
        # Check for sequence stop
        seq_stop_match = re.search(seq_stop_pattern, message)
        if seq_stop_match:
            seq_num = int(seq_stop_match.group(1))
            time_end_seq[seq_num - 1, run_idx] = row['onset']

        # Check for trial onset
        trial_onset_match = re.search(trial_onset_pattern, message)
        if trial_onset_match:
            trial_num = int(trial_onset_match.group(1))
            if last_seq_num is not None:  # Ensure sequence has been identified
                time_start_trial[trial_num - 1, last_seq_num - 1, run_idx] = row['onset']
                # Store ongoing trial in case offset is found later
                ongoing_trials[trial_num] = row['onset']

        # Check for trial offset (ensure it's stored after the onset)
        trial_offset_match = re.search(trial_offset_pattern, message)
        if trial_offset_match:
            trial_num = int(trial_offset_match.group(1))
            # Check if this trial has an ongoing onset recorded
            if trial_num in ongoing_trials:
                if last_seq_num is not None:
                    time_end_trial[trial_num - 1, last_seq_num - 1, run_idx] = row['onset']
                    del ongoing_trials[trial_num]  # Remove from ongoing trials as offset is found
            else:
                # Trial offset found without a matching onset, this means it was out of order
                print(f"Out-of-order trial offset found for trial {trial_num}, but onset wasn't found.")


        # Check for sequence 9 stopped
        if re.search(seq_9_stop_pattern, message):
            time_end_eye[0, run_idx] = row['onset']

    # Handle any missing elements
    for seq_num in range(1, num_seq + 1): 
        num_trials_in_seq = trials_seq[seq_num - 1]  

        print(f"Checking sequence {seq_num} with {num_trials_in_seq} trials.")

        for trial_num in range(num_trials_in_seq):  
            trial_num_in_data = trial_num + 1  # Convert to 1-based index for matching actual trial numbers

            # Check if the element is 0 in the corresponding sequence and trial (onset)
            if time_start_trial[trial_num, seq_num - 1, run_idx] == 0:
                print(f"Missing onset for trial {trial_num_in_data} in sequence {seq_num}. Searching...")

                # Search for the corresponding trial onset pattern (1-based) for the correct sequence
                trial_onset_search_pattern = rf'sequence\s{seq_num}\strial\s{trial_num_in_data}\sonset'

                # Loop through the 'messages' column again to find the onset
                found_onset = False  # Track if the onset is found
                for index, row in df.iterrows():
                    message = row['message']

                    # Ensure that the message is a string (skip if not)
                    if not isinstance(message, str):
                        continue  

                    # Check if the message contains the trial onset pattern for the correct sequence
                    trial_onset_match = re.search(trial_onset_search_pattern, message)
                    if trial_onset_match:
                        # Append or update the time_start_trial array for the found onset
                        time_start_trial[trial_num, seq_num - 1, run_idx] = row['onset']
                        found_onset = True
                        print(f"Found matching onset for trial {trial_num_in_data} in sequence {seq_num} at index {index} with onset {row['onset']}.")
                        break  # Exit loop once the trial onset is found

                if not found_onset:
                    print(f"Could not find onset for trial {trial_num_in_data} in sequence {seq_num}.")

            # Check if the element is 0 for trial offsets in the corresponding sequence and trial
            if time_end_trial[trial_num, seq_num - 1, run_idx] == 0:
                print(f"Missing offset for trial {trial_num_in_data} in sequence {seq_num}. Searching...")

                # Search for the corresponding trial offset pattern (1-based) for the correct sequence
                trial_offset_search_pattern = rf'sequence\s{seq_num}\strial\s{trial_num_in_data}\soffset'

                # Loop through the 'messages' column again to find the offset
                found_offset = False  # Track if the offset is found
                for index, row in df.iterrows():
                    message = row['message']

                    # Ensure that the message is a string (skip if not)
                    if not isinstance(message, str):
                        continue  

                    print(f"Checking message '{message}' at index {index}...")

                    # Check if the message contains the trial offset pattern for the correct sequence
                    trial_offset_match = re.search(trial_offset_search_pattern, message)
                    if trial_offset_match:
                        # Append or update the time_end_trial array for the found offset
                        time_end_trial[trial_num, seq_num - 1, run_idx] = row['onset']
                        found_offset = True
                        print(f"Found matching offset for trial {trial_num_in_data} in sequence {seq_num} at index {index} with offset {row['onset']}.")
                        break  # Exit loop once the trial offset is found

                if not found_offset:
                    print(f"Could not find offset for trial {trial_num_in_data} in sequence {seq_num}.")




data_events = load_event_files(main_dir, subject, ses, task)


# --------------- Save timestampes with event file data as tsv -----------------
# Load the data
for run, path_event_run in enumerate(data_events):
    df_event_run = pd.read_csv(path_event_run, sep="\t")

    # Flatten the time_start_trial array and filter out elements that are not 0
    flattened_time_start = time_start_trial[:,:,run].flatten()
    filtered_time_start = flattened_time_start[flattened_time_start != 0]

    flattened_time_end = time_end_trial[:,:,run].flatten()
    filtered_time_end = flattened_time_end[flattened_time_end != 0]

    # Check if the lengths match 
    if len(filtered_time_start) != len(df_event_run) or len(filtered_time_end) != len(df_event_run):
        print(f"Warning: Mismatch in lengths. Dataframe rows: {len(df_event_run)}, Filtered time_start_trial length: {len(filtered_time_start)}")

    # Add the filtered time_start_trial data to a new column in df_event_run
    df_event_run['trial_time_start'] = filtered_time_start
    df_event_run['trial_time_end'] = filtered_time_end

    #TODO save duration (difference end start) here as well 


    # Print to check the result
    print(df_event_run.head())
    # save as tsv 

    df_event_run.to_csv(f"{file_dir_save}/{subject}_task_{task}_run_0{run+1}_triggers.tsv", index=False, sep="\t")

# ------------------ Save all data needed for saccade analysis ----------------------------
# get amplitude sequence from event files

dfs = []
legend_amp = {1: 4, 2: 6, 3: 8, 4: 10, 5: "none"}
    
for file_path in data_events:
	df = pd.read_csv(file_path, sep='\t')
	dfs.append(df)

appended_df = pd.concat(dfs, ignore_index=True)
amp_sequence_ev = list(appended_df['eyemov_amplitude'])

amp_sequence = [legend_amp[val] if not math.isnan(val) else float('nan') for val in amp_sequence_ev]

# save as h5 
h5_file = '{file_dir}/{sub}_task-{task}_eyedata_sac_stats.h5'.format(

    file_dir=file_dir_save, sub=subject, task=task)

# Remove existing file if it exists
try:
    os.system(f'rm "{h5_file}"')
except:
    pass

# Open a new HDF5 file for this run
with h5py.File(h5_file, "a") as h5file:
    # Create datasets for this run
    h5file.create_dataset(f'time_start_eye', data=time_start_eye, dtype='float32')
    h5file.create_dataset(f'time_end_eye', data=time_end_eye, dtype='float32')
    h5file.create_dataset(f'time_start_seq', data=time_start_seq, dtype='float32')
    h5file.create_dataset(f'time_end_seq', data=time_end_seq, dtype='float64')
    h5file.create_dataset(f'time_start_trial', data=time_start_trial, dtype='float32')
    h5file.create_dataset(f'time_end_trial', data=time_end_trial, dtype='float32')
    h5file.create_dataset(f'amp_sequence', data=amp_sequence_ev, dtype='float32')
        


