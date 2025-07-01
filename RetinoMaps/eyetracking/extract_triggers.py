"""
-----------------------------------------------------------------------------------------
extract_triggers.py
-----------------------------------------------------------------------------------------
Goal of the script:
- extract timestamps of experiment for saccade analysis (PurLoc and SacLoc tasks)
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main directory 
sys.argv[2]: project directory 
sys.argv[3]: subject 
sys.argv[4]: task 
sys.argv[5]: group 
-----------------------------------------------------------------------------------------
Output(s):
Hdf5 file per run with all timestamps
tsv file with events and timestamps
-----------------------------------------------------------------------------------------
To run:
cd ~/projects/pRF_analysis/RetinoMaps/eyetracking/
python extract_triggers.py /scratch/mszinte/data RetinoMaps sub-01 327
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
import re
import os
import sys
import h5py
import json
import numpy as np
import pandas as pd

# Personal imports
sys.path.append("{}/../../analysis_code/utils".format(os.getcwd()))
from eyetrack_utils import load_event_files, extract_data

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
group = sys.argv[4]

# Load general settings
with open('../settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)
tasks = analysis_info['task_intertask'][0]
prf_task_name = analysis_info['prf_task_name']

if prf_task_name in tasks:
    tasks.remove(prf_task_name)


for task in tasks :
    print('Processing {} ...'.format(task))
    # Load inputs and settings
    base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../"))
    settings_path = os.path.join(base_dir, project_dir, '{}_settings.json'.format(task))
    with open(settings_path) as f:
        settings = json.load(f)
    if subject == 'sub-01':
        if task == 'pRF': ses = 'ses-02'
        else: ses = 'ses-01'
    else: ses = settings['session']
    
    
    # Load main experiment settings 
    eye = settings['eye']
    num_run = settings['num_run']
    num_seq = settings['num_seq']
    seq_trs = settings['seq_trs']
    eye_mov_seq = settings['eye_mov_seq']
    trials_seq = settings['trials_seq']
    rads = settings['rads']
    pursuits_tr = np.arange(0,seq_trs,2)
    saccades_tr = np.arange(1,seq_trs,2)
    eyetracking_sampling = settings['eyetrack_sampling']
    screen_size = settings['screen_size']
    ppd = settings['ppd']
    
    file_dir_save = '{}/{}/derivatives/pp_data/{}/eyetracking/stats'.format(
        main_dir, project_dir, subject)
    os.makedirs(file_dir_save, exist_ok=True)
    
    
    # ------------- Trigger extraction -------------------------
    # Extract data from physio and physioevents as dataframes 
    df_event_runs = extract_data(main_dir, project_dir, subject, task, ses, num_run, eye, file_type = "physioevents")
    df_data_runs = extract_data(main_dir, project_dir, subject, task, ses, num_run, eye, file_type = "physio")
    
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
    
    
    
    
    data_events = load_event_files(main_dir, project_dir, subject, ses, task)
    
    
    # --------------- Save timestampes with event file data as tsv -----------------
    # Load the data
    for run, path_event_run in enumerate(data_events):
        df_event_run = pd.read_csv(path_event_run, sep="\t", index_col=0)
    
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
    
        # save as tsv 
        df_event_run.to_csv(f"{file_dir_save}/{subject}_task_{task}_run_0{run+1}_triggers.tsv", index=False, sep="\t")
    
    # ------------------ Save all data needed for saccade analysis ----------------------------
    
    # get amplitude sequence from event files
    df = pd.read_csv(data_events[0], sep='\t')
    
    amp_sequence = list(df['eyemov_amplitude'])
    
    # Get one amplitue per sequence 
    sequence_lengths = [16, 32] 
    sequence_length_index = 0 
    first_elements = [] 
    i = 0  
    while i < len(amp_sequence): 
        first_elements.append(amp_sequence[i]) 
        i += sequence_lengths[sequence_length_index] 
        sequence_length_index = (sequence_length_index + 1) % 2 
        
    
    # save as h5 
    h5_file = '{file_dir}/{sub}_task-{task}_eyedata_sac_stats.h5'.format(file_dir=file_dir_save, sub=subject, task=task)
    
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
        h5file.create_dataset(f'amp_sequence', data=first_elements, dtype='float32')
            
    
# # Define permission cmd
# print('Changing files permissions in {}/{}'.format(main_dir, project_dir))
# os.system("chmod -Rf 771 {}/{}".format(main_dir, project_dir))
# os.system("chgrp -Rf {} {}/{}".format(group, main_dir, project_dir))
    
    