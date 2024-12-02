def extract_data(main_dir, project_dir, subject, task, ses, runs, eye, file_type):
    """
    Load and process eye-tracking data and associated metadata from TSV and JSON files.

    Args:
        main_dir (str): Directory containing the data.
        subject (str): Subject ID.
        task (str): Task name.
        ses (str): Session identifier.
        runs (int): Number of runs.
        eye (str): Eye being tracked (e.g., 'left', 'right').
        file_type (str): Type of the file (e.g., 'recording').

    Returns:
        list: A list of pandas dataframes, each containing the data for one run, with columns defined by the JSON metadata.
    """

    import json 
    import pandas as pd
    df_runs = []
    for run in range(runs):
        json_file_path = f'{main_dir}/{project_dir}/{subject}/{ses}/func/{subject}_{ses}_task-{task}_run-0{run+1}_recording-{eye}_{file_type}.json' #could be eyetrack 
        tsv_file_path = f"{main_dir}/{project_dir}/{subject}/{ses}/func/{subject}_{ses}_task-{task}_run-0{run+1}_recording-{eye}_{file_type}.tsv.gz" #could be eyetrack instead of eyeData
        

        with open(json_file_path, 'r') as file:
            json_data = json.load(file)

        # Extract column names from the JSON
        column_names = json_data['Columns']

        df = pd.read_csv(
            tsv_file_path, 
            compression='gzip', 
            delimiter='\t', 
            header=None,  
            names=column_names,  # Use the column names from JSON
            na_values='n/a'  # Treat 'n/a' as NaN
        )

        df_runs.append(df)

    return df_runs

def extract_eye_data_and_triggers(df_event, df_data, onset_pattern, offset_pattern): 
    """
    Extract eye-tracking data and trial trigger events from event and data frames.

    Args:
        df_event (pd.DataFrame): Event dataframe (physioevents)
        df_data (pd.DataFrame): Eye-tracking data with timestamps (physio)
        onset_pattern (str): Regex pattern to detect the start of trials.
        offset_pattern (str): Regex pattern to detect the end of trials.

    Returns:
        tuple: Numpy array of eye-tracking data within the trial period, start time, and end time of the trial.
    """

    # Extract triggers
    # Initialize arrays to store results
    import re 
    import pandas as pd
    import matplotlib.pyplot as plt
    
    time_start_eye = 0
    time_end_eye = 0
    # Loop through the 'messages' column to extract the patterns

    for index, row in df_event.iterrows():
        message = row['message']
        
        if pd.isna(message):
            continue  # Skip if NaN
        
        # Check for sequence 1 started
        if re.search(onset_pattern, message):
            time_start_eye = row['onset']  # Store by run index
            
    
        # Check for sequence 9 stopped
        if re.search(offset_pattern, message):
            time_end_eye = row['onset']
           

    
    # Filter for only timestamps between first and last trial 
    eye_data_run = df_data[(df_data['timestamp'] >= time_start_eye) & 
                    (df_data['timestamp'] <= time_end_eye)]


    eye_data_run_array = eye_data_run[['timestamp', 'x_coordinate', 'y_coordinate', 'pupil_size']].to_numpy()

    plt_1 = plt.figure(figsize=(15, 6))
    plt.title("Experiment relevant timeseries")
    plt.xlabel('x-coordinate', fontweight='bold')
    plt.plot(eye_data_run_array[:,1])
    plt.show()


    return eye_data_run_array, time_start_eye, time_end_eye

"""
------------------- Preproc functions --------------------------------------

"""

import numpy as np 

def blinkrm_pupil_off(samples, sampling_rate=1000, addms2blink=50, smoothing_duration=200):
    """
    Replace blinks in eye-tracking data (where pupil size is zero) with NaN and extend blink duration for smoothing.

    Args:
        samples (np.array): 4D array of eye-tracking data (time, X, Y, pupil).
        sampling_rate (int): Sampling rate of the data, default is 1000 Hz.

    Returns:
        np.array: Cleaned eye-tracking data with blinks replaced by NaN.
"""
    import numpy as np 
    import matplotlib.pyplot as plt
    print(' - blink replacement with NaN')
    blink_duration_extension = int(sampling_rate / 1000 * addms2blink)
    
    # Detect blinks based on pupil size being 0
    blink_indices = np.where(samples[:, 3] == 0)[0]
    
    blink_bool = np.zeros(len(samples), dtype=bool)
    
    for idx in blink_indices:
        blink_bool[idx] = True
    
    # Adding 50 ms extension to the detected blinks
    for idx in blink_indices:
        start_idx = max(0, idx - blink_duration_extension)
        end_idx = min(len(samples), idx + blink_duration_extension + 1)
        blink_bool[start_idx:end_idx] = True
    
    # Add smoothing around each blink (100 ms before and after)
    smth_kernel = np.ones(int(sampling_rate / 1000 * smoothing_duration)) / (sampling_rate / 1000 * smoothing_duration)
    extended_blink_bool = np.convolve(blink_bool, smth_kernel, mode='same') > 0
    
    # Replace blink points in the samples with NaN
    cleaned_samples = samples.copy()
    cleaned_samples[extended_blink_bool, 1:] = np.nan 

    plt_2 = plt.figure(figsize=(15, 6))
    plt.title("Blink removed timeseries")
    plt.xlabel('x-coordinate', fontweight='bold')
    plt.plot(cleaned_samples[:,1])
    plt.show()
    
    return cleaned_samples

def convert_to_dva(eye_data, center, ppd):
    """
    Convert eye-tracking data to degrees of visual angle (dva), center it, and flip Y-axis.

    Args:
        eye_data (np.array): Eye-tracking data in screen pixel coordinates (n_samples, n_features).
        center (tuple): Screen center (x, y) in pixels.
        ppd (float): Pixels per degree for the conversion.

    Returns:
        np.array: Converted eye-tracking data (centered and converted to dva).
    """
    eye_data[:, 1] = (eye_data[:, 1] - center[0]) / ppd
    eye_data[:, 2] = -1.0 * (eye_data[:, 2] - center[1]) / ppd
    return eye_data


def downsample_to_tr(original_data, eyetracking_rate):
    """
    Downsample eye-tracking data to match the temporal resolution of functional MRI TRs.

    Args:
        original_data (np.array): 1D array of eye-tracking data (e.g., X or Y coordinates).
        eyetracking_rate (int): The sampling rate of the eye-tracking data.

    Returns:
        np.array: Resampled data reshaped to match TRs.
    """
    from scipy.signal import resample
    target_points_per_tr = 10  # 10 data points per 1.2 seconds
    tr_duration = 1.2  # 1.2 sec
    target_rate = target_points_per_tr / tr_duration  # 8.33 Hz

    # Calculate total number of data points in target rate
    eyetracking_in_sec = len(original_data) / eyetracking_rate  # 185 sec
    total_target_points = int(eyetracking_in_sec * target_rate) # 1541

    # Resample the data
    downsampled_data = resample(original_data, total_target_points) # resample into amount of wanted data points

    # Reshape into TRs
    num_trs = int(eyetracking_in_sec / tr_duration) 

    reshaped_data = downsampled_data[:num_trs * target_points_per_tr].reshape(num_trs, target_points_per_tr)

    # Check new shape
    print(reshaped_data.shape)

    return reshaped_data


def downsample_to_targetrate(original_data, eyetracking_rate, target_rate):
    """
    Downsample eye-tracking data to a specified target rate.

    Args:
        original_data (np.array): Eye-tracking data array with columns for timestamp, X, Y, and pupil size.
        eyetracking_rate (int): Sampling rate of the original data.
        target_rate (float): Desired target sampling rate.

    Returns:
        np.array: Downsampled eye-tracking data.
    """
    from scipy.signal import resample
    # Calculate total number of data points in target rate

    eyetracking_in_sec = len(original_data) / eyetracking_rate  
    total_target_points = int(eyetracking_in_sec * target_rate) 
    downsampled_t = resample(original_data[:,0], total_target_points)  # resample into amount of wanted data points
    downsampled_x = resample(original_data[:,1], total_target_points)
    downsampled_y = resample(original_data[:,2], total_target_points)
    downsampled_p = resample(original_data[:,3], total_target_points)

    downsampled_data = np.stack((downsampled_t,
                                 downsampled_x,
                                 downsampled_y, 
                                 downsampled_p))

    return downsampled_data


def moving_average_smoothing(dataframe, eyetracking_rate, window_duration): 
    """
    Apply moving average smoothing to eye-tracking data using a specified window size.

    Args:
        dataframe (pd.DataFrame): Eye-tracking data with 'x' and 'y' columns.
        eyetracking_rate (int): Sampling rate of the data.
        window_duration (int): Duration of the moving average window in milliseconds.

    Returns:
        pd.DataFrame: Smoothed data with rolling averages applied to X and Y coordinates.
    """
    import matplotlib.pyplot as plt
    
    # window duration to sec 
    window_duration = window_duration / 1000
    window_size = int(eyetracking_rate * window_duration)
    

    # Calculate SMA
    dataframe['x'] = dataframe['x'].rolling(window=window_size).mean()
    dataframe['y'] = dataframe['y'].rolling(window=window_size).mean()

    plt_2 = plt.figure(figsize=(15, 6))
    plt.title("Smoothed timeseries")
    plt.xlabel('x-coordinate', fontweight='bold')
    plt.plot(dataframe['x'])
    plt.show()

    return dataframe

def gaussian_smoothing(df, column, sigma):
    """
    Apply Gaussian smoothing to a specified column in a dataframe.

    Args:
        df (pd.DataFrame): Dataframe containing the eye-tracking data.
        column (str): Column name of the data to smooth.
        sigma (float): Standard deviation of the Gaussian kernel.

    Returns:
        np.array: Smoothed data.
    """

    # Apply Gaussian smoothing with convolution
    from scipy.ndimage import gaussian_filter1d
    return gaussian_filter1d(df[column], sigma=sigma)



import numpy as np
from scipy.signal import detrend, resample
import matplotlib.pyplot as plt

def detrending(eyetracking_1D, subject, ses, run, fixation_column, task, design_dir_save): 
    """
    Remove linear trends from eye-tracking data and median-center it during fixation periods for drift correction.

    Args:
        eyetracking_1D (np.array): 1D array of eye-tracking data to detrend.
        task (str): Task type, currently 'pRF' or other.

    Returns:
        np.array: Detrended eye-tracking data with trends removed and median-centered.
    """
    
    # Load and resample fixation data 
    fixation_trials = load_design_matrix_fixations(subject, ses, run, fixation_column, task, design_dir_save)  # Requires design matrix from task (see create_design_matrix.py)
    resampled_fixation_type = resample(fixation_trials, len(eyetracking_1D))
    fixation_bool = resampled_fixation_type > 0.5

    fixation_data = eyetracking_1D[fixation_bool]

    # Fit a linear model for the trend during fixation periods
    fixation_indices = np.where(fixation_bool)[0]
    trend_coefficients = np.polyfit(fixation_indices, fixation_data, deg=1)

    # Apply the linear trend to the entire dataset
    full_indices = np.arange(len(eyetracking_1D))
    linear_trend_full = np.polyval(trend_coefficients, full_indices)

    # Subtract the trend from the full dataset
    detrended_full_data = eyetracking_1D - linear_trend_full

    # Median centering using numpy's median function for consistency with numpy arrays
    fixation_median = np.median(detrended_full_data)
    detrended_full_data -= fixation_median

    # Plot the original and detrended data
    plt.plot(eyetracking_1D, label="Original Data")
    plt.plot(detrended_full_data, label="Detrended Data")
    plt.title("Detrended Full Eye Data")
    plt.xlabel("Time")
    plt.ylabel("Detrended Eye Position")
    plt.legend()
    plt.show()

    return detrended_full_data


    

"""

------------- Other Utils ----------------------------------------------

"""

def interpol_nans(eyetracking_data):
    """
    Interpolate missing (NaN) values in eye-tracking data, filling gaps with the nearest valid data.

    Args:
        eyetracking_data (np.array): Eye-tracking data containing NaN values.

    Returns:
        np.array: Interpolated data with NaNs replaced.
    """

    print("- interpolating data")
    nan_indices = np.isnan(eyetracking_data)

    # Fill NaNs at the start and end with nearest valid values
    if nan_indices[0]:  # If the first value is NaN
        first_valid_idx = np.where(~nan_indices)[0][0]
        eyetracking_data[:first_valid_idx] = eyetracking_data[first_valid_idx]
        
    if nan_indices[-1]:  # If the last value is NaN
        last_valid_idx = np.where(~nan_indices)[0][-1]
        eyetracking_data[last_valid_idx+1:] = eyetracking_data[last_valid_idx]

    # Now interpolate remaining NaNs
    eyetracking_no_nans = np.nan_to_num(eyetracking_data)
    eyetracking_signal_interpolated = np.interp(np.arange(len(eyetracking_data)),
                                               np.where(~nan_indices)[0],
                                               eyetracking_no_nans[~nan_indices])
    
    return eyetracking_signal_interpolated

def load_event_files(main_dir, project_dir, subject, ses, task): 
    """
    Load event files from eye-tracking experiments.

    Args:
        main_dir (str): Main directory containing all experiment data.
        project_dir (str): Main project directory
        subject (str): Subject ID.
        ses (str): Session identifier.
        task (str): Task name.

    Returns:
        list: Sorted list of event file paths.
    """
    import glob
    
    data_events = sorted(glob.glob(r'{main_dir}/{project_dir}/{sub}/{ses}/func/{sub}_{ses}_task-{task}_*_events*.tsv'.format(
        main_dir=main_dir, project_dir=project_dir, sub=subject, ses = ses, task = task)))
    print(f'{main_dir}/{project_dir}/{subject}/{ses}/func/{subject}_{ses}_task-{task}_*_events*.tsv')
    assert len(data_events) > 0, "No event files found"

    return data_events

def load_design_matrix_fixations(subject, ses, run, fixation_column, task, design_dir_save): 
    """
    Load the design matrix and extract fixation trial information.

    Args:
        fixation_column (str): Column name in the design matrix that contains fixation data.

    Returns:
        np.array: Array containing fixation trial information.
    """

    import pandas as pd
   
    design_matrix = pd.read_csv(f"{design_dir_save}/{subject}/{subject}_{ses}_task-{task}_run-0{run+1}_design_matrix.tsv", sep ="\t")
    fixation_trials = np.array(design_matrix[fixation_column])

    return fixation_trials