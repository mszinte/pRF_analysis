def extract_data(main_dir, subject, task, ses, runs, eye, file_type):
    '''read tsv.gsv data as pandas dataframes. Use information from json as columns'''
    import json 
    import pandas as pd
    df_runs = []
    for run in range(runs):
        json_file_path = f'{main_dir}/{subject}/{ses}/func/{subject}_{ses}_task-{task}_run-0{run+1}_eyeData_recording-{eye}_{file_type}.json'
        tsv_file_path = f"{main_dir}/{subject}/{ses}/func/{subject}_{ses}_task-{task}_run-0{run+1}_eyeData_recording-{eye}_{file_type}.tsv.gz"
        

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
# Extract triggers
# Initialize arrays to store results
    import re 
    import pandas as pd
    import matplotlib.pyplot as plt
    

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

def blinkrm_pupil_off(samples, sampling_rate=1000):
    """
    ----------------------------------------------------------------------
    blinkrm_pupil_off(samples, sampling_rate)
    ----------------------------------------------------------------------
    Goal of the function :
    Replace blinks indicated by pupil = 0 with nan. 
    ----------------------------------------------------------------------
    Input(s) :
    samples: 4D data to be adapted (t,X,Y,P)
    sampling_rate: 1000 by default
    ----------------------------------------------------------------------
    Output(s) :
    cleaned_samples
    ----------------------------------------------------------------------
    Function created by Sina Kling
    ----------------------------------------------------------------------
    """
    import numpy as np 
    import matplotlib.pyplot as plt
    print(' - blink replacement with NaN')
    addms2blink = 50  # ms added to end of blink
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
    
    # Add another 100 ms to the beginning and ending of each blink
    smth_kernel = np.ones(int(sampling_rate / 1000 * 200)) / (sampling_rate / 1000 * 200)
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
    ----------------------------------------------------------------------
    downsample_to_tr(original_data)
    ----------------------------------------------------------------------
    Goal of the function :
    Resample eyetracking signal to wanted resolution
    ----------------------------------------------------------------------
    Input(s) :
    original_data : 1D eyetracking data to be resampled (x or y coordinates)
    ----------------------------------------------------------------------
    Output(s) :
    reshaped_data : 1D resampled and reshaped data (x or y coordinates)
    ----------------------------------------------------------------------
    Function created by Sina Kling
    ----------------------------------------------------------------------
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
    
    # window duration to sec 
    window_duration = window_duration / 1000
    window_size = int(eyetracking_rate * window_duration)
    

    # Calculate SMA
    dataframe['x_coordinate'] = dataframe['x_coordinate'].rolling(window=window_size).mean()
    dataframe['y_coordinate'] = dataframe['y_coordinate'].rolling(window=window_size).mean()


    return dataframe

def gaussian_smoothing(df, column, sigma):
    # Apply Gaussian smoothing with convolution
    from scipy.ndimage import gaussian_filter1d
    return gaussian_filter1d(df[column], sigma=sigma)


def linear_detrending(eyetracking_1D): 
    from scipy.signal import detrend
    import numpy as np
    from scipy.signal import resample
    import matplotlib.pyplot as plt

    fixation_trials = load_design_matrix_fixations('trial_type_fixation')
    resampled_fixation_type = resample(fixation_trials, len(eyetracking_1D))

    # Convert resampled fixation type to a boolean mask
    fixation_bool = resampled_fixation_type > 0.5  

    # Filter the eye data based on fixation periods
    fixation_data = eyetracking_1D[fixation_bool]

    # Detrend the fixation data (remove linear trend)
    detrended_fixation_data = detrend(fixation_data)

    # Generate indices for fixation and the full time series
    fixation_indices = np.where(fixation_bool)[0]
    full_indices = np.arange(len(eyetracking_1D))

    # Fit a linear model 
    trend_coefficients = np.polyfit(fixation_indices, fixation_data, deg=1)

    # Use this model to predict the linear trend across the entire time series
    linear_trend_full = np.polyval(trend_coefficients, full_indices)

    # Subtract the linear trend from the entire dataset
    detrended_full_data = eyetracking_1D - linear_trend_full

    # Plot the detrended full dataset
    plt.plot(eyetracking_1D)
    plt.plot(detrended_full_data)
    plt.title("Detrended Full Eye Data")
    plt.xlabel("Time")
    plt.ylabel("Detrended Eye Position")
    plt.show()

    return detrended_full_data

    

"""

------------- Other Utils ----------------------------------------------

"""

def interpol_nans(eyetracking_data):
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

def load_event_files(main_dir, subject, ses, task): 
    import glob
    data_events = sorted(glob.glob(r'{exp_dir}/{sub}/{ses}/func/{sub}_{ses}_task-{task}_*_events*.tsv'.format(exp_dir=main_dir, sub=subject, ses = ses, task = task)))
    if subject == 'sub-01': 
        ses = 'ses-01'
        data_events = sorted(glob.glob(r'{exp_dir}/{sub}/{ses}/func/{sub}_{ses}_task-{task}_*_events*.tsv'.format(exp_dir=main_dir, sub=subject, ses = ses, task = task)))
        
    assert len(data_events) > 0, "No event files found"

    return data_events


def load_design_matrix_fixations(fixation_column): 
    "path: default project/"
    import pandas as pd
    design_matrix = pd.read_csv("/Users/sinakling/Desktop/design_matrix.csv")
    fixation_trials = np.array(design_matrix[fixation_column])

    return fixation_trials