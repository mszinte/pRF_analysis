#%%
import pandas as pd
import numpy as np
import json
import sys
# path of utils folder  
sys.path.insert(0, "/Users/sinakling/projects/pRF_analysis/analysis_code/utils") #TODO make path general
from eyetrack_utils import *

from nilearn.plotting import plot_design_matrix

def create_design_matrix(event_file, json_file):
    # Load the event file into a DataFrame
    event_data = pd.read_csv(event_file,sep='\t')

    # Load the JSON file with event descriptions
    with open(json_file, 'r') as f:
        event_description = json.load(f)

    # Initialize an empty DataFrame for the design matrix
    n_trials = event_data.shape[0]
    design_matrix = pd.DataFrame(index=range(n_trials))

    # Loop over each key in the JSON file to create design matrix columns
    for key, value in event_description.items():
        # Check if the key has "Levels" (categorical variable)
        if 'Levels' in value:
            # Get all possible levels
            levels = value['Levels']
            # Create a column for each level
            for level_key, level_desc in levels.items():
                # Create boolean values: 1 if the event matches the level, else 0
                design_matrix[f"{key}_{level_desc.replace(' ', '_')}"] = (event_data[key] == int(level_key)).astype(int)
        else:
            # If no levels, just copy the values from the event file as-is
            design_matrix[key] = event_data[key]

    # Return the design matrix
    return design_matrix


with open('/Users/sinakling/Desktop/calib_convers/DeepMReyeCalib_behavior_settings.json') as f:
    settings = json.load(f)

# Define file list
main_dir = settings.get('main_dir_mac')

subject = 'sub-01'
task = 'DeepMReyeCalibTraining'
ses = 'ses-01'
eye = 'eye1'

num_run = settings.get('num_run')

file_dir_save = '/Users/sinakling/Desktop/calib_convers'

data_events = load_event_files(main_dir, subject, ses, task)
event_descr_json = f'/Users/sinakling/Desktop/calib_convers/task-DeepMReyeCalibTraining_events.json' 


for run in range(num_run):
    event_file = data_events[run]

    design_matrix = create_design_matrix(event_file, event_descr_json)

    # Save the design matrix
    #print(design_matrix.head())
    design_matrix.to_csv(f"{file_dir_save}/{task}_design_matrix.csv", index=False)

    plot_design_matrix(design_matrix)
# %%

# %%

# %%
