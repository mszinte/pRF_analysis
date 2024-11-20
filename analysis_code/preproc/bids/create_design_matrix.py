"""
-----------------------------------------------------------------------------------------
create_design_matrix.py
-----------------------------------------------------------------------------------------
Goal of the script:
Create design matrix from the task BIDS json files and _events.tsv files
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject name
sys.argv[4]: group of shared data (e.g. 327)
-----------------------------------------------------------------------------------------
Output(s):

-----------------------------------------------------------------------------------------
To run locally: 
cd ~/projects/pRF_analysis/analysis_code/preproc/bids/
python create_design_matrix.py /scratch/mszinte/data RetinoMaps pRF sub-01 327

-----------------------------------------------------------------------------------------
"""

# General imports
import pandas as pd
import numpy as np
import json
import sys
import os
import ipdb
deb = ipdb.set_trace

# path of utils folder  
from nilearn.plotting import plot_design_matrix
sys.path.append("{}/../../utils".format(os.getcwd()))
from eyetrack_utils import load_event_files

# inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
task = sys.argv[4]
group = sys.argv[5]

# Load settings
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../"))
settings_path = os.path.join(base_dir, project_dir, f'{task}_settings.json')
with open(settings_path) as f:
    settings = json.load(f)


def create_design_matrix(event_file, json_file):
    # Load the event file into a DataFrame
    event_data = pd.read_csv(event_file, sep='\t')

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

# Define file list
num_run = settings.get('num_run')
file_dir_save = f'{main_dir}/{project_dir}/derivatives/exp_design/{subject}'
os.makedirs(file_dir_save, exist_ok=True)

# Load events files
data_events = load_event_files(main_dir, subject, ses, task)
event_descr_json = f'{main_dir}/{project_dir}/task-{task}_events.json' 

#### WE ARE HERE #####

for run in range(num_run):
    event_file = data_events[run]
    design_matrix = create_design_matrix(event_file, event_descr_json)

    
    
    # Save the design matrix
    #print(design_matrix.head())
    design_matrix.to_csv(f"{file_dir_save}/{task}_design_matrix.csv", index=False)

    plot_design_matrix(design_matrix)
