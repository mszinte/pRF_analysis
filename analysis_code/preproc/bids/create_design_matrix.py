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
sys.argv[4]: task
sys.argv[5]: group of shared data (e.g. 327)
-----------------------------------------------------------------------------------------
Output(s):

-----------------------------------------------------------------------------------------
To run locally: 
cd ~/projects/pRF_analysis/analysis_code/preproc/bids/
python create_design_matrix.py /scratch/mszinte/data RetinoMaps sub-01 pRF 327
-----------------------------------------------------------------------------------------
"""

# General imports
import pandas as pd
import numpy as np
import json
import sys
import os
import matplotlib.pyplot as plt
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
                # Check for NaN or invalid entries in event_data[key]
                if event_data[key].isnull().any() or (event_data[key].astype(str) == 'n/a').any():
                    design_matrix[f"{key}_{level_desc.replace(' ', '_')}"] = np.nan
                else:
                    # Convert level_key to int safely
                    try:
                        level_key_int = int(level_key)
                        design_matrix[f"{key}_{level_desc.replace(' ', '_')}"] = (event_data[key] == level_key_int).astype(int)
                    except ValueError:
                        # Handle the case where level_key cannot be converted to int
                        design_matrix[f"{key}_{level_desc.replace(' ', '_')}"] = np.nan
        else:
            # If no levels, just copy the values from the event file as-is
            if event_data[key].isnull().any() or (event_data[key].astype(str) == 'n/a').any():
                design_matrix[key] = np.nan
            else:
                design_matrix[key] = event_data[key]
        
    # Return the design matrix
    return design_matrix

# Define file list
num_run = settings.get('num_run')
file_dir_save = f'{main_dir}/{project_dir}/derivatives/exp_design/{subject}'
os.makedirs(file_dir_save, exist_ok=True)

# Load events files
if subject == 'sub-01':
    if task == 'pRF': ses = 'ses-02'
    else: ses = 'ses-01'
else: ses = settings['session']
data_events = load_event_files(main_dir, project_dir, subject, ses, task)
event_descr_json = f'{main_dir}/{project_dir}/task-{task}_events.json' 

for run in range(num_run):

    event_file = data_events[run]
    design_matrix = create_design_matrix(event_file, event_descr_json)

    # draw the design matrix and save it
    design_matrix_img_file = os.path.basename(event_file).replace('_events.tsv', '_design_matrix.png')
    plot_design_matrix(design_matrix)
    plt.savefig(os.path.join(file_dir_save, design_matrix_img_file), format='png')
    
    # Save the design matrix
    design_matrix_file = os.path.basename(event_file).replace('_events.tsv', '_design_matrix.tsv')
    design_matrix.to_csv(os.path.join(file_dir_save, design_matrix_file), sep='\t', index=False)

# Define permission cmd
print('Changing files permissions in {}/{}'.format(main_dir, project_dir))
os.system("chmod -Rf 771 {}/{}".format(main_dir, project_dir))
os.system("chgrp -Rf {} {}/{}".format(group, main_dir, project_dir))