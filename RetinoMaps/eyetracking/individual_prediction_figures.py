"""
-----------------------------------------------------------------------------------------
individual_prediction_figures.py
-----------------------------------------------------------------------------------------
Goal of the script:
Create figures for each subject showing the eyetraces vs the prediction over time 
per sequence
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main directory 
sys.argv[2]: project directory 
sys.argv[3]: subject 
sys.argv[4]: task 
sys.argv[5]: group 
-----------------------------------------------------------------------------------------
Output(s):
{subject}_task-{task}_run-01_1_prediction.pdf
{subject}_task-{task}_run-01_2_prediction.pdf
{subject}_task-{task}_run-01_3_prediction.pdf
{subject}_task-{task}_run-01_4_prediction.pdf
-----------------------------------------------------------------------------------------
To run:
cd /projects/pRF_analysis/RetinoMaps/eyetracking/
python individual_prediction_figures.py /scratch/mszinte/data RetinoMaps "[sub-01,sub-02,sub-03]" pRF 327
-----------------------------------------------------------------------------------------
"""
import pandas as pd
import json
import numpy as np
import re
import matplotlib.pyplot as plt
import glob 
import os
import sys
import math 
import h5py
import scipy.io 
import plotly.graph_objects as go

# path of utils folder  
sys.path.append("{}/../../analysis_code/utils".format(os.getcwd()))
from eyetrack_utils import load_event_files
from sac_utils import plotly_layout_template

# --------------------- Load settings and inputs -------------------------------------

def load_inputs():
    main_dir = sys.argv[1]
    project_dir = sys.argv[2]
    subject = sys.argv[3]
    task = sys.argv[4]
    group = sys.argv[5]

    return main_dir, project_dir, subject, task, group 

def load_events(main_dir, subject, ses, task): 
    data_events = load_event_files(main_dir, subject, ses, task)
    return data_events 


def ensure_save_dir(base_dir, subject):
    save_dir = f"{base_dir}/{subject}/eyetracking"
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    return save_dir

# Load inputs and settings
main_dir, project_dir, subject, task, group = load_inputs()


base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../"))
settings_path = os.path.join(base_dir, project_dir, f'{task}_settings.json')
with open(settings_path) as f:
    settings = json.load(f)


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


#for subject in subjects: 
    
file_dir_save = ensure_save_dir(f'{main_dir}/{project_dir}/derivatives/pp_data', subject)
fig_dir_save = f'{file_dir_save}/figures'
os.makedirs(fig_dir_save, exist_ok=True)

if subject == 'sub-01':
    if task == 'pRF': 
        ses = 'ses-02'
        data_events = load_event_files(main_dir, project_dir, subject, ses, task)
else: 
    if task == 'pRF': 
        ses = 'ses-01'
        data_events = load_event_files(main_dir, project_dir, subject, ses, task)
    else: 
        ses = 'ses-02'
        data_events = load_event_files(main_dir, project_dir, subject, ses, task)


dfs_runs = [pd.read_csv(run, sep="\t") for run in data_events]
all_run_durations = [np.cumsum(dfs['duration'] * 1000) for dfs in dfs_runs]


precision_all_runs = []
precision_one_thrs_list = []

threshold = settings['threshold']


for run in range(num_run):
        
        prediction = pd.read_csv(f"{file_dir_save}/timeseries/{subject}_task-{task}_run_0{run+1}_prediction.tsv.gz", compression='gzip', delimiter='\t')
        pred_x_intpl =  prediction['pred_x']
        pred_y_intpl =  prediction['pred_y']
        
        
        # load eye data 
        eye_data_run_01 = pd.read_csv(f"{file_dir_save}/timeseries/{subject}_task-{task}_run_01_eyedata.tsv.gz", compression='gzip', delimiter='\t')
        eye_data_run_02 = pd.read_csv(f"{file_dir_save}/timeseries/{subject}_task-{task}_run_02_eyedata.tsv.gz", compression='gzip', delimiter='\t')
        eye_data_all_runs = [eye_data_run_01[['timestamp', 'x', 'y', 'pupil_size']].to_numpy(), eye_data_run_02[['timestamp', 'x', 'y', 'pupil_size']].to_numpy()]
        
        # Define the start and end indices for each slice
        slice_indices_mov_seq = [(int(all_run_durations[run][i]), int(all_run_durations[run][i+33])) for i in range(15, 161, 48)]
        for count, (start, end) in enumerate(slice_indices_mov_seq, start=1):
        
            fig = plotly_layout_template(task, run)
            fig.add_trace(go.Scatter(y=eye_data_all_runs[run][start:end][:, 1], showlegend=False, line=dict(color='black', width=2)), row=1, col=1)
            fig.add_trace(go.Scatter(y=pred_x_intpl[start:end], showlegend=False, line=dict(color='blue', width=2)), row=1, col=1)
            fig.add_trace(go.Scatter(y=eye_data_all_runs[run][start:end][:, 2], showlegend=False, line=dict(color='black', width=2)), row=2, col=1)
            fig.add_trace(go.Scatter(y=pred_y_intpl[start:end], showlegend=False, line=dict(color='blue', width=2)), row=2, col=1)
            fig.add_trace(go.Scatter(x=eye_data_all_runs[run][start:end][:, 1], y=eye_data_all_runs[run][start:end][:, 2], showlegend=False, line=dict(color='black', width=2)), row=1, col=2)
            fig.add_trace(go.Scatter(x=pred_x_intpl[start:end], y=pred_y_intpl[start:end], showlegend=False, line=dict(color='blue', width=2)), row=1, col=2)

            fig_fn = f"{fig_dir_save}/{subject}_task-{task}_run-0{run+1}_{count}_prediction.pdf"
            print(f'Saving {fig_fn}')
            #fig.write_image(fig_fn)
            fig.show()


# Define permission cmd
print('Changing files permissions in {}/{}'.format(main_dir, project_dir))
os.system("chmod -Rf 771 {}/{}".format(main_dir, project_dir))
os.system("chgrp -Rf {} {}/{}".format(group, main_dir, project_dir))