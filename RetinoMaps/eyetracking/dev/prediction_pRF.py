import numpy as np
import math
import glob
import pandas as pd
import scipy.io
import json
import matplotlib.pyplot as plt
import h5py
import os
import sys
import plotly.graph_objects as go
import statistics
sys.path.insert(0, "/Users/sinakling/projects/pRF_analysis/RetinoMaps/eyetracking/dev/PurLoc_SacLoc")
from sac_utils import *

# Set path to utils folder
sys.path.insert(0, "/Users/sinakling/projects/pRF_analysis/analysis_code/utils")
from eyetrack_utils import *



def load_inputs():
    subjects_input = sys.argv[1]
    subjects = subjects_input.strip('[]').split(',')
    subjects = [sub.strip() for sub in subjects]
    task = sys.argv[2]
    ses = sys.argv[3]
    return subjects, task, ses

def interp1d(array: np.ndarray, new_len: int) -> np.ndarray:
    la = len(array)
    return np.interp(np.linspace(0, la - 1, num=new_len), np.arange(la), array)

def ensure_save_dir(base_dir, subject):
    save_dir = f"{base_dir}/{subject}/eyetracking"
    os.makedirs(save_dir, exist_ok=True)
    return save_dir

with open('/Users/sinakling/projects/pRF_analysis/RetinoMaps/eyetracking/dev/pRF_behavior_settings.json') as f:
    analysis_info = json.load(f)

main_dir = analysis_info['main_dir_mac']
subjects, task, ses = load_inputs()

def process_subject(subject, task, ses, analysis_info, main_dir):
    file_dir_save = ensure_save_dir(f'{main_dir}/derivatives/pp_data', subject)
    fig_dir_save = f'{file_dir_save}/figures'
    os.makedirs(fig_dir_save, exist_ok=True)
   # h5_filename = f'{file_dir_save}/stats/{subject}_task-{task}_eyedata_sac_stats.h5'
    
   # with h5py.File(h5_filename, 'r') as h5_file:
   #     time_start_trial = np.array(h5_file['time_start_trial'])
   #     time_end_trial = np.array(h5_file['time_end_trial'])
   #     time_start_seq = np.array(h5_file['time_start_seq'])
    #    time_end_seq = np.array(h5_file['time_end_seq'])
  #      time_start_eye = np.array(h5_file['time_start_eye'])
  #      time_end_eye = np.array(h5_file['time_end_eye'])
    
    
    if subject == 'sub-01': 
        ses = 'ses-01'
        data_events = load_event_files(main_dir, subject, ses, task)
        #data_mat = sorted(glob.glob(f'/Users/sinakling/projects/PredictEye/locEMexp/data/{subject}/{ses}/add/*.mat'))
    else: 
        data_events = load_event_files(main_dir, subject, ses, task)
       # data_mat = sorted(glob.glob(f'/Users/sinakling/projects/PredictEye/locEMexp/data/{subject}/{ses}/add/*.mat'))

    dfs_runs = [pd.read_csv(run, sep="\t") for run in data_events]
    
    all_run_durations = [np.cumsum(dfs['duration'] * 1000) for dfs in dfs_runs]
    precision_fraction_list = []
    
    for run in range(5):
        #matfile = scipy.io.loadmat(data_mat[run])
        
        eye_data_run_01 = pd.read_csv(f"{file_dir_save}/timeseries/{subject}_task-{task}_run_01_eyedata.tsv.gz", compression='gzip', delimiter='\t')
        eye_data_run_02 = pd.read_csv(f"{file_dir_save}/timeseries/{subject}_task-{task}_run_02_eyedata.tsv.gz", compression='gzip', delimiter='\t')
        eye_data_run_03 = pd.read_csv(f"{file_dir_save}/timeseries/{subject}_task-{task}_run_03_eyedata.tsv.gz", compression='gzip', delimiter='\t')
        eye_data_run_04 = pd.read_csv(f"{file_dir_save}/timeseries/{subject}_task-{task}_run_04_eyedata.tsv.gz", compression='gzip', delimiter='\t')
        eye_data_run_05 = pd.read_csv(f"{file_dir_save}/timeseries/{subject}_task-{task}_run_05_eyedata.tsv.gz", compression='gzip', delimiter='\t')
        
        eye_data_all_runs = [eye_data_run_01[['timestamp', 'x', 'y', 'pupil_size']].to_numpy(), eye_data_run_02[['timestamp', 'x', 'y', 'pupil_size']].to_numpy(),eye_data_run_03[['timestamp', 'x', 'y', 'pupil_size']].to_numpy(),eye_data_run_04[['timestamp', 'x', 'y', 'pupil_size']].to_numpy(),eye_data_run_05[['timestamp', 'x', 'y', 'pupil_size']].to_numpy()]

        pred_x_intpl = np.zeros(len(eye_data_all_runs[run]))
        pred_y_intpl = np.zeros(len(eye_data_all_runs[run]))
        
        # Define the start and end indices for each slice
        slice_indices_mov_seq = [(int(all_run_durations[run][i]), int(all_run_durations[run][i+33])) for i in range(15, 161, 48)]
        for count, (start, end) in enumerate(slice_indices_mov_seq, start=1):
        
            fig = plotly_layout_template("pRF", 0)
            fig.add_trace(go.Scatter(y=eye_data_all_runs[run][start:end][:, 1], showlegend=False, line=dict(color='black', width=2)), row=1, col=1)
            fig.add_trace(go.Scatter(y=pred_x_intpl[start:end], showlegend=False, line=dict(color='blue', width=2)), row=1, col=1)
            fig.add_trace(go.Scatter(y=eye_data_all_runs[run][start:end][:, 2], showlegend=False, line=dict(color='black', width=2)), row=2, col=1)
            fig.add_trace(go.Scatter(y=pred_y_intpl[start:end], showlegend=False, line=dict(color='blue', width=2)), row=2, col=1)
            fig.add_trace(go.Scatter(x=eye_data_all_runs[run][start:end][:, 1], y=eye_data_all_runs[run][start:end][:, 2], showlegend=False, line=dict(color='black', width=2)), row=1, col=2)
            fig.add_trace(go.Scatter(x=pred_x_intpl[start:end], y=pred_y_intpl[start:end], showlegend=False, line=dict(color='blue', width=2)), row=1, col=2)

            fig_fn = f"{fig_dir_save}/{subject}_task-{task}_run-0{run+1}_{count}_prediction.pdf"
            print(f'Saving {fig_fn}')
            fig.write_image(fig_fn)

        #TODO make function in utils   
        eucl_dist = np.sqrt((eye_data_all_runs[run][:, 1] -  pred_x_intpl) ** 2 +
                            (eye_data_all_runs[run][:, 2] -  pred_y_intpl) ** 2)
        
        precision_fraction = fraction_under_threshold(pred_x_intpl, eucl_dist)
        precision_fraction_list.append(precision_fraction)
        
        precision_file = f"{file_dir_save}/stats/precision_fraction_{subject}_run_0{run+1}.csv"
        np.savetxt(precision_file, precision_fraction, delimiter=",")

    precision_arrays = [np.array(x) for x in precision_fraction_list]
    
    return [np.mean(k) for k in zip(*precision_arrays)]

def process_all_subjects(subjects, task, ses):
    precision_data = {subject: process_subject(subject, task, ses, analysis_info, main_dir) for subject in subjects}
    colormap_subject_dict = {
    'sub-01': '#AA0DFE',
    'sub-02': '#3283FE',
    'sub-03': '#85660D',
    'sub-04': '#782AB6',
    'sub-05': '#565656',
    'sub-06': '#1C8356',
    'sub-07': '#16FF32',
    'sub-08': '#F7E1A0',
    'sub-09': '#E2E2E2',
    'sub-11': '#1CBE4F',
    'sub-12': '#C4451C',
    'sub-13': '#DEA0FD',
    'sub-14': '#FBE426',
    'sub-15': '#325A9B',
    'sub-16': '#FEAF16',
    'sub-17': '#F8A19F',
    'sub-18': '#90AD1C',
    'sub-20': '#F6222E',
    'sub-21': '#1CFFCE',
    'sub-22': '#2ED9FF',
    'sub-23': '#B10DA1',
    'sub-24': '#C075A6',
    'sub-25': '#FC1CBF'}
    generate_final_figure(precision_data, colormap_subject_dict, thresholds=np.linspace(0, 9, 100))

import os
import plotly.graph_objs as go
import plotly.express as px

import plotly.graph_objs as go
import os

def generate_final_figure(precision_data, colormap,thresholds):
    # Generate traces for each subject
    traces = [
        go.Scatter(
            x=thresholds,
            y=precision,
            mode='lines',
            name=f'Subject {subject}',
            line=dict(color=colormap.get(subject, '#000000'))  # Default to black if subject not in colormap
        )
        for subject, precision in precision_data.items()  # Use subject IDs directly for color lookup
    ]

    # Define layout
    layout = go.Layout(
        xaxis=dict(
            title='Euclidean distance error in dva', range=[0, 3], zeroline=True, linecolor='black', showgrid=False, tickmode='linear', dtick=2 
        ),
        yaxis=dict(
            title=r'% ammount of data', range=[0, 1], zeroline=True, linecolor='black', showgrid=False
        ),
        plot_bgcolor='white',
        paper_bgcolor='white',
        font=dict(family="Arial", size=12, color="black"),
        shapes=[dict(
            type="line", x0=1, x1=1, y0=0, y1=1, line=dict(color="black", dash="dash")
        )]
    )

    # Create figure and show it
    fig = go.Figure(data=traces, layout=layout)
    fig.show()

    # Save figure as PDF
    fig_path = f"/Users/sinakling/disks/meso_shared/RetinoMaps/derivatives/pp_data/sub-all"
    fig_fn = f"{fig_path}/pRF_threshold_precision.pdf"
    if not os.path.exists(fig_path):
        os.makedirs(fig_path)
    print(f'Saving {fig_fn}')
    fig.write_image(fig_fn)



process_all_subjects(subjects, task, ses)
