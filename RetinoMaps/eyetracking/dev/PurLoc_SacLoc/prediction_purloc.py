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
sys.path.insert(0, "/Users/sinakling/projects/pRF_analysis/analysis_code/utils")
from sac_utils import *

# Set path to utils folder
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

with open('/Users/sinakling/projects/pRF_analysis/RetinoMaps/PurLoc_settings.json') as f:
    analysis_info = json.load(f)

main_dir = '/Users/sinakling/disks/meso_shared/'
project_dir = "RetinoMaps"
subjects, task, ses = load_inputs()

def process_subject(subject, task, ses, analysis_info, main_dir):
    file_dir_save = ensure_save_dir(f'{main_dir}/{project_dir}/derivatives/pp_data', subject)
    fig_dir_save = f'{file_dir_save}/figures'
    os.makedirs(fig_dir_save, exist_ok=True)
    
    if subject == 'sub-01': 
        ses = 'ses-01'
        data_events = load_event_files(main_dir, project_dir, subject, ses, task)
        data_mat = sorted(glob.glob(f'/Users/sinakling/projects/PredictEye/locEMexp/data/{subject}/{ses}/add/*.mat'))
    else: 
        data_events = load_event_files(main_dir, project_dir, subject, ses, task)
        data_mat = sorted(glob.glob(f'/Users/sinakling/projects/PredictEye/locEMexp/data/{subject}/{ses}/add/*.mat'))

    dfs_runs = [pd.read_csv(run, sep="\t") for run in data_events]
    
    all_run_durations = [np.cumsum(dfs['duration'] * 1000) for dfs in dfs_runs]

    precision_fraction_list = []
    precision_one_thrs_list = []

    threshold = analysis_info['threshold']
    
    for run in range(2):
        matfile = scipy.io.loadmat(data_mat[run])
        pred_x_intpl, pred_y_intpl = predicted_pursuit(
            dfs_runs[run], matfile, analysis_info['center'], analysis_info['ppd']
        )
       
        
        eye_data_run_01 = pd.read_csv(f"{file_dir_save}/timeseries/{subject}_task-{task}_run_01_eyedata.tsv.gz", compression='gzip', delimiter='\t')
        eye_data_run_02 = pd.read_csv(f"{file_dir_save}/timeseries/{subject}_task-{task}_run_02_eyedata.tsv.gz", compression='gzip', delimiter='\t')
        eye_data_all_runs = [eye_data_run_01[['timestamp', 'x', 'y', 'pupil_size']].to_numpy(), eye_data_run_02[['timestamp', 'x', 'y', 'pupil_size']].to_numpy()]
        
        # Define the start and end indices for each slice
       # slice_indices_mov_seq = [(int(all_run_durations[run][i]), int(all_run_durations[run][i+33])) for i in range(15, 161, 48)]
       # for count, (start, end) in enumerate(slice_indices_mov_seq, start=1):
        
       #     fig = plotly_layout_template("PurLoc", 0)
       ##     fig.add_trace(go.Scatter(y=eye_data_all_runs[run][start:end][:, 1], showlegend=False, line=dict(color='black', width=2)), row=1, col=1)
        #    fig.add_trace(go.Scatter(y=pred_x_intpl[start:end], showlegend=False, line=dict(color='blue', width=2)), row=1, col=1)
        #    fig.add_trace(go.Scatter(y=eye_data_all_runs[run][start:end][:, 2], showlegend=False, line=dict(color='black', width=2)), row=2, col=1)
        #    fig.add_trace(go.Scatter(y=pred_y_intpl[start:end], showlegend=False, line=dict(color='blue', width=2)), row=2, col=1)
        #    fig.add_trace(go.Scatter(x=eye_data_all_runs[run][start:end][:, 1], y=eye_data_all_runs[run][start:end][:, 2], showlegend=False, line=dict(color='black', width=2)), row=1, col=2)
        #    fig.add_trace(go.Scatter(x=pred_x_intpl[start:end], y=pred_y_intpl[start:end], showlegend=False, line=dict(color='blue', width=2)), row=1, col=2)

         #   fig_fn = f"{fig_dir_save}/{subject}_task-{task}_run-0{run+1}_{count}_prediction.pdf"
        #    print(f'Saving {fig_fn}')
        #    fig.write_image(fig_fn)
            
        eucl_dist = euclidean_distance_pur(eye_data_all_runs,pred_x_intpl, pred_y_intpl, run)
        
        precision_fraction = fraction_under_threshold(pred_x_intpl, eucl_dist)
        precision_fraction_list.append(precision_fraction)
        
        precision_file = f"{file_dir_save}/stats/precision_fraction_{subject}_run_0{run+1}.csv"
        np.savetxt(precision_file, precision_fraction, delimiter=",")

        precision_one_thrs = fraction_under_one_threshold(pred_x_intpl,eucl_dist,threshold)
        precision_one_thrs_list.append(precision_one_thrs) 

    precision_arrays = [np.array(x) for x in precision_fraction_list]
    precision_one_thrs_mean = np.mean(precision_one_thrs_list)
    
    return {
        "precision_means": [np.mean(k) for k in zip(*precision_arrays)],
        "precision_one_thrs_mean": precision_one_thrs_mean 
    }


def process_all_subjects(subjects, task, ses):

    precision_data = {
        subject: (
            result := process_subject(subject, task, ses, analysis_info, main_dir),
            {
                "precision_means": result["precision_means"],
                "precision_one_thrs_mean": result["precision_one_thrs_mean"]
            }
        )[1]
        for subject in subjects
    }
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
    generate_ranking_figure(precision_data, colormap_subject_dict)

import os
import plotly.graph_objs as go
import plotly.express as px

def generate_final_figure(precision_data, colormap, thresholds):
    traces = [
        go.Scatter(
            x=thresholds,
            y=data["precision_means"],  
            mode='lines',
            name=f'Subject {subject}',
            line=dict(color=colormap.get(subject, '#000000'))  # Default to black if subject not in colormap
        )
        for subject, data in precision_data.items()  # Access subject data dictionaries directly
    ]

    # Define layout
    layout = go.Layout(
        xaxis=dict(
            title='Euclidean distance error in dva', range=[0, 6], zeroline=True, linecolor='black', showgrid=False, tickmode='linear', dtick=2 
        ),
        yaxis=dict(
            title=r'% amount of data', range=[0, 1], zeroline=True, linecolor='black', showgrid=False
        ),
        plot_bgcolor='white',
        paper_bgcolor='white',
        font=dict(family="Arial", size=12, color="black"),
        height=700,
        width=480,  
        shapes=[dict(
            type="line", x0=2, x1=2, y0=0, y1=1, line=dict(color="black", dash="dash")
        )]
    )

    # Create figure and show it
    fig = go.Figure(data=traces, layout=layout)
    fig.show()

    # Save figure as PDF
    fig_path = "/Users/sinakling/disks/meso_shared/RetinoMaps/derivatives/pp_data/group/eyetracking"
    fig_fn = f"{fig_path}/PurLoc_threshold_precision_TEST.pdf"
    if not os.path.exists(fig_path):
        os.makedirs(fig_path)
    print(f'Saving {fig_fn}')
    fig.write_image(fig_fn)


import pandas as pd
import plotly.express as px
import os

def generate_ranking_figure(precision_data, colormap):
    # Sort the precision data by "precision_one_thrs_mean" value
    sorted_precision_data = dict(sorted(precision_data.items(), key=lambda item: item[1]["precision_one_thrs_mean"], reverse=True))
    

    # Prepare data for the plot by creating a DataFrame
    plot_data = {
        "Mean Precision under Threshold": [data["precision_one_thrs_mean"] for data in sorted_precision_data.values()],
        "Category": ["Mean Precision"] * len(sorted_precision_data),  # Single category for all data points
        "Subject": list(sorted_precision_data.keys()),  # Include sorted subject identifiers for coloring
        "Category Position": [0.05 + 0.005 * i for i in range(len(sorted_precision_data))]  # Slightly adjust positions
    }

    df = pd.DataFrame(plot_data)

    # Create a strip plot with a single category
    fig = px.strip(
        df,
        x="Category Position",  # Adjusted positions for x-axis
        y="Mean Precision under Threshold",
        stripmode="overlay",  # Use overlay mode for points
        color="Subject",  # Color by subject
        color_discrete_map=colormap,  # Apply the colormap for each subject
    )

    # Update marker size to adjust spacing
    fig.update_traces(marker=dict(size=8))  # Adjust size as needed

    # Update layout for clarity
    fig.update_layout(
        xaxis=dict(
            title='',  # No title for the x-axis
            range=[0, 0.3],
            showgrid=False,
            showticklabels=False  # Hide x-axis tick labels
        ),
        yaxis=dict(
            title='Mean Precision under Threshold',
            range=[0, 1.05],
            zeroline=True,
            linecolor='black',
            showgrid=True,
        ),
        plot_bgcolor='white',
        paper_bgcolor='white',
        font=dict(family="Arial", size=12, color="black"),
    )

    fig.show()

    # Save figure as PDF
    fig_path = "/Users/sinakling/disks/meso_shared/RetinoMaps/derivatives/pp_data/group/eyetracking"
    fig_fn = f"{fig_path}/PurLoc_threshold_2_ranking_TEST.pdf"
    if not os.path.exists(fig_path):
        os.makedirs(fig_path)
    print(f'Saving {fig_fn}')
    fig.write_image(fig_fn)








    

process_all_subjects(subjects, task, ses)
