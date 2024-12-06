"""
-----------------------------------------------------------------------------------------
stats_figures.py
-----------------------------------------------------------------------------------------
Goal of the script:
Create figures for all subjects together showing the percentage of amount of data of the 
euclidean error under each threshold (precision) as well as under one specific 
threshold 
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main directory 
sys.argv[2]: project directory 
sys.argv[3]: subject 
sys.argv[4]: task 
sys.argv[5]: group 
-----------------------------------------------------------------------------------------
Output(s):
{subject}_{task}_threshold_precision.pdf
group_{task}_threshold_precision.pdf
{subject}_{task}_threshold_ranking.pdf
group_{task}_threshold_ranking.pdf
{subject}_{task}_stats_figure.pdf
group_{task}_stats_figure.pdf
-----------------------------------------------------------------------------------------
To run:
cd ~/projects/pRF_analysis/RetinoMaps/eyetracking/
python stats_figures.py /scratch/mszinte/data RetinoMaps sub-02 pRF 327
python stats_figures.py /scratch/mszinte/data RetinoMaps group pRF 327
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
import plotly.express as px
from plotly.subplots import make_subplots

# path of utils folder  
sys.path.append("{}/../../analysis_code/utils".format(os.getcwd()))
from eyetrack_utils import load_event_files
from sac_utils import plotly_layout_template

# --------------------- Load settings and inputs -------------------------------------

def load_inputs():
    main_dir = sys.argv[1]
    project_dir = sys.argv[2]
    subject = sys.argv[3]
    #subjects = [sub.strip().strip('[]"\'') for sub in subjects_input.split(',')]
    task = sys.argv[4]
    group = sys.argv[5]

    return main_dir, project_dir, subject, task, group 

def load_events(main_dir, subject, ses, task): 
    data_events = load_event_files(main_dir, subject, ses, task)
    return data_events 


def ensure_save_dir(base_dir, subject, is_group=False):
    if is_group:
        save_dir = f"{base_dir}/group/eyetracking"
    else:
        save_dir = f"{base_dir}/{subject}/eyetracking"
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    return save_dir

# Load inputs and settings
main_dir, project_dir, subject_input, task, group = load_inputs()
print(subject_input)

base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../"))
task_settings_path = os.path.join(base_dir, project_dir, f'{task}_settings.json')
with open(task_settings_path) as f:
    settings = json.load(f)

general_settings_path = os.path.join(base_dir, project_dir, 'settings.json')
with open(general_settings_path) as f:
    general_settings = json.load(f)


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

threshold = settings['threshold']


def process_subject(main_dir, project_dir, subject, task): 
    file_dir_save = ensure_save_dir(f'{main_dir}/{project_dir}/derivatives/pp_data', subject)

    precision_summary = pd.read_csv(f"{file_dir_save}/stats/{subject}_task-{task}_precision_summary.tsv", delimiter='\t')
    precision_one_summary = pd.read_csv(f"{file_dir_save}/stats/{subject}_task-{task}_precision_one_threshold_summary.tsv", delimiter='\t')
     
    return {
        "precision_mean": precision_summary["precision_mean"],
        "precision_one_thrs_mean": precision_one_summary["precision_one_thrs_mean"].item()
    }

    
def process_all_subjects(main_dir, project_dir, subjects, subject_input, task, threshold):

    precision_data = {
        subject: (
            result := process_subject(main_dir, project_dir, subject, task),
            {
                "precision_mean": result["precision_mean"],
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

    generate_final_figure(precision_data, colormap_subject_dict, threshold, main_dir, project_dir, subject_input, task, thresholds=np.linspace(0, 9, 100))
    generate_ranking_figure(precision_data, colormap_subject_dict)

    # combined figure function
    generate_combined_figure(precision_data, colormap_subject_dict, threshold, np.linspace(0, 9, 100))



def generate_final_figure(precision_data, colormap, threshold, main_dir, project_dir, subject_input, task, thresholds):
    traces = [
        go.Scatter(
            x=thresholds,
            y=data["precision_mean"],  
            mode='lines',
            name=f'Subject {subject}',
            line=dict(color=colormap.get(subject, '#000000'))  # Default to black if subject not in colormap
        )
        for subject, data in precision_data.items()  
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
            type="line", x0=threshold, x1=threshold, y0=0, y1=1, line=dict(color="black", dash="dash")
        )]
    )

    # Create figure and show it
    fig = go.Figure(data=traces, layout=layout)
    fig.show()

    # Save figure dynamically
    is_group = subject_input == "group"
    save_dir = ensure_save_dir(f"{main_dir}/{project_dir}/derivatives/pp_data", subject_input, is_group)
    fig_fn = f"{save_dir}/{subject_input}_{task}_threshold_precision.pdf" if not is_group else f"{save_dir}/group_{task}_threshold_precision.pdf"

    print(f"Saving {fig_fn}")
    fig.write_image(fig_fn)


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
    is_group = subject_input == "group"
    save_dir = ensure_save_dir(f"{main_dir}/{project_dir}/derivatives/pp_data", subject_input, is_group)
    fig_fn = f"{save_dir}/{subject_input}_{task}_threshold_ranking.pdf" if not is_group else f"{save_dir}/group_{task}_threshold_ranking.pdf"

    print(f"Saving {fig_fn}")
    fig.write_image(fig_fn)



def generate_combined_figure(precision_data, colormap, threshold, thresholds):
    # Create traces for the general figure
    general_traces = [
        go.Scatter(
            x=thresholds,
            y=data["precision_mean"],  
            mode='lines',
            name=f'Subject {subject}',
            line=dict(color=colormap.get(subject, '#000000'))  # Default to black if subject not in colormap
        )
        for subject, data in precision_data.items()
    ]

    # Sort the precision data for the ranking figure
    sorted_precision_data = dict(sorted(precision_data.items(), key=lambda item: item[1]["precision_one_thrs_mean"], reverse=True))
    ranking_traces = [
        go.Scatter(
            x=[0.05 + 0.005 * i],  # Adjusted positions for x-axis
            y=[data["precision_one_thrs_mean"]],
            mode='markers',
            name=f'Subject {subject}',
            marker=dict(size=8, color=colormap.get(subject, '#000000'))  # Default to black if subject not in colormap
        )
        for i, (subject, data) in enumerate(sorted_precision_data.items())
    ]

    # Create a subplot figure
    fig = make_subplots(
        rows=1, cols=2,
        subplot_titles=("General Precision", "Ranking"),
        column_widths=[0.5, 0.5], 
        shared_yaxes=True 
    )

    # Add general figure traces to the first subplot
    for trace in general_traces:
        fig.add_trace(trace, row=1, col=1)

    # Add ranking figure traces to the second subplot
    for trace in ranking_traces:
        fig.add_trace(trace, row=1, col=2)

    # Add vertical dashed line to the general figure subplot
    fig.add_shape(
        type="line",
        x0=threshold, x1=threshold, y0=0, y1=1,  # Line spans the y-axis
        line=dict(color="black", dash="dash"),
        xref="x1",  # Reference to the x-axis of the first subplot
        yref="y1"   # Reference to the y-axis of the first subplot
    )

    # Update layout for both subplots
    fig.update_layout(
       height=700,
        width=960,
        plot_bgcolor='white',
        paper_bgcolor='white',
        font=dict(family="Arial", size=12, color="black"),
        showlegend=False 
    )

    # Update axes for the first subplot
    fig.update_xaxes(
        title='Euclidean distance error in dva', range=[0, 6], row=1, col=1,
        zeroline=True, linecolor='black', showgrid=False, tickmode='linear', dtick=2
    )
    fig.update_yaxes(
        title='% amount of data', range=[0, 1.05], row=1, col=1,
        zeroline=True, linecolor='black', showgrid=False
    )

    # Update axes for the second subplot
    fig.update_xaxes(
        title='', range=[0, 0.3], row=1, col=2,
        showgrid=False, showticklabels=False
    )
    fig.update_yaxes(
        range=[0, 1.05], row=1, col=2,
        zeroline=True, linecolor='black', showgrid=True
    )

    fig.show()

    
    # Save figure dynamically
    is_group = subject_input == "group"
    save_dir = ensure_save_dir(f"{main_dir}/{project_dir}/derivatives/pp_data", subject_input, is_group)
    fig_fn = f"{save_dir}/{subject_input}_{task}_stats_figure.pdf" if not is_group else f"{save_dir}/group_{task}_stats_figure.pdf"

    print(f"Saving {fig_fn}")
    fig.write_image(fig_fn)
    
    


# Check if the input is "group" or a specific subject
if subject_input == "group":
    subjects = general_settings["subjects"]
else:
    subjects = [subject_input]


process_all_subjects(main_dir, project_dir, subjects, subject_input, task, threshold)