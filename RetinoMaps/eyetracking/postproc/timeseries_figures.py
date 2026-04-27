"""
-----------------------------------------------------------------------------------------
timeseries_figures.py
-----------------------------------------------------------------------------------------
Goal of the script:
Create figures for each subject showing the eyetraces vs the prediction over time 
per sequence
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject name
sys.argv[4]: group of shared data (e.g. 327)
-----------------------------------------------------------------------------------------
Output(s):
{subject}_task-{task}_run-01_1_prediction.pdf
{subject}_task-{task}_run-01_2_prediction.pdf
{subject}_task-{task}_run-01_3_prediction.pdf
{subject}_task-{task}_run-01_4_prediction.pdf
-----------------------------------------------------------------------------------------
To run:
cd ~/projects/pRF_analysis/RetinoMaps/eyetracking/postproc
python timeseries_figures.py /scratch/mszinte/data RetinoMaps sub-01 327
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
import os
import sys
import numpy as np
import pandas as pd
import plotly.graph_objects as go

# Personal imports 
sys.path.append("{}/../../../analysis_code/utils".format(os.getcwd()))
from settings_utils import load_settings
from eyetrack_utils import load_event_files
from sac_utils import plotly_layout_template

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
group = sys.argv[4]

# Load general settings
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../"))
eye_tracking_settings_path = os.path.join(base_dir, project_dir, "eyetracking-analysis.yml")
settings_path = os.path.join(base_dir, project_dir, "settings.yml")
settings = load_settings([eye_tracking_settings_path, settings_path])
analysis_info = settings[0]

tasks = analysis_info['eye-tracking_task_names']

# Exception for subject 1 with no pRF eye tracking data
prf_task_name = analysis_info['prf_task_names'][0]
if subject == 'sub-01' and prf_task_name in tasks:
    tasks.remove(prf_task_name)

for task in tasks:
    print('Processing {} ...'.format(task))

    # Load task settings
    base_dir = os.path.abspath(os.path.join(os.getcwd(), "../"))
    task_settings_path = os.path.join(base_dir, "eye-tracking_{}.yml".format(task))
    task_settings = load_settings([task_settings_path, eye_tracking_settings_path, settings_path])[0]

    if subject == 'sub-01':
        if task == 'pRF': ses = 'ses-02'
        else: ses = 'ses-01'
    else:
        ses = task_settings['session']

    # Load main experiment settings
    num_run = task_settings['num_run']
    eyetracking_sampling = task_settings['eyetrack_sampling']

    # Define directories
    eye_tracking_dir = '{}/{}/derivatives/pp_data/{}/eyetracking'.format(main_dir, project_dir, subject)
    fig_dir = '{}/figures'.format(eye_tracking_dir)
    os.makedirs(fig_dir, exist_ok=True)

    # For tasks with event files: compute per-sequence slice indices from event durations.
    # For the rest task: no event files exist — the full run is treated as one sequence.
    is_rest = (task == 'rest')

    if not is_rest:
        try:
            data_events = load_event_files(main_dir, project_dir, subject, ses, task)
            dfs_runs = [pd.read_csv(run, sep="\t") for run in data_events]
        except Exception as e:
            print(f"Error loading or processing event files for {subject}, {ses}, {task}: {e}")
            continue
        all_run_durations = [np.cumsum(dfs['duration'] * 1000) for dfs in dfs_runs]

    for run in range(num_run):

        # Load prediction
        prediction = pd.read_csv(
            f"{eye_tracking_dir}/timeseries/{subject}_task-{task}_run_0{run+1}_prediction.tsv.gz",
            compression='gzip', delimiter='\t'
        )
        pred_x_intpl = prediction['pred_x']
        pred_y_intpl = prediction['pred_y']

        # Load eye data
        eye_data = pd.read_csv(
            f"{eye_tracking_dir}/timeseries/{subject}_task-{task}_run_0{run+1}_eyedata.tsv.gz",
            compression='gzip', delimiter='\t'
        )
        eye_data_np = eye_data[['timestamp', 'x', 'y', 'pupil_size']].to_numpy()

        # Define slice indices:
        # - rest task: single slice covering the full run
        # - other tasks: one slice per movement sequence, derived from event durations
        if is_rest:
            slice_indices = [(0, len(eye_data_np))]
        else:
            slice_indices = [
                (int(all_run_durations[run][i]), int(all_run_durations[run][i + 33]))
                for i in range(15, 161, 48)
            ]

        for count, (start, end) in enumerate(slice_indices, start=1):
            n_samples = end - start
            time_seconds = np.arange(n_samples) / eyetracking_sampling
            duration_seconds = time_seconds[-1]

            fig = plotly_layout_template(task, run, duration_seconds)

            fig.add_trace(go.Scatter(x=time_seconds, y=eye_data_np[start:end][:, 1], showlegend=False, line=dict(color='black', width=2)), row=1, col=1)
            fig.add_trace(go.Scatter(x=time_seconds, y=pred_x_intpl[start:end].values, showlegend=False, line=dict(color='blue', width=2)), row=1, col=1)
            fig.add_trace(go.Scatter(x=time_seconds, y=eye_data_np[start:end][:, 2], showlegend=False, line=dict(color='black', width=2)), row=2, col=1)
            fig.add_trace(go.Scatter(x=time_seconds, y=pred_y_intpl[start:end].values, showlegend=False, line=dict(color='blue', width=2)), row=2, col=1)
            fig.add_trace(go.Scatter(x=eye_data_np[start:end][:, 1], y=eye_data_np[start:end][:, 2], showlegend=False, line=dict(color='black', width=2)), row=1, col=2)
            fig.add_trace(go.Scatter(x=pred_x_intpl[start:end].values, y=pred_y_intpl[start:end].values, showlegend=False, line=dict(color='blue', width=2)), row=1, col=2)

            fig_fn = f"{fig_dir}/{subject}_task-{task}_run-0{run+1}_{count}_prediction.pdf"
            print(f'Saving {fig_fn}')
            fig.write_image(fig_fn)

# # Define permission cmd
# print('Changing files permissions in {}/{}'.format(main_dir, project_dir))
# os.system("chmod -Rf 771 {}/{}".format(main_dir, project_dir))
# os.system("chgrp -Rf {} {}/{}".format(group, main_dir, project_dir))