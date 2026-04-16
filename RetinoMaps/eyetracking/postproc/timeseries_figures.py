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
prf_task_name = analysis_info['prf_task_names'][0]

# Execption for subject 1 with no data for eye tracking
if subject == 'sub-01':
    if prf_task_name in tasks:
        tasks.remove(prf_task_name)

for task in tasks :
    print('Processing {} ...'.format(task))

    # Load task settings
    base_dir = os.path.abspath(os.path.join(os.getcwd(), "../"))
    task_settings_path = os.path.join(base_dir, "eye-tracking_{}.yml".format(task))
    task_settings = load_settings([task_settings_path, eye_tracking_settings_path, settings_path])[0] 
    
    if subject == 'sub-01':
        if task == 'pRF': ses = 'ses-02'
        else: ses = 'ses-01'
    else: ses = task_settings['session']
        
    # Load main experiment settings 
    eye = task_settings['eye']
    num_run = task_settings['num_run']
    eyetracking_sampling = task_settings['eyetrack_sampling']
    ppd = task_settings['ppd']
    
    # Defind directories 
    eye_tracking_dir = '{}/{}/derivatives/pp_data/{}/eyetracking'.format(main_dir, project_dir, subject)
    fig_dir = '{}/figures'.format(eye_tracking_dir)
    os.makedirs(fig_dir, exist_ok=True)

    # Load event files
    try:
        data_events = load_event_files(main_dir, project_dir, subject, ses, task)
        dfs_runs = [pd.read_csv(run, sep="\t") for run in data_events]
    except Exception as e:
        print(f"Error loading or processing event files for {subject}, {ses}, {task}: {e}")
        continue

    dfs_runs = [pd.read_csv(run, sep="\t") for run in data_events]
    all_run_durations = [np.cumsum(dfs['duration'] * 1000) for dfs in dfs_runs]

    precision_all_runs = []
    precision_one_thrs_list = []

    threshold = task_settings['threshold']

    for run in range(num_run):
            
        prediction = pd.read_csv(f"{eye_tracking_dir}/timeseries/{subject}_task-{task}_run_0{run+1}_prediction.tsv.gz", compression='gzip', delimiter='\t')
        pred_x_intpl =  prediction['pred_x']
        pred_y_intpl =  prediction['pred_y']
                
        # load eye data 
        if task != prf_task_name :
            eye_data_run_01 = pd.read_csv(f"{eye_tracking_dir}/timeseries/{subject}_task-{task}_run_01_eyedata.tsv.gz", compression='gzip', delimiter='\t')
            eye_data_run_02 = pd.read_csv(f"{eye_tracking_dir}/timeseries/{subject}_task-{task}_run_02_eyedata.tsv.gz", compression='gzip', delimiter='\t')
            eye_data_all_runs = [eye_data_run_01[['timestamp', 'x', 'y', 'pupil_size']].to_numpy(), eye_data_run_02[['timestamp', 'x', 'y', 'pupil_size']].to_numpy()]
        else : 
            eye_data_run_01 = pd.read_csv(f"{eye_tracking_dir}/timeseries/{subject}_task-{task}_run_01_eyedata.tsv.gz", compression='gzip', delimiter='\t')
            eye_data_run_02 = pd.read_csv(f"{eye_tracking_dir}/timeseries/{subject}_task-{task}_run_02_eyedata.tsv.gz", compression='gzip', delimiter='\t')
            eye_data_run_03 = pd.read_csv(f"{eye_tracking_dir}/timeseries/{subject}_task-{task}_run_03_eyedata.tsv.gz", compression='gzip', delimiter='\t')
            eye_data_run_04 = pd.read_csv(f"{eye_tracking_dir}/timeseries/{subject}_task-{task}_run_04_eyedata.tsv.gz", compression='gzip', delimiter='\t')
            eye_data_run_05 = pd.read_csv(f"{eye_tracking_dir}/timeseries/{subject}_task-{task}_run_05_eyedata.tsv.gz", compression='gzip', delimiter='\t')
            
            eye_data_all_runs = [
                eye_data_run_01[['timestamp', 'x', 'y', 'pupil_size']].to_numpy(),
                eye_data_run_02[['timestamp', 'x', 'y', 'pupil_size']].to_numpy(),
                eye_data_run_03[['timestamp', 'x', 'y', 'pupil_size']].to_numpy(),
                eye_data_run_04[['timestamp', 'x', 'y', 'pupil_size']].to_numpy(),
                eye_data_run_05[['timestamp', 'x', 'y', 'pupil_size']].to_numpy()
            ]
            
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
            fig_fn = f"{fig_dir}/{subject}_task-{task}_run-0{run+1}_{count}_prediction.pdf"
            print(f'Saving {fig_fn}')
            fig.write_image(fig_fn)



# # Define permission cmd
# print('Changing files permissions in {}/{}'.format(main_dir, project_dir))
# os.system("chmod -Rf 771 {}/{}".format(main_dir, project_dir))
# os.system("chgrp -Rf {} {}/{}".format(group, main_dir, project_dir))
    
    