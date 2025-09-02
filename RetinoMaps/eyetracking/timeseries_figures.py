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
cd ~/projects/pRF_analysis/RetinoMaps/eyetracking/
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
import json
import numpy as np
import pandas as pd
import plotly.graph_objects as go

# Personal imports 
sys.path.append("{}/../../analysis_code/utils".format(os.getcwd()))
from eyetrack_utils import load_event_files
from sac_utils import plotly_layout_template

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
group = sys.argv[4]

# Load general settings
with open('../settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)
tasks = analysis_info['task_intertask'][0]
prf_task_name = analysis_info['prf_task_name']

# Execption for subject 1 with no data for eye tracking
if subject == 'sub-01':
    if prf_task_name in tasks:
        tasks.remove(prf_task_name)

for task in tasks :
    print('Processing {} ...'.format(task))

    # Load task settings
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
    
    if subject == 'sub-01':
        ses = 'ses-01'
    else: ses = settings['session']

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

    threshold = settings['threshold']


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
            try :
                fig.add_trace(go.Scatter(y=eye_data_all_runs[run][start:end][:, 1], showlegend=False, line=dict(color='black', width=2)), row=1, col=1)
            except Exception:
                deb()
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
    
    