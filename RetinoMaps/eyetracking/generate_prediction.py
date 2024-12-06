"""
-----------------------------------------------------------------------------------------
generate_prediction.py
-----------------------------------------------------------------------------------------
Goal of the script:
Generate prediction for eyemovements and calculate euclidean distance 
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main directory 
sys.argv[2]: project directory 
sys.argv[3]: subject 
sys.argv[4]: task 
sys.argv[5]: group 
-----------------------------------------------------------------------------------------
Output(s):
tsv of fraction under thresholds
tsv.gz timeseries of Euclidean distance 
tsv.gz timeseries of Prediction
-----------------------------------------------------------------------------------------
To run:
cd /projects/pRF_analysis/RetinoMaps/eyetracking/
python generate_prediction.py /scratch/mszinte/data RetinoMaps sub-01 pRF 327
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
from sac_utils import predicted_pursuit, euclidean_distance_pur, fraction_under_threshold, fraction_under_one_threshold, load_sac_model, euclidean_distance

# --------------------- Load settings and inputs -------------------------------------

def load_settings(settings_file):
    with open(settings_file) as f:
        settings = json.load(f)
    return settings

def load_events(main_dir, subject, ses, task): 
    data_events = load_event_files(main_dir, subject, ses, task)
    return data_events 

def load_inputs():
    return sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5]

def ensure_save_dir(base_dir, subject):
    save_dir = f"{base_dir}/{subject}/eyetracking"
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    return save_dir

# Load inputs and setting
main_dir, project_dir, subject, task, group = load_inputs()

base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../"))
settings_path = os.path.join(base_dir, project_dir, f'{task}_settings.json')
with open(settings_path) as f:
    settings = json.load(f)
if subject == 'sub-01':
    if task == 'pRF': ses = 'ses-02'
    else: ses = 'ses-01'
else: ses = settings['session']


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

precision_all_runs = []
precision_one_thrs_list = []

threshold = settings['threshold']


for run in range(num_run):
    matfile = scipy.io.loadmat(data_mat[run])

    if task == "PurLoc":
        pred_x_intpl, pred_y_intpl = predicted_pursuit(
            dfs_runs[run], settings)
        
        # Save prediction x and y as tsv.gz
        prediction = np.stack((pred_x_intpl, pred_y_intpl), axis=1)
        prediction = pd.DataFrame(prediction, columns=['pred_x', 'pred_y'])
        pred_file_path = f'{file_dir_save}/timeseries/{subject}_task-{task}_run_0{run+1}_prediction.tsv.gz'
        prediction.to_csv(pred_file_path, sep='\t', index=False, compression='gzip')

        eye_data_run_01 = pd.read_csv(f"{file_dir_save}/timeseries/{subject}_task-{task}_run_01_eyedata.tsv.gz", compression='gzip', delimiter='\t')
        eye_data_run_02 = pd.read_csv(f"{file_dir_save}/timeseries/{subject}_task-{task}_run_02_eyedata.tsv.gz", compression='gzip', delimiter='\t')
        eye_data_all_runs = [eye_data_run_01[['timestamp', 'x', 'y', 'pupil_size']].to_numpy(),
                             eye_data_run_02[['timestamp', 'x', 'y', 'pupil_size']].to_numpy()]

        
        eucl_dist = euclidean_distance_pur(eye_data_all_runs,pred_x_intpl, pred_y_intpl, run)
    
    elif task == "pRF": 
        pred_x_intpl = np.zeros(len(eye_data_all_runs[run]))
        pred_y_intpl = np.zeros(len(eye_data_all_runs[run]))

        eye_data_run_01 = pd.read_csv(f"{file_dir_save}/timeseries/{subject}_task-{task}_run_01_eyedata.tsv.gz", compression='gzip', delimiter='\t')
        eye_data_run_02 = pd.read_csv(f"{file_dir_save}/timeseries/{subject}_task-{task}_run_02_eyedata.tsv.gz", compression='gzip', delimiter='\t')
        eye_data_run_03 = pd.read_csv(f"{file_dir_save}/timeseries/{subject}_task-{task}_run_03_eyedata.tsv.gz", compression='gzip', delimiter='\t')
        eye_data_run_04 = pd.read_csv(f"{file_dir_save}/timeseries/{subject}_task-{task}_run_04_eyedata.tsv.gz", compression='gzip', delimiter='\t')
        eye_data_run_05 = pd.read_csv(f"{file_dir_save}/timeseries/{subject}_task-{task}_run_05_eyedata.tsv.gz", compression='gzip', delimiter='\t')
        
        eye_data_all_runs = [eye_data_run_01[['timestamp', 'x', 'y', 'pupil_size']].to_numpy(), 
                             eye_data_run_02[['timestamp', 'x', 'y', 'pupil_size']].to_numpy(),
                             eye_data_run_03[['timestamp', 'x', 'y', 'pupil_size']].to_numpy(),
                             eye_data_run_04[['timestamp', 'x', 'y', 'pupil_size']].to_numpy(),
                             eye_data_run_05[['timestamp', 'x', 'y', 'pupil_size']].to_numpy()]
        
        eucl_dist =  euclidean_distance_pur(eye_data_all_runs,pred_x_intpl, pred_y_intpl, run)


    elif task == "SacLoc": 
        eye_data_run_01 = pd.read_csv(f"{file_dir_save}/timeseries/{subject}_task-{task}_run_01_eyedata.tsv.gz", compression='gzip', delimiter='\t')
        eye_data_run_02 = pd.read_csv(f"{file_dir_save}/timeseries/{subject}_task-{task}_run_02_eyedata.tsv.gz", compression='gzip', delimiter='\t')
        eye_data_all_runs = [eye_data_run_01[['timestamp', 'x', 'y', 'pupil_size']].to_numpy(), 
                            eye_data_run_02[['timestamp', 'x', 'y', 'pupil_size']].to_numpy()]
        
        pred_x_intpl, pred_y_intpl = load_sac_model(file_dir_save, subject, run, eye_data_all_runs[run])

        eucl_dist = euclidean_distance(eye_data_all_runs,pred_x_intpl, pred_y_intpl, run)


    eucl_dist_df = pd.DataFrame(eucl_dist, columns=['ee'])
    # Save eucl_dist as tsv.gz
    ee_file_path = f'{file_dir_save}/timeseries/{subject}_task-{task}_run_0{run+1}_ee.tsv.gz'
    eucl_dist_df.to_csv(ee_file_path, sep='\t', index=False, compression='gzip')


    precision_fraction = fraction_under_threshold(pred_x_intpl, eucl_dist)
    precision_one_thrs = fraction_under_one_threshold(pred_x_intpl,eucl_dist,threshold)
        
    
    # Store precision for this run
    precision_all_runs.append(precision_fraction)
    precision_one_thrs_list.append(precision_one_thrs)



# Combine all precision data into a single DataFrame
precision_df = pd.DataFrame(precision_all_runs).T  # Transpose so each column is a run
precision_one_df = pd.DataFrame(precision_one_thrs_list).T  # Transpose so each column is a run

# Rename columns to match `run_01`, `run_02`, etc.
precision_df.columns = [f"run_{i+1:02d}" for i in range(num_run)]
precision_one_df.columns = [f"run_{i+1:02d}" for i in range(num_run)]

#precision_df["threshold"] = np.linspace(0, 9.0, 100)
# Add a column for the mean across runs
precision_df["precision_mean"] = precision_df.mean(axis=1)
precision_one_df["precision_one_thrs_mean"] = precision_one_df.mean(axis=1)


# Save the DataFrame to a TSV file
output_tsv_file = f"{file_dir_save}/stats/{subject}_task-{task}_precision_summary.tsv"
precision_df.to_csv(output_tsv_file, sep="\t", index=False)

output_one_tsv_file = f"{file_dir_save}/stats/{subject}_task-{task}_precision_one_threshold_summary.tsv"
precision_one_df.to_csv(output_one_tsv_file, sep="\t", index=False)

print(f"Saved precision summary to {output_tsv_file}")
    
# Define permission cmd
print('Changing files permissions in {}/{}'.format(main_dir, project_dir))
os.system("chmod -Rf 771 {}/{}".format(main_dir, project_dir))
os.system("chgrp -Rf {} {}/{}".format(group, main_dir, project_dir))


