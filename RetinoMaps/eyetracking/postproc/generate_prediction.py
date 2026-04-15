"""
-----------------------------------------------------------------------------------------
generate_prediction.py
-----------------------------------------------------------------------------------------
Goal of the script:
Generate prediction for eyemovements and calculate euclidean distance 
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject name
sys.argv[4]: group of shared data (e.g. 327)
-----------------------------------------------------------------------------------------
Output(s):
tsv of fraction under thresholds
tsv.gz timeseries of Euclidean distance 
tsv.gz timeseries of Prediction
-----------------------------------------------------------------------------------------
To run:
cd ~/projects/pRF_analysis/RetinoMaps/eyetracking/postproc
python generate_prediction.py /scratch/mszinte/data RetinoMaps sub-01 327
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
import glob 
import numpy as np
import pandas as pd

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
group = sys.argv[4]

# Personal imports
sys.path.append("{}/../../../analysis_code/utils".format(os.getcwd()))
from settings_utils import load_settings
from sac_utils import predicted_pursuit, euclidean_distance_pur, fraction_under_threshold, fraction_under_one_threshold, load_sac_model, euclidean_distance

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
    
    eye = task_settings['eye']
    num_run = task_settings['num_run']
    eyetracking_sampling = task_settings['eyetrack_sampling']
    ppd = task_settings['ppd']
    
    # Defind directories
    eye_tracking_dir = '{}/{}/derivatives/pp_data/{}/eyetracking'.format(main_dir, project_dir, subject)
    os.makedirs(eye_tracking_dir, exist_ok=True)
    fig_dir = '{}/figures'.format(eye_tracking_dir)
    os.makedirs(fig_dir, exist_ok=True)
          
    events_list = sorted(glob.glob('{}/{}/{}/{}/func/{}_{}_task-{}_*_events*.tsv'.format(
        main_dir, project_dir, subject, ses, subject, ses, task)))
   
    dfs_runs = [pd.read_csv(run, sep="\t") for run in events_list]
    
    precision_all_runs = []
    precision_one_thrs_list = []
    
    threshold = task_settings['threshold']
       
    for run in range(num_run):    
        if task == "PurLoc":
            pred_x_intpl, pred_y_intpl = predicted_pursuit(
                dfs_runs[run], task_settings)
            
            # Save prediction x and y as tsv.gz
            prediction = np.stack((pred_x_intpl, pred_y_intpl), axis=1)
            prediction = pd.DataFrame(prediction, columns=['pred_x', 'pred_y'])
            pred_file_path = f'{eye_tracking_dir}/timeseries/{subject}_task-{task}_run_0{run+1}_prediction.tsv.gz'
            prediction.to_csv(pred_file_path, sep='\t', index=False, compression='gzip')
    
            eye_data_run_01 = pd.read_csv(f"{eye_tracking_dir}/timeseries/{subject}_task-{task}_run_01_eyedata.tsv.gz", compression='gzip', delimiter='\t')
            eye_data_run_02 = pd.read_csv(f"{eye_tracking_dir}/timeseries/{subject}_task-{task}_run_02_eyedata.tsv.gz", compression='gzip', delimiter='\t')
            eye_data_all_runs = [eye_data_run_01[['timestamp', 'x', 'y', 'pupil_size']].to_numpy(),
                                 eye_data_run_02[['timestamp', 'x', 'y', 'pupil_size']].to_numpy()]
    
            
            eucl_dist = euclidean_distance_pur(eye_data_all_runs,pred_x_intpl, pred_y_intpl, run)
        
        elif task == "pRF": 
            eye_data_run_01 = pd.read_csv(f"{eye_tracking_dir}/timeseries/{subject}_task-{task}_run_01_eyedata.tsv.gz", compression='gzip', delimiter='\t')
            eye_data_run_02 = pd.read_csv(f"{eye_tracking_dir}/timeseries/{subject}_task-{task}_run_02_eyedata.tsv.gz", compression='gzip', delimiter='\t')
            eye_data_run_03 = pd.read_csv(f"{eye_tracking_dir}/timeseries/{subject}_task-{task}_run_03_eyedata.tsv.gz", compression='gzip', delimiter='\t')
            eye_data_run_04 = pd.read_csv(f"{eye_tracking_dir}/timeseries/{subject}_task-{task}_run_04_eyedata.tsv.gz", compression='gzip', delimiter='\t')
            eye_data_run_05 = pd.read_csv(f"{eye_tracking_dir}/timeseries/{subject}_task-{task}_run_05_eyedata.tsv.gz", compression='gzip', delimiter='\t')
            
            eye_data_all_runs = [eye_data_run_01[['timestamp', 'x', 'y', 'pupil_size']].to_numpy(), 
                                 eye_data_run_02[['timestamp', 'x', 'y', 'pupil_size']].to_numpy(),
                                 eye_data_run_03[['timestamp', 'x', 'y', 'pupil_size']].to_numpy(),
                                 eye_data_run_04[['timestamp', 'x', 'y', 'pupil_size']].to_numpy(),
                                 eye_data_run_05[['timestamp', 'x', 'y', 'pupil_size']].to_numpy()]
            

            # Save prediction x and y as tsv.gz
            pred_x_intpl = np.zeros(len(eye_data_all_runs[run]))
            pred_y_intpl = np.zeros(len(eye_data_all_runs[run]))
            prediction = np.stack((pred_x_intpl, pred_y_intpl), axis=1)
            prediction = pd.DataFrame(prediction, columns=['pred_x', 'pred_y'])
            pred_file_path = f'{eye_tracking_dir}/timeseries/{subject}_task-{task}_run_0{run+1}_prediction.tsv.gz'
            prediction.to_csv(pred_file_path, sep='\t', index=False, compression='gzip')
            
            eucl_dist =  euclidean_distance(eye_data_all_runs,pred_x_intpl, pred_y_intpl, run)
        
        elif task == "rest": 
            eye_data_run_01 = pd.read_csv(f"{eye_tracking_dir}/timeseries/{subject}_task-{task}_run_01_eyedata.tsv.gz", compression='gzip', delimiter='\t')
            eye_data_run_02 = pd.read_csv(f"{eye_tracking_dir}/timeseries/{subject}_task-{task}_run_02_eyedata.tsv.gz", compression='gzip', delimiter='\t')

            
            eye_data_all_runs = [eye_data_run_01[['timestamp', 'x', 'y', 'pupil_size']].to_numpy(), 
                                 eye_data_run_02[['timestamp', 'x', 'y', 'pupil_size']].to_numpy()]
            

            # Save prediction x and y as tsv.gz
            pred_x_intpl = np.zeros(len(eye_data_all_runs[run]))
            pred_y_intpl = np.zeros(len(eye_data_all_runs[run]))
            prediction = np.stack((pred_x_intpl, pred_y_intpl), axis=1)
            prediction = pd.DataFrame(prediction, columns=['pred_x', 'pred_y'])
            pred_file_path = f'{eye_tracking_dir}/timeseries/{subject}_task-{task}_run_0{run+1}_prediction.tsv.gz'
            prediction.to_csv(pred_file_path, sep='\t', index=False, compression='gzip')
            
            eucl_dist =  euclidean_distance(eye_data_all_runs,pred_x_intpl, pred_y_intpl, run)        
    
    
        elif task == "SacLoc": 
            eye_data_run_01 = pd.read_csv(f"{eye_tracking_dir}/timeseries/{subject}_task-{task}_run_01_eyedata.tsv.gz", compression='gzip', delimiter='\t')
            eye_data_run_02 = pd.read_csv(f"{eye_tracking_dir}/timeseries/{subject}_task-{task}_run_02_eyedata.tsv.gz", compression='gzip', delimiter='\t')
            eye_data_all_runs = [eye_data_run_01[['timestamp', 'x', 'y', 'pupil_size']].to_numpy(), 
                                eye_data_run_02[['timestamp', 'x', 'y', 'pupil_size']].to_numpy()]
            
            pred_x_intpl, pred_y_intpl = load_sac_model(eye_tracking_dir, subject, run, eye_data_all_runs[run])
            # Save prediction x and y as tsv.gz
            prediction = np.stack((pred_x_intpl, pred_y_intpl), axis=1)
            prediction = pd.DataFrame(prediction, columns=['pred_x', 'pred_y'])
            pred_file_path = f'{eye_tracking_dir}/timeseries/{subject}_task-{task}_run_0{run+1}_prediction.tsv.gz'
            prediction.to_csv(pred_file_path, sep='\t', index=False, compression='gzip')
    
            eucl_dist = euclidean_distance(eye_data_all_runs,pred_x_intpl, pred_y_intpl, run)
    
    
        eucl_dist_df = pd.DataFrame(eucl_dist, columns=['ee'])
        # Save eucl_dist as tsv.gz
        ee_file_path = f'{eye_tracking_dir}/timeseries/{subject}_task-{task}_run_0{run+1}_ee.tsv.gz'
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
    output_tsv_file = f"{eye_tracking_dir}/stats/{subject}_task-{task}_precision_summary.tsv"
    precision_df.to_csv(output_tsv_file, sep="\t", index=False)
    
    output_one_tsv_file = f"{eye_tracking_dir}/stats/{subject}_task-{task}_precision_one_threshold_summary.tsv"
    precision_one_df.to_csv(output_one_tsv_file, sep="\t", index=False)
    
    print(f"Saved precision summary to {output_tsv_file}")
    
# # Define permission cmd
# print('Changing files permissions in {}/{}'.format(main_dir, project_dir))
# os.system("chmod -Rf 771 {}/{}".format(main_dir, project_dir))
# os.system("chgrp -Rf {} {}/{}".format(group, main_dir, project_dir))


