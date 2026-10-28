"""
-----------------------------------------------------------------------------------------
saccade_prediction.py
-----------------------------------------------------------------------------------------
Goal of the script:
Generate prediction for SacLoc task (combining target position with saccade 
on- and offsets)
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject name
sys.argv[4]: group of shared data (e.g. 327)
-----------------------------------------------------------------------------------------
Output(s):
tsv.gz timeseries of Prediction
-----------------------------------------------------------------------------------------
To run:
cd ~/projects/pRF_analysis/RetinoMaps/eyetracking/postproc/
python saccade_prediction.py /scratch/mszinte/data RetinoMaps sub-01 327
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
import h5py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Personal imports
sys.path.append("{}/../../../analysis_code/utils".format(os.getcwd()))
from settings_utils import load_settings
from eyetrack_utils import load_event_files
from sac_utils import predicted_saccade, add_missing_sac_rows

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
group = sys.argv[4]
task = 'SacLoc'

# Defind directories 
eye_tracking_dir = '{}/{}/derivatives/pp_data/{}/eyetracking'.format(main_dir, project_dir, subject)

# Load settings 
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../"))
base_dir_general = os.path.abspath(os.path.join(os.getcwd(), "../../../"))
eye_tracking_settings_path = os.path.join(base_dir_general, project_dir, "eyetracking-analysis.yml")
settings_path = os.path.join(base_dir_general, project_dir, "settings.yml")
task_settings_path = os.path.join(base_dir, "eye-tracking_{}.yml".format(task))
task_settings = load_settings([task_settings_path, eye_tracking_settings_path, settings_path])[0] 

if subject == 'sub-01':
    if task == 'pRF': ses = 'ses-02'
    else: ses = 'ses-01'
else: ses = task_settings['session']

trials_out = task_settings["trials_out"]
trials_in = task_settings["trials_in"]

h5_filename = f'{eye_tracking_dir}/stats/{subject}_task-{task}_eyedata_sac_stats.h5'
h5_file = h5py.File(h5_filename,'r')
time_start_trial = np.array(h5_file['time_start_trial'])
time_end_trial = np.array(h5_file['time_end_trial'])
time_start_seq = np.array(h5_file['time_start_seq'])
time_end_seq = np.array(h5_file['time_end_seq'])
time_start_eye = np.array(h5_file['time_start_eye'])
time_end_eye = np.array(h5_file['time_end_eye'])

# Load event files
data_events = load_event_files(main_dir, project_dir, subject, ses, task)
dfs_runs = []

# Read each TSV file 
for i, run in enumerate(data_events):
    df_run = pd.read_csv(run, sep="\t")
    dfs_runs.append(df_run)

# Load eye data 
eye_data_run_01_nan_blink_interpol = pd.read_csv(f"{eye_tracking_dir}/timeseries/{subject}_task-{task}_run_01_eyedata.tsv.gz", compression='gzip', delimiter='\t')
eye_data_run_01_nan_blink_interpol = eye_data_run_01_nan_blink_interpol[['timestamp', 'x', 'y', 'pupil_size']].to_numpy()
eye_data_run_02_nan_blink_interpol = pd.read_csv(f"{eye_tracking_dir}/timeseries/{subject}_task-{task}_run_02_eyedata.tsv.gz", compression='gzip', delimiter='\t')
eye_data_run_02_nan_blink_interpol = eye_data_run_02_nan_blink_interpol[['timestamp', 'x', 'y', 'pupil_size']].to_numpy()

eye_data_all_runs = [eye_data_run_01_nan_blink_interpol,eye_data_run_02_nan_blink_interpol]

# Saccade output 
saccade_output = np.array(h5_file['saccades_output'])

columns = [
    'run', 'sequence', 'trial', 'saccade_num', 'sac_x_onset', 'sac_x_offset', 'sac_y_onset', 'sac_y_offset',
    'sac_t_onset', 'sac_t_offset', 'sac_p_onset', 'sac_p_offset', 'sac_dur', 'sac_vpeak', 'sac_dist', 'sac_amp',
    'sac_dist_ang', 'sac_amp_ang', 'fix_cor', 'sac_cor', 'saccade_task', 'miss_time', 'sac_out_accuracy', 'sac_in_accuracy', 'no_saccade',
    'microsaccade', 'blink_saccade'
]

# Convert vals_all to a DataFrame
df_sacc = pd.DataFrame(saccade_output, columns=columns)

# Filter for correct saccades outwards 
correct_saccades_out = df_sacc[df_sacc['sac_out_accuracy'] == 1]
# Filter for correct saccades inwards 
correct_saccades_in = df_sacc[df_sacc['sac_in_accuracy'] == 1]

# outwards saccades

correct_saccades_out_run_1 = correct_saccades_out[correct_saccades_out['run'] == 0]
correct_saccades_out_run_1_new = add_missing_sac_rows(correct_saccades_out_run_1, 'out')

correct_saccades_out_run_2 = correct_saccades_out[correct_saccades_out['run'] == 1]
correct_saccades_out_run_2_new = add_missing_sac_rows(correct_saccades_out_run_2, 'out')

# Filter the DataFrame to only include rows where the 'trial' column has even numbers
out_filtered_df_run_1 = correct_saccades_out_run_1_new[correct_saccades_out_run_1_new['trial'] % 2 == 0]  #remove any detected saccades from the wrong trials
out_filtered_df_run_2 = correct_saccades_out_run_2_new[correct_saccades_out_run_2_new['trial'] % 2 == 0]

# combine 2 new dataframes to see if it worked
outwards_saccades = [out_filtered_df_run_1,out_filtered_df_run_2]
correct_saccades_out_new = pd.concat(outwards_saccades)

# inwards saccades
correct_saccades_in_run_1 = correct_saccades_in[correct_saccades_in['run'] == 0]
correct_saccades_in_run_1_new = add_missing_sac_rows(correct_saccades_in_run_1, 'in')

correct_saccades_in_run_2 = correct_saccades_in[correct_saccades_in['run'] == 1]
correct_saccades_in_run_2_new = add_missing_sac_rows(correct_saccades_in_run_2, 'in')

# Filter the DataFrame to only include rows where the 'trial' column has even numbers
in_filtered_df_run_1 = correct_saccades_in_run_1_new[correct_saccades_in_run_1_new['trial'] % 2 != 0]  #remove any detected saccades from the wrong trials
in_filtered_df_run_2 = correct_saccades_in_run_2_new[correct_saccades_in_run_2_new['trial'] % 2 != 0]

# combine 2 new dataframes to see if it worked
inwards_saccades = [in_filtered_df_run_1,in_filtered_df_run_2]
correct_saccades_in_new = pd.concat(inwards_saccades)

# Convert timestamps
initial_timestamp_run_1 = eye_data_run_01_nan_blink_interpol[0, 0]
initial_timestamp_run_2 = eye_data_run_02_nan_blink_interpol[0, 0]

# Initialize arrays to store onset and offset times in seconds for each run
onset_all_trials_run_1 = np.zeros_like(time_start_trial[:, :, 0])
offset_all_trials_run_1 = np.zeros_like(time_end_trial[:, :, 0])
onset_all_trials_run_2 = np.zeros_like(time_start_trial[:, :, 1])
offset_all_trials_run_2 = np.zeros_like(time_end_trial[:, :, 1])

# Loop through each trial and apply the conversion for run 1
for i in range(time_start_trial.shape[0]):
    for j in range(time_start_trial.shape[1]):
        # [CORRECTION #65] temps en ms (et non en pas de 100 ms tronques par int() plus bas : decalage de 0 a 100 ms)
        # AVANT :         onset_all_trials_run_1[i, j] = (time_start_trial[i, j, 0] - initial_timestamp_run_1) / 100
        # AVANT :         offset_all_trials_run_1[i, j] = (time_end_trial[i, j, 0] - initial_timestamp_run_1) / 100
        onset_all_trials_run_1[i, j] = (time_start_trial[i, j, 0] - initial_timestamp_run_1)
        offset_all_trials_run_1[i, j] = (time_end_trial[i, j, 0] - initial_timestamp_run_1)

# Loop through each trial and apply the conversion for run 2
for i in range(time_start_trial.shape[0]):
    for j in range(time_start_trial.shape[1]):
        # [CORRECTION #65] temps en ms (et non en pas de 100 ms tronques par int() plus bas : decalage de 0 a 100 ms)
        # AVANT :         onset_all_trials_run_2[i, j] = (time_start_trial[i, j, 1] - initial_timestamp_run_2) / 100
        # AVANT :         offset_all_trials_run_2[i, j] = (time_end_trial[i, j, 1] - initial_timestamp_run_2) / 100
        onset_all_trials_run_2[i, j] = (time_start_trial[i, j, 1] - initial_timestamp_run_2)
        offset_all_trials_run_2[i, j] = (time_end_trial[i, j, 1] - initial_timestamp_run_2)

# Flatten the arrays and filter out non-positive values for run 1
onset_all_trials_flat_run_1 = np.ravel(onset_all_trials_run_1, order='F')
# [CORRECTION #66] cases vides (essais 17-32 des sequences de fixation, = 0) reperees sur les temps bruts : le filtre
# '> 0' apres soustraction retirait aussi l'essai 1 (onset = debut des donnees = 0) et decalait tous les indices d'un essai
# AVANT : onset_all_trials_flat_run_1 = onset_all_trials_flat_run_1[onset_all_trials_flat_run_1 > 0]
onset_all_trials_flat_run_1 = onset_all_trials_flat_run_1[np.ravel(time_start_trial[:, :, 0], order='F') > 0]

offset_all_trials_flat_run_1 = np.ravel(offset_all_trials_run_1, order='F')
# [CORRECTION #66] idem pour les offsets (un indice = un essai, 0-based, 208 essais)
# AVANT : offset_all_trials_flat_run_1 = offset_all_trials_flat_run_1[offset_all_trials_flat_run_1 > 0]
offset_all_trials_flat_run_1 = offset_all_trials_flat_run_1[np.ravel(time_end_trial[:, :, 0], order='F') > 0]

# Flatten the arrays and filter out non-positive values for run 2
onset_all_trials_flat_run_2 = np.ravel(onset_all_trials_run_2, order='F')
# [CORRECTION #66] cases vides (essais 17-32 des sequences de fixation, = 0) reperees sur les temps bruts : le filtre
# '> 0' apres soustraction retirait aussi l'essai 1 (onset = debut des donnees = 0) et decalait tous les indices d'un essai
# AVANT : onset_all_trials_flat_run_2 = onset_all_trials_flat_run_2[onset_all_trials_flat_run_2 > 0]
onset_all_trials_flat_run_2 = onset_all_trials_flat_run_2[np.ravel(time_start_trial[:, :, 1], order='F') > 0]

offset_all_trials_flat_run_2 = np.ravel(offset_all_trials_run_2, order='F')
# [CORRECTION #66] idem pour les offsets (un indice = un essai, 0-based, 208 essais)
# AVANT : offset_all_trials_flat_run_2 = offset_all_trials_flat_run_2[offset_all_trials_flat_run_2 > 0]
offset_all_trials_flat_run_2 = offset_all_trials_flat_run_2[np.ravel(time_end_trial[:, :, 1], order='F') > 0]

# Generate Expected Position 
sac_expected_x, sac_expected_y = predicted_saccade(main_dir=main_dir, project_dir=project_dir, df_run=dfs_runs[0], settings=task_settings)



#----------------------------- MODEL X -------------------------------------------------------
out_filtered_dfs = [out_filtered_df_run_1,out_filtered_df_run_2]
in_filtered_dfs = [in_filtered_df_run_1,in_filtered_df_run_2]


trial_onsets_all = [onset_all_trials_flat_run_1,onset_all_trials_flat_run_2]
trial_offsets_all = [offset_all_trials_flat_run_1,offset_all_trials_flat_run_2]

for run, (out_filtered_df, in_filtered_df) in enumerate(zip(out_filtered_dfs, in_filtered_dfs)):
    eye_data = eye_data_all_runs[run]
    # [CORRECTION #65] temps des saccades en ms depuis le debut des donnees (et non en pas de 100 ms tronques par int())
    # AVANT :     correct_sac_out_seconds_on = (out_filtered_df['sac_t_onset'] - eye_data[0, 0]) / 100
    # AVANT :     correct_sac_out_seconds_off = (out_filtered_df['sac_t_offset'] - eye_data[0, 0]) / 100
    # AVANT :     correct_sac_in_seconds_on = (in_filtered_df['sac_t_onset'] - eye_data[0, 0]) / 100
    # AVANT :     correct_sac_in_seconds_off = (in_filtered_df['sac_t_offset'] - eye_data[0, 0]) / 100
    correct_sac_out_seconds_on = (out_filtered_df['sac_t_onset'] - eye_data[0, 0])
    correct_sac_out_seconds_off = (out_filtered_df['sac_t_offset'] - eye_data[0, 0])

    correct_sac_in_seconds_on = (in_filtered_df['sac_t_onset'] - eye_data[0, 0])
    correct_sac_in_seconds_off = (in_filtered_df['sac_t_offset'] - eye_data[0, 0])

    trial_onsets = trial_onsets_all[run]
    trial_offsets = trial_offsets_all[run]


    total_length = task_settings["total_len_sac_model"]

    # Initialize the model with NaNs
    model_x = np.full(total_length, np.nan)

    current_index = 0

    for idx, (trial_in, trial_out) in enumerate(zip(trials_in, trials_out)):
        # [CORRECTION #66] indices d'essais 0-based coherents (design locEMexp : 16 essais de fixation puis 32 essais
        # alternant out/in). trials_out (yml) = numeros 1-based des essais 'out' (17, 19, ...) ; l'essai 'in' qui suit est
        # g_in = g_out + 1 (trials_in du yml, 16, 18, ..., ne correspond a aucun des deux et n'est plus utilise).
        # Cible (constConfig.m) : saut centre -> peripherie a onset(g_out) + 100 ms, retour au centre a onset(g_in) + 100 ms.
        # NB : sac_expected_x/y[trial_out] = position de l'image 1 de l'essai 0-based trial_out = g_in, c.-a-d. la cible
        # peripherique de la saccade 'out' : indexation correcte, conservee.
        # AVANT : (ligne(s) ajoutee(s))
        g_out = trial_out - 1
        g_in = g_out + 1
        saccade_fix_ms = 100    # const.saccade_fix_dur
        onset_out = int(list(correct_sac_out_seconds_on)[idx])
        offset_out = int(list(correct_sac_out_seconds_off)[idx])
        onset_in = int(list(correct_sac_in_seconds_on)[idx])
        offset_in = int(list(correct_sac_in_seconds_off)[idx])

        
        # Case 4: Both saccade outwards and inwards are missing
        if onset_out < 0 and onset_in < 0:
            print("case 4")
            # [CORRECTION #66] temps de repli = temps des sauts de la cible (avant : trial_onsets[trial_out - 2] = un essai trop
            # tot si l'essai 1 n'etait pas retire, et fin de l'essai trial_in = 16 (0-based) = essai 'out')
            # AVANT :             onset_out = int(trial_onsets[trial_out - 2 ])
            # AVANT :             offset_in = int(trial_offsets[trial_in])
            onset_out = int(trial_onsets[g_out] + saccade_fix_ms)
            offset_in = int(trial_onsets[g_in] + saccade_fix_ms)
            model_x[current_index:onset_out] = 0
            model_x[onset_out:offset_in] = sac_expected_x[trial_out]
            current_index = offset_in
            continue

        # Case 2: Saccade outwards is missing
        if onset_out < 0 or offset_out < 0:
            print("case 2")
            # [CORRECTION #66] saut de la cible vers la peripherie (essai g_out, 0-based)
            # AVANT :             onset_out = int(list(trial_onsets)[trial_out - 2])
            onset_out = int(trial_onsets[g_out] + saccade_fix_ms)
            onset_in = int(list(correct_sac_in_seconds_on)[idx])
            model_x[current_index:onset_out] = 0
            model_x[onset_out:onset_in] = sac_expected_x[trial_out]
            offset_in = int(list(correct_sac_in_seconds_off)[idx])
            model_x[onset_in:offset_in] = np.interp(np.arange(onset_in, offset_in), [onset_in, offset_in], [sac_expected_x[trial_out], 0])
            current_index = offset_in
            continue

        # Case 3: Saccade inwards is missing
        if onset_in < 0 or offset_in < 0:
            print("case 3")
            onset_out = int(list(correct_sac_out_seconds_on)[idx])
            offset_out = int(list(correct_sac_out_seconds_off)[idx])
            # [CORRECTION #66] retour de la cible au centre (essai g_in, 0-based)
            # AVANT :             offset_in = int(list(trial_offsets)[trial_in])
            offset_in = int(trial_onsets[g_in] + saccade_fix_ms)
            model_x[current_index:onset_out] = 0
            model_x[onset_out:offset_out] = np.interp(np.arange(onset_out, offset_out), [onset_out, offset_out], [0, sac_expected_x[trial_out]])
            model_x[offset_out:offset_in] = sac_expected_x[trial_out]
            current_index = offset_in
            continue

        # Case 1: All saccades are present
        print('case 1')
        model_x[current_index:onset_out] = 0
        model_x[onset_out:offset_out] = np.interp(np.arange(onset_out, offset_out), [onset_out, offset_out], [0, sac_expected_x[trial_out]])
        model_x[offset_out:onset_in] = sac_expected_x[trial_out]
        model_x[onset_in:offset_in] = np.interp(np.arange(onset_in, offset_in), [onset_in, offset_in], [sac_expected_x[trial_out], 0])
        
        # Update current_index to the next saccade offset
        current_index = offset_in

    # If there are remaining data points after the last trial, set them to 0
    if current_index < total_length:
        model_x[current_index:] = 0

    plt.figure(figsize=(20, 10))
    plt.plot(model_x, label='Model', color='blue')
    plt.xlabel('Time (arbitrary units)')
    plt.ylabel('X Position')
    plt.title(f'Run {run + 1 } Model Using Saccade Onsets and Offsets with Interpolation')
    plt.grid(True)
    plt.legend()
    # plt.show()

    # ------------------------- Save ---------------------------------------------------- 
    model_dir = '{}/models'.format(eye_tracking_dir)
    os.makedirs(model_dir, exist_ok=True)
    np.save(f"{model_dir}/{subject}_run-0{run+1}_saccade_model_x", model_x) 
  

#---------------------------- MODEL Y ----------------------------------------------------------
# Convert saccade onset times, relative to the start of the recording

out_filtered_dfs = [out_filtered_df_run_1,out_filtered_df_run_2]
in_filtered_dfs = [in_filtered_df_run_1,in_filtered_df_run_2]


trial_onsets_all = [onset_all_trials_flat_run_1,onset_all_trials_flat_run_2]
trial_offsets_all = [offset_all_trials_flat_run_1,offset_all_trials_flat_run_2]

for run, (out_filtered_df, in_filtered_df) in enumerate(zip(out_filtered_dfs, in_filtered_dfs)):
    eye_data = eye_data_all_runs[run]
    # [CORRECTION #65] temps des saccades en ms depuis le debut des donnees (et non en pas de 100 ms tronques par int())
    # AVANT :     correct_sac_out_seconds_on = (out_filtered_df['sac_t_onset'] - eye_data[0, 0]) / 100
    # AVANT :     correct_sac_out_seconds_off = (out_filtered_df['sac_t_offset'] - eye_data[0, 0]) / 100
    # AVANT :     correct_sac_in_seconds_on = (in_filtered_df['sac_t_onset'] - eye_data[0, 0]) / 100
    # AVANT :     correct_sac_in_seconds_off = (in_filtered_df['sac_t_offset'] - eye_data[0, 0]) / 100
    correct_sac_out_seconds_on = (out_filtered_df['sac_t_onset'] - eye_data[0, 0])
    correct_sac_out_seconds_off = (out_filtered_df['sac_t_offset'] - eye_data[0, 0])

    correct_sac_in_seconds_on = (in_filtered_df['sac_t_onset'] - eye_data[0, 0])
    correct_sac_in_seconds_off = (in_filtered_df['sac_t_offset'] - eye_data[0, 0])

    trial_onsets = trial_onsets_all[run]
    trial_offsets = trial_offsets_all[run]

    
    total_length = task_settings["total_len_sac_model"] 

    # Initialize the model with NaNs
    model_y = np.full(total_length, np.nan)

    current_index = 0

    for idx, (trial_in, trial_out) in enumerate(zip(trials_in, trials_out)):
        # [CORRECTION #66] indices d'essais 0-based coherents (design locEMexp : 16 essais de fixation puis 32 essais
        # alternant out/in). trials_out (yml) = numeros 1-based des essais 'out' (17, 19, ...) ; l'essai 'in' qui suit est
        # g_in = g_out + 1 (trials_in du yml, 16, 18, ..., ne correspond a aucun des deux et n'est plus utilise).
        # Cible (constConfig.m) : saut centre -> peripherie a onset(g_out) + 100 ms, retour au centre a onset(g_in) + 100 ms.
        # NB : sac_expected_x/y[trial_out] = position de l'image 1 de l'essai 0-based trial_out = g_in, c.-a-d. la cible
        # peripherique de la saccade 'out' : indexation correcte, conservee.
        # AVANT : (ligne(s) ajoutee(s))
        g_out = trial_out - 1
        g_in = g_out + 1
        saccade_fix_ms = 100    # const.saccade_fix_dur
        onset_out = int(list(correct_sac_out_seconds_on)[idx])
        offset_out = int(list(correct_sac_out_seconds_off)[idx])
        onset_in = int(list(correct_sac_in_seconds_on)[idx])
        offset_in = int(list(correct_sac_in_seconds_off)[idx])
        
        # Case 4: Both saccade outwards and inwards are missing
        if onset_out < 0 and onset_in < 0:
            print("case 4")
            # [CORRECTION #66] temps de repli = temps des sauts de la cible (avant : trial_onsets[trial_out - 2] = un essai trop
            # tot si l'essai 1 n'etait pas retire, et fin de l'essai trial_in = 16 (0-based) = essai 'out')
            # AVANT :             onset_out = int(trial_onsets[trial_out - 2 ])
            # AVANT :             offset_in = int(trial_offsets[trial_in])
            onset_out = int(trial_onsets[g_out] + saccade_fix_ms)
            offset_in = int(trial_onsets[g_in] + saccade_fix_ms)
            # 0 until trial outwards onset
            model_y[current_index:onset_out] = 0
            # Expected position until trial inwards offset
            model_y[onset_out:offset_in] = sac_expected_y[trial_out]
            current_index = offset_in
            continue

        # Case 2: Saccade outwards is missing
        if onset_out < 0 or offset_out < 0:
            print("case 2")
            # [CORRECTION #66] saut de la cible vers la peripherie (essai g_out, 0-based)
            # AVANT :             onset_out = int(list(trial_onsets)[trial_out - 2])
            onset_out = int(trial_onsets[g_out] + saccade_fix_ms)
            onset_in = int(list(correct_sac_in_seconds_on)[idx])
            offset_in = int(list(correct_sac_in_seconds_off)[idx])
            # 0 until trial outwards onset 
            model_y[current_index:onset_out] = 0
            # Expected position until saccade inwards onset
            model_y[onset_out:onset_in] = sac_expected_y[trial_out]
            # Interpolate between saccade inwards onset and offset
            model_y[onset_in:offset_in] = np.interp(np.arange(onset_in, offset_in), [onset_in, offset_in], [sac_expected_y[trial_out], 0])
            current_index = offset_in
            continue

        # Case 3: Saccade inwards is missing
        if onset_in < 0 or offset_in < 0:
            print("case 3")
            onset_out = int(list(correct_sac_out_seconds_on)[idx])
            offset_out = int(list(correct_sac_out_seconds_off)[idx])
            # [CORRECTION #66] retour de la cible au centre (essai g_in, 0-based)
            # AVANT :             offset_in = int(list(trial_offsets)[trial_in])
            offset_in = int(trial_onsets[g_in] + saccade_fix_ms)
            # 0 until saccade outwards onset
            model_y[current_index:onset_out] = 0
            # Interpolate between saccade outwards onset and offset
            model_y[onset_out:offset_out] = np.interp(np.arange(onset_out, offset_out), [onset_out, offset_out], [0, sac_expected_y[trial_out]])
            # Expected position until trial inwards offset
            model_y[offset_out:offset_in] = sac_expected_y[trial_out]
            current_index = offset_in
            continue

        # Case 1: All saccades are present
        print('case 1')
        # 0 until saccade outwards onset 
        model_y[current_index:onset_out] = 0
        # Interpolate between saccade outwards onset and offset
        model_y[onset_out:offset_out] = np.interp(np.arange(onset_out, offset_out), [onset_out, offset_out], [0, sac_expected_y[trial_out]])
        # Expected Position until saccade inwards onset
        model_y[offset_out:onset_in] = sac_expected_y[trial_out]
        # Interpolate between saccade inwards onset and offset 
        model_y[onset_in:offset_in] = np.interp(np.arange(onset_in, offset_in), [onset_in, offset_in], [sac_expected_y[trial_out], 0])
        
        # Update current_index to the next saccade offset
        current_index = offset_in

    # If there are remaining data points after the last trial, set them to 0
    if current_index < total_length:
        model_y[current_index:] = 0


    plt.figure(figsize=(20, 10))
    plt.plot(model_y, label='Model', color='blue')
    plt.xlabel('Time (arbitrary units)')
    plt.ylabel('Y Position')
    plt.title(f'Run {run +1} Model Using Saccade Onsets and Offsets with Interpolation')
    plt.grid(True)
    plt.legend()
    # plt.show()


    # ------------------------- Save ---------------------------------------------------- 
    np.save(f"{eye_tracking_dir}/models/{subject}_run-0{run+1}_saccade_model_y", model_y)



# # Define permission cmd
# print('Changing files permissions in {}/{}'.format(main_dir, project_dir))
# os.system("chmod -Rf 771 {}/{}".format(main_dir, project_dir))
# os.system("chgrp -Rf {} {}/{}".format(group, main_dir, project_dir))
    










