"""
-----------------------------------------------------------------------------------------
sacloc_vdm-builder.py
-----------------------------------------------------------------------------------------
Goal of the script:
Generate vidual design matrix and videos of retinal displacement of SacLoc task.
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject name
sys.argv[4]: run 
sys.argv[5]: group of shared data (e.g. 327)
-----------------------------------------------------------------------------------------
Output(s):
- target + gaze vdm_target_gaze.mp4 in eyetracking frequency
- retinal vdm_preview.mp4 in eyetracking frequency 
- retinal vdm.npy and vdm.mp4 downsampled to TR 
-----------------------------------------------------------------------------------------
To run:
cd ~/projects/pRF_analysis/RetinoMaps/eyetracking/
python sacloc_vdm-builder.py /scratch/mszinte/data RetinoMaps sub-02 01 327
-----------------------------------------------------------------------------------------
Written by Sina Kling (sina.kling@outlook.de)
-----------------------------------------------------------------------------------------
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import h5py
import os
import cv2
import json
import sys 
import yaml 

sys.path.append("{}/../../analysis_code/utils".format(os.getcwd()))
from eyetrack_utils import (
    load_eye_data, load_event_files, load_settings, 
    flatten_time_start_trial, create_visual_frame_target_gaze, 
    create_visual_frame, downsample_vdm_to_tr
)
from sac_utils import predicted_saccade

##################################### CONFIGURATION ##########################################
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
run = int(sys.argv[4][-1]) 
group = sys.argv[5]

# general settings for vdm 
n_pixels = 100
visual_field_size = 40  # degrees 
pixel_size = visual_field_size / n_pixels  # degrees per pixel
# Define target/stimulus size
target_radius = 0.5  
n_TRs = 208

eye_tracking_dir = f"{main_dir}/{project_dir}/derivatives/pp_data/{subject}/eyetracking" 
gaze_dir = f'{main_dir}/{project_dir}/derivatives/pp_data/{subject}/eyetracking/timeseries'
os.makedirs(gaze_dir, exist_ok=True)
vdm_dir = f'{main_dir}/{project_dir}/derivatives/vdm/{subject}'
os.makedirs(vdm_dir, exist_ok=True)

if subject == 'sub-01':
    ses = 'ses-01'
else: 
    ses = 'ses-02'

TASK = "SacLoc"
TIME_SCALE = 1e-3  # Convert ms to seconds
TR = 1.2

# LOAD DATA
event_files = load_event_files(main_dir, project_dir, subject, ses, TASK)
event_runs = []
for filepath in event_files:
    df = pd.read_csv(filepath, sep='\t')
    event_runs.append(df)

# TODO: adapt to new settings

base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../"))


settings_path = os.path.join(base_dir, project_dir, 'eyetracking', 'eye-tracking.yml')
# Loading the YAML settings file
with open(settings_path) as f:
    settings_raw = yaml.safe_load(f)
settings = {key: value['value'] for key, value in settings_raw.items()}

print(f"Loading data from subject {subject} ")
print(f"processing run: {run}")
eye_data_runs = load_eye_data(main_dir, project_dir, subject, TASK)
eye_data = eye_data_runs[run-1] 
gaze_x = eye_data[:,1]
gaze_y = eye_data[:,2]

# LOAD HDF5 FILES TO GET TRIAL ONSET TIMES
h5_filename = f"{eye_tracking_dir}/stats/{subject}_task-{TASK}_eyedata_sac_stats.h5"

print("Loading h5")
with h5py.File(h5_filename, "r") as f:
    time_start = np.array(f["time_start_trial"])  # (max_trials, n_sequences, n_runs)
    time_end   = np.array(f["time_end_trial"])

trial_onsets_ms = flatten_time_start_trial(time_start, run_idx=run-1)
eye_times = (eye_data[:, 0] - eye_data[0, 0]) * TIME_SCALE

# Trial onsets (seconds since first eye sample)
trial_onsets_sec = (trial_onsets_ms - eye_data[0, 0]) * TIME_SCALE

##################################### MAIN ####################################################
# COMPUTE EXPECTED TARGET POSITIONS
target_x, target_y = predicted_saccade(event_runs[run-1], settings)

target_x_expanded = np.full(len(eye_times), np.nan)
target_y_expanded = np.full(len(eye_times), np.nan)

# For each trial, fill from its onset to the next trial's onset
for i in range(len(trial_onsets_sec) - 1):
    # Find start index
    start_idx = np.argmin(np.abs(eye_times - trial_onsets_sec[i]))
    
    # Find end index (next trial)
    end_idx = np.argmin(np.abs(eye_times - trial_onsets_sec[i+1]))
    
    # Fill with NEXT trial's target value
    target_x_expanded[start_idx:end_idx] = target_x[i+1]
    target_y_expanded[start_idx:end_idx] = target_y[i+1]

# Handle last trial
start_idx = np.argmin(np.abs(eye_times - trial_onsets_sec[-1]))
target_x_expanded[start_idx:] = target_x[0]  
target_y_expanded[start_idx:] = target_y[0]


# save target X
np.save(os.path.join(gaze_dir, f"{subject}_SacLoc_expected_target_X"), target_x_expanded)
# save target y
np.save(os.path.join(gaze_dir, f"{subject}_SacLoc_expected_target_Y"), target_y_expanded)
# ------------
# CALCULATE RETINAL POSITION 
retino_x = target_x_expanded - gaze_x
retino_y = target_y_expanded - gaze_y

# save retino X
np.save(os.path.join(gaze_dir, f"{subject}_SacLoc_retino_X"), retino_x)
# save retino y
np.save(os.path.join(gaze_dir, f"{subject}_SacLoc_retino_Y"), retino_y)

# ------------ 
print("Creating target + gaze design matrix...")
# Define visual field parameters


n_timepoints = len(gaze_x)

# Create vdm_target with shape (2, n_pixels, n_pixels, n_timepoints)
vdm_target = np.zeros((2, n_pixels, n_pixels, n_timepoints), dtype=np.float32)

for t in range(n_timepoints):
    result = create_visual_frame_target_gaze(
        target_x_expanded[t], target_y_expanded[t],
        gaze_x[t], gaze_y[t],
        n_pixels=n_pixels,
        visual_field_size=visual_field_size,
        target_radius=target_radius,
        gaze_radius=0.5 
    )
    if result.shape[0] == 2:
        vdm_target[:, :, :, t] = result
    else:
        raise ValueError(f"Unexpected shape from create_visual_frame_target_gaze: {result.shape}")

print(f"VDM target shape: {vdm_target.shape}")  # (2, H, W, T)

# ------------
# Create visual design matrix for all timepoints
print("\nCreating visual design matrix...")
vdm = np.zeros((n_timepoints,n_pixels, n_pixels), dtype=np.float32)

for t in range(n_timepoints):
    vdm[t] = create_visual_frame(
        retino_x[t], 
        retino_y[t],
        n_pixels=n_pixels,
        visual_field_size=visual_field_size,
        target_radius=target_radius
    )

vdm = np.transpose(vdm, (1, 2, 0))
print(f"Visual design matrix shape: {vdm.shape}") # (H, W, T)


vdm_for_downsampling = np.transpose(vdm, (2, 0, 1))  # (T, H, W)


# Downsample original to TR 
print("Downsampling visual design matrix...")
vdm_downsampled = downsample_vdm_to_tr(vdm_for_downsampling, TR=TR, original_sampling_rate=1000)
vdm_downsampled = np.transpose(vdm_downsampled, (1, 2, 0))  # back to (H, W, T_TR)
np.save(os.path.join(vdm_dir, f"{subject}_run_0{run}_task-{TASK}_vdm.npy"), vdm_downsampled)
print(f"VDM downsampled shape: {vdm_downsampled.shape}")

############################# MAKE VIDEOS #################################################
target_duration_seconds = n_TRs * TR  # Total scan duration
preview_fps = 60
total_frames_needed = int(target_duration_seconds * preview_fps)
downsample_factor_preview = max(1, n_timepoints // total_frames_needed)
preview_indices = np.arange(0, n_timepoints, downsample_factor_preview)[:total_frames_needed]
frames_per_tr = int(preview_fps * TR)

print(f"Target duration: {target_duration_seconds:.1f} seconds")
print(f"Preview FPS: {preview_fps}")
print(f"Total frames needed: {total_frames_needed}")

# 1. VIDEO FOR TARGET + GAZE (original resolution, synced to downsampled)
print("Creating target + gaze video (synced)...")
fourcc = cv2.VideoWriter_fourcc(*'mp4v')
vdm_target_video_fn = f'{vdm_dir}/{subject}_run_0{run}_task-{TASK}_vdm_target_gaze.mp4'
out_target = cv2.VideoWriter(vdm_target_video_fn, fourcc, preview_fps, (n_pixels, n_pixels), True)

for idx in preview_indices:
    target_frame = vdm_target[0, :, :, idx]
    gaze_frame = vdm_target[1, :, :, idx]
    rgb = np.zeros((n_pixels, n_pixels, 3), dtype=np.uint8)
    rgb[..., 0] = np.uint8(np.clip(target_frame, 0, 1) * 255)  # Red channel
    rgb[..., 2] = np.uint8(np.clip(gaze_frame, 0, 1) * 255)    # Blue channel
    out_target.write(rgb)

out_target.release()
print(f'✓ Saved target+gaze video: {os.path.getsize(vdm_target_video_fn) / 1e6:.1f} MB')
print(f'  Duration: {len(preview_indices) / preview_fps:.1f} seconds at {preview_fps} fps')


# 2. ORIGINAL VDM VIDEO (original resolution, synced to downsampled)
print("\nCreating original VDM preview video (synced)...")
vdm_video_fn = f'{vdm_dir}/{subject}_run_0{run}_task-{TASK}_vdm_preview.mp4'
out_vdm = cv2.VideoWriter(vdm_video_fn, fourcc, preview_fps, (n_pixels, n_pixels), False)

for idx in preview_indices:
    frame = vdm[:, :, idx]
    out_vdm.write(np.uint8(frame * 255))

out_vdm.release()
print(f'✓ Saved preview video: {os.path.getsize(vdm_video_fn) / 1e6:.1f} MB')
print(f'  Duration: {len(preview_indices) / preview_fps:.1f} seconds at {preview_fps} fps')


# 3. DOWNSAMPLED VDM VIDEO (TR resolution)
print("\nCreating downsampled VDM video (TR resolution)...")
vdm_downsampled_video_fn = f'{vdm_dir}/{subject}_run_0{run}_task-{TASK}_vdm_downsampled.mp4'
out_downsampled = cv2.VideoWriter(vdm_downsampled_video_fn, fourcc, preview_fps, 
                                   (n_pixels, n_pixels), False)

# Normalize to 0-255 range
vdm_down_normalized = (vdm_downsampled - vdm_downsampled.min()) / (vdm_downsampled.max() - vdm_downsampled.min() + 1e-8)

for i in range(n_TRs):
    frame = vdm_down_normalized[:, :, i]
    for _ in range(frames_per_tr):
        out_downsampled.write(np.uint8(frame * 255))

out_downsampled.release()
print(f'✓ Saved downsampled video: {os.path.getsize(vdm_downsampled_video_fn) / 1e6:.1f} MB')
print(f'  Duration: {n_TRs * frames_per_tr / preview_fps:.1f} seconds at {preview_fps} fps')

# Final summary
print("\nVIDEO SUMMARY:")
print(f"  All videos should be ~{target_duration_seconds:.1f} seconds")
print(f"  Preview videos: {len(preview_indices)} frames at {preview_fps} fps")
print(f"  Downsampled video: {n_TRs * frames_per_tr} frames at {preview_fps} fps")

# Define permission cmd
print('Changing files permissions in {}/{}'.format(main_dir, project_dir))
os.system("chmod -Rf 771 {}/{}".format(main_dir, project_dir))
os.system("chgrp -Rf {} {}/{}".format(group, main_dir, project_dir))
print("PROCESSING COMPLETE")

