"""
-----------------------------------------------------------------------------------------
vdm_builder_sacloc.py
-----------------------------------------------------------------------------------------
Goal of the script:
Generate vidual design matrix and videos of retinal displacement of SacLoc task.
Loops over all runs and creates concatenated downsampled VDM.
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject name
sys.argv[4]: group of shared data (e.g. 327)
-----------------------------------------------------------------------------------------
Output(s):
Per run:
- target + gaze vdm_target_gaze.mp4 in eyetracking frequency
- retinal vdm_preview.mp4 in eyetracking frequency 
- retinal vdm_run.npy and vdm.mp4 downsampled to TR 

Concatenated:
- vdm.npy downsampled to TR (concatenated across all runs)
-----------------------------------------------------------------------------------------
To run:
cd ~/projects/pRF_analysis/RetinoMaps/postproc/pmf/fit
python vdm_builder_sacloc.py /scratch/mszinte/data RetinoMaps sub-02 327
-----------------------------------------------------------------------------------------
Written by Sina Kling (sina.kling@outlook.de)
Modified to process all runs and concatenate outputs
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

# Debug
import ipdb
deb = ipdb.set_trace


sys.path.append("{}/../../../../analysis_code/utils".format(os.getcwd()))
from eyetrack_utils import (
    load_eye_data, load_event_files,
    flatten_time_start_trial, create_visual_frame_target_gaze, 
    create_visual_frame, downsample_vdm_to_tr
)
from sac_utils import predicted_saccade

from settings_utils import load_settings

##################################### CONFIGURATION ##########################################
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
group = sys.argv[4]

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

base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../../"))
settings_path = os.path.join(base_dir, project_dir, 'eyetracking', 'eye-tracking.yml')
settings = load_settings([settings_path])
analysis_info = settings[0]

print(f"Loading data from subject {subject}")
print(f"Found {len(event_runs)} runs")

# Load all eye data
eye_data_runs = load_eye_data(main_dir, project_dir, subject, TASK)

# LOAD HDF5 FILES TO GET TRIAL ONSET TIMES
h5_filename = f"{eye_tracking_dir}/stats/{subject}_task-{TASK}_eyedata_sac_stats.h5"

print("Loading h5")
with h5py.File(h5_filename, "r") as f:
    time_start = np.array(f["time_start_trial"])  # (max_trials, n_sequences, n_runs)
    time_end   = np.array(f["time_end_trial"])

# Storage for concatenated VDM
all_vdm_downsampled = []

##################################### MAIN LOOP OVER RUNS ####################################
for run in range(1, len(event_runs) + 1):
    print(f"\n{'='*70}")
    print(f"PROCESSING RUN {run}")
    print(f"{'='*70}")
    
    eye_data = eye_data_runs[run-1] 
    gaze_x = eye_data[:,1]
    gaze_y = eye_data[:,2]

    trial_onsets_ms = flatten_time_start_trial(time_start, run_idx=run-1)
    eye_times = (eye_data[:, 0] - eye_data[0, 0]) * TIME_SCALE

    # Trial onsets (seconds since first eye sample)
    trial_onsets_sec = (trial_onsets_ms - eye_data[0, 0]) * TIME_SCALE

    # COMPUTE EXPECTED TARGET POSITIONS
    target_x, target_y = predicted_saccade(main_dir, project_dir, event_runs[run-1], analysis_info)

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
    np.save(os.path.join(gaze_dir, f"{subject}_SacLoc_run_0{run}_expected_target_X"), target_x_expanded)
    # save target y
    np.save(os.path.join(gaze_dir, f"{subject}_SacLoc_run_0{run}_expected_target_Y"), target_y_expanded)
    
    # CALCULATE RETINAL POSITION 
    retino_x = target_x_expanded - gaze_x
    retino_y = target_y_expanded - gaze_y

    # save retino X
    np.save(os.path.join(gaze_dir, f"{subject}_SacLoc_run_0{run}_retino_X"), retino_x)
    # save retino y
    np.save(os.path.join(gaze_dir, f"{subject}_SacLoc_run_0{run}_retino_Y"), retino_y)

    # ------------ 
    print(f"Creating target + gaze design matrix for run {run}...")

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
    print(f"Creating visual design matrix for run {run}...")
    vdm = np.zeros((n_timepoints, n_pixels, n_pixels), dtype=np.float32)

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
    print(f"Downsampling visual design matrix for run {run}...")
    vdm_downsampled = downsample_vdm_to_tr(vdm_for_downsampling, TR=TR, original_sampling_rate=1000)
    vdm_downsampled = np.transpose(vdm_downsampled, (1, 2, 0))  # back to (H, W, T_TR)
    
    np.save(os.path.join(vdm_dir, f"{subject}_run_0{run}_task-{TASK}_vdm.npy"), vdm_downsampled)
    print(f"VDM downsampled shape: {vdm_downsampled.shape}")
    
    # Store for concatenation
    all_vdm_downsampled.append(vdm_downsampled)

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
    print(f"Creating target + gaze video for run {run} (synced)...")
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
    print(f"Creating original VDM preview video for run {run} (synced)...")
    vdm_video_fn = f'{vdm_dir}/{subject}_run_0{run}_task-{TASK}_vdm_preview.mp4'
    out_vdm = cv2.VideoWriter(vdm_video_fn, fourcc, preview_fps, (n_pixels, n_pixels), False)

    for idx in preview_indices:
        frame = vdm[:, :, idx]
        out_vdm.write(np.uint8(frame * 255))

    out_vdm.release()
    print(f'✓ Saved preview video: {os.path.getsize(vdm_video_fn) / 1e6:.1f} MB')
    print(f'  Duration: {len(preview_indices) / preview_fps:.1f} seconds at {preview_fps} fps')

    # 3. DOWNSAMPLED VDM VIDEO (TR resolution)
    print(f"Creating downsampled VDM video for run {run} (TR resolution)...")
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
    print(f"\nVIDEO SUMMARY FOR RUN {run}:")
    print(f"  All videos should be ~{target_duration_seconds:.1f} seconds")
    print(f"  Preview videos: {len(preview_indices)} frames at {preview_fps} fps")
    print(f"  Downsampled video: {n_TRs * frames_per_tr} frames at {preview_fps} fps")

##################################### CONCATENATE ALL RUNS ##################################
print(f"\n{'='*70}")
print("CONCATENATING ALL RUNS")
print(f"{'='*70}")

# Concatenate along time dimension (last axis)
vdm_all_runs = np.concatenate(all_vdm_downsampled, axis=-1)
concatenated_fn = os.path.join(vdm_dir, f"{subject}_task-{TASK}_vdm.npy")
np.save(concatenated_fn, vdm_all_runs)

print(f"✓ Concatenated VDM shape: {vdm_all_runs.shape}")
print(f"  Individual run shapes: {[vdm.shape for vdm in all_vdm_downsampled]}")
print(f"  Total TRs: {vdm_all_runs.shape[-1]}")
print(f"  Expected: ~{n_TRs * len(event_runs)} TRs")
print(f"  Saved to: {concatenated_fn}")

# Define permission cmd
print(f'\nChanging files permissions in {main_dir}/{project_dir}')
os.system(f"chmod -Rf 771 {main_dir}/{project_dir}")
os.system(f"chgrp -Rf {group} {main_dir}/{project_dir}")
print("PROCESSING COMPLETE")