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
- target-only vdm_target_run.npy downsampled to TR (black and white)

Concatenated:
- vdm.npy downsampled to TR (concatenated across all runs)
- vdm_target.npy downsampled to TR (concatenated across all runs)
-----------------------------------------------------------------------------------------
To run:
cd ~/projects/pRF_analysis/RetinoMaps/postproc/pmf/fit
python vdm_builder_sacloc.py /scratch/mszinte/data RetinoMaps sub-02 327
-----------------------------------------------------------------------------------------
Written by Sina Kling (sina.kling@outlook.de)
-----------------------------------------------------------------------------------------
"""

import pandas as pd
import numpy as np
import os, sys, cv2, ipdb

sys.path.append("{}/../../../../analysis_code/utils".format(os.getcwd()))
from eyetrack_utils import (
    load_eye_data, load_event_files,
    flatten_time_start_trial,
    create_visual_frame_target_gaze,
    create_visual_frame,
)
from sac_utils import predicted_saccade
from settings_utils import load_settings
import h5py

deb = ipdb.set_trace

# ── Config ────────────────────────────────────────────────────────────────────
main_dir    = sys.argv[1]
project_dir = sys.argv[2]
subject     = sys.argv[3]
group       = sys.argv[4]

TASK               = "SacLoc"
TR                 = 1.2
TIME_SCALE         = 1e-3          # ms → seconds
SAMPLES_PER_TR     = int(TR * 1000)  # 1200 samples at 1 kHz
N_PIXELS           = 100
VISUAL_FIELD_SIZE  = 40            # degrees
TARGET_RADIUS      = 0.5           # degrees
CHUNK_SIZE         = 10_000        # timepoints per chunk (controls RAM use)
PREVIEW_FPS        = 60

ses = 'ses-01' if subject == 'sub-01' else 'ses-02'
sub_num = subject[4:]

# ── Directories 
gaze_dir = f"{main_dir}/{project_dir}/derivatives/pp_data/{subject}/eyetracking/timeseries"
vdm_dir  = f"{main_dir}/{project_dir}/derivatives/vdm/{subject}"
os.makedirs(gaze_dir, exist_ok=True)
os.makedirs(vdm_dir,  exist_ok=True)

# ── Load settings + events 
base_dir      = os.path.abspath(os.path.join(os.getcwd(), "../../../../"))
settings_path = os.path.join(base_dir, project_dir, 'eyetracking', 'eye-tracking.yml')
analysis_info = load_settings([settings_path])[0]

event_files = load_event_files(main_dir, project_dir, subject, ses, TASK)
event_runs  = [pd.read_csv(fp, sep='\t') for fp in event_files]

eye_data_runs = load_eye_data(main_dir, project_dir, subject, TASK)

h5_fn = (f"{main_dir}/{project_dir}/derivatives/pp_data/{subject}/eyetracking"
         f"/stats/{subject}_task-{TASK}_eyedata_sac_stats.h5")
with h5py.File(h5_fn, "r") as f:
    time_start = np.array(f["time_start_trial"])

print(f"Subject: {subject}  |  Runs found: {len(event_runs)}")


def build_vdm_chunked(frame_fn, n_timepoints):
    """
    Uses the last sample of each TR window instead of mean,
    avoiding carryover blending from position transitions.
    """
    n_trs = int(np.floor(n_timepoints / SAMPLES_PER_TR))
    out   = np.zeros((n_trs, N_PIXELS, N_PIXELS), dtype=np.float32)

    for tr_idx in range(n_trs):
        # Take the middle sample of the TR window
        t = tr_idx * SAMPLES_PER_TR + SAMPLES_PER_TR // 2
        out[tr_idx] = frame_fn(t)

        if tr_idx % 50 == 0:
            print(f"    TR {tr_idx}/{n_trs}", end='\r')

    print(f"    {n_trs}/{n_trs} TRs processed")
    return np.transpose(out, (1, 2, 0))  # → (H, W, T_TR)




# ── Per-run processing 
all_vdm_downsampled        = []
all_vdm_target_downsampled = []

for run in range(1, len(event_runs) + 1):
    print(f"\n{'='*60}\nRUN {run}\n{'='*60}")

    eye_data = eye_data_runs[run - 1]
    gaze_x   = eye_data[:, 1]
    gaze_y   = eye_data[:, 2]
    eye_times = (eye_data[:, 0] - eye_data[0, 0]) * TIME_SCALE
    n_timepoints = len(gaze_x)

    # Trial onsets
    trial_onsets_ms  = flatten_time_start_trial(time_start, run_idx=run - 1)
    trial_onsets_sec = (trial_onsets_ms - eye_data[0, 0]) * TIME_SCALE

    # Expected target positions
    target_x, target_y = predicted_saccade(main_dir, project_dir,
                                           event_runs[run - 1], analysis_info)
    target_x_exp = np.full(n_timepoints, np.nan)
    target_y_exp = np.full(n_timepoints, np.nan)

    for i in range(len(trial_onsets_sec) - 1):
        s = np.argmin(np.abs(eye_times - trial_onsets_sec[i]))
        e = np.argmin(np.abs(eye_times - trial_onsets_sec[i + 1]))
        target_x_exp[s:e] = target_x[i + 1]
        target_y_exp[s:e] = target_y[i + 1]
    s = np.argmin(np.abs(eye_times - trial_onsets_sec[-1]))
    target_x_exp[s:] = target_x[0]
    target_y_exp[s:] = target_y[0]

    # Retinal position
    retino_x = target_x_exp - gaze_x
    retino_y = target_y_exp - gaze_y

    # Save gaze timeseries
    for tag, arr in [("expected_target_X", target_x_exp),
                     ("expected_target_Y", target_y_exp),
                     ("retino_X", retino_x),
                     ("retino_Y", retino_y)]:
        np.save(os.path.join(gaze_dir, f"{subject}_SacLoc_run_0{run}_{tag}"), arr)

    # ── 1. Retinal VDM (chunked, saved as npy) 
    print("Building retinal VDM...")
    vdm_down = build_vdm_chunked(
        lambda t: create_visual_frame(
            retino_x[t], retino_y[t],
            n_pixels=N_PIXELS,
            visual_field_size=VISUAL_FIELD_SIZE,
            target_radius=TARGET_RADIUS,
        ),
        n_timepoints,
    )
    np.save(os.path.join(vdm_dir, f"{subject}_run_0{run}_task-{TASK}_vdm.npy"), vdm_down)
    print(f"  Retinal VDM shape: {vdm_down.shape}")
    all_vdm_downsampled.append(vdm_down)

    # ── 2. Target-only VDM (chunked, accumulated for concat) ─────────────────
    print("Building target-only VDM...")
    vdm_target_down = build_vdm_chunked(
        lambda t: create_visual_frame(
            target_x_exp[t], target_y_exp[t],
            n_pixels=N_PIXELS,
            visual_field_size=VISUAL_FIELD_SIZE,
            target_radius=TARGET_RADIUS,
        ),
        n_timepoints,
    )
    print(f"  Target VDM shape: {vdm_target_down.shape}")
    all_vdm_target_downsampled.append(vdm_target_down)   # no per-run save

    # ── 3. Target + gaze video
    print("Writing target+gaze video...")
    n_trs_run = vdm_down.shape[-1]
    frames_per_tr    = int(PREVIEW_FPS * TR)
    total_frames     = n_trs_run * frames_per_tr
    # one eye sample per preview frame (at 60 fps ≈ every 16.7 ms → every 17th sample)
    step = max(1, n_timepoints // total_frames)
    preview_indices  = np.arange(0, n_timepoints, step)[:total_frames]

    video_fn = os.path.join(vdm_dir, f"{subject}_run_0{run}_task-{TASK}_vdm_target_gaze.mp4")
    fourcc   = cv2.VideoWriter_fourcc(*'mp4v')
    writer   = cv2.VideoWriter(video_fn, fourcc, PREVIEW_FPS, (N_PIXELS, N_PIXELS), True)

    for idx in preview_indices:
        frame = create_visual_frame_target_gaze(
            target_x_exp[idx], target_y_exp[idx],
            gaze_x[idx], gaze_y[idx],
            n_pixels=N_PIXELS,
            visual_field_size=VISUAL_FIELD_SIZE,
            target_radius=TARGET_RADIUS,
            gaze_radius=0.5,
        )
        rgb = np.zeros((N_PIXELS, N_PIXELS, 3), dtype=np.uint8)
        rgb[..., 0] = np.uint8(np.clip(frame[0], 0, 1) * 255)  # red  = target
        rgb[..., 2] = np.uint8(np.clip(frame[1], 0, 1) * 255)  # blue = gaze
        writer.write(rgb)

    writer.release()
    size_mb = os.path.getsize(video_fn) / 1e6
    print(f"  Saved {video_fn}  ({size_mb:.1f} MB, {len(preview_indices)/PREVIEW_FPS:.1f}s)")


# ── Concatenate retinal VDMs across all runs
print(f"\n{'='*60}\nConcatenating retinal VDMs\n{'='*60}")
vdm_all = np.concatenate(all_vdm_downsampled, axis=-1)  # (H, W, T_total)
concat_fn = os.path.join(vdm_dir, f"{subject}_task-{TASK}_vdm.npy")
np.save(concat_fn, vdm_all)

# Concatenate target-only VDMs across all runs
vdm_target_all = np.concatenate(all_vdm_target_downsampled, axis=-1)  # (H, W, T_total)
concat_target_fn = os.path.join(vdm_dir, f"{subject}_task-{TASK}_vdm_target.npy")
np.save(concat_target_fn, vdm_target_all)

print(f"  Target per-run shapes : {[v.shape for v in all_vdm_target_downsampled]}")
print(f"  Target concatenated   : {vdm_target_all.shape}  →  {concat_target_fn}")
assert vdm_target_all.shape == vdm_all.shape, \
    f"Shape mismatch between retinal {vdm_all.shape} and target {vdm_target_all.shape} VDMs!"

print(f"  Per-run shapes : {[v.shape for v in all_vdm_downsampled]}")
print(f"  Concatenated   : {vdm_all.shape}  →  {concat_fn}")
assert vdm_all.shape[-1] == sum(v.shape[-1] for v in all_vdm_downsampled), \
    "Concatenation length mismatch!"


print(f"\nSetting permissions in {main_dir}/{project_dir}")
os.system(f"chmod -Rf 771 {main_dir}/{project_dir}")
os.system(f"chgrp -Rf {group} {main_dir}/{project_dir}")
print("Done.")