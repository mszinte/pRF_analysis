"""
-----------------------------------------------------------------------------------------
vdm_builder_sacloc.py
-----------------------------------------------------------------------------------------
Goal of the script:
Generate visual design matrix and videos of retinal displacement of SacLoc task.
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
- target + gaze preview:  vdm_target_gaze.mp4    (eyetracking frequency)
- retinal preview:        vdm_preview.mp4         (eyetracking frequency)
- retinal downsampled:    vdm_run_0{run}.mp4      (TR resolution)
Concatenated:
- vdm.npy   downsampled to TR (concatenated across all runs)
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
VISUAL_FIELD_SIZE  = 20            # degrees
TARGET_RADIUS      = 0.5           # degrees
PREVIEW_FPS        = 60

ses = 'ses-01' if subject == 'sub-01' else 'ses-02'
sub_num = subject[4:]

# ── Directories ───────────────────────────────────────────────────────────────
gaze_dir = f"{main_dir}/{project_dir}/derivatives/pp_data/{subject}/eyetracking/timeseries"
vdm_dir  = f"{main_dir}/{project_dir}/derivatives/vdm/{subject}"
os.makedirs(gaze_dir, exist_ok=True)
os.makedirs(vdm_dir,  exist_ok=True)

# ── Load settings + events ────────────────────────────────────────────────────
base_dir      = os.path.abspath(os.path.join(os.getcwd(), "../../../../"))
settings_path = os.path.join(base_dir, project_dir, 'eyetracking-analysis.yml')
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
    Build VDM by averaging all samples within each TR window
    """
    n_trs = int(np.floor(n_timepoints / SAMPLES_PER_TR))
    out   = np.zeros((n_trs, N_PIXELS, N_PIXELS), dtype=np.float32)

    for tr_idx in range(n_trs):
        start = tr_idx * SAMPLES_PER_TR
        end   = min((tr_idx + 1) * SAMPLES_PER_TR, n_timepoints)

        # accumulate 
        acc = np.zeros((N_PIXELS, N_PIXELS), dtype=np.float32)
        count = 0

        for t in range(start, end):
            acc += frame_fn(t)
            count += 1

        out[tr_idx] = acc / max(count, 1)

        if tr_idx % 10 == 0:
            print(f"    TR {tr_idx}/{n_trs}", end='\r')

    print(f"    {n_trs}/{n_trs} TRs processed")

    return np.transpose(out, (1, 2, 0))  # → (H, W, T_TR)


def write_grayscale_video(frames_hw_t, output_path, fps, repeat=1, flip_ud=True):
    H, W, T = frames_hw_t.shape
    fourcc  = cv2.VideoWriter_fourcc(*'mp4v')
    writer  = cv2.VideoWriter(output_path, fourcc, fps, (W, H), False)

    for t in range(T):
        frame = np.uint8(np.clip(frames_hw_t[:, :, t], 0, 1) * 255)
        if flip_ud:
            frame = np.flipud(frame)

        for _ in range(repeat):
            writer.write(frame)

    writer.release()
    print(f"  Saved {output_path}")


# ── Per-run processing ────────────────────────────────────────────────────────
all_vdm_downsampled = []

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

    # ── 1. Retinal VDM downsampled → npy + mp4 per run ───────────────────────
    print("Building retinal VDM (downsampled)...")
    vdm_down = build_vdm_chunked(
        lambda t: create_visual_frame(
            retino_x[t], retino_y[t],
            n_pixels=N_PIXELS,
            visual_field_size=VISUAL_FIELD_SIZE,
            target_radius=TARGET_RADIUS,
        ),
        n_timepoints,
    )
    print(f"  Retinal VDM shape: {vdm_down.shape}")

    # Save npy
    npy_fn = os.path.join(vdm_dir, f"{subject}_run_0{run}_task-{TASK}_vdm.npy")
    np.save(npy_fn, vdm_down)
    print(f"  Saved {npy_fn}")

    # Save mp4 at TR resolution (one frame per TR, played at PREVIEW_FPS)
    mp4_fn = os.path.join(vdm_dir, f"{subject}_run_0{run}_task-{TASK}_vdm.mp4")
    frames_per_tr = 10
    video_fps = frames_per_tr / TR
    write_grayscale_video(vdm_down, mp4_fn, fps=video_fps, repeat=frames_per_tr, flip_ud=False)

    all_vdm_downsampled.append(vdm_down)

    # ── 2. Retinal preview video at eyetracking frequency ─────────────────────
    print("Writing retinal preview video (eyetracking frequency)...")
    n_trs_run    = vdm_down.shape[-1]
    frames_per_tr = int(PREVIEW_FPS * TR)
    total_frames  = n_trs_run * frames_per_tr
    step          = max(1, n_timepoints // total_frames)
    preview_idx   = np.arange(0, n_timepoints, step)[:total_frames]

    preview_fn = os.path.join(vdm_dir, f"{subject}_run_0{run}_task-{TASK}_vdm_preview.mp4")
    fourcc     = cv2.VideoWriter_fourcc(*'mp4v')
    writer     = cv2.VideoWriter(preview_fn, fourcc, PREVIEW_FPS, (N_PIXELS, N_PIXELS), False)
    for idx in preview_idx:
        frame = create_visual_frame(
            retino_x[idx], retino_y[idx],
            n_pixels=N_PIXELS,
            visual_field_size=VISUAL_FIELD_SIZE,
            target_radius=TARGET_RADIUS,
        )
        frame_u8 = np.uint8(np.clip(frame, 0, 1) * 255)
        writer.write(frame_u8)
    writer.release()
    size_mb = os.path.getsize(preview_fn) / 1e6
    print(f"  Saved {preview_fn}  ({len(preview_idx)/PREVIEW_FPS:.1f}s, {size_mb:.1f} MB)")

    # ── 3. Target + gaze video at eyetracking frequency ───────────────────────
    print("Writing target+gaze video (eyetracking frequency)...")
    tg_fn  = os.path.join(vdm_dir, f"{subject}_run_0{run}_task-{TASK}_vdm_target_gaze.mp4")
    fourcc = cv2.VideoWriter_fourcc(*'mp4v')
    writer = cv2.VideoWriter(tg_fn, fourcc, PREVIEW_FPS, (N_PIXELS, N_PIXELS), True)
    for idx in preview_idx:
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
    size_mb = os.path.getsize(tg_fn) / 1e6
    print(f"  Saved {tg_fn}  ({len(preview_idx)/PREVIEW_FPS:.1f}s, {size_mb:.1f} MB)")

# ── Concatenate retinal VDMs across all runs ──────────────────────────────────
print(f"\n{'='*60}\nConcatenating retinal VDMs\n{'='*60}")
vdm_all   = np.concatenate(all_vdm_downsampled, axis=-1)   # (H, W, T_total)
concat_fn = os.path.join(vdm_dir, f"{subject}_task-{TASK}_vdm.npy")
np.save(concat_fn, vdm_all)

assert vdm_all.shape[-1] == sum(v.shape[-1] for v in all_vdm_downsampled), \
    "Concatenation length mismatch!"
print(f"  Per-run shapes : {[v.shape for v in all_vdm_downsampled]}")
print(f"  Concatenated   : {vdm_all.shape}  →  {concat_fn}")

print(f"\nSetting permissions in {main_dir}/{project_dir}")
os.system(f"chmod -Rf 771 {main_dir}/{project_dir}")
os.system(f"chgrp -Rf {group} {main_dir}/{project_dir}")
print("Done.")