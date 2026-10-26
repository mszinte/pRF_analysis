"""
-----------------------------------------------------------------------------------------
odm_builder_purloc.py
-----------------------------------------------------------------------------------------
Goal of the script:
Generate the oculomotor design matrix of the pursuit localizer (PurLoc).
One disc per pursuit trial at the measured gaze position (head-centred):
frame "stop"   -> gaze endpoint (outward ~ periphery, return ~ centre)
frame "vector" -> endpoint - startpoint (outward ~ p, return ~ -p)
Trials where the gaze did not arrive / overshot are masked (left empty).
Loops over all runs and creates the concatenated ODM.
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject name
sys.argv[4]: group of shared data (e.g. 327)
-----------------------------------------------------------------------------------------
Output(s):
Per run:      {subject}_run-0X_task-PurLoc_odm_log.tsv
Concatenated: {subject}_task-PurLoc_odm.npy  (frame "stop"; other frames get a -{frame} suffix)
-----------------------------------------------------------------------------------------
To run:
cd ~/projects/pRF_analysis/RetinoMaps/postproc/pmf/fit
python odm_builder_purloc.py /scratch/mszinte/data RetinoMaps sub-02 meso_0327
-----------------------------------------------------------------------------------------
Written by Sina Kling (sina.kling@outlook.de)
-----------------------------------------------------------------------------------------
"""
 
import sys
import pandas as pd
import numpy as np
import os
import h5py
 
sys.path.append("{}/../../../../analysis_code/utils".format(os.getcwd()))
from settings_utils import load_settings
sys.path.append("{}/../utils".format(os.getcwd()))
from pmf_utils import make_vdm_from_pursuit, save_vdm_video
 
main_dir    = sys.argv[1]
project_dir = sys.argv[2]
subject     = sys.argv[3]
group       = sys.argv[4]
 
# ── Load settings
base_dir      = os.path.abspath(os.path.join(os.getcwd(), "../../../../"))
settings_path = os.path.join(base_dir, project_dir, 'settings.yml')
pmf_settings_path = os.path.join(base_dir, project_dir, 'pmf-analysis.yml')
settings = load_settings([settings_path, pmf_settings_path])[0]
 
n_TRs    = settings['TRs']
TR       = settings['TR']
n_runs   = 2
dm_range = settings['dm_range']
target_radius = settings['target_radius']
 
odm_frame = "sacloc"     # "stop" or "vector" or "sacloc"
canvas_size = 100
 
save_fn = (f"{main_dir}/{project_dir}/derivatives/vdm/{subject}")
os.makedirs(save_fn, exist_ok=True)
odm_suffix = "" if odm_frame == "stop" else f"-{odm_frame}"
 
# ── Load trial timing  ────────────────────────────────────
h5_filename = (
    f"{main_dir}/{project_dir}/derivatives/pp_data/{subject}/eyetracking/stats/"
    f"{subject}_task-PurLoc_eyedata_sac_stats.h5"
)
with h5py.File(h5_filename, 'r') as h5_file:
    time_start_trial = np.array(h5_file['time_start_trial'])    # (trial, sequence, run), eye clock (ms), 0 = no trial
    time_end_trial   = np.array(h5_file['time_end_trial'])
    time_start_eye   = np.array(h5_file['time_start_eye'])      # (1, run)
    amp_sequence     = np.array(h5_file['amp_sequence'])        # ring index per sequence, 5 = fixation
 
# ── Per-run loop ───────────────────────────────────────────────────────────────
vdm_runs = []
 
for run_idx in range(n_runs):
    run_num = run_idx + 1
 
    # Load per-run eyetracking timeseries
    eye_path = (
        f"{main_dir}/{project_dir}/derivatives/pp_data/{subject}/eyetracking/timeseries/"
        f"{subject}_task-PurLoc_run_{run_num:02d}_eyedata.tsv.gz"
    )
    eye_data = pd.read_csv(eye_path, compression='gzip', delimiter='\t')
    initial_timestamp = eye_data['timestamp'].iloc[0]
    print(f"run {run_num}: first eye sample - time_start_eye = {initial_timestamp - time_start_eye[0, run_idx]:.0f} ms")
 
    vdm, log = make_vdm_from_pursuit(
        ts=eye_data,
        t_start=time_start_trial[..., run_idx],
        t_end=time_end_trial[..., run_idx],
        amp=amp_sequence,
        t0=initial_timestamp,
        n_TRs=n_TRs,
        TR=TR,
        canvas_size=canvas_size,
        dva_range=dm_range,
        dot_dva=target_radius,
        frame=odm_frame
    )
    vdm_runs.append(vdm)
 
    # Masking summary
    print(f"run {run_num}: kept {log['valid'].mean():.0%} of pursuit trials")
    print(log.groupby(['ring', 'leg'])['valid'].mean().unstack().round(2))
 
    # Sanity check: disc centroid y in the ODM has the sign of the coded gaze y (row 0 = +y)
    chk = log[log['valid'] & (log['y'].abs() > 1) & log['tr_idx'].between(0, n_TRs - 1)]
    rows = np.arange(canvas_size)[:, None]
    frames = vdm[..., chk['tr_idx'].astype(int).values]
    y_centroid = -((frames.sum(1) * rows).sum(0) / frames.sum((0, 1)) - (canvas_size - 1) / 2)
    print(f"run {run_num}: ODM-y vs gaze-y sign agreement = {np.mean(np.sign(y_centroid) == np.sign(chk['y'].values)):.2f}")
    
    # Save video 
    video_fn = f"{save_fn}/{subject}_run-{run_num:02d}_task-PurLoc_odm.mp4"
    save_vdm_video(vdm, log, output_path=video_fn, dva_range=dm_range, TR=TR)
 
# ── Concatenate and save all runs ──────────────────────────────────────────────
vdm_concat = np.concatenate(vdm_runs, axis=-1)
print(vdm_concat.shape)
 
    

np.save(f"{save_fn}/{subject}_task-PurLoc_odm.npy", vdm_concat)
print(f"Saved concatenated ODM: {save_fn}  shape={vdm_concat.shape}, range= {dm_range} dva, frame = {odm_frame}")