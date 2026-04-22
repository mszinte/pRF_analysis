"""
-----------------------------------------------------------------------------------------
mdm_builder_sac_dir.py
-----------------------------------------------------------------------------------------
Goal of the script:
Generate motion design matrix and videos of saccade direction in retinal space.
Uses saccade analysis from RetinoMaps. 
Loops over all runs and creates concatenated downsampled MDM.
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject name
sys.argv[4]: group of shared data (e.g. 327)
-----------------------------------------------------------------------------------------
Output(s):
Per run:      {subject}_task-SacLoc_saccade_vdm_run_{i}.mp4
Concatenated: {subject}_task-SacLoc_saccade_vdm_concat.npy
-----------------------------------------------------------------------------------------
To run:
cd ~/projects/pRF_analysis/RetinoMaps/postproc/pmf/fit
python mdm_builder_sac_dir.py /scratch/mszinte/data RetinoMaps sub-02 327
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
from sac_utils import add_missing_sac_rows
from pmf_utils import make_vdm_from_saccades, save_vdm_video

main_dir    = sys.argv[1]
project_dir = sys.argv[2]
subject     = sys.argv[3]
group       = sys.argv[4]

# ── Load settings + events 
base_dir      = os.path.abspath(os.path.join(os.getcwd(), "../../../../"))
settings_path = os.path.join(base_dir, project_dir, 'settings.yml')
analysis_info = load_settings([settings_path])[0]

n_TRs   = analysis_info['TRs']
TR      = analysis_info['TR']
n_runs  = 2

save_fn = (f"{main_dir}/{project_dir}/derivatives/vdm/{subject}")

# ── Load saccade stats  ────────────────────────────────────
h5_filename = (
    f"{main_dir}/{project_dir}/derivatives/pp_data/{subject}/eyetracking/stats/"
    f"{subject}_task-SacLoc_eyedata_sac_stats.h5"
)
with h5py.File(h5_filename, 'r') as h5_file:
    saccade_output   = np.array(h5_file['saccades_output'])
    time_start_trial = np.array(h5_file['time_start_trial'])
    time_end_trial   = np.array(h5_file['time_end_trial'])

columns = [
    'run', 'sequence', 'trial', 'saccade_num',
    'sac_x_onset', 'sac_x_offset',
    'sac_y_onset', 'sac_y_offset',
    'sac_t_onset', 'sac_t_offset',
    'sac_p_onset', 'sac_p_offset',
    'sac_dur', 'sac_vpeak', 'sac_dist', 'sac_amp',
    'sac_dist_ang', 'sac_amp_ang',
    'fix_cor', 'sac_cor', 'saccade_task',
    'miss_time', 'sac_out_accuracy', 'sac_in_accuracy',
    'no_saccade', 'microsaccade', 'blink_saccade'
]
df_sacc = pd.DataFrame(saccade_output, columns=columns)

# ── Per-run loop ───────────────────────────────────────────────────────────────
vdm_runs = []

for run_idx in range(n_runs):
    run_num = run_idx + 1   

    # Load per-run eyetracking timeseries
    eye_path = (
        f"{main_dir}/{project_dir}/derivatives/pp_data/{subject}/eyetracking/timeseries/"
        f"{subject}_task-SacLoc_run_{run_num:02d}_eyedata.tsv.gz"
    )
    eye_data = pd.read_csv(eye_path, compression='gzip', delimiter='\t')
    eye_data = eye_data[['timestamp', 'x', 'y', 'pupil_size']].to_numpy()
    initial_timestamp = eye_data[0, 0]

    # Filter saccades for this run (h5 uses 0-based run index)
    sacc_run = df_sacc[df_sacc['run'] == run_idx].copy()

    # Outward saccades: correct, even trials only
    sacc_out = sacc_run[sacc_run['sac_out_accuracy'] == 1].copy()
    sacc_out = add_missing_sac_rows(sacc_out, 'out')
    sacc_out = sacc_out[sacc_out['trial'] % 2 == 0]

    # Inward saccades: correct, odd trials only
    sacc_in = sacc_run[sacc_run['sac_in_accuracy'] == 1].copy()
    sacc_in = add_missing_sac_rows(sacc_in, 'in')
    sacc_in = sacc_in[sacc_in['trial'] % 2 != 0]

    # Convert timestamps to seconds relative to run start
    for df in (sacc_out, sacc_in):
        df['sac_t_onset']  = (df['sac_t_onset']  - initial_timestamp) / 1000
        df['sac_t_offset'] = (df['sac_t_offset'] - initial_timestamp) / 1000

    # Build VDM for this run
    vdm, log = make_vdm_from_saccades(
        sacc_out=sacc_out,
        sacc_in=sacc_in,
        scan_start=0.0,
        n_TRs=n_TRs,
        TR=TR,
        canvas_size=100,
        dva_range=10.0,
        sigma_scale=0.3
    )
    vdm_runs.append(vdm)

    # Save per-run video
    video_fn = f"{save_fn}/{subject}_task-SacLoc_saccade_vdm_run_{run_num:02d}.mp4"

    save_vdm_video(vdm, log, output_path=video_fn)
    print(f"Saved video: {video_fn}")

# ── Concatenate and save all runs ──────────────────────────────────────────────
vdm_concat = np.concatenate(vdm_runs, axis=0)   

np.save(f"{save_fn}/{subject}_task-SacLoc_saccade_mdm.npy", vdm_concat)
print(f"Saved concatenated VDM: {save_fn}  shape={vdm_concat.shape}")