import numpy as np
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import cv2
import os

def make_rotated_gaussian_2d(
    canvas_size=100,
    x_dva=0.0,
    y_dva=0.0,
    sigma_x_dva=1.0,
    sigma_y_dva=0.5,
    theta=0.0,
    dva_range=10.0,
):
    x_axis = np.linspace(-dva_range, dva_range, canvas_size)
    y_axis = np.linspace(-dva_range, dva_range, canvas_size)
    xx, yy = np.meshgrid(x_axis, y_axis, indexing="xy")

    dx = xx - x_dva
    dy = yy - y_dva

    t = np.deg2rad(theta)
    cos_t, sin_t = np.cos(t), np.sin(t)

    dx_rot =  cos_t * dx + sin_t * dy
    dy_rot = -sin_t * dx + cos_t * dy

    canvas = np.exp(
        -(
            dx_rot**2 / (2 * sigma_x_dva**2)
            + dy_rot**2 / (2 * sigma_y_dva**2)
        )
    )
    return canvas




def make_vdm_from_saccades(
    sacc_out,
    sacc_in,
    scan_start,
    n_TRs,
    TR=1.2,
    canvas_size=100,
    dva_range=10.0,
    sigma_scale=0.5
):
    """
    Build VDM from precomputed saccades (from H5 file).

    Parameters
    ----------
    sacc_out : DataFrame
        Correct outward saccades
    sacc_in : DataFrame
        Correct inward saccades
    scan_start : float
        Scan start time 
    n_TRs : int
    TR : float (seconds)
    """
    '''
    Algorithm: 

    For each TR: 
        select saccades whose onset time is within 
        for each saccade: 
            compute direction (dx, dy) for gaussian angle 
            select trial type (outwards, inwards)
            if outwards: 
                define gaussian center at saccade landing position 
            elif inwards: 
                define gaussian center at opposite of saccade landing position 
            compute sigma from amplitude (exponential with scaling factor of 0.3)
            draw gaussian 
    '''

    vdm = np.zeros((canvas_size, canvas_size, n_TRs))
    log = []

    # --- combine and label ---
    sacc_out = sacc_out.copy()
    sacc_out['direction'] = 'out'

    sacc_in = sacc_in.copy()
    sacc_in['direction'] = 'in'

    
    saccades = pd.concat([sacc_out, sacc_in], ignore_index=True)

    # --- main loop ---
    for tr_idx in range(n_TRs):

        tr_start = scan_start + tr_idx * TR
        tr_end   = tr_start + TR

        # select saccades in this TR
        tr_sacc = saccades[
            (saccades['sac_t_onset'] >= tr_start) &
            (saccades['sac_t_onset'] < tr_end)
        ]

        if len(tr_sacc) == 0:
            continue

        # --- select strongest saccade by amplitude ---
        amplitudes = tr_sacc.apply(
            lambda r: np.sqrt(
                (r['sac_x_offset'] - r['sac_x_onset'])**2 +
                (r['sac_y_offset'] - r['sac_y_onset'])**2
            ), axis=1
        )
        sac = tr_sacc.loc[amplitudes.idxmax()]  # single row, a Series

        # --- direction vector ---
        dx = sac['sac_x_offset'] - sac['sac_x_onset']
        dy = sac['sac_y_offset'] - sac['sac_y_onset']
        theta = np.degrees(np.arctan2(dy, dx))

        # --- placement ---
        if sac['direction'] == 'out':
            gx = sac['sac_x_offset']
            gy = sac['sac_y_offset']
        else:
            gx = -sac['sac_x_onset']
            gy = -sac['sac_y_onset']

        # --- amplitude & sigma ---
        amp = np.sqrt(dx**2 + dy**2)
        sigma = max(amp * sigma_scale, 0.5)

        
        # --- draw gaussian once ---
        g = make_rotated_gaussian_2d(
            canvas_size=canvas_size,
            x_dva=gx, y_dva=gy,
            sigma_x_dva=sigma,
            sigma_y_dva=0.5 * sigma,
            theta=theta,
            dva_range=dva_range
        )

        # --- timing  ---
        t_on  = sac['sac_t_onset']
        t_off = sac['sac_t_offset']         
        sac_dur = max(t_off - t_on, 1e-6)   # guard against zero-duration
        
        for overlap_tr in range(n_TRs):
            tr_s = scan_start + overlap_tr * TR
            tr_e = tr_s + TR

            overlap = min(t_off, tr_e) - max(t_on, tr_s)
            if overlap <= 0:
                continue

            weight = overlap / sac_dur      # fraction of saccade in this TR
            vdm[:, :, overlap_tr] = np.maximum(
                vdm[:, :, overlap_tr],
                g * weight
            )

        # --- logging (one row per saccade, not per TR) ---
        log.append({
            'tr': tr_idx,           # TR where onset fell
            'direction': sac['direction'],
            'gx': gx, 'gy': gy,
            'theta': theta,
            'amp': amp, 'sigma': sigma,
            't_onset': t_on,
            't_offset': t_off
        })

    log_df = pd.DataFrame(log)

    return vdm, log_df



def save_vdm_video(vdm, log_df, output_path, scan_start=0.0,
                   dva_range=10.0, preview_fps=60, TR=1.2):
    """
    Save saccade VDM as grayscale mp4 with correct within-TR timing.
    Each video frame corresponds to a real time window; gaussians are shown
    only during frames that overlap [t_onset, t_offset].
    """
    n_trs = vdm.shape[2]
    n_total_frames = int(round(n_trs * TR * preview_fps))
    frame_dur = 1.0 / preview_fps          # seconds per frame

    fourcc = cv2.VideoWriter_fourcc(*'mp4v')
    writer = cv2.VideoWriter(output_path, fourcc, preview_fps,
                             (vdm.shape[1], vdm.shape[0]), False)

    for frame_idx in range(n_total_frames):
        # real time at centre of this frame
        t_frame_start = scan_start + frame_idx * frame_dur
        t_frame_end   = t_frame_start + frame_dur

        # which TR does this frame fall in?
        tr_idx = min(int(t_frame_start / TR - scan_start / TR), n_trs - 1)

        # any saccade whose [t_onset, t_offset] overlaps this frame?
        active = log_df[
            (log_df['t_onset']  < t_frame_end) &
            (log_df['t_offset'] > t_frame_start)
        ]

        if len(active) > 0:
            frame = vdm[:, :, tr_idx]
            frame_uint8 = np.uint8(np.clip(frame, 0, 1) * 255)
        else:
            frame_uint8 = np.zeros((vdm.shape[0], vdm.shape[1]), dtype=np.uint8)

        writer.write(np.flipud(frame_uint8))

    writer.release()
    size_mb = os.path.getsize(output_path) / 1e6
    print(f"Saved {output_path}  ({n_trs * TR:.1f}s, {n_total_frames} frames, {size_mb:.1f} MB)")
