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



    
def make_vdm_from_saccades_gaussian(
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

    sacc_out = sacc_out.copy(); sacc_out['direction'] = 'out'
    sacc_in  = sacc_in.copy();  sacc_in['direction']  = 'in'
    saccades = pd.concat([sacc_out, sacc_in], ignore_index=True)

    for tr_idx in range(n_TRs):
        tr_start = scan_start + tr_idx * TR
        tr_end   = tr_start + TR

        tr_sacc = saccades[
            (saccades['sac_t_onset']  >= tr_start) &
            (saccades['sac_t_offset'] <  tr_end)
        ]
        if len(tr_sacc) == 0:
            continue

        amplitudes = tr_sacc.apply(
            lambda r: np.sqrt(
                (r['sac_x_offset'] - r['sac_x_onset'])**2 +
                (r['sac_y_offset'] - r['sac_y_onset'])**2
            ), axis=1
        )
        sac = tr_sacc.loc[amplitudes.idxmax()]

        dx = sac['sac_x_offset'] - sac['sac_x_onset']
        dy = sac['sac_y_offset'] - sac['sac_y_onset']
        theta = np.degrees(np.arctan2(dy, dx))

        if sac['direction'] == 'out':
            gx = sac['sac_x_offset']
            gy = sac['sac_y_offset']
        else:
            # Mirror of landing position 
            gx = -sac['sac_x_onset']
            gy = -sac['sac_y_onset']

        amp   = np.sqrt(dx**2 + dy**2)
        sigma = max(amp * sigma_scale, 0.5)

        g = make_rotated_gaussian_2d(
            canvas_size=canvas_size,
            x_dva=gx,
            y_dva=gy,
            sigma_x_dva=sigma,
            sigma_y_dva=0.5 * sigma,
            theta=theta,
            dva_range=dva_range
        )
        vdm[:, :, tr_idx] = g

        log.append({
            'tr':        tr_idx,
            'tr_start':  tr_start,       
            'tr_end':    tr_end,          
            'direction': sac['direction'],
            'gx':        gx,
            'gy':        gy,
            'theta':     theta,
            'amp':       amp,
            'sigma':     sigma,
            't_onset':   sac['sac_t_onset'],
            't_offset':  sac['sac_t_offset'],
        })

    return vdm, pd.DataFrame(log)

def make_vdm_from_saccades_blob(
    sacc_out,
    sacc_in,
    scan_start,
    n_TRs,
    TR=1.2,
    canvas_size=100,
    dva_range=15.0,
    dot_dva=0.6,        # radius of the dot in dva
):
    vdm = np.zeros((canvas_size, canvas_size, n_TRs))
    log = []

    # pre-compute pixel-per-dva scale
    px_per_dva = canvas_size / (2 * dva_range)

    sacc_out = sacc_out.copy(); sacc_out['direction'] = 'out'
    sacc_in  = sacc_in.copy();  sacc_in['direction']  = 'in'
    saccades = pd.concat([sacc_out, sacc_in], ignore_index=True)

    for tr_idx in range(n_TRs):
        tr_start = scan_start + tr_idx * TR
        tr_end   = tr_start + TR

        tr_sacc = saccades[
            (saccades['sac_t_onset'] >= tr_start) &
            (saccades['sac_t_onset'] <  tr_end)
        ]
        if len(tr_sacc) == 0:
            continue

        draw_tr = tr_idx + 1
        if draw_tr >= n_TRs:
            continue

        amplitudes = tr_sacc.apply(
            lambda r: np.sqrt(
                (r['sac_x_offset'] - r['sac_x_onset'])**2 +
                (r['sac_y_offset'] - r['sac_y_onset'])**2
            ), axis=1
        )
        sac = tr_sacc.loc[amplitudes.idxmax()]

        dx = sac['sac_x_offset'] - sac['sac_x_onset']
        dy = sac['sac_y_offset'] - sac['sac_y_onset']
        theta = np.degrees(np.arctan2(dy, dx))
        amp   = np.sqrt(dx**2 + dy**2)

        if sac['direction'] == 'out':
            gx = sac['sac_x_offset']
            gy = sac['sac_y_offset']
        else:
            gx = -sac['sac_x_onset']
            gy = -sac['sac_y_onset']

        # --- position of the circle centre on the canvas (convert from dva to pixels) ---
        cx = int(round((gx + dva_range) / (2 * dva_range) * (canvas_size - 1)))
        cy = int(round(((-gy) + dva_range) / (2 * dva_range) * (canvas_size - 1)))

        # --- radius of the circle (in pixels) ---
        r_px = dot_dva * px_per_dva          # dot_dva [dva] × px_per_dva [px/dva] → px

        # --- build a filled circle mask and stamp it into the volume ---
        ys, xs = np.ogrid[:canvas_size, :canvas_size]   # two 1-D index arrays
        mask   = (xs - cx)**2 + (ys - cy)**2 <= r_px**2 # True where inside the circle
        vdm[mask, tr_idx] = 1.0                          # write 1 into those pixels

        r_px = dot_dva * px_per_dva
        ys, xs = np.ogrid[:canvas_size, :canvas_size]
        mask = (xs - cx)**2 + (ys - cy)**2 <= r_px**2
        vdm[mask, tr_idx] = 1.0

        log.append({
            'tr':        tr_idx,
            'tr_start':  tr_start,
            'tr_end':    tr_end,
            'direction': sac['direction'],
            'gx':        gx,
            'gy':        gy,
            'theta':     theta,
            'amp':       amp,
            't_onset':   sac['sac_t_onset'],
            't_offset':  sac['sac_t_offset'],
        })
    return vdm, pd.DataFrame(log)

def make_vdm_from_pursuit(ts, t_start, t_end, amp, t0, n_TRs, TR, canvas_size, dva_range, dot_dva,
                          frame="sacloc", win_ms=100, win_stat="mean", mad_k=3, mad_floor=0.25, back_max_dva=2):
    """
    Build a pursuit oculomotor design matrix for one run: one disc per pursuit trial, from the measured gaze
    (win_stat of x and y over the last win_ms of the trial; the vector frame also uses the first win_ms).

    frame="sacloc": outward at the outward endpoint p, return at -p (p = endpoint of the preceding outward trial)
                    -> same coding as the SacLoc ODM (outward landing, inward = minus inward-saccade onset)
    frame="stop":   disc at each trial's own endpoint (outward ~p, return ~centre)
    frame="vector": disc at endpoint - startpoint (outward ~p, return ~-p)

    ts: eye timeseries (timestamp [ms], x, y [dva]); t_start/t_end: (n_trials, n_seq) for this run, 0 = no trial;
    amp: amp_sequence (5 = fixation block, left empty); t0: run start on the eye clock (ms).
    Masked trials (no gaze, outside the canvas, endpoint far from its ring, return not back at the centre;
    sacloc: also returns whose outward trial is masked) stay empty.
    Discs are drawn as in make_vdm_from_saccades_blob (centre rounded to a pixel, row 0 = +y).
    """
    win_xy = lambda a, b: getattr(ts.loc[(ts.timestamp >= a) & (ts.timestamp < b), ["x", "y"]], win_stat)().values

    # one row per pursuit trial (fixation sequences skipped)
    rows = []
    for seq in np.where(amp < 5)[0]:
        for trial in range(t_end.shape[0]):
            if t_end[trial, seq] == 0:
                continue
            s, e = t_start[trial, seq], t_end[trial, seq]
            end = win_xy(e - win_ms, e)
            x, y = end - win_xy(s, s + win_ms) if frame == "vector" else end
            rows.append(dict(seq=seq, trial=trial, ring=amp[seq],
                             leg="out" if trial % 2 == 0 else "back",           # even = outward, odd = return
                             x=x, y=y, end_ecc=np.hypot(*end),
                             tr_idx=int(((s + e) / 2 - t0) / 1000 // TR)))      # TR holding the trial midpoint
    log = pd.DataFrame(rows)

    # mask outward legs: missing gaze, outside the canvas, endpoint eccentricity far from its ring
    log["ecc"] = np.hypot(log.x, log.y)
    is_out = log.leg == "out"
    med = log.groupby(["ring", "leg"])["ecc"].transform("median")
    mad = log.groupby(["ring", "leg"])["ecc"].transform(lambda v: np.median(np.abs(v - v.median())))
    radial_ok = np.abs(log.ecc - med) <= mad_k * np.maximum(1.4826 * mad, mad_floor)
    finite_ok = log[["x", "y"]].notna().all(axis=1) & (log.ecc <= dva_range)
    back_ok = log.end_ecc <= back_max_dva                                          # return reached the centre
    log["valid"] = finite_ok & np.where(is_out, radial_ok, back_ok & (radial_ok if frame == "vector" else True))

    # sacloc frame: return = minus the endpoint of its outward trial (same sequence, trial - 1)
    if frame == "sacloc":
        out_xy = log[is_out].set_index(["seq", "trial"])[["x", "y", "valid"]]
        pair = pd.MultiIndex.from_arrays([log.loc[~is_out, "seq"], log.loc[~is_out, "trial"] - 1])
        paired = out_xy.reindex(pair)
        log.loc[~is_out, ["x", "y"]] = -paired[["x", "y"]].values
        log.loc[~is_out, "valid"] = log.loc[~is_out, "valid"].values & paired["valid"].fillna(False).values.astype(bool)
        log["ecc"] = np.hypot(log.x, log.y)

    # draw one disc per valid trial in its TR (same rendering as make_vdm_from_saccades_blob)
    yy, xx = np.mgrid[0:canvas_size, 0:canvas_size]
    r_px = dot_dva * canvas_size / (2 * dva_range)
    vdm = np.zeros((canvas_size, canvas_size, n_TRs), dtype=np.float32)
    for _, r in log[log.valid & log.tr_idx.between(0, n_TRs - 1)].iterrows():
        cx = int(round((r.x + dva_range) / (2 * dva_range) * (canvas_size - 1)))
        cy = int(round((-r.y + dva_range) / (2 * dva_range) * (canvas_size - 1)))     # row 0 = +y
        disc = (xx - cx) ** 2 + (yy - cy) ** 2 <= r_px ** 2
        vdm[..., int(r.tr_idx)] = np.maximum(vdm[..., int(r.tr_idx)], disc)

    return vdm, log


def save_vdm_video(vdm, log_df, output_path, scan_start=0.0,
                   dva_range=15.0, preview_fps=60, TR=1.2):
    n_trs = vdm.shape[2]
    n_total_frames = int(round(n_trs * TR * preview_fps))
    frame_dur = 1.0 / preview_fps

    fourcc = cv2.VideoWriter_fourcc(*'mp4v')
    writer = cv2.VideoWriter(output_path, fourcc, preview_fps,
                             (vdm.shape[1], vdm.shape[0]), False)

    for frame_idx in range(n_total_frames):
        t_frame_start = scan_start + frame_idx * frame_dur

        tr_idx = min(int((t_frame_start - scan_start) / TR), n_trs - 1)

        # Show the dot for the entire TR if any saccade was logged in it
        if np.any(vdm[:, :, tr_idx]):
            frame = vdm[:, :, tr_idx]
            frame_uint8 = np.uint8(np.clip(frame, 0, 1) * 255)
        else:
            frame_uint8 = np.zeros((vdm.shape[0], vdm.shape[1]), dtype=np.uint8)

        writer.write(frame_uint8)

    writer.release()
    size_mb = os.path.getsize(output_path) / 1e6
    print(f"Saved {output_path}  ({n_trs * TR:.1f}s, {n_total_frames} frames, {size_mb:.1f} MB)")

