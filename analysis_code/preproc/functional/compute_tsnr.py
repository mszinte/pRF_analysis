"""
-----------------------------------------------------------------------------------------
compute_tsnr.py
-----------------------------------------------------------------------------------------
Goal of the script:
Compute per-run tSNR maps (standard and/or robust) from preprocessed surface timeseries
(post-DCT, pre-z-score), and save an average tSNR map across all runs.
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject name
sys.argv[4]: group of shared data (e.g. 327)
-----------------------------------------------------------------------------------------
Output(s):
Per-run tSNR maps and averaged tSNR map across all runs
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/pRF_analysis/analysis_code/preproc/functional/
2. run python command
python compute_tsnr.py [main directory] [project name] [subject] [group]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/projects/pRF_analysis/analysis_code/preproc/functional/
python compute_tsnr.py /scratch/mszinte/data amblyo7T_prf sub-01 327
-----------------------------------------------------------------------------------------
Written by Martin Szinte (martin.szinte@gmail.com)
-----------------------------------------------------------------------------------------
"""

# Stop warnings
import warnings
warnings.filterwarnings("ignore")

# Debug
import ipdb
deb = ipdb.set_trace

# General imports
import os
import sys
import glob
import datetime
import numpy as np
import nibabel as nb

# Personal imports
sys.path.append("{}/../../utils".format(os.getcwd()))
from surface_utils import load_surface, make_surface_image, compute_tsnr, compute_tsnr_robust
from pycortex_utils import set_pycortex_config_file
from settings_utils import load_settings

# Time
start_time = datetime.datetime.now()

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
group = sys.argv[4]

# Set pycortex db and colormaps
cortex_dir = "{}/{}/derivatives/pp_data/cortex".format(main_dir, project_dir)
set_pycortex_config_file(cortex_dir)

# Load settings
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../"))
settings_path = os.path.join(base_dir, project_dir, "settings.yml")
settings = load_settings([settings_path])
analysis_info = settings[0]

formats = analysis_info['formats']
extensions = analysis_info['extensions']
sessions = analysis_info['sessions']
preproc_prep = analysis_info['preproc_prep']
filtering = analysis_info['filtering']
normalization = analysis_info['normalization']
partial_scan = analysis_info['partial_scan']
starting_TR = analysis_info['starting_TR']
ending_TR = analysis_info['ending_TR']
TR = analysis_info['TR']
tsnr_method = analysis_info['tsnr_method']

# Make output folders
for format_ in formats:
    os.makedirs("{}/{}/derivatives/pp_data/{}/{}/func/{}_{}".format(
        main_dir, project_dir, subject, format_, preproc_prep, filtering), exist_ok=True)

# Dictionary to collect tSNR maps across runs for averaging
tsnr_collect = {'standard': {}, 'robust': {}}

for format_, extension in zip(formats, extensions):
    print('Computing tSNR: {}'.format(format_))
    for session in sessions:
        # Find fMRIprep outputs
        fmriprep_func_fns = glob.glob(
            "{}/{}/derivatives/fmriprep/fmriprep/{}/{}/func/*{}*.{}".format(
                main_dir, project_dir, subject, session, format_, extension))

        if not fmriprep_func_fns:
            print('No files for {}'.format(session))
            continue

        for func_fn in fmriprep_func_fns:
            fn_basename = func_fn.split('/')[-1]
            print('Processing: {}'.format(fn_basename))

            # Load data
            surf_img, surf_data = load_surface(fn=func_fn)

            # Cut TRs
            end_idx = surf_data.shape[0] if ending_TR is None else surf_data.shape[0] - ending_TR
            surf_data = surf_data[starting_TR:end_idx, :]

            # Valid vertices mask
            if partial_scan:
                vertex_means = np.mean(surf_data, axis=0)
                valid_mask = vertex_means > 10000
            else:
                valid_mask = np.ones(surf_data.shape[1], dtype=bool)

            # DCT high-pass filtering
            from nilearn import signal
            from nilearn.glm.first_level import make_first_level_design_matrix
            nb_tr = surf_data.shape[0]
            design_matrix = make_first_level_design_matrix(
                frame_times=np.arange(nb_tr) * TR,
                drift_model='cosine',
                high_pass=analysis_info['high_pass_threshold'])
            cosine_drift = design_matrix.values[:, :-1]
            surf_data = signal.clean(surf_data,
                                     detrend=False,
                                     standardize=False,
                                     confounds=cosine_drift)

            # Build format+hemi key for collecting across runs
            if 'hemi-L' in func_fn:
                fkey = '{}_hemi-L'.format(format_)
            elif 'hemi-R' in func_fn:
                fkey = '{}_hemi-R'.format(format_)
            else:
                fkey = format_

            # Compute and save tSNR maps
            for method in (['standard', 'robust'] if tsnr_method == 'both' else [tsnr_method]):
                fn_label = 'tSNR' if method == 'standard' else 'tSNR-robust'
                func = compute_tsnr if method == 'standard' else compute_tsnr_robust

                tsnr_map = func(surf_data[:, valid_mask])
                tsnr_full = np.full(surf_data.shape[1], np.nan)
                tsnr_full[valid_mask] = tsnr_map

                # Per-run filename: replace normalization label with tSNR label
                out_fn_base = fn_basename.replace(
                    '_bold.{}'.format(extension), '_{}_{}_{}.{}'.format(preproc_prep, filtering, fn_label, extension))
                out_fn = "{}/{}/derivatives/pp_data/{}/{}/func/{}_{}/{}".format(
                    main_dir, project_dir, subject, format_,
                    preproc_prep, filtering, out_fn_base)

                tsnr_img = make_surface_image(data=tsnr_full[np.newaxis, :], source_img=surf_img)
                nb.save(tsnr_img, out_fn)

                # Collect for averaging
                tsnr_collect[method].setdefault(fkey, {'maps': [], 'template': surf_img}).get('maps').append(tsnr_full)
                tsnr_collect[method][fkey]['template'] = surf_img

# Save averaged tSNR maps across all runs
print('Saving averaged tSNR maps...')
for method in (['standard', 'robust'] if tsnr_method == 'both' else [tsnr_method]):
    fn_label = 'tSNR' if method == 'standard' else 'tSNR-robust'
    for fkey, data in tsnr_collect[method].items():
        avg_map = np.nanmedian(np.stack(data['maps'], axis=0), axis=0)
        format_ = fkey.split('_')[0]
        hemi = '_'.join(fkey.split('_')[1:]) if '_' in fkey else None
        hemi_str = '_{}'.format(hemi) if hemi else ''
        space_label = 'space-{}'.format(format_)
        avg_fn = "{}/{}/derivatives/pp_data/{}/{}/func/{}_{}/{}{}_{}_{}_{}_avg.{}".format(
            main_dir, project_dir, subject, format_,
            preproc_prep, filtering,
            subject, hemi_str, space_label, preproc_prep, fn_label,
            'func.gii' if 'fsnative' in format_ else 'dtseries.nii')
        avg_img = make_surface_image(data=avg_map[np.newaxis, :], source_img=data['template'])
        nb.save(avg_img, avg_fn)
        print('Saved avg: {}'.format(avg_fn))

# Time
end_time = datetime.datetime.now()
print("\nStart time:\t{start_time}\nEnd time:\t{end_time}\nDuration:\t{dur}".format(
    start_time=start_time, end_time=end_time, dur=end_time - start_time))

# Permissions
print('Changing files permissions in {}/{}'.format(main_dir, project_dir))
os.system("chmod -Rf 771 {}/{}".format(main_dir, project_dir))
os.system("chgrp -Rf {} {}/{}".format(group, main_dir, project_dir))