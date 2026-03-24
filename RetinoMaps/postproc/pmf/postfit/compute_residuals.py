"""
-----------------------------------------------------------------------------------------
compute_residuals.py
-----------------------------------------------------------------------------------------
Goal of the script:
Compute residuals between BOLD data and PMF gauss predictions.
Finds all necessary files automatically from settings — no separate job submission needed.
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject name (e.g. sub-01)
-----------------------------------------------------------------------------------------
Output(s):
Residual files saved to:
  - {pp_dir}/{subject}/170k/pmf/pmf_residuals/
  - {pp_dir}/{subject}/fsnative/pmf/pmf_residuals/
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/pRF_analysis/RetinoMaps/postproc/pmf/postfit
2. run python command
>> python compute_residuals.py [main directory] [project name] [subject] [n_jobs]
-----------------------------------------------------------------------------------------
Example:
python compute_residuals.py /scratch/mszinte/data RetinoMaps sub-01 
-----------------------------------------------------------------------------------------
Written by Sina Kling (sina.kling@outlook.de)
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
import numpy as np
import nibabel as nb

# Personal imports
script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.abspath(os.path.join(script_dir, "../../../../analysis_code/utils")))
from surface_utils import make_surface_image, load_surface
from pycortex_utils import set_pycortex_config_file
from settings_utils import load_settings


# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]

# Load settings
base_dir = os.path.abspath(os.path.join(script_dir, "../../../../"))
settings_path = os.path.join(base_dir, project_dir, "pmf-settings.yml")
prf_settings_path = os.path.join(base_dir, project_dir, "prf-analysis.yml")
settings = load_settings([settings_path, prf_settings_path])
analysis_info = settings[0]

pmf_task_names = analysis_info['pmf_task_names']
preproc_prep = analysis_info['preproc_prep']
filtering = analysis_info['filtering']
normalization = analysis_info['normalization']
avg_methods = analysis_info['avg_methods']

# Set pycortex config
cortex_dir = "{}/{}/derivatives/pp_data/cortex".format(main_dir, project_dir)
set_pycortex_config_file(cortex_dir)

pp_dir = "{}/{}/derivatives/pp_data".format(main_dir, project_dir)


def find_bold_files(pp_dir, subject, pmf_task_names, preproc_prep, filtering,
                    normalization, avg_methods):
    """Find concat-averaged BOLD files for the pmf task (both 170k and fsnative)."""
    bold_fns = []
    for avg_method in avg_methods:
        if avg_method != 'concat':
            continue

        gii_pattern = (
            "{}/{}/fsnative/func/{}_{}_{}_{}/"
            "*_task-{}*{}*.func.gii"
        ).format(
            pp_dir, subject,
            preproc_prep, filtering, normalization, avg_method,
            pmf_task_names[0], avg_method
        )
        nii_pattern = (
            "{}/{}/170k/func/{}_{}_{}_{}/"
            "*_task-{}*{}*.dtseries.nii"
        ).format(
            pp_dir, subject,
            preproc_prep, filtering, normalization, avg_method,
            pmf_task_names[0], avg_method
        )

        bold_fns.extend(glob.glob(gii_pattern))
        bold_fns.extend(glob.glob(nii_pattern))

    return bold_fns


def find_pmf_gauss_pred(pp_dir, subject, bold_fn, pmf_task_names, preproc_prep,
                        filtering, normalization):
    """
    Find the PMF gauss prediction file matching the format and hemisphere of bold_fn.
    These are the outputs of the pmf fit step (pmf_pred_from_prf.py).
    """
    if bold_fn.endswith('.nii'):
        pattern = (
            "{}/{}/170k/pmf/fit/"
            "{}_task-{}_{}_{}_{}_*pmf-gauss_pred*.nii"
        ).format(
            pp_dir, subject,
            subject, pmf_task_names[0],
            preproc_prep, filtering, normalization
        )

    elif bold_fn.endswith('.gii'):
        fn_basename = os.path.basename(bold_fn)
        if 'hemi-L' in fn_basename:
            hemi = 'hemi-L'
        elif 'hemi-R' in fn_basename:
            hemi = 'hemi-R'
        else:
            raise ValueError(
                "Could not determine hemisphere from filename: {}".format(fn_basename)
            )

        pattern = (
            "{}/{}/fsnative/pmf/fit/"
            "{}_task-{}_{}_{}_{}_{}_*pmf-gauss_pred*.gii"
        ).format(
            pp_dir, subject,
            subject, pmf_task_names[0], hemi,
            preproc_prep, filtering, normalization
        )

    else:
        raise ValueError("Unsupported file extension: {}".format(bold_fn))

    matches = glob.glob(pattern)
    if not matches:
        raise FileNotFoundError(
            "No PMF gauss prediction file found.\n  Pattern: {}".format(pattern)
        )

    return matches[0]


def make_residual_output_path(bold_fn, pp_dir, subject):
    """Derive the output path and filename from the BOLD input file."""
    if bold_fn.endswith('.nii'):
        residuals_dir = "{}/{}/170k/pmf/pmf_residuals".format(pp_dir, subject)
    else:
        residuals_dir = "{}/{}/fsnative/pmf/pmf_residuals".format(pp_dir, subject)

    os.makedirs(residuals_dir, exist_ok=True)

    residual_fn = os.path.basename(bold_fn).replace('bold', 'pmf-residuals')
    return os.path.join(residuals_dir, residual_fn)


bold_fns = find_bold_files(
    pp_dir, subject, pmf_task_names, preproc_prep, filtering, normalization, avg_methods
)

if not bold_fns:
    raise FileNotFoundError(
        "No BOLD input files found for subject '{}', task '{}'.".format(
            subject, pmf_task_names[0]
        )
    )

print("Found {} BOLD file(s) to process.".format(len(bold_fns)))

for bold_fn in bold_fns:
    print("\n" + "=" * 72)
    print("BOLD:       {}".format(bold_fn))

    # Find matching prediction file
    pmf_pred_fn = find_pmf_gauss_pred(
        pp_dir, subject, bold_fn, pmf_task_names, preproc_prep, filtering, normalization
    )
    print("PMF pred:   {}".format(pmf_pred_fn))

    # Derive output path
    residual_fullfn = make_residual_output_path(bold_fn, pp_dir, subject)
    print("Output:     {}".format(residual_fullfn))

    # Load data
    print("Loading PMF prediction data...")
    pmf_pred_img, pmf_pred_data = load_surface(fn=pmf_pred_fn)

    print("Loading BOLD data...")
    bold_img, bold_data = load_surface(fn=bold_fn)

    print("  BOLD shape:     {}".format(bold_data.shape))
    print("  PMF pred shape: {}".format(pmf_pred_data.shape))

    # Sanity check
    if bold_data.shape != pmf_pred_data.shape:
        raise ValueError(
            "Shape mismatch — BOLD {} vs PMF pred {}.".format(
                bold_data.shape, pmf_pred_data.shape
            )
        )

    # NaN mask: exclude voxels that are NaN in either dataset
    nan_mask = (
        np.isnan(bold_data).any(axis=0) |
        np.isnan(pmf_pred_data).any(axis=0)
    )
    clean_vox = np.where(~nan_mask)[0]
    print("  Clean voxels:   {} / {}".format(len(clean_vox), bold_data.shape[1]))

    # Compute residuals on clean voxels only
    residual_data = np.full_like(bold_data, np.nan)
    residual_data[:, clean_vox] = bold_data[:, clean_vox] - pmf_pred_data[:, clean_vox]

    # Save
    residual_img = make_surface_image(data=residual_data, source_img=bold_img)
    nb.save(residual_img, residual_fullfn)
    print("Saved: {}".format(residual_fullfn))

print("\nDone — all residuals computed.")