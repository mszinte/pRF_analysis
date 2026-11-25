"""
-----------------------------------------------------------------------------------------
siam_brainmask.py
-----------------------------------------------------------------------------------------
Goal of the script:
Create a brain mask from SIAM segmentation (without dura mater, vessels, skull and head)
and save it as a pre-computed fMRIPrep derivative (desc-brain_mask). Given to fMRIPrep
with --derivatives, it replaces the ANTs skull-stripping and is injected in FreeSurfer
as brainmask.mgz (replaces the manual edition of the brainmask)
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject (e.g. sub-001)
-----------------------------------------------------------------------------------------
Output(s):
one brain mask per session with a T1w (sub/ses/anat):
derivatives/siam/{subject}/{session}/anat/{subject}_{session}_desc-brain_mask.nii.gz
derivatives/siam/dataset_description.json
-----------------------------------------------------------------------------------------
To run:
0. run after siam_sbatch.py
1. cd to function
>> cd ~/projects/pRF_analysis/analysis_code/preproc/anatomical
2. run python command
python siam_brainmask.py [main] [project] [subject]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/projects/pRF_analysis/analysis_code/preproc/anatomical
python siam_brainmask.py /scratch/ulascombes/data CRASH sub-2100563718
-----------------------------------------------------------------------------------------
Written by Uriel Lascombes (uriel.lascombes@laposte.net)
-----------------------------------------------------------------------------------------
"""
# Stop warnings
import warnings
warnings.filterwarnings("ignore")

# General imports
import os
import sys
import json
import glob
import numpy as np
import nibabel as nb
from scipy import ndimage
import ipdb
deb = ipdb.set_trace

# Personal imports
sys.path.append("{}/../../utils".format(os.getcwd()))
from settings_utils import load_settings

# inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]

# Load settings
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../"))
settings_path = os.path.join(base_dir, project_dir, "settings.yml")
siam_settings_path = os.path.join(base_dir, project_dir, "siam-settings.yml")
settings = load_settings([settings_path, siam_settings_path])
analysis_info = settings[0]
brainmask_labels = analysis_info['siam_brainmask_labels']
brainmask_fill_holes = analysis_info['siam_brainmask_fill_holes']
tentorium_gap = analysis_info['siam_tentorium_gap']
gm_label, cergm_label = 1, 5

# Folders
bids_dir = f"{main_dir}/{project_dir}"
siam_dir = f"{main_dir}/{project_dir}/derivatives/siam"
siam_subject_dir = f"{siam_dir}/{subject}"

# Sessions with a T1w (sub/ses/anat), each session uses its own T1w
session_dirs = sorted(glob.glob(f"{bids_dir}/{subject}/ses-*/anat"))
sessions = [os.path.basename(os.path.dirname(session_dir)) for session_dir in session_dirs]
if len(sessions) == 0: sessions = [None]

for session in sessions:
    if session: 
        t1w_fns = sorted(glob.glob(f"{bids_dir}/{subject}/{session}/anat/{subject}_{session}_*T1w.nii.gz"))
        anat_dir = f"{siam_subject_dir}/{session}/anat"
        mask_fn = f"{anat_dir}/{subject}_{session}_desc-brain_mask.nii.gz"
    else: 
        t1w_fns = sorted(glob.glob(f"{bids_dir}/{subject}/anat/{subject}_*T1w.nii.gz"))
        anat_dir = f"{siam_subject_dir}/anat"
        mask_fn = f"{anat_dir}/{subject}_desc-brain_mask.nii.gz"
    if len(t1w_fns) == 0: 
        print(f"No T1w for {subject} {session}: skipped")
        continue
    if len(t1w_fns) > 1: 
        sys.exit(f"Expected 1 T1w for {subject} {session}, found {len(t1w_fns)}: {t1w_fns}")
    t1w_fn = t1w_fns[0]
    seg_fn = glob.glob(f"{siam_subject_dir}/siamV*_{os.path.basename(t1w_fn)}")
    if len(seg_fn) != 1: 
        sys.exit(f"Expected 1 SIAM segmentation of {os.path.basename(t1w_fn)}, found {len(seg_fn)}: run siam_sbatch.py")
    seg_fn = seg_fn[0]
    print(f"T1w: {t1w_fn}\nSIAM segmentation: {seg_fn}")
    os.makedirs(anat_dir, exist_ok=True)

    # Check that the segmentation is on the T1w grid (required by fMRIPrep)
    t1w_img = nb.load(t1w_fn)
    seg_img = nb.load(seg_fn)
    if t1w_img.shape[:3] != seg_img.shape[:3] or not np.allclose(t1w_img.affine, seg_img.affine, atol=1e-4):
        sys.exit("SIAM segmentation and T1w grids differ: rerun siam_sbatch.py on this T1w")

    # Brain mask from SIAM labels
    seg_data = np.asarray(seg_img.dataobj).astype(int)
    mask_data = np.isin(seg_data, brainmask_labels)

    # Gap between cerebral cortex and cerebellum (tentorium) to avoid pial leaking
    if tentorium_gap > 0:
        struct = ndimage.generate_binary_structure(3, 1)
        near_cortex = ndimage.binary_dilation(seg_data == gm_label, struct, iterations=tentorium_gap)
        mask_data &= ~(near_cortex & (seg_data == cergm_label))

    # Keep the largest connected component
    labeled, nb_components = ndimage.label(mask_data)
    if nb_components > 1:
        sizes = ndimage.sum(mask_data, labeled, range(1, nb_components + 1))
        mask_data = labeled == (np.argmax(sizes) + 1)

    # Fill holes
    if brainmask_fill_holes:
        mask_data = ndimage.binary_fill_holes(mask_data)

    # Save mask on the T1w grid (RAS)
    mask_img = nb.Nifti1Image(mask_data.astype(np.uint8), t1w_img.affine, t1w_img.header)
    mask_img.set_data_dtype(np.uint8)
    # fMRIPrep reorients the T1w to RAS (niworkflows Conform) before applying the mask
    mask_img = nb.as_closest_canonical(mask_img)
    nb.save(mask_img, mask_fn)
    vox_vol = np.prod(t1w_img.header.get_zooms()[:3])
    print(f"Saved {mask_fn} ({mask_data.sum() * vox_vol / 1000:.1f} cm3)")

# Dataset description (derivatives folder read by fMRIPrep --derivatives)
dd_fn = f"{siam_dir}/dataset_description.json"
if not os.path.exists(dd_fn):
    dd = {"Name": "SIAM brain masks", "BIDSVersion": "1.10.0", "DatasetType": "derivative",
          "GeneratedBy": [{"Name": "SIAM", "CodeURL": "https://github.com/romainVala/SIAM"}]}
    with open(dd_fn, 'w') as f:
        json.dump(dd, f, indent=4)
