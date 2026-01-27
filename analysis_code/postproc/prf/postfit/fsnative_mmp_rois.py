"""
-----------------------------------------------------------------------------------------
fsnative_mmp_rois.py
-----------------------------------------------------------------------------------------
Goal of the script:
Load freesurfer and use it to project MMP from fsaverage to fsnative 
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject name
sys.argv[4]: group of shared data (e.g. 327)
-----------------------------------------------------------------------------------------
Output(s):
None
-----------------------------------------------------------------------------------------
To run:
1. cd to function
cd ~/disks/meso_H/projects/pRF_analysis/analysis_code/preproc/anatomical/
2. run python command
python fsnative_mmp_rois.py [main directory] [project name] [subject name] [group]
-----------------------------------------------------------------------------------------
Executions:
RUN LOCALY ! 
cd ~/disks/meso_H/projects/pRF_analysis/analysis_code/preproc/anatomical/
python fsnative_mmp_rois.py ~/disks/meso_shared RetinoMaps sub-01 327
-----------------------------------------------------------------------------------------
Written by Uriel Lascombes (uriel.lascombes@laposte.net)
-----------------------------------------------------------------------------------------
"""

# stop warnings
import warnings
warnings.filterwarnings("ignore")

# Debug 
import ipdb
deb = ipdb.set_trace

# General imports
import os
import sys
import shutil
import numpy as np
import nibabel as nb

# Personal imports
sys.path.append("{}/../../utils".format(os.getcwd()))
from settings_utils import load_settings
from pycortex_utils import set_pycortex_config_file
from pycortex_rois_utils import *

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
group = sys.argv[4]

# Load settings
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../"))
settings_path = os.path.join(base_dir, project_dir, "settings.yml")
prf_settings_path = os.path.join(base_dir, project_dir, "prf-analysis.yml")
settings = load_settings([settings_path, prf_settings_path])
analysis_info = settings[0]
rois_group_mmp = analysis_info['rois-group-mmp']

# Set pycortex db and colormaps
cortex_dir = "{}/{}/derivatives/pp_data/cortex".format(main_dir, project_dir)
set_pycortex_config_file(cortex_dir)

# Load fsnative MMP labels
mmp_fsnative_rois_fn = '{}/{}/derivatives/fmriprep/freesurfer/{}/label'.format(main_dir, project_dir, subject)
labels_lh, ctab_lh, names_lh = nb.freesurfer.read_annot("{}/lh.HCPMMP1.annot".format(mmp_fsnative_rois_fn))
labels_rh, ctab_rh, names_rh = nb.freesurfer.read_annot("{}/rh.HCPMMP1.annot".format(mmp_fsnative_rois_fn))

atlas_fn = '{}/{}/derivatives/pp_data/{}/fsnative/rois/atlas'.format(main_dir, project_dir, subject)
os.makedirs(atlas_fn, exist_ok=True)

# left hemi
id_to_name_lh = {i: names_lh[i].decode() for i in range(len(names_lh))}
rois_mmp_masks_lh = {}

# Clean ROIs labels (delete L_ _ROI to keep only ROI name)
for label_id, roi_name in id_to_name_lh.items():
    if roi_name == "???":
        continue
    if roi_name.startswith("L_"):
        roi_clean = roi_name[len("L_"):]
    else:
        roi_clean = roi_name
    if roi_clean.endswith("_ROI"):
        roi_clean = roi_clean[:-len("_ROI")]

    rois_mmp_masks_lh[roi_clean] = (labels_lh == label_id)
    rois_mmp_masks_lh[roi_clean] = rois_mmp_masks_lh[roi_clean].squeeze()

# Export masks as npz
rois_mmp_lh_fn = '{}/{}_rois-mmp_hemi-L.npz'.format(atlas_fn, subject)
print('saving {}'.format(rois_mmp_lh_fn))
np.savez(rois_mmp_lh_fn, **rois_mmp_masks_lh)

# right hemi
id_to_name_rh = {i: names_rh[i].decode() for i in range(len(names_rh))}
rois_mmp_masks_rh = {}

# Clean ROIs labels (delete L_ _ROI to keep only ROI name)
for label_id, roi_name in id_to_name_rh.items():
    if roi_name == "???":
        continue
    if roi_name.startswith("R_"):
        roi_clean = roi_name[len("R_"):]
    else:
        roi_clean = roi_name
    if roi_clean.endswith("_ROI"):
        roi_clean = roi_clean[:-len("_ROI")]

    rois_mmp_masks_rh[roi_clean] = (labels_rh == label_id)
    rois_mmp_masks_rh[roi_clean] = rois_mmp_masks_rh[roi_clean].squeeze()

# Export masks as npz
rois_mmp_rh_fn = '{}/{}_rois-mmp_hemi-R.npz'.format(atlas_fn, subject)
print('saving {}'.format(rois_mmp_rh_fn))
np.savez(rois_mmp_rh_fn, **rois_mmp_masks_rh)

# brain
rois_mmp_masks_brain = {}
for key in rois_mmp_masks_lh.keys():
    rois_mmp_masks_brain[key] = np.concatenate([rois_mmp_masks_lh[key], rois_mmp_masks_rh[key]], axis=0).squeeze()

# Export masks as npz
rois_mmp_brain_fn = '{}/{}_rois-mmp.npz'.format(atlas_fn, subject)
print('saving {}'.format(rois_mmp_brain_fn))
np.savez(rois_mmp_brain_fn, **rois_mmp_masks_brain)

# make rois-mmp-group dict 
rois_group_mmp_masks_lh = {}
rois_group_mmp_masks_rh = {}
rois_group_mmp_masks_brain = {}

for group_name, roi_list in rois_group_mmp.items():
    masks_lh = [rois_mmp_masks_lh[roi] for roi in roi_list if roi in rois_mmp_masks_lh]
    masks_rh = [rois_mmp_masks_rh[roi] for roi in roi_list if roi in rois_mmp_masks_rh]
    masks_brain = [rois_mmp_masks_brain[roi] for roi in roi_list if roi in rois_mmp_masks_brain]

    if len(masks_lh) > 0:
        rois_group_mmp_masks_lh[group_name] = np.logical_or.reduce(masks_lh)

    if len(masks_rh) > 0:
        rois_group_mmp_masks_rh[group_name] = np.logical_or.reduce(masks_rh)
        
    if len(masks_brain) > 0:
        rois_group_mmp_masks_brain[group_name] = np.logical_or.reduce(masks_brain)
        
# Export masks as npz
rois_group_mmp_lh_fn = '{}/{}_rois-group-mmp_hemi-L.npz'.format(atlas_fn, subject)
print('saving {}'.format(rois_group_mmp_lh_fn))
np.savez(rois_group_mmp_lh_fn, **rois_group_mmp_masks_lh)

rois_group_mmp_rh_fn = '{}/{}_rois-group-mmp_hemi-R.npz'.format(atlas_fn, subject)
print('saving {}'.format(rois_group_mmp_rh_fn))
np.savez(rois_group_mmp_rh_fn, **rois_group_mmp_masks_rh)

rois_group_mmp_brain_fn = '{}/{}_rois-group-mmp.npz'.format(atlas_fn, subject)
print('saving {}'.format(rois_group_mmp_brain_fn))
np.savez(rois_group_mmp_brain_fn, **rois_group_mmp_masks_brain)

# Project borders of rois on overlays
overlays_fn = '{}/db/{}/overlays.svg'.format(cortex_dir, subject)
overlays_rois_mmp_fn = '{}/db/{}/overlays_roi-mmp.svg'.format(cortex_dir, subject)
overlays_rois_group_mmp_fn = '{}/db/{}/overlays_roi-group-mmp.svg'.format(cortex_dir, subject)

# Copy overlays 
shutil.copy(overlays_fn, overlays_rois_mmp_fn)
shutil.copy(overlays_fn, overlays_rois_group_mmp_fn)

# ROIs MMP
rp = ROIpack(subject, rois_mmp_brain_fn)
rp.to_svg(filename=overlays_rois_mmp_fn)

# ROIs group MMP 
rp = ROIpack(subject, rois_group_mmp_brain_fn)
rp.to_svg(filename=overlays_rois_group_mmp_fn)

# # Define permission cmd
# print('Changing files permissions in {}/{}'.format(main_dir, project_dir))
# os.system("chmod -Rf 771 {}/{}".format(main_dir, project_dir))
# os.system("chgrp -Rf {} {}/{}".format(group, main_dir, project_dir))