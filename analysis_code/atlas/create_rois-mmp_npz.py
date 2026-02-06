"""
-----------------------------------------------------------------------------------------
create_rois-mmp_npz.py
-----------------------------------------------------------------------------------------
Goal of the script:
Create 170k npz mmp atlas files and a tsv of rois numbers and names
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: group of shared data (e.g. 327)
-----------------------------------------------------------------------------------------
Output(s):
.npz masks for roi-mmp and roi-group-mmp for both 91k and 170k template 
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/disks/meso_H/projects/[PROJECT]/analysis_code/atlas/
2. run python command
python create_rois-mmp_npz.py [main directory] [project name] [group]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/disks/meso_H/projects/pRF_analysis/analysis_code/atlas
python create_rois-mmp_npz.py /Users/uriel/disks/meso_shared RetinoMaps 327
-----------------------------------------------------------------------------------------
Written by Uriel Lascombes (uriel.lascombes@laposte.net)
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
import shutil
import numpy as np

# Personal imports
sys.path.append("{}/../../analysis_code/utils".format(os.getcwd()))
from pycortex_rois_utils import *
from settings_utils import load_settings
from pycortex_utils import set_pycortex_config_file

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
group = sys.argv[3]

# Load settings
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../"))
settings_path = os.path.join(base_dir, project_dir, "settings.yml")
prf_settings_path = os.path.join(base_dir, project_dir, "prf-analysis.yml")
settings = load_settings([settings_path, prf_settings_path])
analysis_info = settings[0]

rois_group_mmp = analysis_info['rois-group-mmp']
formats = analysis_info['formats_conversion']

# Set pycortex db and colormaps
cortex_dir = "{}/{}/derivatives/pp_data/cortex".format(main_dir, project_dir)
set_pycortex_config_file(cortex_dir)


for full_brain_format, cortex_format in formats.items():
    subject = 'sub-hcp{}'.format(cortex_format)
    
    # Load rois mmp mask
    roi_dir = '{}/{}/derivatives/pp_data/cortex/db/{}/rois'.format(
        main_dir, project_dir, subject)
    roi_mmp_brain_npz_fn = '{}/{}_{}_rois-mmp.npz'.format(
        roi_dir, subject, cortex_format)
    roi_mmp_lh_npz_fn = '{}/{}_{}_hemi-L_rois-mmp.npz'.format(
        roi_dir, subject, cortex_format)
    roi_mmp_rh_npz_fn = '{}/{}_{}_hemi-R_rois-mmp.npz'.format(
        roi_dir, subject, cortex_format)
    
    rois_mmp_masks_brain = dict(np.load(roi_mmp_brain_npz_fn))
    rois_mmp_masks_lh = dict(np.load(roi_mmp_lh_npz_fn))
    rois_mmp_masks_rh = dict(np.load(roi_mmp_rh_npz_fn))
    
    # Load ful brain mask to extend rois masks 
    full_brain_to_cortex_mask_fn = '{}/{}/derivatives/pp_data/cortex/db/{}/masks/{}_cortex_mask.npz'.format(
        main_dir, project_dir, subject, full_brain_format)
    full_brain_to_cortex_mask = dict(np.load(full_brain_to_cortex_mask_fn))['brain_mask']
    
    # make full brain masks
    rois_mmp_masks_full_lh_dict = {}
    rois_mmp_masks_full_rh_dict = {}
    rois_mmp_masks_full_brain_dict = {}
    
    cortex_idx = np.where(full_brain_to_cortex_mask)[0]
    
    n_vert_lh = next(iter(rois_mmp_masks_lh.values())).shape[0]
    n_vert_rh = next(iter(rois_mmp_masks_rh.values())).shape[0]
    
    lh_cortex_idx = cortex_idx[:n_vert_lh]
    rh_cortex_idx = cortex_idx[n_vert_lh:n_vert_lh + n_vert_rh]
    
    # --- Left hemisphere ---
    for roi_name, roi_mask in rois_mmp_masks_lh.items():
        roi_mask_full_lh = np.zeros_like(full_brain_to_cortex_mask, dtype=bool)
        roi_mask_full_lh[lh_cortex_idx] = roi_mask
        rois_mmp_masks_full_lh_dict[roi_name] = roi_mask_full_lh
    
    # --- Right hemisphere ---
    for roi_name, roi_mask in rois_mmp_masks_rh.items():
        roi_mask_full_rh = np.zeros_like(full_brain_to_cortex_mask, dtype=bool)
        roi_mask_full_rh[rh_cortex_idx] = roi_mask
        rois_mmp_masks_full_rh_dict[roi_name] = roi_mask_full_rh
    
    # --- Full brain ROIs ---
    for roi_name, roi_mask in rois_mmp_masks_brain.items():
        roi_mask_full_brain = np.zeros_like(full_brain_to_cortex_mask, dtype=bool)
        roi_mask_full_brain[cortex_idx] = roi_mask
        rois_mmp_masks_full_brain_dict[roi_name] = roi_mask_full_brain
    
    # Export masks as npz
    print('Saving roi mmp npz')
    np.savez('{}/{}_{}_hemi-L_rois-mmp.npz'.format(
        roi_dir, subject, full_brain_format), **rois_mmp_masks_full_lh_dict)
    np.savez('{}/{}_{}_hemi-R_rois-mmp.npz'.format(
        roi_dir, subject, full_brain_format), **rois_mmp_masks_full_rh_dict)
    np.savez('{}/{}_{}_rois-mmp.npz'.format(
        roi_dir, subject, full_brain_format), **rois_mmp_masks_full_brain_dict)
    
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
    print('Saving roi group mmp npz')
    np.savez('{}/{}_{}_hemi-L_rois-group-mmp.npz'.format(
        roi_dir, subject, cortex_format), **rois_group_mmp_masks_lh)
    np.savez('{}/{}_{}_hemi-R_rois-group-mmp.npz'.format(
        roi_dir, subject, cortex_format), **rois_group_mmp_masks_rh)
    np.savez('{}/{}_{}_rois-group-mmp.npz'.format(
        roi_dir, subject, cortex_format), **rois_group_mmp_masks_brain)
    
    # make full brain masks
    rois_group_mmp_masks_full_lh_dict = {}
    rois_group_mmp_masks_full_rh_dict = {}
    rois_group_mmp_masks_full_brain_dict = {}
    
    cortex_idx = np.where(full_brain_to_cortex_mask)[0]
    
    n_vert_lh = next(iter(rois_mmp_masks_lh.values())).shape[0]
    n_vert_rh = next(iter(rois_mmp_masks_rh.values())).shape[0]
    
    lh_cortex_idx = cortex_idx[:n_vert_lh]
    rh_cortex_idx = cortex_idx[n_vert_lh:n_vert_lh + n_vert_rh]
    
    # --- Left hemisphere ---
    for roi_name, roi_mask in rois_group_mmp_masks_lh.items():
        roi_mask_full_lh = np.zeros_like(full_brain_to_cortex_mask, dtype=bool)
        roi_mask_full_lh[lh_cortex_idx] = roi_mask
        rois_group_mmp_masks_full_lh_dict[roi_name] = roi_mask_full_lh
    
    # --- Right hemisphere ---
    for roi_name, roi_mask in rois_group_mmp_masks_rh.items():
        roi_mask_full_rh = np.zeros_like(full_brain_to_cortex_mask, dtype=bool)
        roi_mask_full_rh[rh_cortex_idx] = roi_mask
        rois_group_mmp_masks_full_rh_dict[roi_name] = roi_mask_full_rh
    
    # --- Full brain ROIs ---
    for roi_name, roi_mask in rois_group_mmp_masks_brain.items():
        roi_mask_full_brain = np.zeros_like(full_brain_to_cortex_mask, dtype=bool)
        roi_mask_full_brain[cortex_idx] = roi_mask
        rois_group_mmp_masks_full_brain_dict[roi_name] = roi_mask_full_brain
    
    # Export masks as npz
    print('Saving roi group mmp npz')
    np.savez('{}/{}_{}_hemi-L_rois-group-mmp.npz'.format(
        roi_dir, subject, full_brain_format), **rois_group_mmp_masks_full_lh_dict)
    np.savez('{}/{}_{}_hemi-R_rois-group-mmp.npz'.format(
        roi_dir, subject, full_brain_format), **rois_group_mmp_masks_full_rh_dict)
    np.savez('{}/{}_{}_rois-group-mmp.npz'.format(
        roi_dir, subject, full_brain_format), **rois_group_mmp_masks_full_brain_dict)
    
    # Project borders of rois on overlays
    overlays_fn = '{}/db/{}/overlays.svg'.format(cortex_dir, subject)
    overlays_rois_group_mmp_fn = '{}/db/{}/overlays_roi-group-mmp.svg'.format(cortex_dir, subject)

    # Copy overlays 
    shutil.copy(overlays_fn, overlays_rois_group_mmp_fn)

    # ROIs group MMP 
    rp = ROIpack(subject, '{}/{}_{}_rois-group-mmp.npz'.format(
        roi_dir, subject, cortex_format))
    rp.to_svg(filename=overlays_rois_group_mmp_fn, res=1024)
    
# Define permission cmd
print('Changing files permissions in {}/{}'.format(main_dir, project_dir))
os.system("chmod -Rf 771 {}/{}".format(main_dir, project_dir))
os.system("chgrp -Rf {} {}/{}".format(group, main_dir, project_dir))       
    
    
    
    