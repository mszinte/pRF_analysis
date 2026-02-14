"""
-----------------------------------------------------------------------------------------
create_rois-drawn_npz.py
-----------------------------------------------------------------------------------------
Goal of the script:
Create 170k npz mmp atlas files and a tsv of rois numbers and names
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject name (e.g. sub-01)
sys.argv[4]: group (e.g. 327)
-----------------------------------------------------------------------------------------
Output(s):
.npz masks for of rois-drawn on fsnative
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/pRF_analysis/analysis_code/postproc/prf/postfit
2. run python command
python create_rois-drawn_npz.py [main directory] [project name] [subject] [group]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/projects/pRF_analysis/analysis_code/postproc/prf/postfit
python create_rois-drawn_npz.py /scratch/mszinte/data RetinoMaps sub-01 327
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
import cortex
import numpy as np

# Personal imports
sys.path.append("{}/../../../utils".format(os.getcwd()))
from pycortex_utils import set_pycortex_config_file

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
group = sys.argv[4]

# Set pycortex db and colormaps
cortex_dir = "{}/{}/derivatives/pp_data/cortex".format(main_dir, project_dir)
set_pycortex_config_file(cortex_dir)

overlay_roi_drawn_fn = '{}/db/{}/overlays_rois-drawn.svg'.format(cortex_dir, subject)

rois_drawn_brain_dict = cortex.get_roi_verts(subject=subject, mask=True, overlay_file=overlay_roi_drawn_fn)

# make hemispheres dict
surfs = [cortex.polyutils.Surface(*d) for d in cortex.db.get_surf(subject, "flat")]
surf_lh, surf_rh = surfs[0], surfs[1]
lh_vert_num, rh_vert_num = surf_lh.pts.shape[0], surf_rh.pts.shape[0]

rois_drawn_lh_dict = {}
rois_drawn_rh_dict = {}

for roi in rois_drawn_brain_dict.keys():
    rois_drawn_lh_dict[roi] =  rois_drawn_brain_dict[roi][:lh_vert_num]
    rois_drawn_rh_dict[roi] =  rois_drawn_brain_dict[roi][lh_vert_num:]

# Export masks as npz
print('Saving ROIs mask npz ...')
roi_dir = '{}/{}/derivatives/pp_data/cortex/db/{}/rois'.format(main_dir, project_dir, subject)
os.makedirs(roi_dir, exist_ok=True)
np.savez('{}/{}_hemi-L_fsnative_rois-drawn.npz'.format(roi_dir, subject), **rois_drawn_lh_dict)
np.savez('{}/{}_hemi-R_fsnative_rois-drawn.npz'.format(roi_dir, subject), **rois_drawn_rh_dict)
np.savez('{}/{}_fsnative_rois-drawn.npz'.format(roi_dir, subject), **rois_drawn_brain_dict)

# Define permission cmd
print('Changing files permissions in {}/{}'.format(main_dir, project_dir))
os.system("chmod -Rf 771 {}/{}".format(main_dir, project_dir))
os.system("chgrp -Rf {} {}/{}".format(group, main_dir, project_dir))




