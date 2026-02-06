"""
-----------------------------------------------------------------------------------------
create_mmp_rois_atlas.py
-----------------------------------------------------------------------------------------
Goal of the script:
Create 170k npz mmp atlas files and a tsv of rois numbers and names
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: data project directory
sys.argv[1]: codee project directory
sys.argv[2]: project name (correspond to directory)
-----------------------------------------------------------------------------------------
Output(s):
.npz masks and .tsv
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/[PROJECT]/analysis_code/atlas/
2. run python command
python prf_submit_gridfit_jobs.py [main directory] [project name] [subject] 
                                  [group] [server project]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/projects/pRF_analysis/analysis_code/atlas/
python create_mmp_rois_atlas.py /scratch/mszinte/data /home/mszinte/projects MotConf
python create_mmp_rois_atlas.py /scratch/mszinte/data /home/mszinte/projects RetinoMaps
python create_mmp_rois_atlas.py /scratch/mszinte/data /home/mszinte/projects amblyo_prf
python create_mmp_rois_atlas.py /scratch/mszinte/data /home/mszinte/projects centbids
-----------------------------------------------------------------------------------------
Written by Uriel Lascombes (uriel.lascombes@laposte.net)
Edited by Martin Szinte (mail@martinszinte.net)
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
import yaml
import cortex
import numpy as np
import pandas as pd


# Personals Import 
sys.path.append("{}/../utils".format(os.getcwd()))
from cifti_utils import from_170k_to_59k
from surface_utils import load_surface
from pycortex_utils import set_pycortex_config_file
from settings_utils import load_settings

# Get input
main_dir = sys.argv[1]
code_dir = sys.argv[2]
project_dir = sys.argv[3]

# Load settings
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../"))
settings_path = os.path.join(base_dir, project_dir, "settings.yml")
prf_settings_path = os.path.join(base_dir, project_dir, "prf-analysis.yml")
settings = load_settings([settings_path, prf_settings_path])
analysis_info = settings[0]

formats = analysis_info['formats']
extensions = analysis_info['extensions']
rois = analysis_info['rois']
rois_groups = analysis_info['rois_groups']

# Set pycortex db and colormaps
cortex_dir = "{}/{}/derivatives/pp_data/cortex".format(main_dir, project_dir)
set_pycortex_config_file(cortex_dir)

# Make mmp roi number correspondance TSV
#  Load csv 
sub_170k_cortex_dir = '{}/db/sub-170k'.format(cortex_dir)
mmp_csv_fn = '{}/HCP-MMP1_UniqueRegionList.csv'.format(sub_170k_cortex_dir)
mmp_csv = pd.read_csv(mmp_csv_fn)
mmp_csv['index_col'] = (mmp_csv.index + 1).astype('int32')
mmp_rois_nums_df = mmp_csv.loc[mmp_csv['LR'] == 'L', ['region', 'index_col']].rename(columns={'region': 'roi_name', 'index_col': 'roi_num'})
mmp_rois_nums_df['roi_name'] = mmp_rois_nums_df['roi_name'].astype(str)
mmp_tsv_fn = '{}/mmp_rois_numbers.tsv'.format(sub_170k_cortex_dir)
mmp_rois_nums_df.to_csv(mmp_tsv_fn, sep="\t", na_rep='NaN', index=False)

# Make MMP masks
_170k_dir_ ='{}/{}/derivatives/pp_data/cortex/db/sub-170k'.format(
    main_dir, project_dir)
_170k_fn = 'sub-170k_template.dtseries.nii'

mmp_rois_dir = '{}/{}/derivatives/pp_data/cortex/db/sub-170k/'.format(
    main_dir, project_dir)
mmp_rois_fn = 'mmp_atlas_rois_59k.npz'

# Load any 170k data to acces 59k mask 
img_170k, data_170k = load_surface('{}/{}'.format(_170k_dir_, _170k_fn))
rois_dict_59k = dict(np.load('{}/{}'.format(mmp_rois_dir, mmp_rois_fn)))

# Load surface 
surfs = [cortex.polyutils.Surface(*d) for d in cortex.db.get_surf('sub-170k', "flat")]
surf_lh, surf_rh = surfs[0], surfs[1]
lh_vert_num, rh_vert_num = surf_lh.pts.shape[0], surf_rh.pts.shape[0]

# Decompose the 170k datas 
results = from_170k_to_59k(img=img_170k, data=data_170k, return_concat_hemis=False, return_59k_mask=True)
data_59k_L = results['data_L']
data_59k_R = results['data_R']
mask_59k = results['mask_59k']

# Make hemi masks 59k
mask_59k_L = mask_59k[:lh_vert_num]
mask_59k_R = mask_59k[-rh_vert_num:]

# Make rois dict 
rois_dict_59k_L = {roi: data[:lh_vert_num] for roi, data in rois_dict_59k.items()}
rois_dict_59k_R = {roi: data[-rh_vert_num:] for roi, data in rois_dict_59k.items()}

# Make the final 170k Rois dict 
n_vertrx_170k = data_170k.shape[1]
rois_dict_170k = {}
for key, value in rois_dict_59k.items():
    rois_dict_170k[key] = value[mask_59k]
    vertex_to_add = n_vertrx_170k - len(value[mask_59k])
    rois_dict_170k[key] = np.concatenate((value[mask_59k], np.full(vertex_to_add, False)))
        
rois_dict_170k_L = {}
for key, value in rois_dict_59k_L.items():
    rois_dict_170k_L[key] = value[mask_59k_L]
    vertex_to_add = n_vertrx_170k - len(value[mask_59k_L])
    rois_dict_170k_L[key] = np.concatenate((value[mask_59k_L], np.full(vertex_to_add, False)))

rois_dict_170k_R = {}
for key, value in rois_dict_59k_R.items():
    rois_dict_170k_R[key] = value[mask_59k_R]
    vertex_to_add = n_vertrx_170k - (len(value[mask_59k_R])+len(value[mask_59k_L]))
    rois_dict_170k_R[key] = np.concatenate((np.full(len(value[mask_59k_L]), False),value[mask_59k_R], np.full(vertex_to_add, False)))

# Define groups and their corresponding ROIs in a dictionary
groups = dict(zip(rois, rois_groups))

# Brain 59k group
rois_dict_59k_group = {}
for group_name, rois_keys in groups.items():
    # Create an array of False for the group
    group_array = np.array([False] * len(next(iter(rois_dict_59k.values()))))
    for key in rois_keys:
        if key in rois_dict_59k:  # Check for key existence
            # Combine arrays with logical OR
            group_array = np.logical_or(group_array, rois_dict_59k[key])
    # Add the result under the group key
    rois_dict_59k_group[group_name] = group_array
    
# Left hemisphere 59k group
rois_dict_59k_group_L = {}
for group_name, rois_keys in groups.items():
    # Create an array of False for the group
    group_array = np.array([False] * len(next(iter(rois_dict_59k_L.values()))))
    for key in rois_keys:
        if key in rois_dict_59k_L:  # Check for key existence
            # Combine arrays with logical OR
            group_array = np.logical_or(group_array, rois_dict_59k_L[key])
    # Add the result under the group key
    rois_dict_59k_group_L[group_name] = group_array

# Right hemisphere 59k group 
rois_dict_59k_group_R = {}
for group_name, rois_keys in groups.items():
    # Create an array of False for the group
    group_array = np.array([False] * len(next(iter(rois_dict_59k_R.values()))))
    for key in rois_keys:
        if key in rois_dict_59k_R:  # Check for key existence
            # Combine arrays with logical OR
            group_array = np.logical_or(group_array, rois_dict_59k_R[key])
    # Add the result under the group key
    rois_dict_59k_group_R[group_name] = group_array
    
# Brain 170k group
rois_dict_170k_group = {}
for group_name, rois_keys in groups.items():
    # Create an array of False for the group
    group_array = np.array([False] * len(next(iter(rois_dict_170k.values()))))
    for key in rois_keys:
        if key in rois_dict_170k:  # Check for key existence
            # Combine arrays with logical OR
            group_array = np.logical_or(group_array, rois_dict_170k[key])
    # Add the result under the group key
    rois_dict_170k_group[group_name] = group_array

# Left hemisphere 170k group
rois_dict_170k_group_L = {}
for group_name, rois_keys in groups.items():
    # Create an array of False for the group
    group_array = np.array([False] * len(next(iter(rois_dict_170k_L.values()))))
    for key in rois_keys:
        if key in rois_dict_170k_L:  # Check for key existence
            # Combine arrays with logical OR
            group_array = np.logical_or(group_array, rois_dict_170k_L[key])
    # Add the result under the group key
    rois_dict_170k_group_L[group_name] = group_array

# Right hemisphere 170k group 
rois_dict_170k_group_R = {}
for group_name, rois_keys in groups.items():
    # Create an array of False for the group
    group_array = np.array([False] * len(next(iter(rois_dict_170k_R.values()))))
    for key in rois_keys:
        if key in rois_dict_170k_R:  # Check for key existence
            # Combine arrays with logical OR
            group_array = np.logical_or(group_array, rois_dict_170k_R[key])
    # Add the result under the group key
    rois_dict_170k_group_R[group_name] = group_array

# Save atlases
atlas_dir = '{}/{}/analysis_code/atlas'.format(code_dir, project_dir)
print(atlas_dir)
os.makedirs(atlas_dir, exist_ok=True)
print("Saving atlases in: {}...".format(atlas_dir))

# Brain 170k 
rois_170k_fn = 'mmp_atlas_rois_170k.npz'
np.savez('{}/{}'.format(atlas_dir, rois_170k_fn), **rois_dict_170k)

# Left hemisphere 170k 
rois_170k_fn_L = 'mmp_atlas_rois_170k_hemi-L.npz'
np.savez('{}/{}'.format(atlas_dir,rois_170k_fn_L), **rois_dict_170k_L)

# Right hemisphere 170k 
rois_170k_fn_R = 'mmp_atlas_rois_170k_hemi-R.npz'
np.savez('{}/{}'.format(atlas_dir,rois_170k_fn_R), **rois_dict_170k_R)

# Brain 59k 
rois_59k_fn = 'mmp_atlas_rois_59k.npz'
np.savez('{}/{}'.format(atlas_dir,rois_59k_fn), **rois_dict_59k)

# Left hemisphere 59k 
rois_59k_fn_L = 'mmp_atlas_rois_59k_hemi-L.npz'
np.savez('{}/{}'.format(atlas_dir,rois_59k_fn_L), **rois_dict_59k_L)

# right hemisphere 59k 
rois_59k_fn_R = 'mmp_atlas_rois_59k_hemi-R.npz'
np.savez('{}/{}'.format(atlas_dir,rois_59k_fn_R), **rois_dict_59k_R)

# Brain 59k group
rois_59k_group_fn = 'mmp_group_atlas_rois_59k.npz'
np.savez('{}/{}'.format(atlas_dir,rois_59k_group_fn), **rois_dict_59k_group)

# Brain 170k group
rois_170k_group_fn = 'mmp_group_atlas_rois_170k.npz'
np.savez('{}/{}'.format(atlas_dir,rois_170k_group_fn), **rois_dict_170k_group)

# Left hemisphere 170k group 
rois_170k_group_fn_L = 'mmp_group_atlas_rois_170k_hemi-L.npz'
np.savez('{}/{}'.format(atlas_dir,rois_170k_group_fn_L), **rois_dict_170k_group_L)

# Right hemisphere 170k group
rois_170k_group_fn_R = 'mmp_group_atlas_rois_170k_hemi-R.npz'
np.savez('{}/{}'.format(atlas_dir,rois_170k_group_fn_R), **rois_dict_170k_group_R)

# left hemisphere 59k group 
rois_59k_group_fn_L = 'mmp_group_atlas_rois_59k_hemi-L.npz'
np.savez('{}/{}'.format(atlas_dir,rois_59k_group_fn_L), **rois_dict_59k_group_L)

# Right hemisphere 59k group 
rois_59k_group_fn_R = 'mmp_group_atlas_rois_59k_hemi-R.npz'
np.savez('{}/{}'.format(atlas_dir,rois_59k_group_fn_R), **rois_dict_59k_group_R)

print('done')