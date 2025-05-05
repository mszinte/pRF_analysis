"""
-----------------------------------------------------------------------------------------
stats_corr.py
-----------------------------------------------------------------------------------------
Goal of the script:
Compute the zfisher file of resting state corrected for fdr_alpha level
Compute the winner index across roi (1:V1 (not 0 !!!); 2: V2; ...; 12: mPCS)
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject
-----------------------------------------------------------------------------------------
Output(s):
gifti files
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/[PROJECT]/analysis_code/postproc/rest/
2. run python command
python stats_corr.py [main directory] [project name] [subject name] 
-----------------------------------------------------------------------------------------
Exemple:
cd ~/projects/pRF_analysis/analysis_code/postproc/rest/
python stats_corr.py /scratch/mszinte/data RetinoMaps sub-01
-----------------------------------------------------------------------------------------
Written by Martin Szinte (martin.szinte@gmail.com)
Edited by Marco Bedini (marco.bedini@univ-amu.fr)
-----------------------------------------------------------------------------------------
"""

# Stop warnings
import warnings
warnings.filterwarnings("ignore")

# Debug
import ipdb
deb = ipdb.set_trace

import json
import numpy as np
import os, sys
import nibabel as nb

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
use_fisher = False
seed_debug = False
roi_debug = True
avg_method = 'median' # 'mean' or 'median'

sys.path.append("{}/../../utils".format(os.getcwd()))
from surface_utils import load_surface, make_surface_image
from maths_utils import linear_regression_surf

# Functions
def find_winner(pvalue_fdr_alpha):
    winner_indices = []

    # Iterate over each column to find the winner index
    for col in range(pvalue_fdr_alpha.shape[1]):
        column_data = pvalue_fdr_alpha[:, col]
        
        # Check if all values are NaN
        if np.all(np.isnan(column_data)):
            winner_indices.append(-1)  # Append -1 for all-NaN columns
        else:
            winner_index = np.nanargmax(column_data) + 1
            winner_indices.append(winner_index)

    # Convert winner_indices to a numpy array
    winner_indices = np.array(winner_indices)

    # Get the corresponding values, handling NaNs
    winner_values = np.array([pvalue_fdr_alpha[idx, col] if idx != -1 else np.nan 
                              for col, idx in enumerate(winner_indices)])

    # Change -1 in winner_indices to NaN
    winner_indices = np.where(winner_indices == -1, np.nan, winner_indices)

    # Stack the results
    data_with_winners = np.vstack((pvalue_fdr_alpha, winner_indices, winner_values))

    return data_with_winners

def save_data(data_array, template_img, data_fn):
    data_array_img = make_surface_image(data=data_array, source_img=template_img)

    wta_darray = nb.gifti.GiftiDataArray(data_array.astype(np.float32))
    nb.save(nb.GiftiImage(darrays=[wta_darray]), data_fn)

# Define settings
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../"))
settings_path = os.path.join(base_dir, project_dir, "settings.json")

with open(settings_path) as f:
    json_s = f.read()
    settings = json.loads(json_s)

task_names = settings['task_names']
task_name = task_names[0]
rois = settings['rois']
fdr_alpha = settings['fdr_alpha']

# Define finename
project = settings['project_name']
pp_data_folder = 'derivatives/pp_data'
atlas_folder = f'{base_dir}/analysis_code/atlas'
vertex_num_91k = 64984

# Create output folder
coor_dir = f'{main_dir}/{project_dir}/{pp_data_folder}/{subject}/91k/{task_name}/corr'
os.makedirs(coor_dir, exist_ok=True)

# Load resting-state dense timeseries
timeseries_fn = f'{main_dir}/{project_dir}/{pp_data_folder}/{subject}/91k/{task_name}/timeseries/{subject}_ses-01_task-{task_name}_space-fsLR_den-91k_desc-denoised_bold.dtseries.nii'
timeseries_img, timeseries_data = load_surface(timeseries_fn)
timeseries_data = timeseries_data[:, 0:vertex_num_91k]
print(f'Timeseries data shape: {timeseries_data.shape}')

# Load seed mask file as array (12,64984) from V1 to mPCS
for roi_num, roi in enumerate(rois):
    seed_mask_lh_fn = f'{main_dir}/{project_dir}/{pp_data_folder}/{subject}/91k/{task_name}/seed/{subject}_91k_intertask_Sac_Pur_vision-pursuit-saccade_lh_{roi}.shape.gii'
    seed_mask_lh_img, seed_mask_lh_data = load_surface(seed_mask_lh_fn)
    seed_mask_rh_fn = f'{main_dir}/{project_dir}/{pp_data_folder}/{subject}/91k/{task_name}/seed/{subject}_91k_intertask_Sac_Pur_vision-pursuit-saccade_rh_{roi}.shape.gii'
    seed_mask_rh_img, seed_mask_rh_data = load_surface(seed_mask_rh_fn)

    roi_seed_mask_data = np.hstack((seed_mask_lh_data, seed_mask_rh_data))

    if roi_num == 0:
        seed_mask_data = roi_seed_mask_data
    else:
        seed_mask_data = np.vstack((seed_mask_data, roi_seed_mask_data))

print(f'Seed mask data shape: {seed_mask_data.shape}')

# Load leaveout macro-region target masks as array (12,64984) from V1 to mPCS
for roi_num, roi in enumerate(rois):
    leaveout_target_mask_lh_fn = f'{atlas_folder}/macro_regions/leaveout/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_discarded_leaveout_bin_lh_{roi}.shape.gii'
    leaveout_target_mask_lh_img, leaveout_target_mask_lh_data = load_surface(leaveout_target_mask_lh_fn)
    leaveout_target_mask_rh_fn = f'{atlas_folder}/macro_regions/leaveout/atlas-Glasser_space-fsLR_den-32k_filtered_ROIs_discarded_leaveout_bin_rh_{roi}.shape.gii'
    leaveout_target_mask_rh_img, leaveout_target_mask_rh_data = load_surface(leaveout_target_mask_rh_fn)

    roi_leaveout_target_mask_data = np.hstack((leaveout_target_mask_lh_data, leaveout_target_mask_rh_data))

    if roi_num == 0:
        leaveout_target_mask_data = roi_leaveout_target_mask_data
    else:
        leaveout_target_mask_data = np.vstack((leaveout_target_mask_data, roi_leaveout_target_mask_data))

print(f'Leave-out target mask data shape: {leaveout_target_mask_data.shape}')


# ROWS : slope / intercept / rvalue / fisher-z / pvalue / stderr / trs / pvalue adjusted for alpha1 / pvalue adjusted for alpha2`
rvalue_row = 2
pvalue_fdr_alpha1_row = -2
pvalue_fdr_alpha2_row = -1
pvalue_fdr_alpha1 = np.zeros((len(rois), vertex_num_91k)) * np.nan
pvalue_fdr_alpha2 = np.zeros((len(rois), vertex_num_91k)) * np.nan

# Correlation
for roi_num, roi in enumerate(rois):

    
    # Get the seed data as before
    seed_data = timeseries_data[:, seed_mask_data[roi_num,:].astype(bool)]
    
    # Create a copy of the original timeseries_data
    target_data = np.full(timeseries_data.shape, np.nan)  # Initialize with NaN
    
    # Apply the mask to the target_data
    target_data[:, leaveout_target_mask_data[roi_num, :].astype(bool)] = timeseries_data[:, leaveout_target_mask_data[roi_num,:].astype(bool)]
    

    # debugging: to accelerate debugging take first 10
    if seed_debug == True:
        seed_data = seed_data[:,:3]

    # debugging: to accelerate debugging take first 2 rois
    if roi_debug and roi_num >= 2:
        break

    print(f'{roi} Seed timeseries: {seed_data.shape}')
    print(f'{roi} Target timeseries: {target_data.shape}')

    for target_col in range(target_data.shape[1]):
        target_column_data = target_data[:, target_col]
        
        target_column_repeated = np.tile(target_column_data[:, np.newaxis], (1, seed_data.shape[1]))
        
        results = linear_regression_surf(bold_signal=target_column_repeated,
                                         model_prediction=seed_data,
                                         correction='fdr_tsbh', 
                                         alpha=fdr_alpha,
                                         use_fisher=use_fisher
                                        )
        
        # median of r values for significant correlations as a function of level of correction (fdr_alpha1, fdr_alpha2)
        if avg_method == 'median':
            pvalue_fdr_alpha1[roi_num, target_col] = np.median(results[rvalue_row, :][results[pvalue_fdr_alpha1_row,:]<=fdr_alpha[0]])
            pvalue_fdr_alpha2[roi_num, target_col] = np.median(results[rvalue_row, :][results[pvalue_fdr_alpha2_row,:]<=fdr_alpha[0]])
        elif avg_method == 'mean':
            pvalue_fdr_alpha1[roi_num, target_col] = np.mean(results[rvalue_row, :][results[pvalue_fdr_alpha1_row,:]<=fdr_alpha[0]])
            pvalue_fdr_alpha2[roi_num, target_col] = np.mean(results[rvalue_row, :][results[pvalue_fdr_alpha2_row,:]<=fdr_alpha[0]])

# Make the files with all 12 ROIS seeds, the winner index (13) and the winner r value (14)
pvalue_fdr_alpha1_with_winners = find_winner(pvalue_fdr_alpha1)
pvalue_fdr_alpha2_with_winners = find_winner(pvalue_fdr_alpha2)

# FDR-ALPHA1 : 0.05
if use_fisher: stats_txt = f'{avg_method}_fisher-z'
else: stats_txt = f'{avg_method}_pearson-r'

pvalue_fdr_alpha1_with_winners_lh = pvalue_fdr_alpha1_with_winners[:,:32492].T
pvalue_fdr_alpha1_with_winners_rh = pvalue_fdr_alpha1_with_winners[:,-32492:].T

pvalue_fdr_alpha1_with_winners_lh_fn = f'{coor_dir}/{subject}_ses-01_task-{task_name}_space-fsLR_den-91k_desc-full_corr_{stats_txt}_fdr_alpha-05_L.shape.gii'
pvalue_fdr_alpha1_with_winners_rh_fn = f'{coor_dir}/{subject}_ses-01_task-{task_name}_space-fsLR_den-91k_desc-full_corr_{stats_txt}_fdr_alpha-05_R.shape.gii'
save_data(pvalue_fdr_alpha1_with_winners_lh, seed_mask_lh_img, pvalue_fdr_alpha1_with_winners_lh_fn)
save_data(pvalue_fdr_alpha1_with_winners_rh, seed_mask_rh_img, pvalue_fdr_alpha1_with_winners_rh_fn)

# FDR-ALPHA2 : 0.01
pvalue_fdr_alpha2_with_winners_lh = pvalue_fdr_alpha2_with_winners[:,:32492].T
pvalue_fdr_alpha2_with_winners_rh = pvalue_fdr_alpha2_with_winners[:,-32492:].T

pvalue_fdr_alpha2_with_winners_lh_fn = f'{coor_dir}/{subject}_ses-01_task-{task_name}_space-fsLR_den-91k_desc-full_corr_{stats_txt}_fdr_alpha-01_L.shape.gii'
pvalue_fdr_alpha2_with_winners_rh_fn = f'{coor_dir}/{subject}_ses-01_task-{task_name}_space-fsLR_den-91k_desc-full_corr_{stats_txt}_fdr_alpha-01_R.shape.gii'
save_data(pvalue_fdr_alpha2_with_winners_lh, seed_mask_lh_img, pvalue_fdr_alpha2_with_winners_lh_fn)
save_data(pvalue_fdr_alpha2_with_winners_rh, seed_mask_rh_img, pvalue_fdr_alpha2_with_winners_rh_fn)
