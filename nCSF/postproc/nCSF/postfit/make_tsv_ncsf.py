"""
-----------------------------------------------------------------------------------------
make_tsv_ncsf.py
-----------------------------------------------------------------------------------------
Goal of the script:
Create TSV file with all css analysis output
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject name (e.g. sub-01)
sys.argv[4]: server group (e.g. 327)
-----------------------------------------------------------------------------------------
Output(s):
TSV file
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/pRF_analysis/nCSF/postproc/nCSF/postfit/
2. run python command
>> python make_tsv_ncsf.py [main directory] [project name] 
                          [subject num] [group] [analysis folder - optional]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/projects/pRF_analysis/nCSF/postproc/nCSF/postfit/
python make_tsv_ncsf.py /scratch/mszinte/data nCSF sub-01 327
-----------------------------------------------------------------------------------------
Written by Uriel Lascombes (uriel.lascombes@laposte.net) 
and Martin Szinte (martin.szinte@gmail.com)
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
import numpy as np
import pandas as pd

# Personal import
sys.path.append("{}/../../../../analysis_code/utils".format(os.getcwd()))
from surface_utils import load_surface
from settings_utils import load_settings
from pycortex_utils import get_rois, set_pycortex_config_file

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
group = sys.argv[4]

# Load settings
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../../"))
settings_path = os.path.join(base_dir, project_dir, "settings.yml")
nCSF_settings_path = os.path.join(base_dir, project_dir, "nCSF-analysis.yml")
prf_settings_path = os.path.join(base_dir, project_dir, "prf-analysis.yml")
settings = load_settings([settings_path, nCSF_settings_path, prf_settings_path])
analysis_info = settings[0]

formats = analysis_info['formats']
extensions = analysis_info['extensions']
maps_names_ncsf = analysis_info['maps_names_ncsf']
maps_names_ncsf_stats = analysis_info['maps_names_ncsf_stats']
preproc_prep = analysis_info['preproc_prep']
filtering = analysis_info['filtering']
normalization = analysis_info['normalization']
avg_methods = analysis_info['avg_methods']
rois_methods = analysis_info['rois_methods']
pycortex_subject_template = analysis_info['pycortex_subject_template']
ncsf_task_name = analysis_info['nCSF_task_name']

# Set pycortex db and colormaps
cortex_dir = "{}/{}/derivatives/pp_data/cortex".format(main_dir, project_dir)
set_pycortex_config_file(cortex_dir)

# Loop over format
for avg_method in avg_methods:
    if 'loo' in avg_method:
        maps_names_ncsf = maps_names_ncsf + ['ncsf_loo_rsq']
    else: 
        maps_names_ncsf = maps_names_ncsf

    maps_names = maps_names_ncsf + maps_names_ncsf_stats
    
    for format_, extension in zip(formats, extensions):
        
        # define list of rois for each format
        rois_methods_format = rois_methods[format_]
        for rois_method_format in rois_methods_format:

            if rois_method_format == 'rois-drawn':
                rois = analysis_info[rois_method_format]
            elif rois_method_format == 'rois-group-mmp':
                rois = list(analysis_info[rois_method_format].keys())

            ncsf_dir = "{}/{}/derivatives/pp_data/{}/{}/ncsf".format(
                main_dir, project_dir, subject, format_)
            
            if not os.path.isdir(ncsf_dir):
                print(f"[SKIP] ncsf_dir not found for format={format_}: {ncsf_dir}")
                continue
            
            ncsf_deriv_dir = "{}/ncsf_derivatives".format(ncsf_dir)
    
            tsv_dir = "{}/tsv".format(ncsf_dir)
            os.makedirs(tsv_dir, exist_ok=True)
                
            print(f'{ncsf_task_name} - {avg_method} - {format_}')
            tsv_fn = '{}/{}_task-{}_{}_{}_{}_{}_{}_ncsf_deriv.tsv'.format(
                tsv_dir, subject, ncsf_task_name,
                preproc_prep, filtering, normalization, avg_method, rois_method_format)
    
            # Load all data
            df_rois = pd.DataFrame()
            if format_ == 'fsnative':
                pycortex_subject = subject
                for hemi in ['hemi-L', 'hemi-R']:

                    # Derivatives
                    deriv_fn = '{}/{}_task-{}_{}_{}_{}_{}_{}_ncsf_deriv.func.gii'.format(
                        ncsf_deriv_dir, subject, ncsf_task_name, hemi,
                        preproc_prep, filtering, normalization, avg_method)
                    print(f'loading {deriv_fn}')
                    deriv_img, deriv_mat = load_surface(deriv_fn)
                    
                    # Stats
                    stats_fn = '{}/{}_task-{}_{}_{}_{}_{}_{}_ncsf_stats.func.gii'.format(
                        ncsf_deriv_dir, subject, ncsf_task_name, hemi,
                        preproc_prep, filtering, normalization, avg_method)
                    print(f'loading {stats_fn}')
                    stats_img, stats_mat = load_surface(stats_fn)
                                        
                    # Combine all derivatives
                    all_deriv_mat = np.concatenate((deriv_mat, stats_mat))
                    
                    # get MMP rois img
                    roi_dir = '{}/{}/derivatives/pp_data/{}/{}/rois'.format(main_dir, project_dir, subject, format_)
                    roi_fn = '{}/{}_{}_{}_{}_{}_rois-mmp.func.gii'.format(roi_dir, subject, hemi,
                                                                           preproc_prep, filtering, normalization)
                    roi_img, roi_mat = load_surface(roi_fn)
                
                    # get MMP rois numbers tsv
                    mmp_rois_numbers_tsv_fn = os.path.join(base_dir, "analysis_code", "atlas", "mmp_rois_numbers.tsv")
                    mmp_rois_numbers_df = pd.read_table(mmp_rois_numbers_tsv_fn, sep="\t")
                    
                    # Replace rois nums by names
                    roi_mat_names = np.vectorize(lambda x: dict(zip(mmp_rois_numbers_df['roi_num'], 
                                                                    mmp_rois_numbers_df['roi_name'])).get(x, x))(roi_mat)
        
                    roi_verts = get_rois(subject=pycortex_subject, 
                                          surf_format=format_, 
                                          rois_type=rois_method_format, 
                                          mask=True,
                                          rois=rois, 
                                          hemis=hemi)
                    
                    # Create and combine pandas df for each roi and brain hemisphere
                    print('Creating dataframe...')
                    for roi in roi_verts.keys():
                        data_dict = {col: all_deriv_mat[col_idx, roi_verts[roi]] for col_idx, col in enumerate(maps_names)}
                        data_dict['roi'] = np.array([roi] * all_deriv_mat[:, roi_verts[roi]].shape[1])
                        data_dict['roi_mmp'] = roi_mat_names[0, roi_verts[roi]]                                
                        data_dict['subject'] = np.array([subject] * all_deriv_mat[:, roi_verts[roi]].shape[1])
                        data_dict['hemi'] = np.array([hemi] * all_deriv_mat[:, roi_verts[roi]].shape[1])
                        data_dict['num_vert'] = np.where(roi_verts[roi])[0]
                        df_rois = pd.concat([df_rois, pd.DataFrame(data_dict)], ignore_index=True)
                
            elif format_ == '170k':
                pycortex_subject = pycortex_subject_template

                # Derivatives
                deriv_fn = '{}/{}_task-{}_{}_{}_{}_{}_ncsf_deriv.dtseries.nii'.format(
                    ncsf_deriv_dir, subject, ncsf_task_name,
                    preproc_prep, filtering, normalization, avg_method)
                print(f'loading {deriv_fn}')
                deriv_img, deriv_mat = load_surface(deriv_fn)
                
                # Stats
                stats_fn = '{}/{}_task-{}_{}_{}_{}_{}_ncsf_stats.dtseries.nii'.format(
                    ncsf_deriv_dir, subject, ncsf_task_name,
                    preproc_prep, filtering, normalization, avg_method)
                print(f'loading {stats_fn}')
                stats_img, stats_mat = load_surface(stats_fn)
                

                # Combine all derivatives
                all_deriv_mat = np.concatenate((deriv_mat, stats_mat))
                 
                # get MMP rois img
                roi_dir = '{}/{}/derivatives/pp_data/{}/{}/rois'.format(main_dir, project_dir, subject, format_)
                roi_fn = '{}/{}_{}_{}_{}_rois-mmp.dtseries.nii'.format(roi_dir, subject, 
                                                                       preproc_prep, filtering, normalization)
                roi_img, roi_mat = load_surface(roi_fn)
            
                # get MMP rois numbers tsv
                mmp_rois_numbers_tsv_fn = os.path.join(base_dir, "analysis_code", "atlas", "mmp_rois_numbers.tsv")
                mmp_rois_numbers_df = pd.read_table(mmp_rois_numbers_tsv_fn, sep="\t")
                
                # Replace rois nums by names
                roi_mat_names = np.vectorize(lambda x: dict(zip(mmp_rois_numbers_df['roi_num'], 
                                                                mmp_rois_numbers_df['roi_name'])).get(x, x))(roi_mat)
        
                # Create and combine pandas df for each roi and brain hemisphere
                print('Creating dataframe...')
                for hemi in ['hemi-L', 'hemi-R']:
                    # Get roi masks
                    roi_verts = get_rois(subject=pycortex_subject, 
                                          surf_format=format_, 
                                          rois_type=rois_method_format, 
                                          mask=True,
                                          rois=rois, 
                                          hemis=hemi)
        
                    for roi in roi_verts.keys():
                        data_dict = {col: all_deriv_mat[col_idx, roi_verts[roi]] for col_idx, col in enumerate(maps_names)}
                        data_dict['roi'] = np.array([roi] * all_deriv_mat[:, roi_verts[roi]].shape[1])
                        data_dict['roi_mmp'] = roi_mat_names[0, roi_verts[roi]]                                
                        data_dict['subject'] = np.array([subject] * all_deriv_mat[:, roi_verts[roi]].shape[1])
                        data_dict['hemi'] = np.array([hemi] * all_deriv_mat[:, roi_verts[roi]].shape[1])
                        data_dict['num_vert'] = np.where(roi_verts[roi])[0]
                        df_rois = pd.concat([df_rois, pd.DataFrame(data_dict)], ignore_index=True)
    
            df_rois[['hrf_1', 'hrf_2']] = df_rois[['hrf_1', 'hrf_2']].fillna('non_computed')
            print('Saving tsv: {}'.format(tsv_fn))
            df_rois.to_csv(tsv_fn, sep="\t", na_rep='NaN', index=False)

# Define permission cmd
print('Changing files permissions in {}/{}'.format(main_dir, project_dir))
os.system("chmod -Rf 771 {}/{}".format(main_dir, project_dir))
os.system("chgrp -Rf {} {}/{}".format(group, main_dir, project_dir))