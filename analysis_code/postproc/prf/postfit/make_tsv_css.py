"""
-----------------------------------------------------------------------------------------
make_tsv_css.py
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
>> cd ~/projects/pRF_analysis/analysis_code/postproc/prf/postfit/
2. run python command
>> python make_tsv_css.py [main directory] [project name] 
                          [subject num] [group] [analysis folder - optional]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/projects/pRF_analysis/analysis_code/postproc/prf/postfit/
python make_tsv_css.py /scratch/mszinte/data RetinoMaps sub-01 327
python make_tsv_css.py /scratch/mszinte/data RetinoMaps sub-170k 327
-----------------------------------------------------------------------------------------
Written by Martin Szinte (martin.szinte@gmail.com)
and Uriel Lascombes (uriel.lascombes@laposte.net)
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
import json
import numpy as np
import pandas as pd

# Personal import
sys.path.append("{}/../../../utils".format(os.getcwd()))
from pycortex_utils import get_rois, set_pycortex_config_file
from surface_utils import load_surface

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
group = sys.argv[4]

# Define analysis parameters
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../../"))
settings_path = os.path.join(base_dir, project_dir, "settings.json")

with open(settings_path) as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)
if subject == 'sub-170k': formats = ['170k']
else: formats = analysis_info['formats']
extensions = analysis_info['extensions']
rois = analysis_info['rois']
maps_names_pcm = analysis_info['maps_names_pcm']
maps_names_css_stats = analysis_info['maps_names_css_stats']
prf_task_names = analysis_info['prf_task_names']
preproc_prep = analysis_info['preproc_prep']
filtering = analysis_info['filtering']
normalization = analysis_info['normalization']
avg_methods = analysis_info['avg_methods']

# Set pycortex db and colormaps
cortex_dir = "{}/{}/derivatives/pp_data/cortex".format(main_dir, project_dir)
set_pycortex_config_file(cortex_dir)

# Loop over format

for avg_method in avg_methods:
    if 'loo' in avg_method:
        maps_names_css = analysis_info['maps_names_css_loo']
    else: 
        maps_names_css = analysis_info['maps_names_css']

    maps_names = maps_names_css + maps_names_css_stats + maps_names_pcm
    
    for format_, extension in zip(formats, extensions):

        prf_dir = "{}/{}/derivatives/pp_data/{}/{}/prf".format(
            main_dir, project_dir, subject, format_)
        prf_deriv_dir = "{}/prf_derivatives".format(prf_dir)

        tsv_dir = "{}/tsv".format(prf_dir)
        os.makedirs(tsv_dir, exist_ok=True)
        for prf_task_name in prf_task_names:
            
            print(f'{prf_task_name} - {avg_method} - {format_}')
            tsv_fn = '{}/{}_task-{}_{}_{}_{}_{}_prf-css-deriv.tsv'.format(
                tsv_dir, subject, prf_task_name,
                preproc_prep, filtering, normalization, avg_method)
    
            # Load all data
            df_rois = pd.DataFrame()
            if format_ == 'fsnative':
                pycortex_subject = subject
                
                for hemi in ['hemi-L', 'hemi-R']:

                    # Derivatives
                    deriv_fn = '{}/{}_task-{}_{}_{}_{}_{}_{}_prf-css_deriv.func.gii'.format(
                        prf_deriv_dir, subject, prf_task_name, hemi,
                        preproc_prep, filtering, normalization, avg_method)
                    print(f'loading {deriv_fn}')
                    deriv_img, deriv_mat = load_surface(deriv_fn)
                    
                    # Stats
                    stats_fn = '{}/{}_task-{}_{}_{}_{}_{}_{}_prf-css_stats.func.gii'.format(
                        prf_deriv_dir, subject, prf_task_name, hemi,
                        preproc_prep, filtering, normalization, avg_method)
                    print(f'loading {stats_fn}')
                    stats_img, stats_mat = load_surface(stats_fn)
                    

                    # pCM
                    pcm_fn = '{}/{}_task-{}_{}_{}_{}_{}_{}_prf-css_pcm.func.gii'.format(
                        prf_deriv_dir, subject, prf_task_name, hemi,
                        preproc_prep, filtering, normalization, avg_method)
                    print(f'loading {pcm_fn}')
                    pcm_img, pcm_mat = load_surface(pcm_fn)
                    

                    # Combine all derivatives
                    all_deriv_mat = np.concatenate((deriv_mat, stats_mat, pcm_mat))
        
                    # Get roi masks
                    roi_verts = get_rois(subject=pycortex_subject, 
                                          return_concat_hemis=False, 
                                          return_hemi=hemi, 
                                          rois=rois,
                                          mask=True, 
                                          atlas_name=None, 
                                          surf_size=None)
                
                    # Create and combine pandas df for each roi and brain hemisphere
                    print('Creating dataframe...')
                    for roi in roi_verts.keys():
                        data_dict = {col: all_deriv_mat[col_idx, roi_verts[roi]] for col_idx, col in enumerate(maps_names)}
                        data_dict['roi'] = np.array([roi] * all_deriv_mat[:, roi_verts[roi]].shape[1])
                        data_dict['subject'] = np.array([subject] * all_deriv_mat[:, roi_verts[roi]].shape[1])
                        data_dict['hemi'] = np.array([hemi] * all_deriv_mat[:, roi_verts[roi]].shape[1])
                        data_dict['num_vert'] = np.where(roi_verts[roi])[0]
                        df_rois = pd.concat([df_rois, pd.DataFrame(data_dict)], ignore_index=True)
                
            elif format_ == '170k':
                pycortex_subject = 'sub-170k'

                # Derivatives
                deriv_fn = '{}/{}_task-{}_{}_{}_{}_{}_prf-css_deriv.dtseries.nii'.format(
                    prf_deriv_dir, subject, prf_task_name,
                    preproc_prep, filtering, normalization, avg_method)
                print(f'loading {deriv_fn}')
                deriv_img, deriv_mat = load_surface(deriv_fn)
                
                # Stats
                stats_fn = '{}/{}_task-{}_{}_{}_{}_{}_prf-css_stats.dtseries.nii'.format(
                    prf_deriv_dir, subject, prf_task_name,
                    preproc_prep, filtering, normalization, avg_method)
                print(f'loading {stats_fn}')
                stats_img, stats_mat = load_surface(stats_fn)
                
                # pRF CM
                pcm_fn = '{}/{}_task-{}_{}_{}_{}_{}_prf-css_pcm.dtseries.nii'.format(
                    prf_deriv_dir, subject, prf_task_name,
                    preproc_prep, filtering, normalization, avg_method)
                print(f'loading {pcm_fn}')
                pcm_img, pcm_mat = load_surface(pcm_fn)
                
                # Combine all derivatives
                all_deriv_mat = np.concatenate((deriv_mat, stats_mat, pcm_mat))
            
                # Get roi masks
                roi_verts_L, roi_verts_R = get_rois(subject=subject,
                                                    return_concat_hemis=False,
                                                    return_hemi=None,
                                                    rois=rois,
                                                    mask=True,
                                                    atlas_name='mmp_group',
                                                    surf_size='170k')
                
                # get MMP rois img
                roi_dir = '{}/{}/derivatives/pp_data/{}/{}/rois'.format(main_dir, project_dir, subject, format_)
                roi_fn = '{}/{}_rois_mmp.dtseries.nii'.format(roi_dir, subject)
                roi_img, roi_mat = load_surface(roi_fn)
            
                # get MMP rois numbers tsv
                mmp_rois_numbers_tsv_fn = '{}/db/sub-170k/mmp_rois_numbers.tsv'.format(cortex_dir)
                mmp_rois_numbers_df = pd.read_table(mmp_rois_numbers_tsv_fn, sep="\t")
                
                # Replace rois nums by names
                roi_mat_names = np.vectorize(lambda x: dict(zip(mmp_rois_numbers_df['roi_num'], 
                                                                mmp_rois_numbers_df['roi_name'])).get(x, x))(roi_mat)
        
                # Create and combine pandas df for each roi and brain hemisphere
                print('Creating dataframe...')
                for hemi in ['hemi-L', 'hemi-R']:
                    if hemi == 'hemi-L': roi_verts = roi_verts_L
                    elif hemi == 'hemi-R': roi_verts = roi_verts_R
        
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