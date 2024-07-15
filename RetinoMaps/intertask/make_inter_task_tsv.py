"""
-----------------------------------------------------------------------------------------
make_inter_task_tsv.py
-----------------------------------------------------------------------------------------
Goal of the script:
Create TSV file with all css analysis output and glm outputs
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
>> cd ~/projects/[PROJECT]/analysis_code/postproc/prf/postfit/
2. run python command
>> python make_tsv_css.py [main directory] [project name] [subject num] [group]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/projects/RetinoMaps/analysis_code/postproc/inter_task_analysis
python make_inter_task_tsv.py /scratch/mszinte/data RetinoMaps sub-01 327
python make_inter_task_tsv.py /scratch/mszinte/data RetinoMaps sub-170k 327
-----------------------------------------------------------------------------------------
Written by Martin Szinte (martin.szinte@gmail.com)
Edited by Uriel Lascombes (uriel.lascombes@laposte.net)
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
import pandas as pd
import numpy as np

# Personal import
sys.path.append("{}/../../utils".format(os.getcwd()))
from pycortex_utils import get_rois, set_pycortex_config_file
from surface_utils import load_surface

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
group = sys.argv[4]

# Define analysis parameters
with open('../../settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)
if subject == 'sub-170k': formats = ['170k']
else: formats = analysis_info['formats']
extensions = analysis_info['extensions']
prf_task_name = analysis_info['prf_task_name']
rois = analysis_info['rois']
maps_names_css = analysis_info['maps_names_css']
maps_names_pcm = analysis_info['maps_names_pcm']
maps_names_css_stats = analysis_info['maps_names_css_stats']
maps_names_vert_area = analysis_info["maps_names_vert_area"]
maps_names_inter_task = analysis_info["maps_names_inter_task"]
glm_code_names = analysis_info['glm_code_names']
maps_names = maps_names_css + maps_names_pcm + maps_names_css_stats + maps_names_vert_area + maps_names_inter_task

# Set pycortex db and colormaps
cortex_dir = "{}/{}/derivatives/pp_data/cortex".format(main_dir, project_dir)
set_pycortex_config_file(cortex_dir)

# Loop over format
for format_, pycortex_subject in zip(formats, [subject, 'sub-170k']):
    
    # Define directories and fn
    prf_dir = "{}/{}/derivatives/pp_data/{}/{}/prf".format(main_dir, project_dir, subject, format_)
    prf_deriv_dir = "{}/prf_derivatives".format(prf_dir)
    vert_area_dir = "{}/{}/derivatives/pp_data/{}/{}/vertex_area".format(main_dir, project_dir, subject, format_)
    inter_task_dir = "{}/{}/derivatives/pp_data/{}/{}/inter_task".format(main_dir, project_dir, subject, format_)
    
    tsv_dir = "{}/{}/derivatives/pp_data/{}/{}/inter_task/tsv/".format(main_dir, project_dir, subject, format_)
    os.makedirs(tsv_dir, exist_ok=True)
    tsv_fn = '{}/{}_inter_task-all_derivatives.tsv'.format(tsv_dir, subject)

    # Load all data
    df_rois = pd.DataFrame()
    if format_ == 'fsnative':
        for hemi in ['hemi-L', 'hemi-R']:
            
            # Derivatives
            deriv_avg_fn = '{}/{}_task-{}_{}_fmriprep_dct_prf-deriv-loo-avg_css.func.gii'.format(
                prf_deriv_dir, subject, prf_task_name, hemi)
            deriv_img, deriv_mat = load_surface(deriv_avg_fn)
            
            # pcm
            pcm_avg_fn = '{}/{}_task-{}_{}_fmriprep_dct_prf-pcm-loo-avg_css.func.gii'.format(
                prf_deriv_dir, subject, prf_task_name, hemi)
            pcm_img, pcm_mat = load_surface(pcm_avg_fn)

            # Stats
            stats_avg_fn = '{}/{}_task-{}_{}_fmriprep_dct_loo-avg_prf-stats.func.gii'.format(
                prf_deriv_dir, subject, prf_task_name, hemi)
            stats_img, stats_mat = load_surface(stats_avg_fn)
            
            # Vertex area
            vertex_area_fn = '{}/{}_{}_vertex_area.func.gii'.format(vert_area_dir, subject, hemi)
            vertex_area_img, vertex_area_mat = load_surface(vertex_area_fn)
            
            # Inter task
            inter_task_fn = '{}/{}_{}_inter_task.func.gii'.format(inter_task_dir, subject, hemi)
            inter_task_img, inter_task_mat = load_surface(inter_task_fn)

            # Combine all derivatives
            all_deriv_mat = np.concatenate((deriv_mat, pcm_mat, stats_mat, vertex_area_mat, inter_task_mat))

            # Get roi mask
            roi_verts = get_rois(subject=subject, 
                                 return_concat_hemis=False, 
                                 return_hemi=hemi, 
                                 rois=rois,
                                 mask=True, 
                                 atlas_name=None, 
                                 surf_size=None)

            # Create and combine pandas df for each roi and brain hemisphere
            for roi in roi_verts.keys():
                data_dict = {col: all_deriv_mat[col_idx, roi_verts[roi]] for col_idx, col in enumerate(maps_names)}
                data_dict['roi'] = np.array([roi] * all_deriv_mat[:, roi_verts[roi]].shape[1])
                data_dict['subject'] = np.array([subject] * all_deriv_mat[:, roi_verts[roi]].shape[1])
                data_dict['hemi'] = np.array([hemi] * all_deriv_mat[:, roi_verts[roi]].shape[1])
                data_dict['num_vert'] = np.where(roi_verts[roi])[0]
                df_rois = pd.concat([df_rois, pd.DataFrame(data_dict)], ignore_index=True)
                

                
    elif format_ == '170k':
        

        # Derivatives
        deriv_avg_fn = '{}/{}_task-{}_fmriprep_dct_prf-deriv-loo-avg_css.dtseries.nii'.format(
            prf_deriv_dir, subject, prf_task_name)
        deriv_img, deriv_mat = load_surface(deriv_avg_fn)
        
        # pcm
        pcm_avg_fn = '{}/{}_task-{}_fmriprep_dct_prf-pcm-loo-avg_css.dtseries.nii'.format(
            prf_deriv_dir, subject, prf_task_name)
        pcm_img, pcm_mat = load_surface(pcm_avg_fn)

        # Stats
        stats_avg_fn = '{}/{}_task-{}_fmriprep_dct_loo-avg_prf-stats.dtseries.nii'.format(
            prf_deriv_dir, subject, prf_task_name)
        stats_img, stats_mat = load_surface(stats_avg_fn)
        
        # Vertex area
        vertex_area_fn = '{}/{}_vertex_area.dtseries.nii'.format(vert_area_dir, subject)
        vertex_area_img, vertex_area_mat = load_surface(vertex_area_fn)

        # Inter task
        inter_task_fn = '{}/{}_inter_task.dtseries.nii'.format(inter_task_dir, subject)
        inter_task_img, inter_task_mat = load_surface(inter_task_fn)

        # Combine all derivatives
        all_deriv_mat = np.concatenate((deriv_mat, pcm_mat, stats_mat, vertex_area_mat, inter_task_mat))

        # Get roi mask
        roi_verts_L, roi_verts_R = get_rois(subject=subject,
                                            return_concat_hemis=False,
                                            return_hemi=None,
                                            rois=rois,
                                            mask=True,
                                            atlas_name='mmp_group',
                                            surf_size='170k')

        # Create and combine pandas df for each roi and brain hemisphere
        for hemi in ['hemi-L', 'hemi-R']:
            if hemi == 'hemi-L': roi_verts = roi_verts_L
            elif hemi == 'hemi-R': roi_verts = roi_verts_R

            for roi in roi_verts.keys():
                data_dict = {col: all_deriv_mat[col_idx, roi_verts[roi]] for col_idx, col in enumerate(maps_names)}
                data_dict['roi'] = np.array([roi] * all_deriv_mat[:, roi_verts[roi]].shape[1])
                data_dict['subject'] = np.array([subject] * all_deriv_mat[:, roi_verts[roi]].shape[1])
                data_dict['hemi'] = np.array([hemi] * all_deriv_mat[:, roi_verts[roi]].shape[1])
                data_dict['num_vert'] = np.where(roi_verts[roi])[0]
                df_rois = pd.concat([df_rois, pd.DataFrame(data_dict)], ignore_index=True)
      
    # Replace numbers by code value 
    code_to_name = {v: k for k, v in glm_code_names.items()}
    for column in maps_names_inter_task:
        if column in df_rois.columns:
            df_rois[column] = df_rois[column].map(code_to_name)
            
    print('Saving tsv: {}'.format(tsv_fn))
    df_rois.to_csv(tsv_fn, sep="\t", na_rep='NaN', index=False)

# Define permission cmd
print('Changing files permissions in {}/{}'.format(main_dir, project_dir))
os.system("chmod -Rf 771 {}/{}".format(main_dir, project_dir))
os.system("chgrp -Rf {} {}/{}".format(group, main_dir, project_dir))