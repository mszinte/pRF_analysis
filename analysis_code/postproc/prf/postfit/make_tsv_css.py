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
sys.argv[5]: OPTIONAL main analysis folder (e.g. prf_em_ctrl)
-----------------------------------------------------------------------------------------
Output(s):
TSV file
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/[PROJECT]/analysis_code/postproc/prf/postfit/
2. run python command
>> python make_tsv_css.py [main directory] [project name] 
                          [subject num] [group] [analysis folder - optional]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/projects/pRF_analysis/analysis_code/postproc/prf/postfit/

python make_tsv_css.py /scratch/mszinte/data MotConf sub-01 327
python make_tsv_css.py /scratch/mszinte/data MotConf sub-170k 327

python make_tsv_css.py /scratch/mszinte/data RetinoMaps sub-01 327
python make_tsv_css.py /scratch/mszinte/data RetinoMaps sub-170k 327

python make_tsv_css.py /scratch/mszinte/data amblyo_prf sub-01 327
python make_tsv_css.py /scratch/mszinte/data amblyo_prf sub-170k 327
python make_tsv_css.py /scratch/mszinte/data centbids sub-21002475231 327 
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
if len(sys.argv) > 5: output_folder = sys.argv[5]
else: output_folder = "prf"

# Define analysis parameters
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../../"))
settings_path = os.path.join(base_dir, project_dir, "settings.json")

with open(settings_path) as f:
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
maps_names = maps_names_css + maps_names_css_stats + maps_names_pcm

# Set pycortex db and colormaps
cortex_dir = "{}/{}/derivatives/pp_data/cortex".format(main_dir, project_dir)
set_pycortex_config_file(cortex_dir)

# Loop over format
for format_ in formats:
    # Define directories and fn
    prf_dir = "{}/{}/derivatives/pp_data/{}/{}/{}".format(
        main_dir, project_dir, subject, format_, output_folder)
    prf_deriv_dir = "{}/prf_derivatives".format(prf_dir)
    tsv_dir = "{}/tsv".format(prf_dir)
    os.makedirs(tsv_dir, exist_ok=True)
    tsv_fn = '{}/{}_css-all_derivatives.tsv'.format(tsv_dir, subject)

    # Load all data
    df_rois = pd.DataFrame()
    if format_ == 'fsnative':
        pycortex_subject = subject
        for hemi in ['hemi-L', 'hemi-R']:
            
            # Derivatives
            deriv_median_fn = '{}/{}_task-{}_{}_fmriprep_dct_avg_prf-deriv_css_loo-median.func.gii'.format(
                prf_deriv_dir, subject, prf_task_name, hemi)
            deriv_img, deriv_mat = load_surface(deriv_median_fn)
            
            # CM
            pcm_median_fn = '{}/{}_task-{}_{}_fmriprep_dct_avg_prf-pcm_css_loo-median.func.gii'.format(
                prf_deriv_dir, subject, prf_task_name, hemi)
            pcm_img, pcm_mat = load_surface(pcm_median_fn)

            # Stats
            stats_median_fn = '{}/{}_task-{}_{}_fmriprep_dct_avg_prf-stats_loo-median.func.gii'.format(
                prf_deriv_dir, subject, prf_task_name, hemi)
            stats_img, stats_mat = load_surface(stats_median_fn)
            
            # Combine all derivatives
            all_deriv_mat = np.concatenate((deriv_mat, stats_mat, pcm_mat))

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
        pycortex_subject = 'sub-170k'
        # Derivatives
        deriv_median_fn = '{}/{}_task-{}_fmriprep_dct_avg_prf-deriv_css_loo-median.dtseries.nii'.format(
            prf_deriv_dir, subject, prf_task_name)
        deriv_img, deriv_mat = load_surface(deriv_median_fn)
        
        # CM
        pcm_median_fn = '{}/{}_task-{}_fmriprep_dct_avg_prf-pcm_css_loo-median.dtseries.nii'.format(
            prf_deriv_dir, subject, prf_task_name)
        pcm_img, pcm_mat = load_surface(pcm_median_fn)

        # Stats
        stats_median_fn = '{}/{}_task-{}_fmriprep_dct_avg_prf-stats_loo-median.dtseries.nii'.format(
            prf_deriv_dir, subject, prf_task_name)
        stats_img, stats_mat = load_surface(stats_median_fn)
        
        # Combine all derivatives
        all_deriv_mat = np.concatenate((deriv_mat, stats_mat, pcm_mat))
        
        # Get roi mask
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

    print('Saving tsv: {}'.format(tsv_fn))
    df_rois.to_csv(tsv_fn, sep="\t", na_rep='NaN', index=False)

## Define permission cmd
#print('Changing files permissions in {}/{}'.format(main_dir, project_dir))
#os.system("chmod -Rf 771 {}/{}".format(main_dir, project_dir))
#os.system("chgrp -Rf {} {}/{}".format(group, main_dir, project_dir))
