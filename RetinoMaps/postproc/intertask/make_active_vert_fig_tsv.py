"""
-----------------------------------------------------------------------------------------
make_active_vert_fig_tsv.py
-----------------------------------------------------------------------------------------
Goal of the script:
Make tsv for active vertex figures
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject name (e.g. sub-01)
sys.argv[4]: server group (e.g. 327)
-----------------------------------------------------------------------------------------
Output(s):
TSV for active vertex figures
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/pRF_analysis/RetinoMaps/postproc/intertask/
2. run python command
python make_rois_fig.py [main directory] [project name] [subject] [group]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/projects/pRF_analysis/RetinoMaps/postproc/intertask/
python make_active_vert_fig_tsv.py /scratch/mszinte/data RetinoMaps sub-01 327
python make_active_vert_fig_tsv.py /scratch/mszinte/data RetinoMaps sub-170k 327
python make_active_vert_fig_tsv.py /scratch/mszinte/data RetinoMaps group 327
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

# Imports
import os
import sys
import numpy as np
import pandas as pd

# Personal import
sys.path.append("{}/../../../analysis_code/utils".format(os.getcwd()))
from settings_utils import load_settings

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
group = sys.argv[4]

# load settings
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../"))
settings_path = os.path.join(base_dir, project_dir, "settings.yml")
prf_settings_path = os.path.join(base_dir, project_dir, "prf-analysis.yml")
glm_settings_path = os.path.join(base_dir, project_dir, "glm-analysis.yml")
settings = load_settings([settings_path, prf_settings_path, glm_settings_path])
analysis_info = settings[0]  


formats = analysis_info['formats']
extensions = analysis_info['extensions']
rois_methods = analysis_info['rois_methods']
TRs = analysis_info['TRs']
group_tasks = analysis_info['task_intertask']
subjects_to_group = analysis_info['subjects']
preproc_prep = analysis_info['preproc_prep']
filtering = analysis_info['filtering']
normalization = analysis_info['normalization']

# Threshold settings
ecc_threshold = analysis_info['ecc_th']
size_threshold = analysis_info['size_th']
rsqr_threshold = analysis_info['rsqr_th']
amplitude_threshold = analysis_info['prf_amp_th']
stats_threshold = analysis_info['stats_th']
n_threshold = analysis_info['n_th']

rsq2use = 'prf_loo_rsq'
avg_method = 'loo-avg'

for tasks in group_tasks : 
    if 'SacVELoc' in tasks: intertask_group = 'SacVE-PurVE-pRF'
    else : intertask_group = 'Sac-Pur-pRF'
    # Format loop
    for format_, extension in zip(formats, extensions):
        print(format_)
        # define list of rois for each format
        rois_methods_format = rois_methods[format_]
        for rois_method_format in rois_methods_format:

            if rois_method_format == 'rois-drawn':
                rois = analysis_info[rois_method_format]
            elif rois_method_format == 'rois-group-mmp':
                rois = list(analysis_info[rois_method_format].keys())

        
            # Individual subject analysis
            if 'group' not in subject:
                # Directories
                fn_spec = "task-{}_{}_{}_{}_{}_{}".format(
                    intertask_group, preproc_prep, filtering,
                    normalization, avg_method, rois_method_format)
                
                tsv_dir = "{}/{}/derivatives/pp_data/{}/{}/intertask/tsv".format(main_dir, project_dir, subject, format_)
                tsv_fn = '{}/{}_{}_intertask-deriv.tsv'.format(
                    tsv_dir, subject, fn_spec)
                
                # Exception if no data for one format (e.g template subject)
                if not os.path.isdir(tsv_dir):
                    print(f"[SKIP] tsv_dir not found for format={format_}: {tsv_dir}")
                    continue
                
                data = pd.read_table(tsv_fn, sep="\t")
                raw_data = data.copy()
                
                # replace nan in pcm by non_computed 
                data.loc[data['all'] != 'vision', 
    ['n_neighbor', 'pcm_median', 'vert_geo_dist_median', 'vert_prf_dist_median']] = data.loc[data['all'] != 'vision', 
    ['n_neighbor', 'pcm_median', 'vert_geo_dist_median', 'vert_prf_dist_median']].fillna('non_computed')
    
                # Threshold data (replace by nan)
                if stats_threshold == 0.05: stats_col = 'corr_pvalue_5pt'
                elif stats_threshold == 0.01: stats_col = 'corr_pvalue_1pt'
                data.loc[(data.amplitude < amplitude_threshold[0]) |
                         (data.prf_ecc < ecc_threshold[0]) | (data.prf_ecc > ecc_threshold[1]) |
                         (data.prf_size < size_threshold[0]) | (data.prf_size > size_threshold[1]) | 
                         (data.prf_n < n_threshold[0]) | (data.prf_n > n_threshold[1]) | 
                         (data[rsq2use] < rsqr_threshold) |
                         (data[stats_col] > stats_threshold)] = np.nan
                data = data.dropna()
    
                for rois_type in ['roi', 'roi_mmp']:
                    if rois_type == 'roi':
                        rois = analysis_info[rois_method_format]
                    elif rois_type == 'roi_mmp':
                        mmp_rois_numbers_tsv_fn = os.path.join(base_dir, "analysis_code", "atlas", "mmp_rois_numbers.tsv")
                        mmp_rois_numbers_df = pd.read_table(mmp_rois_numbers_tsv_fn, sep="\t")
                        rois = mmp_rois_numbers_df['roi_name'].tolist()
                    
                    # Active vertex ROI 
                    subject_rois_area_categorie_df = pd.DataFrame()
                    for roi in rois : 
                        # Compute categorie proportions 
                        n_vert_roi = raw_data.loc[raw_data[rois_type] == roi].shape[0]
                        n_vert_roi_saccade = data.loc[(data[rois_type] == roi) & (data['saccade'] == 'saccade')].shape[0]
                        n_vert_roi_pursuit = data.loc[(data[rois_type] == roi) & (data['pursuit'] == 'pursuit')].shape[0]
                        n_vert_roi_vision = data.loc[(data[rois_type] == roi) & (data['vision'] == 'vision')].shape[0]
                        n_vert_roi_vision_and_pursuit_and_saccade = data.loc[(data[rois_type] == roi) & 
                                                                             (data['saccade'] == 'saccade') & 
                                                                             (data['pursuit'] == 'pursuit') & 
                                                                             (data['vision'] == 'vision')].shape[0]
                        # Compute percentage active vertex in rois
                        if n_vert_roi == 0 :
                            percent_saccade = percent_pursuit = percent_vision = percent_vision_and_pursuit_and_saccade = 0
                        else : 
                            percent_saccade = ((n_vert_roi_saccade * 100)/n_vert_roi)
                            percent_pursuit = ((n_vert_roi_pursuit * 100)/n_vert_roi)
                            percent_vision = ((n_vert_roi_vision * 100)/n_vert_roi)
                            percent_vision_and_pursuit_and_saccade = ((n_vert_roi_vision_and_pursuit_and_saccade * 100)/n_vert_roi)
        
                        active_vertex_roi_categorie_df = pd.DataFrame({'subject':[subject], 
                                                                       rois_type: [roi], 
                                                                       'saccade': [percent_saccade], 
                                                                       'pursuit': [percent_pursuit], 
                                                                       'vision': [percent_vision], 
                                                                       'all': [percent_vision_and_pursuit_and_saccade]})
                
                
                        subject_rois_area_categorie_df = pd.concat([subject_rois_area_categorie_df, active_vertex_roi_categorie_df], ignore_index=True)
                      
                      
                    subject_rois_area_categorie_melt_df = subject_rois_area_categorie_df.melt(id_vars=['subject', rois_type], 
                                                                                              value_vars=['saccade', 'pursuit', 'vision', 'all'], 
                                                                                              var_name='categorie', 
                                                                                              value_name='percentage_active')
        
                    # Export subject DF 
                    tsv_active_vertex_roi_fn = "{}/{}_{}_active-vertex-{}.tsv".format(
                        tsv_dir, subject, fn_spec, rois_type)
                    print('Saving tsv: {}'.format(tsv_active_vertex_roi_fn))
                    subject_rois_area_categorie_melt_df.to_csv(tsv_active_vertex_roi_fn, sep="\t", na_rep='NaN', index=False)
            
            # Group analysis
            else : 
                for rois_type in ['roi', 'roi_mmp']:
                    for i, subject_to_group in enumerate(subjects_to_group):
                        # Load and concat subjects
                        fn_spec = "task-{}_{}_{}_{}_{}_{}".format(
                            intertask_group, preproc_prep, filtering,
                            normalization, avg_method, rois_method_format)
                        
                        tsv_dir = "{}/{}/derivatives/pp_data/{}/{}/intertask/tsv".format(main_dir, project_dir, subject_to_group, format_)
                        
                        # Active Vertex roi
                        # -----------------
                        tsv_active_vertex_roi_fn = "{}/{}_{}_active-vertex-{}.tsv".format(
                            tsv_dir, subject_to_group, fn_spec, rois_type)
                        df_active_vertex_roi_indiv = pd.read_table(tsv_active_vertex_roi_fn, sep="\t")
                        if i == 0: df_active_vertex_roi = df_active_vertex_roi_indiv.copy()
                        else: df_active_vertex_roi = pd.concat([df_active_vertex_roi, df_active_vertex_roi_indiv])
                        
    
                    # Median across subjects
                    # Active Vertex roi
                    # -----------------
                    df_active_vertex_roi_median_df = df_active_vertex_roi.groupby([rois_type,'categorie'], sort=False).median(numeric_only=True).reset_index().rename(columns={'percentage_active': 'median'})
                    df_active_vertex_roi_ci_low_df = df_active_vertex_roi.groupby([rois_type, 'categorie'], sort=False).quantile(0.025, numeric_only=True).reset_index().rename(columns={'percentage_active': 'ci_low'})
                    df_active_vertex_roi_ci_high_df = df_active_vertex_roi.groupby([rois_type, 'categorie'], sort=False).quantile(0.975, numeric_only=True).reset_index().rename(columns={'percentage_active': 'ci_high'})
                    
                    group_active_vertex_roi_melt_df = pd.concat([df_active_vertex_roi_median_df, 
                                                                 df_active_vertex_roi_ci_low_df.drop(columns=['roi', 'categorie']), 
                                                                 df_active_vertex_roi_ci_high_df.drop(columns=['roi', 'categorie'])], 
                                                                axis=1).reset_index(drop=True)
                    # Export DF 
                    fn_spec = "task-{}_{}_{}_{}_{}_{}".format(
                        intertask_group, preproc_prep, filtering,
                        normalization, avg_method, rois_method_format)
                    
                    tsv_dir = "{}/{}/derivatives/pp_data/{}/{}/intertask/tsv".format(main_dir, project_dir, subject_to_group, format_)
                    tsv_active_vertex_roi_fn = "{}/{}_{}_active-vertex-{}.tsv".format(
                        tsv_dir, subject_to_group, fn_spec, rois_type)
                    print('Saving tsv: {}'.format(tsv_active_vertex_roi_fn))
                group_active_vertex_roi_melt_df.to_csv(tsv_active_vertex_roi_fn, sep="\t", na_rep='NaN', index=False)
                
                   
# # Define permission cmd
# print('Changing files permissions in {}/{}'.format(main_dir, project_dir))
# os.system("chmod -Rf 771 {}/{}".format(main_dir, project_dir))
# os.system("chgrp -Rf {} {}/{}".format(group, main_dir, project_dir))             
            
            
            
