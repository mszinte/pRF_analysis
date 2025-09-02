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
>> cd ~/projects/pRF_analysis/RetinoMaps/intertask/
2. run python command
python make_rois_fig.py [main directory] [project name] [subject] [group]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/projects/pRF_analysis/RetinoMaps/intertask/
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
import json
import numpy as np
import pandas as pd

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
group = sys.argv[4]

# Load settings
with open('../settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)
if subject == 'sub-170k': 
    formats = ['170k']
    extensions = ['dtseries.nii']
else: 
    formats = analysis_info['formats']
    extensions = analysis_info['extensions']
TRs = analysis_info['TRs']
rois = analysis_info['rois']
rois_groups = analysis_info['rois_groups']
group_tasks = analysis_info['task_intertask']
subjects_to_group = analysis_info['subjects']

# Threshold settings
ecc_threshold = analysis_info['ecc_th']
size_threshold = analysis_info['size_th']
rsqr_threshold = analysis_info['rsqr_th']
amplitude_threshold = analysis_info['amplitude_th']
stats_threshold = analysis_info['stats_th']
n_threshold = analysis_info['n_th']

for tasks in group_tasks : 
    if 'SacVELoc' in tasks: suffix = 'SacVE_PurVE'
    else : suffix = 'Sac_Pur'
    for format_ in formats :
        # Individual subject analysis
        if 'group' not in subject:
            # Directories
            intertask_tsv_dir ='{}/{}/derivatives/pp_data/{}/{}/intertask/tsv'.format(
            main_dir, project_dir, subject, format_)
            os.makedirs(intertask_tsv_dir, exist_ok=True)
            
            # Load subject TSV
            data = pd.read_table('{}/{}_intertask-all_derivatives_{}.tsv'.format(intertask_tsv_dir, subject, suffix))
            
            data.loc[(data.amplitude < amplitude_threshold) |
                      (data.prf_ecc < ecc_threshold[0]) | (data.prf_ecc > ecc_threshold[1]) |
                      (data.prf_size < size_threshold[0]) | (data.prf_size > size_threshold[1]) | 
                      (data.prf_n < n_threshold[0]) | (data.prf_n > n_threshold[1]) | 
                      (data.prf_loo_r2 < rsqr_threshold)] = np.nan
            data = data.dropna()

            # Active vertex ROI 
            subject_rois_area_categorie_df = pd.DataFrame()
            for roi in rois : 
                # Compute categorie proportions 
                n_vert_roi = data.loc[data['roi'] == roi].shape[0]
                n_vert_roi_saccade = data.loc[(data['roi'] == roi) & (data['saccade'] == 'saccade')].shape[0]
                n_vert_roi_pursuit = data.loc[(data['roi'] == roi) & (data['pursuit'] == 'pursuit')].shape[0]
                n_vert_roi_vision = data.loc[(data['roi'] == roi) & (data['vision'] == 'vision')].shape[0]
                n_vert_roi_vision_and_pursuit_and_saccade = data.loc[(data['roi'] == roi) & 
                                                                     (data['saccade'] == 'saccade') & 
                                                                     (data['pursuit'] == 'pursuit') & 
                                                                     (data['vision'] == 'vision'), 'vert_area'].shape[0]
                # Compute percentage active vertex in rois
                if n_vert_roi == 0 :
                    percent_saccade = percent_pursuit = percent_vision = percent_vision_and_pursuit_and_saccade = 0
                else : 
                    percent_saccade = ((n_vert_roi_saccade * 100)/n_vert_roi)
                    percent_pursuit = ((n_vert_roi_pursuit * 100)/n_vert_roi)
                    percent_vision = ((n_vert_roi_vision * 100)/n_vert_roi)
                    percent_vision_and_pursuit_and_saccade = ((n_vert_roi_vision_and_pursuit_and_saccade * 100)/n_vert_roi)

        
                active_vertex_roi_categorie_df = pd.DataFrame({'subject':[subject], 
                                                               'roi': [roi], 
                                                               'saccade': [percent_saccade], 
                                                               'pursuit': [percent_pursuit], 
                                                               'vision': [percent_vision], 
                                                               'all': [percent_vision_and_pursuit_and_saccade]})
        
        
                subject_rois_area_categorie_df = pd.concat([subject_rois_area_categorie_df, active_vertex_roi_categorie_df], ignore_index=True)
              
              
            subject_rois_area_categorie_melt_df = subject_rois_area_categorie_df.melt(id_vars=['subject', 'roi'], 
                                                                                      value_vars=['saccade', 'pursuit', 'vision', 'all'], 
                                                                                      var_name='categorie', 
                                                                                      value_name='percentage_active')
            
            # Export subject DF 
            tsv_active_vertex_roi_fn = "{}/{}_active_vertex_roi_{}.tsv".format(intertask_tsv_dir, subject, suffix)
            print('Saving tsv: {}'.format(tsv_active_vertex_roi_fn))
            subject_rois_area_categorie_melt_df.to_csv(tsv_active_vertex_roi_fn, sep="\t", na_rep='NaN', index=False)
        
            if format_ == '170k':
                # Active vertex MMP ROI 
                subject_active_percent_rois_categorie_df = pd.DataFrame()
                for n_roi, roi in enumerate(rois):
                    for n_roi_mmp, roi_mmp in enumerate(rois_groups[n_roi]):
                        # Roi DF
                        roi_area = data.loc[(data['roi'] == roi) & (data['roi_mmp'] == roi_mmp)]
                
                        # Compute amount of vertex in MMP roi
                        n_vert_roi = data.loc[data['roi_mmp'] == roi_mmp].shape[0]
                        n_vert_roi_saccade = data.loc[(data['roi_mmp'] == roi_mmp) & (data['saccade'] == 'saccade')].shape[0]
                        n_vert_roi_pursuit = data.loc[(data['roi_mmp'] == roi_mmp) & (data['pursuit'] == 'pursuit')].shape[0]
                        n_vert_roi_vision = data.loc[(data['roi_mmp'] == roi_mmp) & (data['vision'] == 'vision')].shape[0]
                        n_vert_roi_vision_and_pursuit_and_saccade = data.loc[(data['roi_mmp'] == roi_mmp) & 
                                                                             (data['saccade'] == 'saccade') & 
                                                                             (data['pursuit'] == 'pursuit') & 
                                                                             (data['vision'] == 'vision')].shape[0]
                
                        # Compute percentage active vertex in mmp rois
                        if n_vert_roi == 0 :
                            percent_saccade = percent_pursuit = percent_vision = percent_vision_and_pursuit_and_saccade = 0
                        else : 
                            percent_saccade = ((n_vert_roi_saccade * 100) / n_vert_roi)
                            percent_pursuit = ((n_vert_roi_pursuit * 100) / n_vert_roi)
                            percent_vision = ((n_vert_roi_vision * 100) / n_vert_roi)
                            percent_vision_and_pursuit_and_saccade = ((n_vert_roi_vision_and_pursuit_and_saccade * 100) / n_vert_roi)
                    
                
                        active_percent_roi_categorie_df = pd.DataFrame({'subject':[subject], 
                                                                        'roi': [roi], 
                                                                        'roi_mmp': [roi_mmp], 
                                                                        'saccade': [percent_saccade], 
                                                                        'pursuit': [percent_pursuit], 
                                                                        'vision': [percent_vision], 
                                                                        'all': [percent_vision_and_pursuit_and_saccade]})
                
                        subject_active_percent_rois_categorie_df = pd.concat([subject_active_percent_rois_categorie_df, 
                                                                              active_percent_roi_categorie_df], 
                                                                             ignore_index=True)
                        
                subject_active_percent_rois_categorie_melt_df = subject_active_percent_rois_categorie_df.melt(id_vars=['subject', 'roi', 'roi_mmp'], 
                                                                                                              value_vars=['saccade', 'pursuit', 'vision', 'all'], 
                                                                                                              var_name='categorie', 
                                                                                                              value_name='percentage_active')
                
                # Export subject DF 
                tsv_active_vertex_roi_mmp_fn = "{}/{}_active_vertex_roi_mmp_{}.tsv".format(intertask_tsv_dir, subject, suffix)
                print('Saving tsv: {}'.format(tsv_active_vertex_roi_mmp_fn))
                subject_active_percent_rois_categorie_melt_df.to_csv(tsv_active_vertex_roi_mmp_fn, sep="\t", na_rep='NaN', index=False)
        
        # Group analysis
        else : 
            for i, subject_to_group in enumerate(subjects_to_group):
                
                # Load and concat subjects
                intertask_tsv_dir = '{}/{}/derivatives/pp_data/{}/{}/intertask/tsv'.format(
                    main_dir, project_dir, subject_to_group, format_)
                
                # Active Vertex roi
                # -----------------
                tsv_active_vertex_roi_fn = "{}/{}_active_vertex_roi_{}.tsv".format(intertask_tsv_dir, subject_to_group, suffix)
                df_active_vertex_roi_indiv = pd.read_table(tsv_active_vertex_roi_fn, sep="\t")
                if i == 0: df_active_vertex_roi = df_active_vertex_roi_indiv.copy()
                else: df_active_vertex_roi = pd.concat([df_active_vertex_roi, df_active_vertex_roi_indiv])
                
                if format_ == '170k':
                    # Active Vertex mmp roi
                    # ---------------------
                    tsv_active_vertex_roi_mmp_fn = "{}/{}_active_vertex_roi_mmp_{}.tsv".format(intertask_tsv_dir, subject_to_group, suffix)
                    df_active_vertex_roi_mmp_indiv = pd.read_table(tsv_active_vertex_roi_mmp_fn, sep="\t")
                    if i == 0: df_active_vertex_roi_mmp = df_active_vertex_roi_mmp_indiv.copy()
                    else: df_active_vertex_roi_mmp = pd.concat([df_active_vertex_roi_mmp, df_active_vertex_roi_mmp_indiv])

                
            
            # Median across subjects
            # Active Vertex roi
            # -----------------
            df_active_vertex_roi_median_df = df_active_vertex_roi.groupby(['roi','categorie'], sort=False).median(numeric_only=True).reset_index().rename(columns={'percentage_active': 'median'})
            df_active_vertex_roi_ci_low_df = df_active_vertex_roi.groupby(['roi', 'categorie'], sort=False).quantile(0.025, numeric_only=True).reset_index().rename(columns={'percentage_active': 'ci_low'})
            df_active_vertex_roi_ci_high_df = df_active_vertex_roi.groupby(['roi', 'categorie'], sort=False).quantile(0.975, numeric_only=True).reset_index().rename(columns={'percentage_active': 'ci_high'})
            
            group_active_vertex_roi_melt_df = pd.concat([df_active_vertex_roi_median_df, 
                                                         df_active_vertex_roi_ci_low_df.drop(columns=['roi', 'categorie']), 
                                                         df_active_vertex_roi_ci_high_df.drop(columns=['roi', 'categorie'])], 
                                                        axis=1).reset_index(drop=True)
            # Export DF 
            intertask_tsv_dir = '{}/{}/derivatives/pp_data/{}/{}/intertask/tsv'.format(
                main_dir, project_dir, subject, format_)
            tsv_active_vertex_roi_fn = "{}/{}_active_vertex_roi_{}.tsv".format(intertask_tsv_dir, subject, suffix)
            print('Saving tsv: {}'.format(tsv_active_vertex_roi_fn))
            group_active_vertex_roi_melt_df.to_csv(tsv_active_vertex_roi_fn, sep="\t", na_rep='NaN', index=False)
            
            if format_ == '170k':
                # Active Vertex MMP roi
                # ---------------------
                df_active_vertex_mmp_roi_median_df = df_active_vertex_roi_mmp.groupby(['roi','roi_mmp', 'categorie'], sort=False).median(numeric_only=True).reset_index().rename(columns={'percentage_active': 'median'})
                df_active_vertex_mmp_roi_ci_low_df = df_active_vertex_roi_mmp.groupby(['roi','roi_mmp', 'categorie'], sort=False).quantile(0.025, numeric_only=True).reset_index().rename(columns={'percentage_active': 'ci_low'})
                df_active_vertex_mmp_roi_ci_high_df = df_active_vertex_roi_mmp.groupby(['roi','roi_mmp', 'categorie'], sort=False).quantile(0.975, numeric_only=True).reset_index().rename(columns={'percentage_active': 'ci_high'})
                
                group_active_vertex_mmp_roi_melt_df = pd.concat([df_active_vertex_mmp_roi_median_df, 
                                                                 df_active_vertex_mmp_roi_ci_low_df.drop(columns=['roi', 'roi_mmp', 'categorie']), 
                                                                 df_active_vertex_mmp_roi_ci_high_df.drop(columns=['roi', 'roi_mmp', 'categorie'])], 
                                                                axis=1).reset_index(drop=True)
                    
                # Export DF
                tsv_active_vertex_roi_mmp_fn = "{}/{}_active_vertex_roi_mmp_{}.tsv".format(intertask_tsv_dir, subject, suffix)
                print('Saving tsv: {}'.format(tsv_active_vertex_roi_mmp_fn))
                group_active_vertex_mmp_roi_melt_df.to_csv(tsv_active_vertex_roi_mmp_fn, sep="\t", na_rep='NaN', index=False)
               
# # Define permission cmd
# print('Changing files permissions in {}/{}'.format(main_dir, project_dir))
# os.system("chmod -Rf 771 {}/{}".format(main_dir, project_dir))
# os.system("chgrp -Rf {} {}/{}".format(group, main_dir, project_dir))             
            
            
            
