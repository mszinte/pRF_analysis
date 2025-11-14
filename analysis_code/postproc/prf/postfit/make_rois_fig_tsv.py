"""
-----------------------------------------------------------------------------------------
make_rois_fig_tsv.py
-----------------------------------------------------------------------------------------
Goal of the script:
Make ROIs figure specific TSV
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject name (e.g. sub-01)
sys.argv[4]: server group (e.g. 327)
sys.argv[5]: OPTIONAL main analysis folder (e.g. prf_em_ctrl)
-----------------------------------------------------------------------------------------
Output(s):
CSS analysis tsv
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/[PROJECT]/analysis_code/postproc/prf/postfit/
2. run python command
python make_rois_fig.py [main directory] [project name] [subject] 
                        [group] [analysis folder - optional]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/projects/pRF_analysis/analysis_code/postproc/prf/postfit/

python make_rois_fig_tsv.py /scratch/mszinte/data MotConf sub-01 327
python make_rois_fig_tsv.py /scratch/mszinte/data MotConf sub-170k 327
python make_rois_fig_tsv.py /scratch/mszinte/data MotConf group 327

python make_rois_fig_tsv.py /scratch/mszinte/data RetinoMaps sub-01 327
python make_rois_fig_tsv.py /scratch/mszinte/data RetinoMaps sub-170k 327
python make_rois_fig_tsv.py /scratch/mszinte/data RetinoMaps group 327

python make_rois_fig_tsv.py /scratch/mszinte/data amblyo_prf sub-01 327
python make_rois_fig_tsv.py /scratch/mszinte/data amblyo_prf sub-170k 327
python make_rois_fig_tsv.py /scratch/mszinte/data amblyo_prf group 327
python make_rois_fig_tsv.py /scratch/mszinte/data amblyo_prf group_patient 327
python make_rois_fig_tsv.py /scratch/mszinte/data amblyo_prf group_aniso 327
python make_rois_fig_tsv.py /scratch/mszinte/data amblyo_prf group_mixed 327
python make_rois_fig_tsv.py /scratch/mszinte/data amblyo_prf group_strab 327
python make_rois_fig_tsv.py /scratch/mszinte/data amblyo_prf group_control 327
python make_rois_fig_tsv.py /scratch/mszinte/data amblyo_prf group_excluded 327
-----------------------------------------------------------------------------------------
Written by Uriel Lascombes (uriel.lascombes@laposte.net)
Edited by Martin Szinte (martin.szinte@gmail.com)
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
from maths_utils import make_prf_distribution_df, weighted_nan_median, weighted_nan_percentile

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
group = sys.argv[4]
if len(sys.argv) > 5: output_folder = sys.argv[5]
else: output_folder = "prf"

# Load settings
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../../"))
settings_path = os.path.join(base_dir, project_dir, "settings.json")

with open(settings_path) as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)
if subject == 'sub-170k': 
    formats = ['170k']
    extensions = ['dtseries.nii']
else: 
    formats = analysis_info['formats']
    extensions = analysis_info['extensions']
rois = analysis_info['rois']
ecc_threshold = analysis_info['ecc_th']
size_threshold = analysis_info['size_th']
rsqr_threshold = analysis_info['rsqr_th']
amplitude_threshold = analysis_info['amplitude_th']
stats_threshold = analysis_info['stats_th']
n_threshold = analysis_info['n_th']
subjects_to_group = analysis_info['subjects']
if subject == 'sub-170k': 
    formats = ['170k']
    extensions = ['dtseries.nii']
else: 
    formats = analysis_info['formats']
    extensions = analysis_info['extensions']
rois = analysis_info['rois']

fig_settings_path = os.path.join(base_dir, project_dir, "figure_settings.json")
with open(fig_settings_path) as f:
    json_s = f.read()
    figure_info = json.loads(json_s)
num_ecc_size_bins = figure_info['num_ecc_size_bins']
num_ecc_pcm_bins = figure_info['num_ecc_pcm_bins']
num_polar_angle_bins = figure_info['num_polar_angle_bins']
max_ecc = figure_info['max_ecc']
screen_side = figure_info['screen_side']
gaussian_mesh_grain = figure_info['gaussian_mesh_grain']
hot_zone_percent = figure_info['hot_zone_percent']

# Format loop
for format_, extension in zip(formats, extensions):
    print(format_)
    
    # Individual subject analysis
    if 'group' not in subject:
        print('Subject {} is processed'.format(subject))
        tsv_dir = '{}/{}/derivatives/pp_data/{}/{}/{}/tsv'.format(
            main_dir, project_dir, subject, format_, output_folder)
        os.makedirs(tsv_dir, exist_ok=True)
        
        tsv_fn = '{}/{}_css-all_derivatives.tsv'.format(tsv_dir, subject)
        data = pd.read_table(tsv_fn, sep="\t")
       
        # Keep a raw data df 
        data_raw = data.copy()
        
        # Threshold data (replace by nan)
        if stats_threshold == 0.05: stats_col = 'corr_pvalue_5pt'
        elif stats_threshold == 0.01: stats_col = 'corr_pvalue_1pt'
        data.loc[(data.amplitude < amplitude_threshold) |
                 (data.prf_ecc < ecc_threshold[0]) | (data.prf_ecc > ecc_threshold[1]) |
                 (data.prf_size < size_threshold[0]) | (data.prf_size > size_threshold[1]) | 
                 (data.prf_n < n_threshold[0]) | (data.prf_n > n_threshold[1]) | 
                 (data.prf_loo_r2 < rsqr_threshold) |
                 (data[stats_col] > stats_threshold)] = np.nan
        data = data.dropna()
        
        # ROI active proportion
        # ---------------------
        
        # Number of vert per roi
        n_vert_tot_roi = (data_raw.groupby('roi', sort=False)
                          .size()
                          .reset_index(name='n_vert_tot'))
        
        # Number of significant vert for 0.05 p vale correction
        n_vert_corr_5pt = (data_raw[data_raw['corr_pvalue_5pt'] < 0.05]
                           .groupby('roi', sort=False)
                           .size()
                           .reset_index(name='n_vert_corr_pvalue_5pt'))
        
        # Number of significant vert for 0.01 p vale correction
        n_vert_corr_1pt = (data_raw[data_raw['corr_pvalue_1pt'] < 0.01]
                           .groupby('roi', sort=False)
                           .size()
                           .reset_index(name='n_vert_corr_pvalue_1pt'))
        
        
        tmp_df = n_vert_tot_roi.merge(n_vert_corr_5pt, on='roi', how='left') \
                               .merge(n_vert_corr_1pt, on='roi', how='left')
        
        tmp_df = tmp_df.fillna(0)
        
        ratio_5pt = tmp_df[['roi']].copy()
        ratio_1pt = tmp_df[['roi']].copy()
        
        ratio_5pt['ratio_5pt'] = tmp_df['n_vert_corr_pvalue_5pt'] / tmp_df['n_vert_tot']
        ratio_1pt['ratio_1pt'] = tmp_df['n_vert_corr_pvalue_1pt'] / tmp_df['n_vert_tot']
        
        df_roi_active_vert = (n_vert_tot_roi
                              .merge(n_vert_corr_5pt, on='roi', how='left')
                              .merge(n_vert_corr_1pt, on='roi', how='left')
                              .merge(ratio_1pt, on='roi', how='left')
                              .merge(ratio_5pt, on='roi', how='left'))
        
        df_roi_active_vert = df_roi_active_vert.fillna(0)
           
        # Export tsv
        tsv_roi_active_vert_fn = "{}/{}_prf_active_vert.tsv".format(tsv_dir, subject)
        print('Saving tsv: {}'.format(tsv_roi_active_vert_fn))
        df_roi_active_vert.to_csv(tsv_roi_active_vert_fn, sep="\t", na_rep='NaN', index=False)
    
        # Violins
        # -------
        df_violins = data
        
        tsv_violins_fn = "{}/{}_prf_violins.tsv".format(tsv_dir, subject)
        print('Saving tsv: {}'.format(tsv_violins_fn))
        df_violins.to_csv(tsv_violins_fn, sep="\t", na_rep='NaN', index=False)
        
        # Parameters median
        # ------------------
        for num_roi, roi in enumerate(rois):
            df_roi = data.loc[(data.roi == roi)]

            df_params_median_roi = pd.DataFrame()
            df_params_median_roi['roi'] = [roi]
        
            df_params_median_roi['prf_loo_r2_weighted_median'] = weighted_nan_median(df_roi.prf_loo_r2, weights=df_roi.prf_loo_r2)
            df_params_median_roi['prf_loo_r2_ci_down'] = weighted_nan_percentile(df_roi.prf_loo_r2, df_roi.prf_loo_r2, 25)
            df_params_median_roi['prf_loo_r2_ci_up'] = weighted_nan_percentile(df_roi.prf_loo_r2, df_roi.prf_loo_r2, 75)
            
            df_params_median_roi['prf_size_weighted_median'] = weighted_nan_median(df_roi.prf_size, weights=df_roi.prf_loo_r2)
            df_params_median_roi['prf_size_ci_down'] = weighted_nan_percentile(df_roi.prf_size, df_roi.prf_loo_r2, 25)
            df_params_median_roi['prf_size_ci_up'] = weighted_nan_percentile(df_roi.prf_size, df_roi.prf_loo_r2, 75)
            
            df_params_median_roi['prf_ecc_weighted_median'] = weighted_nan_median(df_roi.prf_ecc, weights=df_roi.prf_loo_r2)
            df_params_median_roi['prf_ecc_ci_down'] = weighted_nan_percentile(df_roi.prf_ecc, df_roi.prf_loo_r2, 25)
            df_params_median_roi['prf_ecc_ci_up'] = weighted_nan_percentile(df_roi.prf_ecc, df_roi.prf_loo_r2, 75)
            
            df_params_median_roi['prf_n_weighted_median'] = weighted_nan_median(df_roi.prf_n, weights=df_roi.prf_loo_r2)
            df_params_median_roi['prf_n_ci_down'] = weighted_nan_percentile(df_roi.prf_n, df_roi.prf_loo_r2, 25)
            df_params_median_roi['prf_n_ci_up'] = weighted_nan_percentile(df_roi.prf_n, df_roi.prf_loo_r2, 75)
             
            df_params_median_roi['pcm_median_weighted_median'] = weighted_nan_median(df_roi.pcm_median, weights=df_roi.prf_loo_r2)
            df_params_median_roi['pcm_median_ci_down'] = weighted_nan_percentile(df_roi.pcm_median, df_roi.prf_loo_r2, 25)
            df_params_median_roi['pcm_median_ci_up'] = weighted_nan_percentile(df_roi.pcm_median, df_roi.prf_loo_r2, 75)
    
            if num_roi == 0: df_params_median = df_params_median_roi
            else: df_params_median = pd.concat([df_params_median, df_params_median_roi])

        tsv_params_median_fn = "{}/{}_prf_params_median.tsv".format(tsv_dir, subject)
        print('Saving tsv: {}'.format(tsv_params_median_fn))
        df_params_median.to_csv(tsv_params_median_fn, sep="\t", na_rep='NaN', index=False)
        
        # Ecc.size
        # --------
        # ecc_bins = np.concatenate(([0],np.linspace(0.26, 1, num_ecc_size_bins)**2 * max_ecc))
        # ecc_bins = np.concatenate(([0],np.linspace(0.3, 1, num_ecc_size_bins)**2 * max_ecc))
        ecc_bins = np.linspace(0, max_ecc, num_ecc_size_bins+1) 
        for num_roi, roi in enumerate(rois):
            df_roi = data.loc[(data.roi == roi)]
            df_bins = df_roi.groupby(pd.cut(df_roi['prf_ecc'], bins=ecc_bins))
            df_ecc_size_bin = pd.DataFrame()
            df_ecc_size_bin['roi'] = [roi]*num_ecc_size_bins
            df_ecc_size_bin['num_bins'] = np.arange(num_ecc_size_bins)
            df_ecc_size_bin['prf_ecc_bins'] = df_bins.apply(lambda x: weighted_nan_median(x['prf_ecc'].values, x['prf_loo_r2'].values)).values
            df_ecc_size_bin['prf_size_bins_median'] = df_bins.apply(lambda x: weighted_nan_median(x['prf_size'].values, x['prf_loo_r2'].values)).values
            df_ecc_size_bin['prf_loo_r2_bins_median'] = np.array(df_bins['prf_loo_r2'].median())
            df_ecc_size_bin['prf_size_bins_ci_upper_bound'] = df_bins.apply(lambda x: weighted_nan_percentile(x['prf_size'].values, x['prf_loo_r2'].values, 75)).values
            df_ecc_size_bin['prf_size_bins_ci_lower_bound'] = df_bins.apply(lambda x: weighted_nan_percentile(x['prf_size'].values, x['prf_loo_r2'].values, 25)).values
            if num_roi == 0: df_ecc_size_bins = df_ecc_size_bin
            else: df_ecc_size_bins = pd.concat([df_ecc_size_bins, df_ecc_size_bin]) 
        
        df_ecc_size = df_ecc_size_bins
        
        tsv_ecc_size_fn = "{}/{}_prf_ecc_size.tsv".format(tsv_dir, subject)
        print('Saving tsv: {}'.format(tsv_ecc_size_fn))
        df_ecc_size.to_csv(tsv_ecc_size_fn, sep="\t", na_rep='NaN', index=False)
        
        # Ecc.pCM
        # --------
        data_pcm = data
        # ecc_bins = np.concatenate(([0],np.linspace(0.26, 1, num_ecc_pcm_bins)**2 * max_ecc))
        # ecc_bins = np.concatenate(([0],np.linspace(0.3, 1, num_ecc_size_bins)**2 * max_ecc))
        ecc_bins = np.linspace(0, max_ecc, num_ecc_size_bins+1) 
        for num_roi, roi in enumerate(rois):
            df_roi = data_pcm.loc[(data.roi == roi)]
            df_bins = df_roi.groupby(pd.cut(df_roi['prf_ecc'], bins=ecc_bins))
            df_ecc_pcm_bin = pd.DataFrame()
            df_ecc_pcm_bin['roi'] = [roi]*num_ecc_pcm_bins
            df_ecc_pcm_bin['num_bins'] = np.arange(num_ecc_pcm_bins)
            df_ecc_pcm_bin['prf_ecc_bins'] = df_bins.apply(lambda x: weighted_nan_median(x['prf_ecc'].values, x['prf_loo_r2'].values)).values
            df_ecc_pcm_bin['prf_pcm_bins_median'] = df_bins.apply(lambda x: weighted_nan_median(x['pcm_median'].values, x['prf_loo_r2'].values)).values
            df_ecc_pcm_bin['prf_loo_r2_bins_median'] = np.array(df_bins['prf_loo_r2'].median())
            df_ecc_pcm_bin['prf_pcm_bins_ci_upper_bound'] = df_bins.apply(lambda x: weighted_nan_percentile(x['pcm_median'].values, x['prf_loo_r2'].values, 75)).values
            df_ecc_pcm_bin['prf_pcm_bins_ci_lower_bound'] = df_bins.apply(lambda x: weighted_nan_percentile(x['pcm_median'].values, x['prf_loo_r2'].values, 25)).values
            if num_roi == 0: df_ecc_pcm_bins = df_ecc_pcm_bin
            else: df_ecc_pcm_bins = pd.concat([df_ecc_pcm_bins, df_ecc_pcm_bin])  
    
        df_ecc_pcm = df_ecc_pcm_bins

        tsv_ecc_pcm_fn = "{}/{}_prf_ecc_pcm.tsv".format(tsv_dir, subject)
        print('Saving tsv: {}'.format(tsv_ecc_pcm_fn))
        df_ecc_pcm.to_csv(tsv_ecc_pcm_fn, sep="\t", na_rep='NaN', index=False)
        
        # Polar angle
        # -----------
        theta_slices = np.linspace(0, 360, num_polar_angle_bins, endpoint=False)
        data['prf_polar_angle'] = np.mod(np.degrees(np.angle(data.polar_real + 1j * data.polar_imag)), 360) 
        hemis = ['hemi-L', 'hemi-R', 'hemi-LR']
        for i, hemi in enumerate(hemis):
            hemi_values = ['hemi-L', 'hemi-R'] if hemi == 'hemi-LR' else [hemi]
            for j, roi in enumerate(rois): #
                df = data.loc[(data.roi==roi) & (data.hemi.isin(hemi_values))]
                if len(df): 
                    df_bins = df.groupby(pd.cut(df['prf_polar_angle'], bins=num_polar_angle_bins))
                    loo_r2_sum = df_bins['prf_loo_r2'].sum()
                else: loo_r2_sum = [np.nan]*num_polar_angle_bins
    
                df_polar_angle_bin = pd.DataFrame()
                df_polar_angle_bin['roi'] = [roi]*(num_polar_angle_bins)
                df_polar_angle_bin['hemi'] = [hemi]*(num_polar_angle_bins)
                df_polar_angle_bin['num_bins'] = np.arange((num_polar_angle_bins))
                df_polar_angle_bin['theta_slices'] = np.array(theta_slices)
                df_polar_angle_bin['loo_r2_sum'] = np.array(loo_r2_sum)
                
                if j == 0 and i == 0: df_polar_angle_bins = df_polar_angle_bin
                else: df_polar_angle_bins = pd.concat([df_polar_angle_bins, df_polar_angle_bin])  
                    
        df_polar_angle = df_polar_angle_bins

        tsv_polar_angle_fn = "{}/{}_prf_polar_angle.tsv".format(tsv_dir, subject)
        print('Saving tsv: {}'.format(tsv_polar_angle_fn))
        df_polar_angle.to_csv(tsv_polar_angle_fn, sep="\t", na_rep='NaN', index=False)
        
        # Contralaterality
        # ----------------         
        for j, roi in enumerate(rois):
            df_rh = data.loc[(data.roi == roi) & (data.hemi == 'hemi-R')]
            df_lh = data.loc[(data.roi == roi) & (data.hemi == 'hemi-L')]
            try: contralaterality_prct = (sum(df_rh.loc[df_rh.prf_x < 0].prf_loo_r2) + \
                                         sum(df_lh.loc[df_lh.prf_x > 0].prf_loo_r2)) / \
                                        (sum(df_rh.prf_loo_r2) + sum(df_lh.prf_loo_r2))
            except: contralaterality_prct = np.nan
            
            df_contralaterality_roi = pd.DataFrame()
            df_contralaterality_roi['roi'] = [roi]
            df_contralaterality_roi['contralaterality_prct'] = np.array(contralaterality_prct)
    
            if j == 0: df_contralaterality = df_contralaterality_roi
            else: df_contralaterality = pd.concat([df_contralaterality, df_contralaterality_roi]) 

        tsv_contralaterality_fn = "{}/{}_prf_contralaterality.tsv".format(tsv_dir, subject)
        print('Saving tsv: {}'.format(tsv_contralaterality_fn))
        df_contralaterality.to_csv(tsv_contralaterality_fn, sep="\t", na_rep='NaN', index=False)
            
        # Spatial distribution 
        # --------------------  
        hemis = ['hemi-L', 'hemi-R', 'hemi-LR']
        for i, hemi in enumerate(hemis):
            hemi_values = ['hemi-L', 'hemi-R'] if hemi == 'hemi-LR' else [hemi]
            data_hemi = data.loc[data.hemi.isin(hemi_values)]
            df_distribution_hemi = make_prf_distribution_df(
                data_hemi, rois, screen_side, gaussian_mesh_grain)
    
            df_distribution_hemi['hemi'] = [hemi] * len(df_distribution_hemi)
            if i == 0: df_distribution = df_distribution_hemi
            else: df_distribution = pd.concat([df_distribution, df_distribution_hemi])

        tsv_distribution_fn = "{}/{}_prf_distribution.tsv".format(tsv_dir, subject)
        print('Saving tsv: {}'.format(tsv_distribution_fn))
        df_distribution.to_csv(tsv_distribution_fn, sep="\t", na_rep='NaN', index=False)
        
        # # Spatial distribution hot zone barycentre
        # # ----------------------------------------
        # hemis = ['hemi-L', 'hemi-R', 'hemi-LR']
        # for i, hemi in enumerate(hemis):
        #     hemi_values = ['hemi-L', 'hemi-R'] if hemi == 'hemi-LR' else [hemi]
        #     df_distribution_hemi = df_distribution.loc[df_distribution.hemi.isin(hemi_values)]
        #     df_barycentre_hemi = make_prf_barycentre_df(
        #         df_distribution_hemi, rois, screen_side, gaussian_mesh_grain, hot_zone_percent=hot_zone_percent)
            
        #     df_barycentre_hemi['hemi'] = [hemi] * len(df_barycentre_hemi)
        #     if i == 0: df_barycentre = df_barycentre_hemi
        #     else: df_barycentre = pd.concat([df_barycentre, df_barycentre_hemi])
        
        # tsv_barycentre_fn = "{}/{}_prf_barycentre.tsv".format(tsv_dir, subject)
        # print('Saving tsv: {}'.format(tsv_barycentre_fn))
        # df_barycentre.to_csv(tsv_barycentre_fn, sep="\t", na_rep='NaN', index=False)
        
    # Group Analysis    
    else :
        print('group')
        for i, subject_to_group in enumerate(subjects_to_group):
            tsv_dir = '{}/{}/derivatives/pp_data/{}/{}/{}/tsv'.format(
                main_dir, project_dir, subject_to_group, format_, output_folder)
    
            # ROI surface areas 
            # -----------------
            tsv_roi_area_fn = "{}/{}_prf_active_vert.tsv".format(tsv_dir, subject_to_group)
            df_roi_area_indiv = pd.read_table(tsv_roi_area_fn, sep="\t")
            if i == 0: df_roi_area = df_roi_area_indiv.copy()
            else: df_roi_area = pd.concat([df_roi_area, df_roi_area_indiv])
    
            # Violins
            # -------
            tsv_violins_fn = "{}/{}_prf_violins.tsv".format(tsv_dir, subject_to_group)
            df_violins_indiv = pd.read_table(tsv_violins_fn, sep="\t")
            if i == 0: df_violins = df_violins_indiv.copy()
            else: df_violins = pd.concat([df_violins, df_violins_indiv])
    
            # Parameters median
            # ------------------
            # use df_violins
    
            # Ecc.size
            # --------
            tsv_ecc_size_fn = "{}/{}_prf_ecc_size.tsv".format(tsv_dir, subject_to_group)
            df_ecc_size_indiv = pd.read_table(tsv_ecc_size_fn, sep="\t")
            if i == 0: df_ecc_size = df_ecc_size_indiv.copy()
            else: df_ecc_size = pd.concat([df_ecc_size, df_ecc_size_indiv])
    
            # Ecc.pCM
            # -------
            tsv_ecc_pcm_fn = "{}/{}_prf_ecc_pcm.tsv".format(tsv_dir, subject_to_group)
            df_ecc_pcm_indiv = pd.read_table(tsv_ecc_pcm_fn, sep="\t")
            if i == 0: df_ecc_pcm = df_ecc_pcm_indiv.copy()
            else: df_ecc_pcm = pd.concat([df_ecc_pcm, df_ecc_pcm_indiv])
    
            # Polar angle
            # -----------
            tsv_polar_angle_fn = "{}/{}_prf_polar_angle.tsv".format(tsv_dir, subject_to_group)
            df_polar_angle_indiv = pd.read_table(tsv_polar_angle_fn, sep="\t")
            if i == 0: df_polar_angle = df_polar_angle_indiv.copy()
            else: df_polar_angle = pd.concat([df_polar_angle, df_polar_angle_indiv])
    
            # Contralaterality
            # ----------------
            tsv_contralaterality_fn = "{}/{}_prf_contralaterality.tsv".format(tsv_dir, subject_to_group)
            df_contralaterality_indiv = pd.read_table(tsv_contralaterality_fn, sep="\t")
            if i == 0: df_contralaterality = df_contralaterality_indiv.copy()
            else: df_contralaterality = pd.concat([df_contralaterality, df_contralaterality_indiv])
            
            # Spatial distribution 
            # -------------------
            tsv_distribution_fn = "{}/{}_prf_distribution.tsv".format(tsv_dir, subject_to_group)
            df_distribution_indiv = pd.read_table(tsv_distribution_fn, sep="\t")
            mesh_indiv = df_distribution_indiv.drop(columns=['roi', 'x', 'y', 'hemi']).values
            others_columns = df_distribution_indiv[['roi', 'x', 'y', 'hemi']]

            if i == 0: mesh_group = np.expand_dims(mesh_indiv, axis=0)
            else: mesh_group = np.vstack((mesh_group, np.expand_dims(mesh_indiv, axis=0)))
           
        # Median and saving tsv
        tsv_dir = '{}/{}/derivatives/pp_data/{}/{}/{}/tsv'.format(
            main_dir, project_dir, subject, format_, output_folder)
        os.makedirs(tsv_dir, exist_ok=True)
        
        # ROI surface areas 
        # -----------------
        df_roi_area = df_roi_area.groupby(['roi'], sort=False).median().reset_index()
        tsv_roi_area_fn = "{}/{}_prf_active_vert.tsv".format(tsv_dir, subject)
        print('Saving tsv: {}'.format(tsv_roi_area_fn))
        df_roi_area.to_csv(tsv_roi_area_fn, sep="\t", na_rep='NaN', index=False)
        
        # Violins
        # -------
        df_violins = df_violins # no averaging
        tsv_violins_fn = "{}/{}_prf_violins.tsv".format(tsv_dir, subject)
        print('Saving tsv: {}'.format(tsv_violins_fn))
        df_violins.to_csv(tsv_violins_fn, sep="\t", na_rep='NaN', index=False)
    
        # Parameters median
        # ------------------
        df_params_median = df_violins
        
        # compute median 
        colnames = ['prf_loo_r2', 'prf_size', 'prf_ecc', 'prf_n', 'pcm_median']
        df_params_median_indiv = df_params_median.groupby(['roi', 'subject'])[['prf_loo_r2']].apply(
            lambda x: weighted_nan_median(x['prf_loo_r2'], df_params_median.loc[x.index, 'prf_loo_r2'])).reset_index(name='prf_loo_r2_weighted_median')
        
        for colname in colnames[1:]:
            df_params_median_indiv['{}_weighted_median'.format(colname)] = df_params_median.groupby(['roi', 'subject'])[[colname, 'prf_loo_r2']].apply(
                lambda x: weighted_nan_median(x[colname], df_params_median.loc[x.index, 'prf_loo_r2'])).reset_index()[0]
        df_params_med_median = df_params_median_indiv.groupby(['roi'])[[colname + '_weighted_median' for colname in colnames]].median()

        # compute Ci
        df_params_median_ci = pd.DataFrame()
        for colname in colnames:
            df_params_median_ci['{}_ci_down'.format(colname)] = df_params_median_indiv.groupby(['roi']).apply(
                lambda x: weighted_nan_percentile(x['{}_weighted_median'.format(colname)], x['prf_loo_r2_weighted_median'], 2.5)) 
            df_params_median_ci['{}_ci_up'.format(colname)] = df_params_median_indiv.groupby(['roi']).apply(
                lambda x: weighted_nan_percentile(x['{}_weighted_median'.format(colname)], x['prf_loo_r2_weighted_median'], 97.5)) 

        df_params_median = pd.concat([df_params_med_median, df_params_median_ci], axis=1).reset_index()
        tsv_params_median_fn = "{}/{}_prf_params_median.tsv".format(tsv_dir, subject)
        print('Saving tsv: {}'.format(tsv_params_median_fn))
        df_params_median.to_csv(tsv_params_median_fn, sep="\t", na_rep='NaN', index=False)
        
        # Ecc.size
        # --------
        df_ecc_size = df_ecc_size.groupby(['roi', 'num_bins'], sort=False).median().reset_index()
        tsv_ecc_size_fn = "{}/{}_prf_ecc_size.tsv".format(tsv_dir, subject)
        print('Saving tsv: {}'.format(tsv_ecc_size_fn))
        df_ecc_size.to_csv(tsv_ecc_size_fn, sep="\t", na_rep='NaN', index=False)
    
        # Ecc.pCM
        # -------
        df_ecc_pcm = df_ecc_pcm.groupby(['roi', 'num_bins'], sort=False).median().reset_index()
        tsv_ecc_pcm_fn = "{}/{}_prf_ecc_pcm.tsv".format(tsv_dir, subject)
        print('Saving tsv: {}'.format(tsv_ecc_pcm_fn))
        df_ecc_pcm.to_csv(tsv_ecc_pcm_fn, sep="\t", na_rep='NaN', index=False)
    
        # Polar angle
        # -----------
        df_polar_angle = df_polar_angle.groupby(['roi', 'hemi', 'num_bins'], sort=False).median().reset_index()
        tsv_polar_angle_fn = "{}/{}_prf_polar_angle.tsv".format(tsv_dir, subject)
        print('Saving tsv: {}'.format(tsv_polar_angle_fn))
        df_polar_angle.to_csv(tsv_polar_angle_fn, sep="\t", na_rep='NaN', index=False)
    
        # Contralaterality
        # ----------------
        df_contralaterality = df_contralaterality.groupby(['roi'], sort=False).median().reset_index()
        tsv_contralaterality_fn = "{}/{}_prf_contralaterality.tsv".format(tsv_dir, subject)
        print('Saving tsv: {}'.format(tsv_contralaterality_fn))
        df_contralaterality.to_csv(tsv_contralaterality_fn, sep="\t", na_rep='NaN', index=False)
        
        # Spatial distribution 
        # -------------------
        tsv_distribution_fn = "{}/{}_prf_distribution.tsv".format(tsv_dir, subject)
        print('Saving tsv: {}'.format(tsv_distribution_fn))
        median_mesh = np.median(mesh_group, axis=0)
        df_distribution = pd.DataFrame(median_mesh)

        # Concatenating non-numeric columns back to the dataframe
        df_distribution = pd.concat([others_columns, df_distribution], axis=1)
        df_distribution.to_csv(tsv_distribution_fn, sep="\t", na_rep='NaN', index=False)
        
#        # Spatial distribution hot zone barycentre 
#        # ----------------------------------------
#        hemis = ['hemi-L', 'hemi-R', 'hemi-LR']
#        for j, hemi in enumerate(hemis):
#            hemi_values = ['hemi-L', 'hemi-R'] if hemi == 'hemi-LR' else [hemi]
#            df_distribution_hemi = df_distribution.loc[df_distribution.hemi.isin(hemi_values)]
#            df_barycentre_hemi = make_prf_barycentre_df(
#                df_distribution_hemi, rois, screen_side, gaussian_mesh_grain, hot_zone_percent=hot_zone_percent)
#            
#            df_barycentre_hemi['hemi'] = [hemi] * len(df_barycentre_hemi)
#            if j == 0: df_barycentre = df_barycentre_hemi
#            else: df_barycentre = pd.concat([df_barycentre, df_barycentre_hemi])
#            
#        tsv_barycentre_fn = "{}/{}_prf_barycentre.tsv".format(tsv_dir, subject)
#        print('Saving tsv: {}'.format(tsv_barycentre_fn))
#        df_barycentre.to_csv(tsv_barycentre_fn, sep="\t", na_rep='NaN', index=False)
    
# # Define permission cmd
# print('Changing files permissions in {}/{}'.format(main_dir, project_dir))
# os.system("chmod -Rf 771 {}/{}".format(main_dir, project_dir))
# os.system("chgrp -Rf {} {}/{}".format(group, main_dir, project_dir))    
