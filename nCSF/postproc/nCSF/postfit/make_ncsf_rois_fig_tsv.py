"""
-----------------------------------------------------------------------------------------
make_ncsf_rois_fig_tsv.py
-----------------------------------------------------------------------------------------
Goal of the script:
Make ROIs figure specific TSV
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject name (e.g. sub-01)
sys.argv[4]: server group (e.g. 327)
-----------------------------------------------------------------------------------------
Output(s):
nCSF figures tsv
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/pRF_analysis/nCSF/postproc/nCSF/postfit/
2. run python command
python make_ncsf_rois_fig_tsv.py [main directory] [project name] [subject] [group]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/projects/pRF_analysis/nCSF/postproc/nCSF/postfit/
python make_ncsf_rois_fig_tsv.py /scratch/mszinte/data nCSF sub-01 327
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
import numpy as np
import pandas as pd

# Personal import
sys.path.append("{}/../../../../analysis_code/utils".format(os.getcwd()))
from settings_utils import load_settings
from maths_utils import weighted_nan_median, weighted_nan_percentile

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
group = sys.argv[4]

# Load settings
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../../"))
settings_path = os.path.join(base_dir, project_dir, "settings.yml")
prf_settings_path = os.path.join(base_dir, project_dir, "prf-analysis.yml")
ncsf_settings_path = os.path.join(base_dir, project_dir, "nCSF-analysis.yml")
figure_settings_path = os.path.join(base_dir, project_dir, "figure-settings.yml")
settings = load_settings([settings_path, prf_settings_path, figure_settings_path, ncsf_settings_path])
analysis_info = settings[0]

# formats = analysis_info['formats']
# extensions = analysis_info['extensions']


formats = [analysis_info['formats'][0]] #######################################
extensions = [analysis_info['extensions'][0]]

rois_methods = analysis_info['rois_methods']
subjects_to_group = analysis_info['subjects']
preproc_prep = analysis_info['preproc_prep']
filtering = analysis_info['filtering']
normalization = analysis_info['normalization']
avg_methods = analysis_info['avg_methods']
stats_threshold = analysis_info['stats_th']
ncsf_task_name = analysis_info['nCSF_task_name']
ecc_SFp_num_bins = analysis_info['ecc_SFp_num_bins']
ecc_SFp_max = analysis_info['ecc_size_max'][0]
num_polar_SFp_bins = 8


# Main loop
for avg_method in avg_methods:
    if 'loo' in avg_method: 
        rsq2use_task = 'ncsf_loo_rsq'
        rsq2use = 'loo_rsq'
        
    else: 
        rsq2use_task = 'ncsf_rsq'
        rsq2use = 'rsq'
    rsq2use_median_task = 'ncsf_prf_median_{}'.format(rsq2use)
    
    for format_, extension in zip(formats, extensions):        
        # define list of rois for each format
        rois_methods_format = rois_methods[format_]
        for rois_method_format in rois_methods_format:

            if rois_method_format == 'rois-drawn':
                rois = analysis_info[rois_method_format]
            elif rois_method_format == 'rois-group-mmp':
                rois = list(analysis_info[rois_method_format].keys())
            
            # Individual subject analysis
            if 'group' not in subject:
                tsv_dir = '{}/{}/derivatives/pp_data/{}/{}/ncsf/tsv'.format(
                    main_dir, project_dir, subject, format_)
                
                # Exception if no data for one format (e.g template subject)
                if not os.path.isdir(tsv_dir):
                    print(f"[SKIP] tsv_dir not found for format={format_}: {tsv_dir}")
                    continue
                
                fn_spec = "task-{}_{}_{}_{}_{}_{}".format(
                    ncsf_task_name, preproc_prep, filtering,
                    normalization, avg_method, rois_method_format)
                tsv_fn = '{}/{}_{}_ncsf_deriv.tsv'.format(
                    tsv_dir, subject, fn_spec)
                data = pd.read_table(tsv_fn, sep="\t")
               
                # Keep a raw data df 
                data_raw = data.copy()
                
                # Threshold data (replace by nan)
                if stats_threshold == 0.05: stats_col = 'nCSF_corr_pvalue_5pt'
                elif stats_threshold == 0.01: stats_col = 'nCSF_corr_pvalue_1pt'
                data.loc[(data[stats_col] > stats_threshold)] = np.nan
                
                data = data.dropna()
                
                # median rsq between nCSf and pRF
                data[rsq2use_median_task] = (data[['ncsf_{}'.format(rsq2use), 'prf_{}'.format(rsq2use)]].median(axis=1))
                
                # Compute polar angle in degree
                data['prf_polar_angle'] = np.mod(np.degrees(np.angle(data.polar_real + 1j * data.polar_imag)), 360)

                # ROI active proportion
                # ---------------------
                
                # Number of vert per roi
                n_vert_tot_roi = (data_raw.groupby('roi', sort=False)
                                  .size()
                                  .reset_index(name='n_vert_tot'))
                
                # Number of significant vert for 0.05 p vale correction
                n_vert_corr_5pt = (data_raw[data_raw['nCSF_corr_pvalue_5pt'] < 0.05]
                                   .groupby('roi', sort=False)
                                   .size()
                                   .reset_index(name='n_vert_corr_pvalue_5pt'))
                
                # Number of significant vert for 0.01 p vale correction
                n_vert_corr_1pt = (data_raw[data_raw['nCSF_corr_pvalue_1pt'] < 0.01]
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
                tsv_roi_active_vert_fn = "{}/{}_{}_ncsf_active-vert.tsv".format(
                    tsv_dir, subject, fn_spec)
                print('Saving tsv: {}'.format(tsv_roi_active_vert_fn))
                df_roi_active_vert.to_csv(tsv_roi_active_vert_fn, sep="\t", na_rep='NaN', index=False)
            
                # Violins
                # -------
                df_violins = data
                tsv_violins_fn = "{}/{}_{}_ncsf_violins.tsv".format(
                    tsv_dir, subject, fn_spec)
                print('Saving tsv: {}'.format(tsv_violins_fn))
                df_violins.to_csv(tsv_violins_fn, sep="\t", na_rep='NaN', index=False)
                
                # Histogrammes
                # ------------
                df_histogrammes = data
                tsv_histogrammes_fn = "{}/{}_{}_ncsf_histogrammes.tsv".format(
                    tsv_dir, subject, fn_spec)
                print('Saving tsv: {}'.format(tsv_histogrammes_fn))
                df_histogrammes.to_csv(tsv_histogrammes_fn, sep="\t", na_rep='NaN', index=False)
                
                
                # Parameters median
                # ------------------
                for num_roi, roi in enumerate(rois):
                    df_roi = data.loc[(data.roi == roi)]
        
                    df_params_median_roi = pd.DataFrame()
                    df_params_median_roi['roi'] = [roi]

                    df_params_median_roi[f'{rsq2use_task}_weighted_median'] = weighted_nan_median(df_roi[rsq2use_task], weights=df_roi[rsq2use_task])
                    df_params_median_roi[f'{rsq2use_task}_ci_down'] = weighted_nan_percentile(df_roi[rsq2use_task], df_roi[rsq2use_task], 25)
                    df_params_median_roi[f'{rsq2use_task}_ci_up'] = weighted_nan_percentile(df_roi[rsq2use_task], df_roi[rsq2use_task], 75)
                    
                    df_params_median_roi['width_r_weighted_median'] = weighted_nan_median(df_roi.width_r, weights=df_roi[rsq2use_task])
                    df_params_median_roi['width_r_ci_down'] = weighted_nan_percentile(df_roi.width_r, df_roi[rsq2use_task], 25)
                    df_params_median_roi['width_r_ci_up'] = weighted_nan_percentile(df_roi.width_r, df_roi[rsq2use_task], 75)
                
                    df_params_median_roi['SFp_weighted_median'] = weighted_nan_median(df_roi.SFp, weights=df_roi[rsq2use_task])
                    df_params_median_roi['SFp_ci_down'] = weighted_nan_percentile(df_roi.SFp, df_roi[rsq2use_task], 25)
                    df_params_median_roi['SFp_ci_up'] = weighted_nan_percentile(df_roi.SFp, df_roi[rsq2use_task], 75)
                    
                    df_params_median_roi['CSp_median'] = weighted_nan_median(df_roi.CSp, weights=df_roi[rsq2use_task])
                    df_params_median_roi['CSp_ci_down'] = weighted_nan_percentile(df_roi.CSp, df_roi[rsq2use_task], 25)
                    df_params_median_roi['CSp_ci_up'] = weighted_nan_percentile(df_roi.CSp, df_roi[rsq2use_task], 75)
                     
                    df_params_median_roi['width_l_weighted_median'] = weighted_nan_median(df_roi.width_l, weights=df_roi[rsq2use_task])
                    df_params_median_roi['width_l_ci_down'] = weighted_nan_percentile(df_roi.width_l, df_roi[rsq2use_task], 25)
                    df_params_median_roi['width_l_ci_up'] = weighted_nan_percentile(df_roi.width_l, df_roi[rsq2use_task], 75)
                
                    df_params_median_roi['crf_exp_weighted_median'] = weighted_nan_median(df_roi.crf_exp, weights=df_roi[rsq2use_task])
                    df_params_median_roi['crf_exp_ci_down'] = weighted_nan_percentile(df_roi.crf_exp, df_roi[rsq2use_task], 25)
                    df_params_median_roi['crf_exp_ci_up'] = weighted_nan_percentile(df_roi.crf_exp, df_roi[rsq2use_task], 75)
                
                    df_params_median_roi['auc_weighted_median'] = weighted_nan_median(df_roi.auc, weights=df_roi[rsq2use_task])
                    df_params_median_roi['auc_ci_down'] = weighted_nan_percentile(df_roi.auc, df_roi[rsq2use_task], 25)
                    df_params_median_roi['auc_ci_up'] = weighted_nan_percentile(df_roi.auc, df_roi[rsq2use_task], 75)
                
                    df_params_median_roi['normalize_auc_weighted_median'] = weighted_nan_median(df_roi.normalize_auc, weights=df_roi[rsq2use_task])
                    df_params_median_roi['normalize_auc_ci_down'] = weighted_nan_percentile(df_roi.normalize_auc, df_roi[rsq2use_task], 25)
                    df_params_median_roi['normalize_auc_ci_up'] = weighted_nan_percentile(df_roi.normalize_auc, df_roi[rsq2use_task], 75)
                    
                    df_params_median_roi['SFmax_weighted_median'] = weighted_nan_median(df_roi.SFmax, weights=df_roi[rsq2use_task])
                    df_params_median_roi['SFmax_ci_down'] = weighted_nan_percentile(df_roi.SFmax, df_roi[rsq2use_task], 25)
                    df_params_median_roi['SFmax_ci_up'] = weighted_nan_percentile(df_roi.SFmax, df_roi[rsq2use_task], 75)
                            
                    if num_roi == 0: df_params_median = df_params_median_roi
                    else: df_params_median = pd.concat([df_params_median, df_params_median_roi])
                
                tsv_params_median_fn = "{}/{}_{}_ncsf_params-median.tsv".format(
                    tsv_dir, subject, fn_spec)
                print('Saving tsv: {}'.format(tsv_params_median_fn))
                df_params_median.to_csv(tsv_params_median_fn, sep="\t", na_rep='NaN', index=False)
                
                # Ecc.SFp
                # --------
                x = np.linspace(0, 1, ecc_SFp_num_bins + 1)
                ecc_bins = (x**2) * ecc_SFp_max
                
                for num_roi, roi in enumerate(rois):
                    df_roi = data.loc[(data.roi == roi)]
                    df_bins = df_roi.groupby(pd.cut(df_roi['prf_ecc'], bins=ecc_bins, include_lowest=True), observed=False)
                    df_ecc_SFp_bin = pd.DataFrame()
                    df_ecc_SFp_bin['roi'] = [roi]*ecc_SFp_num_bins
                    df_ecc_SFp_bin['num_bins'] = np.arange(ecc_SFp_num_bins)  
                    df_ecc_SFp_bin['prf_ecc_bins'] = df_bins.apply(lambda x: weighted_nan_median(x['prf_ecc'].values, x[rsq2use_median_task].values)).values
                    df_ecc_SFp_bin['SFp_bins_median'] = df_bins.apply(lambda x: weighted_nan_median(x['SFp'].values, x[rsq2use_median_task].values)).values
                    df_ecc_SFp_bin[f'{rsq2use_median_task}_bins_median'] = np.array(df_bins[rsq2use_median_task].median())
                    df_ecc_SFp_bin['SFp_bins_ci_upper_bound'] = df_bins.apply(lambda x: weighted_nan_percentile(x['SFp'].values, x[rsq2use_median_task].values, 75)).values
                    df_ecc_SFp_bin['SFp_bins_ci_lower_bound'] = df_bins.apply(lambda x: weighted_nan_percentile(x['SFp'].values, x[rsq2use_median_task].values, 25)).values
                    if num_roi == 0: df_ecc_SFp_bins = df_ecc_SFp_bin
                    else: df_ecc_SFp_bins = pd.concat([df_ecc_SFp_bins, df_ecc_SFp_bin]) 
                
                df_ecc_SFp = df_ecc_SFp_bins
                
                tsv_ecc_SFp_fn = "{}/{}_{}_ncsf_ecc-SFp.tsv".format(
                    tsv_dir, subject, fn_spec)
                print('Saving tsv: {}'.format(tsv_ecc_SFp_fn))
                df_ecc_SFp.to_csv(tsv_ecc_SFp_fn, sep="\t", na_rep='NaN', index=False)
                
                # Ecc.auc
                # --------
                x = np.linspace(0, 1, ecc_SFp_num_bins + 1)
                ecc_bins = (x**2) * ecc_SFp_max
                
                for num_roi, roi in enumerate(rois):
                    df_roi = data.loc[(data.roi == roi)]
                    df_bins = df_roi.groupby(pd.cut(df_roi['prf_ecc'], bins=ecc_bins, include_lowest=True), observed=False)
                    df_ecc_auc_bin = pd.DataFrame()
                    df_ecc_auc_bin['roi'] = [roi]*ecc_SFp_num_bins
                    df_ecc_auc_bin['num_bins'] = np.arange(ecc_SFp_num_bins)  
                    df_ecc_auc_bin['prf_ecc_bins'] = df_bins.apply(lambda x: weighted_nan_median(x['prf_ecc'].values, x[rsq2use_median_task].values)).values
                    df_ecc_auc_bin['auc_bins_median'] = df_bins.apply(lambda x: weighted_nan_median(x['auc'].values, x[rsq2use_median_task].values)).values
                    df_ecc_auc_bin[f'{rsq2use_median_task}_bins_median'] = np.array(df_bins[rsq2use_median_task].median())
                    df_ecc_auc_bin['auc_bins_ci_upper_bound'] = df_bins.apply(lambda x: weighted_nan_percentile(x['auc'].values, x[rsq2use_median_task].values, 75)).values
                    df_ecc_auc_bin['auc_bins_ci_lower_bound'] = df_bins.apply(lambda x: weighted_nan_percentile(x['auc'].values, x[rsq2use_median_task].values, 25)).values
                    if num_roi == 0: df_ecc_auc_bins = df_ecc_auc_bin
                    else: df_ecc_auc_bins = pd.concat([df_ecc_auc_bins, df_ecc_auc_bin]) 
                
                df_ecc_auc = df_ecc_auc_bins
                
                tsv_ecc_auc_fn = "{}/{}_{}_ncsf_ecc-auc.tsv".format(
                    tsv_dir, subject, fn_spec)
                print('Saving tsv: {}'.format(tsv_ecc_auc_fn))
                df_ecc_auc.to_csv(tsv_ecc_auc_fn, sep="\t", na_rep='NaN', index=False)
                
                
                # Size.SFp
                # --------
                x = np.linspace(0, 1, ecc_SFp_num_bins + 1)
                size_bins = (x**2) * ecc_SFp_max
                
                for num_roi, roi in enumerate(rois):
                    df_roi = data.loc[(data.roi == roi)]
                    df_bins = df_roi.groupby(pd.cut(df_roi['prf_size'], bins=size_bins, include_lowest=True), observed=False)
                    df_size_sfp_bin = pd.DataFrame()
                    df_size_sfp_bin['roi'] = [roi]*ecc_SFp_num_bins
                    df_size_sfp_bin['num_bins'] = np.arange(ecc_SFp_num_bins)  
                    df_size_sfp_bin['prf_size_bins'] = df_bins.apply(lambda x: weighted_nan_median(x['prf_size'].values, x[rsq2use_median_task].values)).values
                    df_size_sfp_bin['sfp_bins_median'] = df_bins.apply(lambda x: weighted_nan_median(x['SFp'].values, x[rsq2use_median_task].values)).values
                    df_size_sfp_bin[f'{rsq2use_median_task}_bins_median'] = np.array(df_bins[rsq2use_median_task].median())
                    df_size_sfp_bin['sfp_bins_ci_upper_bound'] = df_bins.apply(lambda x: weighted_nan_percentile(x['SFp'].values, x[rsq2use_median_task].values, 75)).values
                    df_size_sfp_bin['sfp_bins_ci_lower_bound'] = df_bins.apply(lambda x: weighted_nan_percentile(x['SFp'].values, x[rsq2use_median_task].values, 25)).values
                    if num_roi == 0: df_size_sfp_bins = df_size_sfp_bin
                    else: df_size_sfp_bins = pd.concat([df_size_sfp_bins, df_size_sfp_bin]) 
                
                df_size_sfp = df_size_sfp_bins
                
                tsv_size_sfp_fn = "{}/{}_{}_ncsf_size-SFp.tsv".format(
                    tsv_dir, subject, fn_spec)
                print('Saving tsv: {}'.format(tsv_size_sfp_fn))
                df_size_sfp.to_csv(tsv_size_sfp_fn, sep="\t", na_rep='NaN', index=False)
                
                
                # Size.auc
                # --------
                x = np.linspace(0, 1, ecc_SFp_num_bins + 1)
                size_bins = (x**2) * ecc_SFp_max

                for num_roi, roi in enumerate(rois):
                    df_roi = data.loc[(data.roi == roi)]
                    df_bins = df_roi.groupby(pd.cut(df_roi['prf_size'], bins=size_bins, include_lowest=True), observed=False)
                    df_size_auc_bin = pd.DataFrame()
                    df_size_auc_bin['roi'] = [roi]*ecc_SFp_num_bins
                    df_size_auc_bin['num_bins'] = np.arange(ecc_SFp_num_bins)  
                    df_size_auc_bin['prf_size_bins'] = df_bins.apply(lambda x: weighted_nan_median(x['prf_size'].values, x[rsq2use_median_task].values)).values
                    df_size_auc_bin['auc_bins_median'] = df_bins.apply(lambda x: weighted_nan_median(x['auc'].values, x[rsq2use_median_task].values)).values
                    df_size_auc_bin[f'{rsq2use_median_task}_bins_median'] = np.array(df_bins[rsq2use_median_task].median())
                    df_size_auc_bin['auc_bins_ci_upper_bound'] = df_bins.apply(lambda x: weighted_nan_percentile(x['auc'].values, x[rsq2use_median_task].values, 75)).values
                    df_size_auc_bin['auc_bins_ci_lower_bound'] = df_bins.apply(lambda x: weighted_nan_percentile(x['auc'].values, x[rsq2use_median_task].values, 25)).values
                    if num_roi == 0: df_size_auc_bins = df_size_auc_bin
                    else: df_size_auc_bins = pd.concat([df_size_auc_bins, df_size_auc_bin]) 
                
                df_size_auc = df_size_auc_bins
                
                tsv_size_auc_fn = "{}/{}_{}_ncsf_size-auc.tsv".format(
                    tsv_dir, subject, fn_spec)
                print('Saving tsv: {}'.format(tsv_size_auc_fn))
                df_size_auc.to_csv(tsv_size_auc_fn, sep="\t", na_rep='NaN', index=False)
                
                # Polar.SFp
                # ---------
                for num_roi, roi in enumerate(rois):
                    df_roi = data.loc[(data.roi == roi)]
                    df_bins = df_roi.groupby(pd.cut(df_roi['prf_polar_angle'], bins=num_polar_SFp_bins))
                    df_polar_SFp_bin = pd.DataFrame()
                    df_polar_SFp_bin['roi'] = [roi]*num_polar_SFp_bins
                    df_polar_SFp_bin['num_bins'] = np.arange(num_polar_SFp_bins)
                    df_polar_SFp_bin['prf_polar_bins'] = df_bins.apply(lambda x: weighted_nan_median(x['prf_polar_angle'].values, x[rsq2use_median_task].values)).values
                    df_polar_SFp_bin['SFp_bins_median'] = df_bins.apply(lambda x: weighted_nan_median(x['SFp'].values, x[rsq2use_median_task].values)).values
                    df_polar_SFp_bin[f'{rsq2use_median_task}_bins_median'] = np.array(df_bins[rsq2use_median_task].median())
                    df_polar_SFp_bin['SFp_bins_ci_upper_bound'] = df_bins.apply(lambda x: weighted_nan_percentile(x['SFp'].values, x[rsq2use_median_task].values, 75)).values
                    df_polar_SFp_bin['SFp_bins_ci_lower_bound'] = df_bins.apply(lambda x: weighted_nan_percentile(x['SFp'].values, x[rsq2use_median_task].values, 25)).values
                
                    if num_roi == 0: df_polar_SFp_bins = df_polar_SFp_bin
                    else: df_polar_SFp_bins = pd.concat([df_polar_SFp_bins, df_polar_SFp_bin]) 
                
                df_polar_sfp = df_polar_SFp_bins
                
                tsv_polar_sfp_fn = "{}/{}_{}_ncsf_polar-SFp.tsv".format(
                    tsv_dir, subject, fn_spec)
                print('Saving tsv: {}'.format(tsv_polar_sfp_fn))
                df_polar_sfp.to_csv(tsv_polar_sfp_fn, sep="\t", na_rep='NaN', index=False)
                             
            # Group Analysis    
            else :
                print('group')
                for i, subject_to_group in enumerate(subjects_to_group):
                    tsv_dir = '{}/{}/derivatives/pp_data/{}/{}/ncsf/tsv'.format(
                        main_dir, project_dir, subject_to_group, format_)

                    fn_spec = "task-{}_{}_{}_{}_{}_{}".format(
                        ncsf_task_name, preproc_prep, filtering, normalization, avg_method, rois_method_format)
            
                    # ROI active vertices
                    # -------------------
                    tsv_roi_area_fn = "{}/{}_{}_ncsf_active-vert.tsv".format(
                        tsv_dir, subject_to_group, fn_spec)
                    df_roi_area_indiv = pd.read_table(tsv_roi_area_fn, sep="\t")
                    if i == 0: df_roi_area = df_roi_area_indiv.copy()
                    else: df_roi_area = pd.concat([df_roi_area, df_roi_area_indiv])
            
                    # Violins
                    # -------
                    tsv_violins_fn = "{}/{}_{}_ncsf_violins.tsv".format(
                        tsv_dir, subject_to_group, fn_spec)
                    df_violins_indiv = pd.read_table(tsv_violins_fn, sep="\t")
                    if i == 0: df_violins = df_violins_indiv.copy()
                    else: df_violins = pd.concat([df_violins, df_violins_indiv])
            
                    # Parameters median
                    # ------------------
                    # use df_violins see below
            
                    # Ecc.Sfp
                    # --------
                    tsv_ecc_sfp_fn = "{}/{}_{}_ncsf_ecc-SFp.tsv".format(
                        tsv_dir, subject_to_group, fn_spec)
                    df_ecc_sfp_indiv = pd.read_table(tsv_ecc_sfp_fn, sep="\t")
                    if i == 0: df_ecc_sfp = df_ecc_sfp_indiv.copy()
                    else: df_ecc_sfp = pd.concat([df_ecc_sfp, df_ecc_sfp_indiv])

                    # Ecc.auc
                    # -------
                    tsv_ecc_auc_fn = "{}/{}_{}_ncsf_ecc-auc.tsv".format(
                        tsv_dir, subject_to_group, fn_spec)
                    df_ecc_auc_indiv = pd.read_table(tsv_ecc_auc_fn, sep="\t")
                    if i == 0: df_ecc_auc = df_ecc_auc_indiv.copy()
                    else: df_ecc_auc = pd.concat([df_ecc_auc, df_ecc_auc_indiv])
              
                # Median and saving tsv
                tsv_dir = '{}/{}/derivatives/pp_data/{}/{}/ncsf/tsv'.format(
                    main_dir, project_dir, subject, format_)
                os.makedirs(tsv_dir, exist_ok=True)
                
                # ROI surface areas 
                # -----------------
                df_roi_area = df_roi_area.groupby(['roi'], sort=False).median().reset_index()
                tsv_roi_area_fn =  "{}/{}_{}_ncsf_active-vert.tsv".format(
                    tsv_dir, subject, fn_spec)
                print('Saving tsv: {}'.format(tsv_roi_area_fn))
                df_roi_area.to_csv(tsv_roi_area_fn, sep="\t", na_rep='NaN', index=False)
                
                # Violins
                # -------
                df_violins = df_violins # no averaging
                tsv_violins_fn = "{}/{}_{}_ncsf_violins.tsv".format(
                        tsv_dir, subject, fn_spec)
                print('Saving tsv: {}'.format(tsv_violins_fn))
                df_violins.to_csv(tsv_violins_fn, sep="\t", na_rep='NaN', index=False)

                # Parameters median
                # ------------------
                df_params_median = df_violins
                
                # compute median 
                colnames = [rsq2use_task, 'width_r', 'SFp', 'CSp', 'width_l', 'crf_exp', 'auc', 'normalize_auc']
                df_params_median_indiv = df_params_median.groupby(['roi', 'subject'])[[rsq2use_task]].apply(
                    lambda x: weighted_nan_median(x[rsq2use_task], df_params_median.loc[x.index, rsq2use_task])).reset_index(name=f'{rsq2use_task}_weighted_median')
                
                for colname in colnames[1:]:
                    df_params_median_indiv['{}_weighted_median'.format(colname)] = df_params_median.groupby(['roi', 'subject'])[[colname, rsq2use_task]].apply(
                        lambda x: weighted_nan_median(x[colname], df_params_median.loc[x.index, rsq2use_task])).reset_index()[0]
                df_params_med_median = df_params_median_indiv.groupby(['roi'])[[colname + '_weighted_median' for colname in colnames]].median()
        
                # compute Ci
                df_params_median_ci = pd.DataFrame()
                for colname in colnames:
                    df_params_median_ci['{}_ci_down'.format(colname)] = df_params_median_indiv.groupby(['roi']).apply(
                        lambda x: weighted_nan_percentile(x['{}_weighted_median'.format(colname)], x[f'{rsq2use_task}_weighted_median'], 2.5)) 
                    df_params_median_ci['{}_ci_up'.format(colname)] = df_params_median_indiv.groupby(['roi']).apply(
                        lambda x: weighted_nan_percentile(x['{}_weighted_median'.format(colname)], x[f'{rsq2use_task}_weighted_median'], 97.5)) 
        
                df_params_median = pd.concat([df_params_med_median, df_params_median_ci], axis=1).reset_index()
                tsv_params_median_fn = "{}/{}_{}_ncsf_params-median.tsv".format(
                        tsv_dir, subject, fn_spec)
                print('Saving tsv: {}'.format(tsv_params_median_fn))
                df_params_median.to_csv(tsv_params_median_fn, sep="\t", na_rep='NaN', index=False)
                
                # Ecc.SFp
                # --------
                df_ecc_sfp = df_ecc_sfp.groupby(['roi', 'num_bins'], sort=False).median().reset_index()
                tsv_ecc_sfp_fn = "{}/{}_{}_ncsf_ecc-SFp.tsv".format(
                        tsv_dir, subject, fn_spec)
                print('Saving tsv: {}'.format(tsv_ecc_sfp_fn))
                df_ecc_sfp.to_csv(tsv_ecc_sfp_fn, sep="\t", na_rep='NaN', index=False)
            
                # Ecc.pCM
                # -------
                df_ecc_auc = df_ecc_auc.groupby(['roi', 'num_bins'], sort=False).median().reset_index()
                tsv_ecc_auc_fn = "{}/{}_{}_ncsf_ecc-pcm.tsv".format(
                        tsv_dir, subject, fn_spec)
                print('Saving tsv: {}'.format(tsv_ecc_auc_fn))
                df_ecc_auc.to_csv(tsv_ecc_sfp_fn, sep="\t", na_rep='NaN', index=False)
            
# Define permission cmd
print('Changing files permissions in {}/{}'.format(main_dir, project_dir))
os.system("chmod -Rf 771 {}/{}".format(main_dir, project_dir))
os.system("chgrp -Rf {} {}/{}".format(group, main_dir, project_dir))    