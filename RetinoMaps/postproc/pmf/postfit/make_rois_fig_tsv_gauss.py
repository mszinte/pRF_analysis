"""
-----------------------------------------------------------------------------------------
make_rois_fig_tsv_gauss.py
-----------------------------------------------------------------------------------------
Goal of the script:
Make ROIs figure-specific TSVs for Gaussian pRF fit
(no stats file, no pCM — thresholding on rsq, ecc, size, amplitude only)
Produces: active-vert, violins, params-median, ecc-size,
          polar-angle, contralaterality, distribution, barycentre
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject name (e.g. sub-01 or group)
sys.argv[4]: server group (e.g. 327)
-----------------------------------------------------------------------------------------
Output(s):
Multiple TSV files per subject / format / task
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/pRF_analysis/RetinoMaps/postproc/pmf/postfit/
2. run python command
python make_rois_fig_tsv_gauss.py [main directory] [project name] [subject] [group]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/projects/pRF_analysis/RetinoMaps/postproc/pmf/postfit/
python make_rois_fig_tsv_gauss.py /scratch/mszinte/data RetinoMaps sub-01 327 bold 
python make_rois_fig_tsv_gauss.py /scratch/mszinte/data RetinoMaps sub-01 327 residuals
python make_rois_fig_tsv_gauss.py /scratch/mszinte/data RetinoMaps sub-hcp1.6mm 327
python make_rois_fig_tsv_gauss.py /scratch/mszinte/data RetinoMaps group 327
-----------------------------------------------------------------------------------------
Written by Uriel Lascombes (uriel.lascombes@laposte.net)
Edited by Martin Szinte (martin.szinte@gmail.com)
and Sina Kling (sina.kling@outlook.de)
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
script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.abspath(os.path.join(script_dir, "../../../../analysis_code/utils")))
from settings_utils import load_settings
from maths_utils import (make_gauss_prf_distribution_df, weighted_nan_median,
                         weighted_nan_percentile, make_prf_barycentre_df)

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
group = sys.argv[4]
deriv_condition = sys.argv[5]

# Load settings
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../../"))
settings_path = os.path.join(base_dir, project_dir, "settings.yml")
prf_settings_path = os.path.join(base_dir, project_dir, "pmf-settings.yml")
figure_settings_path = os.path.join(base_dir, project_dir, "figure-settings.yml")
settings = load_settings([settings_path, prf_settings_path, figure_settings_path])
analysis_info = settings[0]

formats = analysis_info['formats']
extensions = analysis_info['extensions']
rois_methods = analysis_info['rois_methods']
ecc_threshold = analysis_info['ecc_th']
size_threshold = analysis_info['size_th']
rsqr_threshold = analysis_info['rsqr_th']
stats_threshold = analysis_info['stats_th']
amplitude_threshold = analysis_info['prf_amp_th']
subjects_to_group = analysis_info['subjects']
preproc_prep = analysis_info['preproc_prep']
filtering = analysis_info['filtering']
normalization = analysis_info['normalization']
avg_methods = ['concat']
prf_task_names = ['SacLoc']
ecc_size_num_bins = analysis_info['ecc_size_num_bins']
ecc_size_max = analysis_info['ecc_size_max']
polar_angle_num_bins = analysis_info['polar_angle_num_bins']
distribution_max_ecc = analysis_info['distribution_max_ecc']
distribution_mesh_grain = analysis_info['distribution_mesh_grain']
hot_zone_percent = analysis_info['hot_zone_percent']

# Define fit_type-specific filename patterns
fit_type_config = {
    'residuals': {
        'deriv_tag':  'pmf-residuals-gauss_deriv',
        'output_tag': 'pmf-residuals',
    },
    'bold': {
        'deriv_tag':  'pmf2-gauss_deriv',
        'output_tag': 'pmf-bold',
    },
    'prf': {
        'deriv_tag':  'prf-gauss_deriv',
        'output_tag': 'prf-control',
    },
}

cfg        = fit_type_config[deriv_condition]
deriv_tag  = cfg['deriv_tag']
output_tag = cfg['output_tag']

# NOTE: No stats threshold, no pCM — Gaussian fit only uses rsq, ecc, size, amplitude

# rsq column is always 'prf_rsq' for Gaussian (no loo variant stored separately)
rsq2use = 'prf_rsq'

# Main loop
for avg_method in avg_methods:

    for format_, extension in zip(formats, extensions):

        # Define list of rois for each format
        rois_methods_format = rois_methods[format_]
        for rois_method_format in rois_methods_format:

            if rois_method_format == 'rois-drawn':
                rois = analysis_info[rois_method_format]
            elif rois_method_format == 'rois-group-mmp':
                rois = list(analysis_info[rois_method_format].keys())

            for prf_task_name in prf_task_names:

                print(f'{prf_task_name} - {avg_method} - {format_}')

                # ------------------------------------------------------------------
                # Individual subject analysis
                # ------------------------------------------------------------------
                if 'group' not in subject:
                    tsv_dir = '{}/{}/derivatives/pp_data/{}/{}/pmf/tsv'.format(
                        main_dir, project_dir, subject, format_)

                    if not os.path.isdir(tsv_dir):
                        print(f"[SKIP] tsv_dir not found for format={format_}: {tsv_dir}")
                        continue

                    fn_spec = "task-{}_{}_{}_{}_{}_{}".format(
                        prf_task_name, preproc_prep, filtering,
                        normalization, avg_method, rois_method_format)

                    tsv_fn = '{}/{}_{}_{}.tsv'.format(
                        tsv_dir, subject, fn_spec, deriv_tag)
                    data = pd.read_table(tsv_fn, sep="\t")

                    # Keep a raw data df (for active-vert counting)
                    data_raw = data.copy()

                    # ------------------------------------------------------------------
                    # Threshold data (replace rows violating any criterion with NaN)
                    # ------------------------------------------------------------------
                    if stats_threshold == 0.05: stats_col = 'corr_pvalue_5pt'
                    elif stats_threshold == 0.01: stats_col = 'corr_pvalue_1pt'
                    data.loc[
                        (data['amplitude'] < amplitude_threshold[0]) |
                        (data['prf_ecc'] < ecc_threshold[0]) |
                        (data['prf_ecc'] > ecc_threshold[1]) |
                        (data['prf_size'] < size_threshold[0]) |
                        (data['prf_size'] > size_threshold[1]) |
                        (data[rsq2use] < rsqr_threshold)
                    ] = np.nan

                    data = data.dropna()

                    # ------------------------------------------------------------------
                    # ROI active proportion
                    # ------------------------------------------------------------------
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
                    tsv_roi_active_vert_fn = "{}/{}_{}_{}-gauss_active-vert.tsv".format(
                        tsv_dir, subject, fn_spec, output_tag)
                    print('Saving tsv: {}'.format(tsv_roi_active_vert_fn))
                    df_roi_active_vert.to_csv(tsv_roi_active_vert_fn, sep="\t", na_rep='NaN', index=False)
                
                    # ------------------------------------------------------------------
                    # Violins (thresholded data — one row per vertex)
                    # ------------------------------------------------------------------
                    df_violins = data
                    tsv_violins_fn = "{}/{}_{}_{}-gauss_violins.tsv".format(
                        tsv_dir, subject, fn_spec, output_tag)
                    print('Saving tsv: {}'.format(tsv_violins_fn))
                    df_violins.to_csv(
                        tsv_violins_fn, sep="\t", na_rep='NaN', index=False)

                    # ------------------------------------------------------------------
                    # Parameters median (rsq, size, ecc — no pCM, no n)
                    # ------------------------------------------------------------------
                    for num_roi, roi in enumerate(rois):
                        df_roi = data.loc[(data.roi == roi)]

                        df_params_median_roi = pd.DataFrame()
                        df_params_median_roi['roi'] = [roi]

                        df_params_median_roi[f'{rsq2use}_weighted_median'] = weighted_nan_median(
                            df_roi[rsq2use], weights=df_roi[rsq2use])
                        df_params_median_roi[f'{rsq2use}_ci_down'] = weighted_nan_percentile(
                            df_roi[rsq2use], df_roi[rsq2use], 25)
                        df_params_median_roi[f'{rsq2use}_ci_up'] = weighted_nan_percentile(
                            df_roi[rsq2use], df_roi[rsq2use], 75)

                        df_params_median_roi['prf_size_weighted_median'] = weighted_nan_median(
                            df_roi['prf_size'], weights=df_roi[rsq2use])
                        df_params_median_roi['prf_size_ci_down'] = weighted_nan_percentile(
                            df_roi['prf_size'], df_roi[rsq2use], 25)
                        df_params_median_roi['prf_size_ci_up'] = weighted_nan_percentile(
                            df_roi['prf_size'], df_roi[rsq2use], 75)

                        df_params_median_roi['prf_ecc_weighted_median'] = weighted_nan_median(
                            df_roi['prf_ecc'], weights=df_roi[rsq2use])
                        df_params_median_roi['prf_ecc_ci_down'] = weighted_nan_percentile(
                            df_roi['prf_ecc'], df_roi[rsq2use], 25)
                        df_params_median_roi['prf_ecc_ci_up'] = weighted_nan_percentile(
                            df_roi['prf_ecc'], df_roi[rsq2use], 75)

                        if num_roi == 0:
                            df_params_median = df_params_median_roi
                        else:
                            df_params_median = pd.concat(
                                [df_params_median, df_params_median_roi])

                    tsv_params_median_fn = "{}/{}_{}_{}-gauss_params-median.tsv".format(
                        tsv_dir, subject, fn_spec, output_tag)
                    print('Saving tsv: {}'.format(tsv_params_median_fn))
                    df_params_median.to_csv(
                        tsv_params_median_fn, sep="\t", na_rep='NaN', index=False)

                    # ------------------------------------------------------------------
                    # Ecc. × Size
                    # ------------------------------------------------------------------
                    ecc_bins = np.linspace(0, ecc_size_max[1], ecc_size_num_bins + 1)
                    for num_roi, roi in enumerate(rois):
                        df_roi = data.loc[(data.roi == roi)]
                        df_bins = df_roi.groupby(pd.cut(df_roi['prf_ecc'], bins=ecc_bins))

                        df_ecc_size_bin = pd.DataFrame()
                        df_ecc_size_bin['roi'] = [roi] * ecc_size_num_bins
                        df_ecc_size_bin['num_bins'] = np.arange(ecc_size_num_bins)
                        df_ecc_size_bin['prf_ecc_bins'] = df_bins.apply(
                            lambda x: weighted_nan_median(
                                x['prf_ecc'].values, x[rsq2use].values)).values
                        df_ecc_size_bin['prf_size_bins_median'] = df_bins.apply(
                            lambda x: weighted_nan_median(
                                x['prf_size'].values, x[rsq2use].values)).values
                        df_ecc_size_bin[f'{rsq2use}_bins_median'] = np.array(
                            df_bins[rsq2use].median())
                        df_ecc_size_bin['prf_size_bins_ci_upper_bound'] = df_bins.apply(
                            lambda x: weighted_nan_percentile(
                                x['prf_size'].values, x[rsq2use].values, 75)).values
                        df_ecc_size_bin['prf_size_bins_ci_lower_bound'] = df_bins.apply(
                            lambda x: weighted_nan_percentile(
                                x['prf_size'].values, x[rsq2use].values, 25)).values

                        if num_roi == 0:
                            df_ecc_size_bins = df_ecc_size_bin
                        else:
                            df_ecc_size_bins = pd.concat([df_ecc_size_bins, df_ecc_size_bin])

                    tsv_ecc_size_fn = "{}/{}_{}_{}-gauss_ecc-size.tsv".format(
                        tsv_dir, subject, fn_spec, output_tag)
                    print('Saving tsv: {}'.format(tsv_ecc_size_fn))
                    df_ecc_size_bins.to_csv(
                        tsv_ecc_size_fn, sep="\t", na_rep='NaN', index=False)

                    # ------------------------------------------------------------------
                    # Polar angle
                    # ------------------------------------------------------------------
                    theta_slices = np.linspace(0, 360, polar_angle_num_bins, endpoint=False)
                    data['prf_polar_angle'] = np.mod(
                        np.degrees(np.angle(
                            data['polar_real'] + 1j * data['polar_imag'])), 360)

                    hemis = ['hemi-L', 'hemi-R', 'hemi-LR']
                    for i, hemi in enumerate(hemis):
                        hemi_values = ['hemi-L', 'hemi-R'] if hemi == 'hemi-LR' else [hemi]
                        for j, roi in enumerate(rois):
                            df = data.loc[
                                (data.roi == roi) & (data.hemi.isin(hemi_values))]
                            if len(df):
                                df_bins = df.groupby(
                                    pd.cut(df['prf_polar_angle'], bins=polar_angle_num_bins))
                                rsq_sum = df_bins[rsq2use].sum()
                            else:
                                rsq_sum = [np.nan] * polar_angle_num_bins

                            df_polar_angle_bin = pd.DataFrame()
                            df_polar_angle_bin['roi'] = [roi] * polar_angle_num_bins
                            df_polar_angle_bin['hemi'] = [hemi] * polar_angle_num_bins
                            df_polar_angle_bin['num_bins'] = np.arange(polar_angle_num_bins)
                            df_polar_angle_bin['theta_slices'] = np.array(theta_slices)
                            df_polar_angle_bin['rsq_sum'] = np.array(rsq_sum)

                            if j == 0 and i == 0:
                                df_polar_angle_bins = df_polar_angle_bin
                            else:
                                df_polar_angle_bins = pd.concat(
                                    [df_polar_angle_bins, df_polar_angle_bin])

                    tsv_polar_angle_fn = "{}/{}_{}_{}-gauss_polar-angle.tsv".format(
                        tsv_dir, subject, fn_spec, output_tag)
                    print('Saving tsv: {}'.format(tsv_polar_angle_fn))
                    df_polar_angle_bins.to_csv(
                        tsv_polar_angle_fn, sep="\t", na_rep='NaN', index=False)

                    # ------------------------------------------------------------------
                    # Contralaterality
                    # ------------------------------------------------------------------
                    for j, roi in enumerate(rois):
                        df_rh = data.loc[(data.roi == roi) & (data.hemi == 'hemi-R')]
                        df_lh = data.loc[(data.roi == roi) & (data.hemi == 'hemi-L')]
                        try:
                            contralaterality_prct = (
                                sum(df_rh.loc[df_rh['prf_x'] < 0][rsq2use]) +
                                sum(df_lh.loc[df_lh['prf_x'] > 0][rsq2use])
                            ) / (sum(df_rh[rsq2use]) + sum(df_lh[rsq2use]))
                        except Exception:
                            contralaterality_prct = np.nan

                        df_contralaterality_roi = pd.DataFrame()
                        df_contralaterality_roi['roi'] = [roi]
                        df_contralaterality_roi['contralaterality_prct'] = np.array(
                            contralaterality_prct)

                        if j == 0:
                            df_contralaterality = df_contralaterality_roi
                        else:
                            df_contralaterality = pd.concat(
                                [df_contralaterality, df_contralaterality_roi])

                    tsv_contralaterality_fn = "{}/{}_{}_{}-gauss_contralaterality.tsv".format(
                        tsv_dir, subject, fn_spec, output_tag)
                    print('Saving tsv: {}'.format(tsv_contralaterality_fn))
                    df_contralaterality.to_csv(
                        tsv_contralaterality_fn, sep="\t", na_rep='NaN', index=False)

                    # ------------------------------------------------------------------
                    # Spatial distribution
                    # ------------------------------------------------------------------
                    hemis = ['hemi-L', 'hemi-R', 'hemi-LR']
                    for i, hemi in enumerate(hemis):
                        hemi_values = ['hemi-L', 'hemi-R'] if hemi == 'hemi-LR' else [hemi]
                        data_hemi = data.loc[data.hemi.isin(hemi_values)]
                        df_distribution_hemi = make_gauss_prf_distribution_df(
                            data_hemi, rois, distribution_max_ecc,
                            distribution_mesh_grain, rsq2use)

                        df_distribution_hemi['hemi'] = [hemi] * len(df_distribution_hemi)
                        if i == 0:
                            df_distribution = df_distribution_hemi
                        else:
                            df_distribution = pd.concat(
                                [df_distribution, df_distribution_hemi])

                    tsv_distribution_fn = "{}/{}_{}_{}-gauss_distribution.tsv".format(
                        tsv_dir, subject, fn_spec, output_tag)
                    print('Saving tsv: {}'.format(tsv_distribution_fn))
                    df_distribution.to_csv(
                        tsv_distribution_fn, sep="\t", na_rep='NaN', index=False)

                    # ------------------------------------------------------------------
                    # Spatial distribution hot-zone barycentre
                    # ------------------------------------------------------------------
                    hemis = ['hemi-L', 'hemi-R', 'hemi-LR']
                    for i, hemi in enumerate(hemis):
                        hemi_values = ['hemi-L', 'hemi-R'] if hemi == 'hemi-LR' else [hemi]
                        df_distribution_hemi = df_distribution.loc[
                            df_distribution.hemi.isin(hemi_values)]
                        df_barycentre_hemi = make_prf_barycentre_df(
                            df_distribution_hemi, rois, distribution_max_ecc,
                            distribution_mesh_grain, hot_zone_percent=hot_zone_percent)

                        df_barycentre_hemi['hemi'] = [hemi] * len(df_barycentre_hemi)
                        if i == 0:
                            df_barycentre = df_barycentre_hemi
                        else:
                            df_barycentre = pd.concat([df_barycentre, df_barycentre_hemi])

                    tsv_barycentre_fn = "{}/{}_{}_{}-gauss_barycentre.tsv".format(
                        tsv_dir, subject, fn_spec, output_tag)
                    print('Saving tsv: {}'.format(tsv_barycentre_fn))
                    df_barycentre.to_csv(
                        tsv_barycentre_fn, sep="\t", na_rep='NaN', index=False)

                # ------------------------------------------------------------------
                # Group analysis
                # ------------------------------------------------------------------
                else:
                    print('group')
                    for i, subject_to_group in enumerate(subjects_to_group):
                        tsv_dir = '{}/{}/derivatives/pp_data/{}/{}/pmf/tsv'.format(
                            main_dir, project_dir, subject_to_group, format_)

                        fn_spec = "task-{}_{}_{}_{}_{}_{}".format(
                            prf_task_name, preproc_prep, filtering,
                            normalization, avg_method, rois_method_format)

                        # ROI active vertices
                        tsv_roi_area_fn = "{}/{}_{}_{}-gauss_active-vert.tsv".format(
                            tsv_dir, subject_to_group, fn_spec, output_tag)
                        df_roi_area_indiv = pd.read_table(tsv_roi_area_fn, sep="\t")
                        if i == 0:
                            df_roi_area = df_roi_area_indiv.copy()
                        else:
                            df_roi_area = pd.concat([df_roi_area, df_roi_area_indiv])

                        # Violins
                        tsv_violins_fn = "{}/{}_{}_{}-gauss_violins.tsv".format(
                            tsv_dir, subject_to_group, fn_spec, output_tag)
                        df_violins_indiv = pd.read_table(tsv_violins_fn, sep="\t")
                        if i == 0:
                            df_violins = df_violins_indiv.copy()
                        else:
                            df_violins = pd.concat([df_violins, df_violins_indiv])

                        # Ecc. × Size
                        tsv_ecc_size_fn = "{}/{}_{}_{}-gauss_ecc-size.tsv".format(
                            tsv_dir, subject_to_group, fn_spec, output_tag)
                        df_ecc_size_indiv = pd.read_table(tsv_ecc_size_fn, sep="\t")
                        if i == 0:
                            df_ecc_size = df_ecc_size_indiv.copy()
                        else:
                            df_ecc_size = pd.concat([df_ecc_size, df_ecc_size_indiv])

                        # Polar angle
                        tsv_polar_angle_fn = "{}/{}_{}_{}-gauss_polar-angle.tsv".format(
                            tsv_dir, subject_to_group, fn_spec, output_tag)
                        df_polar_angle_indiv = pd.read_table(tsv_polar_angle_fn, sep="\t")
                        if i == 0:
                            df_polar_angle = df_polar_angle_indiv.copy()
                        else:
                            df_polar_angle = pd.concat([df_polar_angle, df_polar_angle_indiv])

                        # Contralaterality
                        tsv_contralaterality_fn = "{}/{}_{}_{}-gauss_contralaterality.tsv".format(
                            tsv_dir, subject_to_group, fn_spec, output_tag)
                        df_contralaterality_indiv = pd.read_table(
                            tsv_contralaterality_fn, sep="\t")
                        if i == 0:
                            df_contralaterality = df_contralaterality_indiv.copy()
                        else:
                            df_contralaterality = pd.concat(
                                [df_contralaterality, df_contralaterality_indiv])

                        # Spatial distribution
                        tsv_distribution_fn = "{}/{}_{}_{}-gauss_distribution.tsv".format(
                            tsv_dir, subject_to_group, fn_spec, output_tag)
                        df_distribution_indiv = pd.read_table(tsv_distribution_fn, sep="\t")
                        mesh_indiv = df_distribution_indiv.drop(
                            columns=['roi', 'x', 'y', 'hemi']).values
                        others_columns = df_distribution_indiv[['roi', 'x', 'y', 'hemi']]

                        if i == 0:
                            mesh_group = np.expand_dims(mesh_indiv, axis=0)
                        else:
                            mesh_group = np.vstack(
                                (mesh_group, np.expand_dims(mesh_indiv, axis=0)))

                    # Group output directory
                    tsv_dir = '{}/{}/derivatives/pp_data/{}/{}/pmf/tsv'.format(
                        main_dir, project_dir, subject, format_)
                    os.makedirs(tsv_dir, exist_ok=True)

                    fn_spec = "task-{}_{}_{}_{}_{}_{}".format(
                        prf_task_name, preproc_prep, filtering,
                        normalization, avg_method, rois_method_format)

                    # ROI active vertices (median across subjects)
                    df_roi_area = (df_roi_area
                                   .groupby(['roi'], sort=False)
                                   .median()
                                   .reset_index())
                    tsv_roi_area_fn = "{}/{}_{}_{}-gauss_active-vert.tsv".format(
                        tsv_dir, subject, fn_spec, output_tag)
                    print('Saving tsv: {}'.format(tsv_roi_area_fn))
                    df_roi_area.to_csv(
                        tsv_roi_area_fn, sep="\t", na_rep='NaN', index=False)

                    # Violins (no averaging — all individual vertices pooled)
                    tsv_violins_fn = "{}/{}_{}_{}-gauss_violins.tsv".format(
                        tsv_dir, subject, fn_spec, output_tag)
                    print('Saving tsv: {}'.format(tsv_violins_fn))
                    df_violins.to_csv(
                        tsv_violins_fn, sep="\t", na_rep='NaN', index=False)

                    # Parameters median (computed from pooled violins df)
                    df_params_median = df_violins
                    colnames = [rsq2use, 'prf_size', 'prf_ecc']  # no pCM, no n

                    df_params_median_indiv = df_params_median.groupby(
                        ['roi', 'subject'])[[rsq2use]].apply(
                        lambda x: weighted_nan_median(
                            x[rsq2use], df_params_median.loc[x.index, rsq2use])
                    ).reset_index(name=f'{rsq2use}_weighted_median')

                    for colname in colnames[1:]:
                        df_params_median_indiv[f'{colname}_weighted_median'] = (
                            df_params_median.groupby(['roi', 'subject'])[[colname, rsq2use]]
                            .apply(lambda x: weighted_nan_median(
                                x[colname], df_params_median.loc[x.index, rsq2use]))
                            .reset_index()[0])

                    df_params_med_median = df_params_median_indiv.groupby(['roi'])[
                        [f'{c}_weighted_median' for c in colnames]].median()

                    df_params_median_ci = pd.DataFrame()
                    for colname in colnames:
                        df_params_median_ci[f'{colname}_ci_down'] = (
                            df_params_median_indiv.groupby(['roi']).apply(
                                lambda x: weighted_nan_percentile(
                                    x[f'{colname}_weighted_median'],
                                    x[f'{rsq2use}_weighted_median'], 2.5)))
                        df_params_median_ci[f'{colname}_ci_up'] = (
                            df_params_median_indiv.groupby(['roi']).apply(
                                lambda x: weighted_nan_percentile(
                                    x[f'{colname}_weighted_median'],
                                    x[f'{rsq2use}_weighted_median'], 97.5)))

                    df_params_median = pd.concat(
                        [df_params_med_median, df_params_median_ci], axis=1).reset_index()
                    tsv_params_median_fn = "{}/{}_{}_{}-gauss_params-median.tsv".format(
                        tsv_dir, subject, fn_spec, output_tag)
                    print('Saving tsv: {}'.format(tsv_params_median_fn))
                    df_params_median.to_csv(
                        tsv_params_median_fn, sep="\t", na_rep='NaN', index=False)

                    # Ecc. × Size (median across subjects per bin)
                    df_ecc_size = (df_ecc_size
                                   .groupby(['roi', 'num_bins'], sort=False)
                                   .median()
                                   .reset_index())
                    tsv_ecc_size_fn = "{}/{}_{}_{}-gauss_ecc-size.tsv".format(
                        tsv_dir, subject, fn_spec, output_tag)
                    print('Saving tsv: {}'.format(tsv_ecc_size_fn))
                    df_ecc_size.to_csv(
                        tsv_ecc_size_fn, sep="\t", na_rep='NaN', index=False)

                    # Polar angle (median across subjects per bin)
                    df_polar_angle = (df_polar_angle
                                      .groupby(['roi', 'hemi', 'num_bins'], sort=False)
                                      .median()
                                      .reset_index())
                    tsv_polar_angle_fn = "{}/{}_{}_{}-gauss_polar-angle.tsv".format(
                        tsv_dir, subject, fn_spec, output_tag)
                    print('Saving tsv: {}'.format(tsv_polar_angle_fn))
                    df_polar_angle.to_csv(
                        tsv_polar_angle_fn, sep="\t", na_rep='NaN', index=False)

                    # Contralaterality (median across subjects)
                    df_contralaterality = (df_contralaterality
                                           .groupby(['roi'], sort=False)
                                           .median()
                                           .reset_index())
                    tsv_contralaterality_fn = "{}/{}_{}_{}-gauss_contralaterality.tsv".format(
                        tsv_dir, subject, fn_spec, output_tag)
                    print('Saving tsv: {}'.format(tsv_contralaterality_fn))
                    df_contralaterality.to_csv(
                        tsv_contralaterality_fn, sep="\t", na_rep='NaN', index=False)

                    # Spatial distribution (voxel-wise median across subjects)
                    tsv_distribution_fn = "{}/{}_{}_prf-gauss_distribution.tsv".format(
                        tsv_dir, subject, fn_spec)
                    print('Saving tsv: {}'.format(tsv_distribution_fn))
                    median_mesh = np.median(mesh_group, axis=0)
                    df_distribution = pd.concat(
                        [others_columns.reset_index(drop=True),
                         pd.DataFrame(median_mesh)], axis=1)
                    df_distribution.to_csv(
                        tsv_distribution_fn, sep="\t", na_rep='NaN', index=False)

                    # Spatial distribution hot-zone barycentre
                    hemis = ['hemi-L', 'hemi-R', 'hemi-LR']
                    for j, hemi in enumerate(hemis):
                        hemi_values = ['hemi-L', 'hemi-R'] if hemi == 'hemi-LR' else [hemi]
                        df_distribution_hemi = df_distribution.loc[
                            df_distribution.hemi.isin(hemi_values)]
                        df_barycentre_hemi = make_prf_barycentre_df(
                            df_distribution_hemi, rois, distribution_max_ecc,
                            distribution_mesh_grain, hot_zone_percent=hot_zone_percent)

                        df_barycentre_hemi['hemi'] = [hemi] * len(df_barycentre_hemi)
                        if j == 0:
                            df_barycentre = df_barycentre_hemi
                        else:
                            df_barycentre = pd.concat([df_barycentre, df_barycentre_hemi])

                    tsv_barycentre_fn = "{}/{}_{}_{}-gauss_barycentre.tsv".format(
                        tsv_dir, subject, fn_spec, output_tag)
                    print('Saving tsv: {}'.format(tsv_barycentre_fn))
                    df_barycentre.to_csv(
                        tsv_barycentre_fn, sep="\t", na_rep='NaN', index=False)

# Define permission cmd
print('Changing files permissions in {}/{}'.format(main_dir, project_dir))
os.system("chmod -Rf 771 {}/{}".format(main_dir, project_dir))
os.system("chgrp -Rf {} {}/{}".format(group, main_dir, project_dir))