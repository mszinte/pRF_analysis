"""
-----------------------------------------------------------------------------------------
make_corr_tsv.py
-----------------------------------------------------------------------------------------
Goal of the script:
Make per-subject TSVs of pRF parameter correlations between eyes (AE/RE vs FE/LE).
For each task, format, ROI method, and parameter, data from the two eyes are merged,
quantile-binned by AE/RE, and weighted medians (+ 95% CI) are computed per bin.
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory (e.g. /home/mszinte/disks/meso_S/data)
sys.argv[2]: project name (e.g. amblyo7T_prf)
sys.argv[3]: subject name (e.g. sub-17)
sys.argv[4]: server group (e.g. 327)
-----------------------------------------------------------------------------------------
Output(s):
- Per-subject, per-task, per-parameter TSVs with binned correlations:
  {subject}_{fn_spec}_{corr_param}-corr.tsv
- Per-subject, per-task merged two-eye TSV:
  {subject}_{fn_spec_combined}_prf-css-deriv.tsv
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/pRF_analysis/analysis_code/postproc/prf/postfit/
2. run python command
python make_correlations_param_eyes_tsv.py [main directory] [project name] [subject] [group]
-----------------------------------------------------------------------------------------
Example:
cd ~/projects/pRF_analysis/analysis_code/postproc/prf/postfit/
python make_corr_tsv.py /home/mszinte/disks/meso_S/data amblyo7T_prf sub-17 327
-----------------------------------------------------------------------------------------
Written by Martin Szinte (martin.szinte@gmail.com)
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
sys.path.append("{}/../../../utils".format(os.getcwd()))
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
figure_settings_path = os.path.join(base_dir, project_dir, "figure-settings.yml")
settings = load_settings([settings_path, prf_settings_path, figure_settings_path])
analysis_info = settings[0]

formats = analysis_info['formats']
extensions = analysis_info['extensions']
rois_methods = analysis_info['rois_methods']
ecc_threshold = analysis_info['ecc_th']
size_threshold = analysis_info['size_th']
rsqr_threshold = analysis_info['rsqr_th']
amplitude_threshold = analysis_info['prf_amp_th']
stats_threshold = analysis_info['stats_th']
n_threshold = analysis_info['n_th']
preproc_prep = analysis_info['preproc_prep']
filtering = analysis_info['filtering']
normalization = analysis_info['normalization']
avg_methods = analysis_info['avg_methods']

# Parameters specific to eye correlation analysis
corr_bin_number = analysis_info['corr_bin_number']
corr_params = analysis_info['corr_params']
prf_tasks_eyes_names = analysis_info['prf_tasks_eyes_names']

# Load participant info
participants_path = os.path.join(main_dir, project_dir, 'participants.tsv')
participants_df = pd.read_table(participants_path)
amblyopic_eyes = dict(zip(participants_df['participant_id'], participants_df['amblyopic_eye']))
groups = dict(zip(participants_df['participant_id'], participants_df['group']))
amblyopic_eye = amblyopic_eyes[subject]
subject_group = groups[subject]

# Main loop
for avg_method in avg_methods:
    if 'loo' in avg_method: rsq2use = 'prf_loo_rsq'
    else: rsq2use = 'prf_rsq'

    for format_, extension in zip(formats, extensions):
        rois_methods_format = rois_methods[format_]
        for rois_method_format in rois_methods_format:

            if rois_method_format == 'rois-drawn':
                rois = analysis_info[rois_method_format]
            elif rois_method_format == 'rois-group-mmp':
                rois = list(analysis_info[rois_method_format].keys())

            for prf_task_eyes_names in prf_tasks_eyes_names:

                data_task_eye = {}
                for prf_task_eye_name in prf_task_eyes_names:

                    print(f'Loading: {prf_task_eye_name} - {avg_method} - {format_}')

                    tsv_dir = '{}/{}/derivatives/pp_data/{}/{}/prf/tsv'.format(
                        main_dir, project_dir, subject, format_)

                    if not os.path.isdir(tsv_dir):
                        print(f"[SKIP] tsv_dir not found for format={format_}: {tsv_dir}")
                        continue

                    fn_spec = "task-{}_{}_{}_{}_{}_{}" .format(
                        prf_task_eye_name, preproc_prep, filtering,
                        normalization, avg_method, rois_method_format)
                    tsv_fn = '{}/{}_{}_prf-css-deriv.tsv'.format(tsv_dir, subject, fn_spec)
                    data = pd.read_table(tsv_fn, sep="\t")

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

                    # Get eye_val
                    if 'RightEye' in prf_task_eye_name:
                        eye_val = 'right eye'
                    elif 'LeftEye' in prf_task_eye_name:
                        eye_val = 'left eye'
                    data['eye_val'] = eye_val

                    # Get eye_type
                    if subject_group == 'patient':
                        amblyopic_side = amblyopic_eye[0]
                        if eye_val == 'right eye' and amblyopic_side == 'R':
                            eye_type = 'amblyopic eye'
                        elif eye_val == 'right eye' and amblyopic_side == 'L':
                            eye_type = 'fellow eye'
                        elif eye_val == 'left eye' and amblyopic_side == 'L':
                            eye_type = 'amblyopic eye'
                        elif eye_val == 'left eye' and amblyopic_side == 'R':
                            eye_type = 'fellow eye'
                    elif subject_group == 'control':
                        eye_type = eye_val
                    data['eye_type'] = eye_type

                    # Merge two eyes
                    if len(data_task_eye) == 0:
                        data_task_eye['first'] = data
                        first_eye_type = eye_type
                    else:
                        second_eye_type = eye_type

                        if first_eye_type in ['amblyopic eye', 'right eye']:
                            suffix1, suffix2 = '_AE-RE', '_FE-LE'
                        else:
                            suffix1, suffix2 = '_FE-LE', '_AE-RE'
                            data_task_eye['first'], data = data, data_task_eye['first']

                        data = pd.merge(data_task_eye['first'], data,
                                        on=['num_vert', 'roi', 'roi_mmp', 'subject', 'hemi', 'trs'],
                                        suffixes=(suffix1, suffix2))

                # Save the merged two-eye TSV
                fn_spec_combined = "task-{}_{}_{}_{}_{}_{}" .format(
                    prf_task_eyes_names[0].replace('RightEye', '').replace('LeftEye', ''),
                    preproc_prep, filtering, normalization, avg_method, rois_method_format)

                output_fn = '{}/{}_{}_prf-css-deriv.tsv'.format(tsv_dir, subject, fn_spec_combined)
                data.to_csv(output_fn, sep='\t', index=False)
                print(f"Saving: {output_fn}")

                # PARAMETER CORRELATIONS AE/RE vs. FE/LE
                # ----------------------------------------
                for corr_param in corr_params:

                    for num_roi, roi in enumerate(rois):

                        df_roi = data.loc[(data.roi == roi)]

                        # Bin by AE-RE (X-axis / independent variable)
                        df_bin_AR = df_roi.groupby(pd.qcut(df_roi[f'{corr_param}_AE-RE'],
                                    q=corr_bin_number,
                                    duplicates='drop'))

                        df_bin = pd.DataFrame()
                        df_bin['roi'] = [roi] * len(df_bin_AR)
                        df_bin['num_bins'] = np.arange(len(df_bin_AR))
                        df_bin['n_vertex'] = df_bin_AR.size().values

                        # AE-RE (X-axis)
                        df_bin[f'{corr_param}_AE-RE_median'] = df_bin_AR.apply(
                            lambda x: weighted_nan_median(x[f'{corr_param}_AE-RE'].values,
                                                          x[f'{rsq2use}_AE-RE'].values)).values
                        df_bin[f'{corr_param}_AE-RE_ci_upper_bound'] = df_bin_AR.apply(
                            lambda x: weighted_nan_percentile(x[f'{corr_param}_AE-RE'].values,
                                                              x[f'{rsq2use}_AE-RE'].values, 97.5)).values
                        df_bin[f'{corr_param}_AE-RE_ci_lower_bound'] = df_bin_AR.apply(
                            lambda x: weighted_nan_percentile(x[f'{corr_param}_AE-RE'].values,
                                                              x[f'{rsq2use}_AE-RE'].values, 2.5)).values

                        # FE-LE (Y-axis) — from same AE-RE bins
                        df_bin[f'{corr_param}_FE-LE_median'] = df_bin_AR.apply(
                            lambda x: weighted_nan_median(x[f'{corr_param}_FE-LE'].values,
                                                          x[f'{rsq2use}_FE-LE'].values)).values
                        df_bin[f'{corr_param}_FE-LE_ci_upper_bound'] = df_bin_AR.apply(
                            lambda x: weighted_nan_percentile(x[f'{corr_param}_FE-LE'].values,
                                                              x[f'{rsq2use}_FE-LE'].values, 97.5)).values
                        df_bin[f'{corr_param}_FE-LE_ci_lower_bound'] = df_bin_AR.apply(
                            lambda x: weighted_nan_percentile(x[f'{corr_param}_FE-LE'].values,
                                                              x[f'{rsq2use}_FE-LE'].values, 2.5)).values

                        # Weight by AE-RE R²
                        df_bin[f'{rsq2use}_median'] = np.array(df_bin_AR[f'{rsq2use}_AE-RE'].median())

                        if num_roi == 0: df_bins = df_bin
                        else: df_bins = pd.concat([df_bins, df_bin])

                    tsv_fn = "{}/{}_{}_{}-corr.tsv".format(tsv_dir, subject, fn_spec_combined, corr_param)
                    print('Saving tsv: {}'.format(tsv_fn))
                    df_bins.to_csv(tsv_fn, sep="\t", na_rep='NaN', index=False)

# Define permission cmd
print('Changing files permissions in {}/{}'.format(main_dir, project_dir))
os.system("chmod -Rf 771 {}/{}".format(main_dir, project_dir))
os.system("chgrp -Rf {} {}/{}".format(group, main_dir, project_dir))