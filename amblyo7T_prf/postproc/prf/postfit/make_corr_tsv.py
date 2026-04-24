"""
-----------------------------------------------------------------------------------------
make_corr_tsv.py
-----------------------------------------------------------------------------------------
Goal of the script:
Make per-subject or per-group TSVs of pRF parameter correlations between eyes
(AE/RE vs FE/LE). For individual subjects, data from the two eyes are merged,
quantile-binned by AE/RE, and weighted medians (+ 95% CI) are computed per bin.
For group subjects (group-patient, group-control), individual subject TSVs are loaded
and aggregated: bin medians are averaged across subjects, and CI bounds reflect
across-subject variability (2.5/97.5 percentiles of per-subject medians).
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory (e.g. /home/mszinte/disks/meso_S/data)
sys.argv[2]: project name (e.g. amblyo7T_prf)
sys.argv[3]: subject name (e.g. sub-17, group-patient, group-control)
sys.argv[4]: server group (e.g. 327)
-----------------------------------------------------------------------------------------
Output(s):
- Per-subject/group, per-task, per-parameter TSVs with binned correlations:
  {subject}_{fn_spec}_{corr_param}-corr.tsv
- Per-subject, per-task merged two-eye TSV:
  {subject}_{fn_spec_combined}_prf-css-deriv.tsv
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/pRF_analysis/amblyo7T_prf/postproc/prf/postfit/
2. run python command
python make_corr_tsv.py [main directory] [project name] [subject] [group]
-----------------------------------------------------------------------------------------
Example:
cd ~/projects/pRF_analysis/amblyo7T_prf/postproc/prf/postfit/
python make_corr_tsv.py /scratch/mszinte/data amblyo7T_prf sub-17 327
python make_corr_tsv.py /scratch/mszinte/data amblyo7T_prf group-patient 327
python make_corr_tsv.py /scratch/mszinte/data amblyo7T_prf group-control 327
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
figure_settings_path = os.path.join(base_dir, project_dir, "figure-settings.yml")
settings = load_settings([settings_path, prf_settings_path, figure_settings_path])
analysis_info = settings[0]

formats = analysis_info['formats']
extensions = analysis_info['extensions']
rois_methods = analysis_info['rois_methods']


ecc_threshold = analysis_info['corr_ecc_th']
size_threshold = analysis_info['corr_size_th']
rsqr_threshold = analysis_info['corr_rsqr_th']
amplitude_threshold = analysis_info['corr_prf_amp_th']
stats_threshold = analysis_info['corr_stats_th']
n_threshold = analysis_info['corr_n_th']

preproc_prep = analysis_info['preproc_prep']
filtering = analysis_info['filtering']
normalization = analysis_info['normalization']
avg_methods = analysis_info['avg_methods']

# Parameters specific to eye correlation analysis
corr_bin_number = analysis_info['corr_bin_number']
corr_bin_eye = analysis_info['corr_bin_eye']
corr_other_eye = 'FE-LE' if corr_bin_eye == 'AE-RE' else 'AE-RE'
corr_params = analysis_info['corr_params']
prf_tasks_eyes_names = analysis_info['prf_tasks_eyes_names']

# Load participant info
participants_path = os.path.join(main_dir, project_dir, 'participants.tsv')
participants_df = pd.read_table(participants_path)
amblyopic_eyes = dict(zip(participants_df['participant_id'], participants_df['amblyopic_eye']))
subject_groups = dict(zip(participants_df['participant_id'], participants_df['group']))

# Determine if individual or group run
if 'group' not in subject:
    amblyopic_eye = amblyopic_eyes[subject]
    subject_group = subject_groups[subject]
else:
    # Identify which group and load the corresponding subject list
    if 'patient' in subject:
        subjects_to_group = analysis_info['group_patient']
    elif 'control' in subject:
        subjects_to_group = analysis_info['group_control']

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

                fn_spec_combined = "task-{}_{}_{}_{}_{}_{}" .format(
                    prf_task_eyes_names[0].replace('RightEye', '').replace('LeftEye', ''),
                    preproc_prep, filtering, normalization, avg_method, rois_method_format)

                # Individual subject analysis
                # ---------------------------
                if 'group' not in subject:

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
                        data = data.reset_index(drop=True)

                        # Store original task name for traceability
                        data['task_eye_name'] = prf_task_eye_name

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
                        print(f"[DEBUG2] {prf_task_eye_name} | eye_type assigned={eye_type} | unique in col={data['eye_type'].unique()}")

                        # Merge two eyes
                        if len(data_task_eye) == 0:
                            data_task_eye['first'] = data.copy()
                            first_eye_type = eye_type
                        else:
                            second_eye_type = eye_type

                            if first_eye_type in ['amblyopic eye', 'right eye']:
                                suffix1, suffix2 = '_AE-RE', '_FE-LE'
                            else:
                                data_task_eye['first'], data = data, data_task_eye['first']
                                suffix1, suffix2 = '_AE-RE', '_FE-LE'

                            # Reassign eye_val and eye_type after potential swap
                            # using .loc to ensure proper in-place assignment
                            data_task_eye['first'].loc[:, 'eye_val'] = 'amblyopic eye' if subject_group == 'patient' else 'right eye'
                            data_task_eye['first'].loc[:, 'eye_type'] = 'amblyopic eye' if subject_group == 'patient' else 'right eye'
                            data.loc[:, 'eye_val'] = 'fellow eye' if subject_group == 'patient' else 'left eye'
                            data.loc[:, 'eye_type'] = 'fellow eye' if subject_group == 'patient' else 'left eye'

                            print(f"[DEBUG] subject={subject}, first_eye_type={first_eye_type}, suffix1={suffix1}, task_eye_name_first={data_task_eye['first']['task_eye_name'].iloc[0]}, task_eye_name_second={data['task_eye_name'].iloc[0]}")

                            data = pd.merge(data_task_eye['first'], data,
                                            on=['num_vert', 'roi', 'roi_mmp', 'subject', 'hemi', 'trs'],
                                            suffixes=(suffix1, suffix2))

                    # Save the merged two-eye TSV
                    output_fn = '{}/{}_{}_prf-css-deriv.tsv'.format(tsv_dir, subject, fn_spec_combined)
                    data.to_csv(output_fn, sep='\t', index=False)
                    print(f"Saving: {output_fn}")

                    # Sanity check: verify correct TSV files ended up in AE-RE and FE-LE columns
                    ae_tasks = data['task_eye_name_AE-RE'].unique()
                    fe_tasks = data['task_eye_name_FE-LE'].unique()
                    if subject_group == 'patient':
                        amblyopic_side = amblyopic_eye[0]
                        ae_ok = all(('RightEye' in t if amblyopic_side == 'R' else 'LeftEye' in t) for t in ae_tasks)
                        fe_ok = all(('LeftEye' in t if amblyopic_side == 'R' else 'RightEye' in t) for t in fe_tasks)
                    else:
                        ae_ok = all('RightEye' in t for t in ae_tasks)
                        fe_ok = all('LeftEye' in t for t in fe_tasks)
                    status = '✓' if ae_ok and fe_ok else '✗ MISMATCH'
                    print(f"[CHECK] {subject} | AE-RE from={ae_tasks} | FE-LE from={fe_tasks} | {status}")

                    # PARAMETER CORRELATIONS AE/RE vs. FE/LE
                    # ----------------------------------------
                    for corr_param in corr_params:

                        for num_roi, roi in enumerate(rois):

                            df_roi = data.loc[(data.roi == roi)]

                            # Bin by corr_bin_eye (X-axis / independent variable)
                            df_bin_x = df_roi.groupby(pd.qcut(df_roi[f'{corr_param}_{corr_bin_eye}'],
                                        q=corr_bin_number,
                                        duplicates='drop'))

                            df_bin = pd.DataFrame()
                            df_bin['roi'] = [roi] * len(df_bin_x)
                            df_bin['num_bins'] = np.arange(len(df_bin_x))
                            df_bin['n_vertex'] = df_bin_x.size().values

                            # X-axis eye (binning eye = corr_bin_eye)
                            df_bin[f'{corr_param}_{corr_bin_eye}_median'] = df_bin_x.apply(
                                lambda x: weighted_nan_median(x[f'{corr_param}_{corr_bin_eye}'].values,
                                                              x[f'{rsq2use}_{corr_bin_eye}'].values)).values
                            df_bin[f'{corr_param}_{corr_bin_eye}_ci_upper_bound'] = df_bin_x.apply(
                                lambda x: weighted_nan_percentile(x[f'{corr_param}_{corr_bin_eye}'].values,
                                                                  x[f'{rsq2use}_{corr_bin_eye}'].values, 97.5)).values
                            df_bin[f'{corr_param}_{corr_bin_eye}_ci_lower_bound'] = df_bin_x.apply(
                                lambda x: weighted_nan_percentile(x[f'{corr_param}_{corr_bin_eye}'].values,
                                                                  x[f'{rsq2use}_{corr_bin_eye}'].values, 2.5)).values

                            # Y-axis eye (other eye = corr_other_eye) — from same corr_bin_eye bins
                            df_bin[f'{corr_param}_{corr_other_eye}_median'] = df_bin_x.apply(
                                lambda x: weighted_nan_median(x[f'{corr_param}_{corr_other_eye}'].values,
                                                              x[f'{rsq2use}_{corr_other_eye}'].values)).values
                            df_bin[f'{corr_param}_{corr_other_eye}_ci_upper_bound'] = df_bin_x.apply(
                                lambda x: weighted_nan_percentile(x[f'{corr_param}_{corr_other_eye}'].values,
                                                                  x[f'{rsq2use}_{corr_other_eye}'].values, 97.5)).values
                            df_bin[f'{corr_param}_{corr_other_eye}_ci_lower_bound'] = df_bin_x.apply(
                                lambda x: weighted_nan_percentile(x[f'{corr_param}_{corr_other_eye}'].values,
                                                                  x[f'{rsq2use}_{corr_other_eye}'].values, 2.5)).values

                            # Weight by corr_bin_eye R²
                            df_bin[f'{rsq2use}_median'] = np.array(df_bin_x[f'{rsq2use}_{corr_bin_eye}'].median())

                            if num_roi == 0: df_bins = df_bin
                            else: df_bins = pd.concat([df_bins, df_bin])

                        tsv_fn = "{}/{}_{}_{}-corr.tsv".format(tsv_dir, subject, fn_spec_combined, corr_param)
                        print('Saving tsv: {}'.format(tsv_fn))
                        df_bins.to_csv(tsv_fn, sep="\t", na_rep='NaN', index=False)

                # Group analysis
                # --------------
                else:
                    print(f'Group: {subject} - {avg_method} - {format_}')

                    for corr_param in corr_params:

                        # Load and concatenate individual subject TSVs
                        for i, subject_to_group in enumerate(subjects_to_group):
                            tsv_dir_indiv = '{}/{}/derivatives/pp_data/{}/{}/prf/tsv'.format(
                                main_dir, project_dir, subject_to_group, format_)
                            tsv_fn = "{}/{}_{}_{}-corr.tsv".format(
                                tsv_dir_indiv, subject_to_group, fn_spec_combined, corr_param)
                            df_indiv = pd.read_table(tsv_fn, sep="\t")
                            if i == 0: df_all = df_indiv.copy()
                            else: df_all = pd.concat([df_all, df_indiv])

                        # Median across subjects per roi/bin
                        median_cols = [f'{corr_param}_{corr_bin_eye}_median',
                                       f'{corr_param}_{corr_other_eye}_median',
                                       f'{rsq2use}_median',
                                       'n_vertex']
                        df_group = df_all.groupby(['roi', 'num_bins'], sort=False)[median_cols].median().reset_index()

                        # CI across subjects (2.5/97.5 percentiles of per-subject medians)
                        for eye in [corr_bin_eye, corr_other_eye]:
                            df_group[f'{corr_param}_{eye}_ci_upper_bound'] = (
                                df_all.groupby(['roi', 'num_bins'], sort=False)
                                [f'{corr_param}_{eye}_median']
                                .quantile(0.975)
                                .reset_index()[f'{corr_param}_{eye}_median'])
                            df_group[f'{corr_param}_{eye}_ci_lower_bound'] = (
                                df_all.groupby(['roi', 'num_bins'], sort=False)
                                [f'{corr_param}_{eye}_median']
                                .quantile(0.025)
                                .reset_index()[f'{corr_param}_{eye}_median'])

                        # Save group TSV
                        tsv_dir_group = '{}/{}/derivatives/pp_data/{}/{}/prf/tsv'.format(
                            main_dir, project_dir, subject, format_)
                        os.makedirs(tsv_dir_group, exist_ok=True)
                        tsv_fn = "{}/{}_{}_{}-corr.tsv".format(
                            tsv_dir_group, subject, fn_spec_combined, corr_param)
                        print('Saving tsv: {}'.format(tsv_fn))
                        df_group.to_csv(tsv_fn, sep="\t", na_rep='NaN', index=False)

# Define permission cmd
# print('Changing files permissions in {}/{}'.format(main_dir, project_dir))
# os.system("chmod -Rf 771 {}/{}".format(main_dir, project_dir))
# os.system("chgrp -Rf {} {}/{}".format(group, main_dir, project_dir))