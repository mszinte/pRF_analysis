"""
-----------------------------------------------------------------------------------------
make_ecc_size_pcm_tsv.py
-----------------------------------------------------------------------------------------
Goal of the script:
Make per-subject or per-group TSVs of eccentricity vs pRF size and eccentricity vs pCM,
separately for each eye condition (AE-RE and FE-LE). Data is loaded from per-eye
prf-css-deriv.tsv files, thresholded, and binned by eccentricity. The eye condition
labeling follows the same logic as make_active_vert_tsv.py: for left-amblyopic patients
the LeftEye file is labeled AE-RE and for right-amblyopic the RightEye file is labeled
AE-RE. Controls always have RightEye=AE-RE (RE), LeftEye=FE-LE (LE).

For individual subjects:
- Each eye condition is processed independently (no vertex merging)
- Eccentricity bins are computed using equal-width bins from 0 to ecc_size_max
- Weighted median and 25/75 CI for pRF size and pCM per bin

For group subjects (group-patient, group-control):
- Individual subject TSVs are loaded and median across subjects per roi/eye_condition/bin
- 2.5/97.5 CI bounds across subjects
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory (e.g. /home/mszinte/disks/meso_S/data)
sys.argv[2]: project name (e.g. amblyo7T_prf)
sys.argv[3]: subject name (e.g. sub-17, group-patient, group-control)
sys.argv[4]: server group (e.g. 327)
-----------------------------------------------------------------------------------------
Output(s):
- Per-subject/group, per-task TSV with eccentricity-binned size and pCM per eye condition:
  {subject}_{fn_spec_combined}_ecc-size-pcm.tsv
  columns: roi, eye_condition, num_bins, prf_ecc_bins,
           prf_size_bins_median, prf_size_bins_ci_upper_bound, prf_size_bins_ci_lower_bound,
           prf_pcm_bins_median, prf_pcm_bins_ci_upper_bound, prf_pcm_bins_ci_lower_bound,
           prf_rsq_median
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/pRF_analysis/amblyo7T_prf/postproc/prf/postfit/
2. run python command
python make_ecc_size_pcm_tsv.py [main directory] [project name] [subject] [group]
-----------------------------------------------------------------------------------------
Example:
cd ~/projects/pRF_analysis/amblyo7T_prf/postproc/prf/postfit/
python make_ecc_size_pcm_tsv.py /scratch/mszinte/data amblyo7T_prf sub-17 327
python make_ecc_size_pcm_tsv.py /scratch/mszinte/data amblyo7T_prf group-patient 327
python make_ecc_size_pcm_tsv.py /scratch/mszinte/data amblyo7T_prf group-control 327
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
preproc_prep = analysis_info['preproc_prep']
filtering = analysis_info['filtering']
normalization = analysis_info['normalization']
avg_methods = analysis_info['avg_methods']
prf_tasks_eyes_names = analysis_info['prf_tasks_eyes_names']

ecc_threshold = analysis_info['ecc_th']
size_threshold = analysis_info['size_th']
rsqr_threshold = analysis_info['rsqr_th']
amplitude_threshold = analysis_info['prf_amp_th']
stats_threshold = analysis_info['stats_th']
n_threshold = analysis_info['n_th']

ecc_size_num_bins = analysis_info['ecc_size_num_bins']
ecc_pcm_num_bins = analysis_info['ecc_pcm_num_bins']
ecc_size_max = analysis_info['ecc_size_max']
ecc_pcm_max = analysis_info['ecc_pcm_max']
ecc_bin_power = analysis_info['ecc_bin_power']

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

                    print(f'Processing make_ecc_size_pcm_tsv.py for: {subject}')

                    tsv_dir = '{}/{}/derivatives/pp_data/{}/{}/prf/tsv'.format(
                        main_dir, project_dir, subject, format_)

                    if not os.path.isdir(tsv_dir):
                        print(f"[SKIP] tsv_dir not found for format={format_}: {tsv_dir}")
                        continue

                    df_eyes = []

                    for prf_task_eye_name in prf_task_eyes_names:

                        print(f'Loading: {prf_task_eye_name} - {avg_method} - {format_}')

                        fn_spec = "task-{}_{}_{}_{}_{}_{}" .format(
                            prf_task_eye_name, preproc_prep, filtering,
                            normalization, avg_method, rois_method_format)
                        tsv_fn = '{}/{}_{}_prf-css-deriv.tsv'.format(tsv_dir, subject, fn_spec)

                        data = pd.read_table(tsv_fn, sep="\t")

                        # Threshold data
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

                        # Get eye_val and eye_condition
                        if 'RightEye' in prf_task_eye_name:
                            eye_val = 'right eye'
                        elif 'LeftEye' in prf_task_eye_name:
                            eye_val = 'left eye'

                        if subject_group == 'patient':
                            amblyopic_side = amblyopic_eye[0]
                            if eye_val == 'right eye' and amblyopic_side == 'R': eye_type = 'amblyopic eye'
                            elif eye_val == 'right eye' and amblyopic_side == 'L': eye_type = 'fellow eye'
                            elif eye_val == 'left eye' and amblyopic_side == 'L': eye_type = 'amblyopic eye'
                            elif eye_val == 'left eye' and amblyopic_side == 'R': eye_type = 'fellow eye'
                        elif subject_group == 'control':
                            eye_type = eye_val

                        if eye_type in ['amblyopic eye', 'right eye']: eye_condition = 'AE-RE'
                        else: eye_condition = 'FE-LE'

                        # Ecc bins for size: power function spacing, finer near fovea
                        ecc_size_bins = ecc_size_max[0] * (np.linspace(0, 1, ecc_size_num_bins + 1) ** ecc_bin_power)
                        # Ecc bins for pCM: power function spacing, finer near fovea
                        ecc_pcm_bins = ecc_pcm_max[0] * (np.linspace(0, 1, ecc_pcm_num_bins + 1) ** ecc_bin_power)

                        for num_roi, roi in enumerate(rois):
                            df_roi = data.loc[(data.roi == roi)]

                            # Size bins
                            df_size_bins = df_roi.groupby(pd.cut(df_roi['prf_ecc'], bins=ecc_size_bins))
                            # pCM bins
                            df_pcm_bins = df_roi.groupby(pd.cut(df_roi['prf_ecc'], bins=ecc_pcm_bins))

                            df_bin = pd.DataFrame()
                            df_bin['roi'] = [roi] * ecc_size_num_bins
                            df_bin['eye_condition'] = eye_condition
                            df_bin['num_bins'] = np.arange(ecc_size_num_bins)
                            df_bin['subject'] = subject

                            df_bin['prf_ecc_bins'] = df_size_bins.apply(
                                lambda x: weighted_nan_median(x['prf_ecc'].values, x[rsq2use].values)).values
                            df_bin['prf_size_bins_median'] = df_size_bins.apply(
                                lambda x: weighted_nan_median(x['prf_size'].values, x[rsq2use].values)).values
                            df_bin['prf_size_bins_ci_upper_bound'] = df_size_bins.apply(
                                lambda x: weighted_nan_percentile(x['prf_size'].values, x[rsq2use].values, 75)).values
                            df_bin['prf_size_bins_ci_lower_bound'] = df_size_bins.apply(
                                lambda x: weighted_nan_percentile(x['prf_size'].values, x[rsq2use].values, 25)).values
                            df_bin[f'{rsq2use}_bins_median'] = df_size_bins.apply(
                                lambda x: weighted_nan_median(x[rsq2use].values, x[rsq2use].values)).values
                            df_bin[f'{rsq2use}_bins_ci_upper_bound'] = df_size_bins.apply(
                                lambda x: weighted_nan_percentile(x[rsq2use].values, x[rsq2use].values, 75)).values
                            df_bin[f'{rsq2use}_bins_ci_lower_bound'] = df_size_bins.apply(
                                lambda x: weighted_nan_percentile(x[rsq2use].values, x[rsq2use].values, 25)).values
                            df_bin['n_vert_bins'] = df_size_bins.size().values
                            df_bin['prf_pcm_bins_median'] = df_pcm_bins.apply(
                                lambda x: weighted_nan_median(x['pcm_median'].values, x[rsq2use].values)).values
                            df_bin['prf_pcm_bins_ci_upper_bound'] = df_pcm_bins.apply(
                                lambda x: weighted_nan_percentile(x['pcm_median'].values, x[rsq2use].values, 75)).values
                            df_bin['prf_pcm_bins_ci_lower_bound'] = df_pcm_bins.apply(
                                lambda x: weighted_nan_percentile(x['pcm_median'].values, x[rsq2use].values, 25)).values

                            if num_roi == 0: df_bins = df_bin
                            else: df_bins = pd.concat([df_bins, df_bin])

                        df_eyes.append(df_bins)

                    # Combine both eye conditions
                    df_ecc_size_pcm = pd.concat(df_eyes, ignore_index=True)

                    # Save TSV
                    tsv_fn = '{}/{}_{}_ecc-size-pcm.tsv'.format(tsv_dir, subject, fn_spec_combined)
                    print('Saving tsv: {}'.format(tsv_fn))
                    df_ecc_size_pcm.to_csv(tsv_fn, sep='\t', na_rep='NaN', index=False)

                # Group analysis
                # --------------
                else:
                    print(f'Group: {subject} - {avg_method} - {format_}')

                    for i, subject_to_group in enumerate(subjects_to_group):
                        tsv_dir_indiv = '{}/{}/derivatives/pp_data/{}/{}/prf/tsv'.format(
                            main_dir, project_dir, subject_to_group, format_)
                        tsv_fn = '{}/{}_{}_ecc-size-pcm.tsv'.format(
                            tsv_dir_indiv, subject_to_group, fn_spec_combined)
                        df_indiv = pd.read_table(tsv_fn, sep="\t")
                        if i == 0: df_all = df_indiv.copy()
                        else: df_all = pd.concat([df_all, df_indiv])

                    # Median across subjects per roi/eye_condition/bin
                    median_cols = ['prf_ecc_bins',
                                   'prf_size_bins_median', 'prf_size_bins_ci_upper_bound', 'prf_size_bins_ci_lower_bound',
                                   'prf_pcm_bins_median', 'prf_pcm_bins_ci_upper_bound', 'prf_pcm_bins_ci_lower_bound',
                                   f'{rsq2use}_bins_median', f'{rsq2use}_bins_ci_upper_bound', f'{rsq2use}_bins_ci_lower_bound',
                                   'n_vert_bins']
                    df_group = df_all.groupby(['roi', 'eye_condition', 'num_bins'],
                                              sort=False)[median_cols].median().reset_index()

                    # CI across subjects (2.5/97.5 percentiles)
                    for col in ['prf_size_bins_median', 'prf_pcm_bins_median',
                                f'{rsq2use}_bins_median', 'n_vert_bins']:
                        df_group[f'{col}_ci_lower'] = (
                            df_all.groupby(['roi', 'eye_condition', 'num_bins'], sort=False)[col]
                            .quantile(0.025).reset_index()[col])
                        df_group[f'{col}_ci_upper'] = (
                            df_all.groupby(['roi', 'eye_condition', 'num_bins'], sort=False)[col]
                            .quantile(0.975).reset_index()[col])

                    # Save group TSV
                    tsv_dir_group = '{}/{}/derivatives/pp_data/{}/{}/prf/tsv'.format(
                        main_dir, project_dir, subject, format_)
                    os.makedirs(tsv_dir_group, exist_ok=True)
                    tsv_fn = '{}/{}_{}_ecc-size-pcm.tsv'.format(tsv_dir_group, subject, fn_spec_combined)
                    print('Saving tsv: {}'.format(tsv_fn))
                    df_group.to_csv(tsv_fn, sep='\t', na_rep='NaN', index=False)

# Define permission cmd
# print('Changing files permissions in {}/{}'.format(main_dir, project_dir))
# os.system("chmod -Rf 771 {}/{}".format(main_dir, project_dir))
# os.system("chgrp -Rf {} {}/{}".format(group, main_dir, project_dir))