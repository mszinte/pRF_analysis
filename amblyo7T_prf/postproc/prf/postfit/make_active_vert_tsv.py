"""
-----------------------------------------------------------------------------------------
make_active_vert_tsv.py
-----------------------------------------------------------------------------------------
Goal of the script:
Make per-subject or per-group TSVs of significant vertex counts per ROI and eye condition
(AE-RE and FE-LE). For individual subjects, data from each eye condition TSV is loaded
without thresholding, and significant vertices at FDR 0.05 and 0.01 are counted per ROI.
For group subjects (group-patient, group-control), individual subject TSVs are loaded
and aggregated: median across subjects, with 2.5/97.5 percentile CI bounds.
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory (e.g. /home/mszinte/disks/meso_S/data)
sys.argv[2]: project name (e.g. amblyo7T_prf)
sys.argv[3]: subject name (e.g. sub-17, group-patient, group-control)
sys.argv[4]: server group (e.g. 327)
-----------------------------------------------------------------------------------------
Output(s):
- Per-subject/group, per-task TSV with active vertex counts:
  {subject}_{fn_spec_combined}_active-vert.tsv
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/pRF_analysis/amblyo7T_prf/postproc/prf/postfit/
2. run python command
python make_active_vert_tsv.py [main directory] [project name] [subject] [group]
-----------------------------------------------------------------------------------------
Example:
cd ~/projects/pRF_analysis/amblyo7T_prf/postproc/prf/postfit/
python make_active_vert_tsv.py /scratch/mszinte/data amblyo7T_prf sub-17 327
python make_active_vert_tsv.py /scratch/mszinte/data amblyo7T_prf group-patient 327
python make_active_vert_tsv.py /scratch/mszinte/data amblyo7T_prf group-control 327
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

                    df_active_vert_eyes = []

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

                        # Keep raw data — no thresholding
                        data_raw = pd.read_table(tsv_fn, sep="\t")

                        # Get eye_val
                        if 'RightEye' in prf_task_eye_name:
                            eye_val = 'right eye'
                        elif 'LeftEye' in prf_task_eye_name:
                            eye_val = 'left eye'

                        # Get eye_type
                        if subject_group == 'patient':
                            amblyopic_side = amblyopic_eye[0]
                            if eye_val == 'right eye' and amblyopic_side == 'R': eye_type = 'amblyopic eye'
                            elif eye_val == 'right eye' and amblyopic_side == 'L': eye_type = 'fellow eye'
                            elif eye_val == 'left eye' and amblyopic_side == 'L': eye_type = 'amblyopic eye'
                            elif eye_val == 'left eye' and amblyopic_side == 'R': eye_type = 'fellow eye'
                        elif subject_group == 'control':
                            eye_type = eye_val

                        # Map to AE-RE / FE-LE label
                        if eye_type in ['amblyopic eye', 'right eye']: eye_condition = 'AE-RE'
                        else: eye_condition = 'FE-LE'

                        # Count vertices per ROI
                        n_vert_tot_roi = (data_raw.groupby('roi', sort=False)
                                          .size().reset_index(name='n_vert_tot'))

                        n_vert_corr_5pt = (data_raw[data_raw['corr_pvalue_5pt'] < 0.05]
                                           .groupby('roi', sort=False)
                                           .size().reset_index(name='n_vert_corr_pvalue_5pt'))

                        n_vert_corr_1pt = (data_raw[data_raw['corr_pvalue_1pt'] < 0.01]
                                           .groupby('roi', sort=False)
                                           .size().reset_index(name='n_vert_corr_pvalue_1pt'))

                        df_eye = (n_vert_tot_roi
                                  .merge(n_vert_corr_5pt, on='roi', how='left')
                                  .merge(n_vert_corr_1pt, on='roi', how='left')
                                  .fillna(0))

                        df_eye['ratio_5pt'] = df_eye['n_vert_corr_pvalue_5pt'] / df_eye['n_vert_tot']
                        df_eye['ratio_1pt'] = df_eye['n_vert_corr_pvalue_1pt'] / df_eye['n_vert_tot']
                        df_eye['eye_condition'] = eye_condition
                        df_eye['subject'] = subject

                        df_active_vert_eyes.append(df_eye)

                    # Combine both eye conditions
                    df_active_vert = pd.concat(df_active_vert_eyes, ignore_index=True)

                    # Save TSV
                    tsv_active_vert_fn = '{}/{}_{}_active-vert.tsv'.format(
                        tsv_dir, subject, fn_spec_combined)
                    print('Saving tsv: {}'.format(tsv_active_vert_fn))
                    df_active_vert.to_csv(tsv_active_vert_fn, sep='\t', na_rep='NaN', index=False)

                # Group analysis
                # --------------
                else:
                    print(f'Group: {subject} - {avg_method} - {format_}')

                    # Load and concatenate individual subject TSVs
                    for i, subject_to_group in enumerate(subjects_to_group):
                        tsv_dir_indiv = '{}/{}/derivatives/pp_data/{}/{}/prf/tsv'.format(
                            main_dir, project_dir, subject_to_group, format_)
                        tsv_fn = '{}/{}_{}_active-vert.tsv'.format(
                            tsv_dir_indiv, subject_to_group, fn_spec_combined)
                        df_indiv = pd.read_table(tsv_fn, sep="\t")
                        if i == 0: df_all = df_indiv.copy()
                        else: df_all = pd.concat([df_all, df_indiv])

                    # Median across subjects per roi/eye_condition
                    median_cols = ['n_vert_tot', 'n_vert_corr_pvalue_5pt',
                                   'n_vert_corr_pvalue_1pt', 'ratio_5pt', 'ratio_1pt']
                    df_group = df_all.groupby(['roi', 'eye_condition'], sort=False)[median_cols].median().reset_index()

                    # CI across subjects (2.5/97.5 percentiles)
                    for col in ['ratio_5pt', 'ratio_1pt']:
                        df_group[f'{col}_ci_upper'] = (
                            df_all.groupby(['roi', 'eye_condition'], sort=False)[col]
                            .quantile(0.975).reset_index()[col])
                        df_group[f'{col}_ci_lower'] = (
                            df_all.groupby(['roi', 'eye_condition'], sort=False)[col]
                            .quantile(0.025).reset_index()[col])

                    # Save group TSV
                    tsv_dir_group = '{}/{}/derivatives/pp_data/{}/{}/prf/tsv'.format(
                        main_dir, project_dir, subject, format_)
                    os.makedirs(tsv_dir_group, exist_ok=True)
                    tsv_active_vert_fn = '{}/{}_{}_active-vert.tsv'.format(
                        tsv_dir_group, subject, fn_spec_combined)
                    print('Saving tsv: {}'.format(tsv_active_vert_fn))
                    df_group.to_csv(tsv_active_vert_fn, sep='\t', na_rep='NaN', index=False)

# Define permission cmd
# print('Changing files permissions in {}/{}'.format(main_dir, project_dir))
# os.system("chmod -Rf 771 {}/{}".format(main_dir, project_dir))
# os.system("chgrp -Rf {} {}/{}".format(group, main_dir, project_dir))
