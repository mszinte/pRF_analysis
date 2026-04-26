"""
-----------------------------------------------------------------------------------------
make_corr_tsv.py
-----------------------------------------------------------------------------------------
Goal of the script:
Make per-subject or per-group TSVs of pRF parameter correlations between eyes
(FE/LE x-axis vs AE/RE y-axis).

For individual subjects:
- Data from RightEye and LeftEye runs are loaded separately, thresholded, and merged
  on common vertices (inner join on num_vert). The AE/FE labeling accounts for which
  physical eye is amblyopic (from participants.tsv): for left-amblyopic patients the
  LeftEye file is assigned _AE-RE and RightEye is assigned _FE-LE (swap), and vice
  versa for right-amblyopic. Controls always have RightEye=_AE-RE, LeftEye=_FE-LE.
- For each parameter and ROI, a weighted 2D KDE is computed on a fixed
  corr_grid_size x corr_grid_size grid (x = FE/LE, y = AE/RE), weighted by the
  mean R² of both eyes.
- A weighted Deming (orthogonal) regression is computed per ROI, weighted by mean R²,
  along with Pearson r², p-value, n_vertex, and mean R² weight.

For group subjects (group-patient, group-control):
- KDE density grids are averaged cell by cell across subjects (each subject contributes
  equally regardless of vertex count).
- Regression stats (slope, intercept, r2_pearson) are averaged across subjects with
  2.5/97.5 percentile CI bounds.
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory (e.g. /home/mszinte/disks/meso_S/data)
sys.argv[2]: project name (e.g. amblyo7T_prf)
sys.argv[3]: subject name (e.g. sub-17, group-patient, group-control)
sys.argv[4]: server group (e.g. 327)
-----------------------------------------------------------------------------------------
Output(s):
- Per-subject/group, per-task, per-parameter KDE density TSVs:
  {subject}_{fn_spec_combined}_{corr_param}-corr.tsv
  columns: roi, x_grid, y_grid, density
- Per-subject/group, per-task, per-parameter regression stats TSVs:
  {subject}_{fn_spec_combined}_{corr_param}-regression.tsv
  columns: roi, slope, intercept, r2_pearson, p_value, n_vertex, mean_r2_weight
  (+ ci bounds for group)
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
from scipy.stats import gaussian_kde

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
corr_params = analysis_info['corr_params']
corr_grid_size = analysis_info['corr_grid_size']
corr_param_ranges = {
    'prf_rsq':    analysis_info['prf_rsq_param_range'],
    'prf_x':      analysis_info['prf_x_param_range'],
    'prf_y':      analysis_info['prf_y_param_range'],
    'prf_size':   analysis_info['prf_size_param_range'],
    'prf_ecc':    analysis_info['prf_ecc_param_range'],
    'pcm_median': analysis_info['pcm_median_param_range']
}
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
                            data_task_eye['first'].loc[:, 'eye_val'] = 'amblyopic eye' if subject_group == 'patient' else 'right eye'
                            data_task_eye['first'].loc[:, 'eye_type'] = 'amblyopic eye' if subject_group == 'patient' else 'right eye'
                            data.loc[:, 'eye_val'] = 'fellow eye' if subject_group == 'patient' else 'left eye'
                            data.loc[:, 'eye_type'] = 'fellow eye' if subject_group == 'patient' else 'left eye'

                            data = pd.merge(data_task_eye['first'], data,
                                            on=['num_vert', 'roi', 'roi_mmp', 'subject', 'hemi', 'trs'],
                                            suffixes=(suffix1, suffix2))

                    # Save the merged two-eye TSV
                    output_fn = '{}/{}_{}_prf-css-deriv.tsv'.format(tsv_dir, subject, fn_spec_combined)
                    data.to_csv(output_fn, sep='\t', index=False)
                    print(f"Saving: {output_fn}")

                    # PARAMETER CORRELATIONS — WEIGHTED 2D KDE + REGRESSION STATS
                    # -------------------------------------------------------------
                    from scipy.stats import pearsonr
                    for corr_param in corr_params:

                        param_range = corr_param_ranges[corr_param]
                        x_grid = np.linspace(param_range[0], param_range[1], corr_grid_size)
                        y_grid = np.linspace(param_range[0], param_range[1], corr_grid_size)
                        xx, yy = np.meshgrid(x_grid, y_grid)
                        grid_points = np.vstack([xx.ravel(), yy.ravel()])

                        df_reg_rows = []
                        for num_roi, roi in enumerate(rois):

                            df_roi = data.loc[(data.roi == roi)].copy()

                            x_vals = df_roi[f'{corr_param}_FE-LE'].values  # x = FE (reference)
                            y_vals = df_roi[f'{corr_param}_AE-RE'].values  # y = AE (outcome)
                            r2_ae = df_roi[f'{rsq2use}_AE-RE'].values
                            r2_fe = df_roi[f'{rsq2use}_FE-LE'].values
                            r2_combined = (r2_ae + r2_fe) / 2

                            # Remove NaNs
                            mask = ~(np.isnan(x_vals) | np.isnan(y_vals) | np.isnan(r2_combined))
                            x_vals = x_vals[mask]
                            y_vals = y_vals[mask]
                            r2_combined = r2_combined[mask]

                            n_vertex = len(x_vals)
                            mean_r2_weight = np.nanmean(r2_combined)

                            # Normalize weights
                            r2_norm = r2_combined / r2_combined.sum()

                            # Compute weighted 2D KDE
                            kde = gaussian_kde(np.vstack([x_vals, y_vals]), weights=r2_norm)
                            density = kde(grid_points).reshape(corr_grid_size, corr_grid_size)

                            # Save KDE as long-format dataframe
                            df_kde_roi = pd.DataFrame({
                                'roi': roi,
                                'x_grid': xx.ravel(),
                                'y_grid': yy.ravel(),
                                'density': density.ravel()
                            })

                            if num_roi == 0: df_kde = df_kde_roi
                            else: df_kde = pd.concat([df_kde, df_kde_roi])

                            # Deming regression + stats
                            if n_vertex > 2:
                                weights = r2_norm
                                x_mean = np.sum(weights * x_vals)
                                y_mean = np.sum(weights * y_vals)
                                sxx = np.sum(weights * (x_vals - x_mean) ** 2)
                                syy = np.sum(weights * (y_vals - y_mean) ** 2)
                                sxy = np.sum(weights * (x_vals - x_mean) * (y_vals - y_mean))
                                slope = (syy - sxx + np.sqrt((syy - sxx) ** 2 + 4 * sxy ** 2)) / (2 * sxy)
                                intercept = y_mean - slope * x_mean
                                r_val, p_val = pearsonr(x_vals, y_vals)
                                r2_pearson = r_val ** 2
                            else:
                                slope = np.nan
                                intercept = np.nan
                                r2_pearson = np.nan
                                p_val = np.nan

                            df_reg_rows.append({
                                'roi': roi,
                                'slope': slope,
                                'intercept': intercept,
                                'r2_pearson': r2_pearson,
                                'p_value': p_val,
                                'n_vertex': n_vertex,
                                'mean_r2_weight': mean_r2_weight
                            })

                        # Save KDE TSV
                        tsv_fn = "{}/{}_{}_{}-corr.tsv".format(tsv_dir, subject, fn_spec_combined, corr_param)
                        print('Saving tsv: {}'.format(tsv_fn))
                        df_kde.to_csv(tsv_fn, sep="\t", na_rep='NaN', index=False)

                        # Save regression TSV
                        df_reg = pd.DataFrame(df_reg_rows)
                        reg_tsv_fn = "{}/{}_{}_{}-regression.tsv".format(tsv_dir, subject, fn_spec_combined, corr_param)
                        print('Saving tsv: {}'.format(reg_tsv_fn))
                        df_reg.to_csv(reg_tsv_fn, sep="\t", na_rep='NaN', index=False)

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

                        # Average density across subjects per roi/grid cell
                        df_group = df_all.groupby(['roi', 'x_grid', 'y_grid'], sort=False)[['density']].mean().reset_index()

                        # Save group KDE TSV
                        tsv_dir_group = '{}/{}/derivatives/pp_data/{}/{}/prf/tsv'.format(
                            main_dir, project_dir, subject, format_)
                        os.makedirs(tsv_dir_group, exist_ok=True)
                        tsv_fn = "{}/{}_{}_{}-corr.tsv".format(
                            tsv_dir_group, subject, fn_spec_combined, corr_param)
                        print('Saving tsv: {}'.format(tsv_fn))
                        df_group.to_csv(tsv_fn, sep="\t", na_rep='NaN', index=False)

                        # Group regression — load and average per-subject regression TSVs
                        for i, subject_to_group in enumerate(subjects_to_group):
                            tsv_dir_indiv = '{}/{}/derivatives/pp_data/{}/{}/prf/tsv'.format(
                                main_dir, project_dir, subject_to_group, format_)
                            reg_tsv_fn = "{}/{}_{}_{}-regression.tsv".format(
                                tsv_dir_indiv, subject_to_group, fn_spec_combined, corr_param)
                            df_reg_indiv = pd.read_table(reg_tsv_fn, sep="\t")
                            if i == 0: df_reg_all = df_reg_indiv.copy()
                            else: df_reg_all = pd.concat([df_reg_all, df_reg_indiv])

                        # Mean slope and intercept across subjects
                        df_reg_group = df_reg_all.groupby(['roi'], sort=False)[['slope', 'intercept', 'r2_pearson', 'n_vertex', 'mean_r2_weight']].mean().reset_index()

                        # CI across subjects (2.5/97.5 percentiles)
                        for col in ['slope', 'intercept', 'r2_pearson']:
                            df_reg_group[f'{col}_ci_lower'] = df_reg_all.groupby(['roi'], sort=False)[col].quantile(0.025).reset_index()[col]
                            df_reg_group[f'{col}_ci_upper'] = df_reg_all.groupby(['roi'], sort=False)[col].quantile(0.975).reset_index()[col]

                        reg_tsv_fn = "{}/{}_{}_{}-regression.tsv".format(
                            tsv_dir_group, subject, fn_spec_combined, corr_param)
                        print('Saving tsv: {}'.format(reg_tsv_fn))
                        df_reg_group.to_csv(reg_tsv_fn, sep="\t", na_rep='NaN', index=False)

# Define permission cmd
# print('Changing files permissions in {}/{}'.format(main_dir, project_dir))
# os.system("chmod -Rf 771 {}/{}".format(main_dir, project_dir))
# os.system("chgrp -Rf {} {}/{}".format(group, main_dir, project_dir))