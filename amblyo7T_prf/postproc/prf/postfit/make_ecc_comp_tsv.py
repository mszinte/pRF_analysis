"""
-----------------------------------------------------------------------------------------
make_ecc_comp_tsv.py
-----------------------------------------------------------------------------------------
Goal of the script:
For each subject, load the thresholded prf-css-deriv.tsv (same thresholds as
make_ecc_size_pcm_tsv.py), compute R²-weighted median of pRF R², pRF size, pRF
eccentricity, pRF CM (pcm_median) and n_vert per roi × eye_condition × ecc_category.

ecc_category:
  - 'all'        : all thresholded vertices (no extra eccentricity restriction)
  - 'foveal'     : prf_ecc <= foveal_ecc_bound
  - 'peripheral' : prf_ecc >  foveal_ecc_bound

eye_condition follows the same labeling logic as make_ecc_size_pcm_tsv.py:
  patients  → AE-RE / FE-LE
  controls  → AE-RE (RE) / FE-LE (LE) / CTRL (median of RE and LE)

Then aggregate across subjects and run statistics:
  group-patient: AE vs FE   → paired Wilcoxon
                 AE vs CTRL → Mann-Whitney U (CTRL from control subjects)
                 FE vs CTRL → Mann-Whitney U
  group-control: RE vs LE   → paired Wilcoxon
  All FDR-corrected with Benjamini-Hochberg (padjust='fdr_bh'), two-sided.

Output TSVs saved under:
  per-subject   : {main_dir}/{project_dir}/derivatives/pp_data/{subject}/{format_}/prf/tsv/
  group-patient : {main_dir}/{project_dir}/derivatives/pp_data/group-patient/{format_}/prf/tsv/
  group-control : {main_dir}/{project_dir}/derivatives/pp_data/group-control/{format_}/prf/tsv/
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory (e.g. /scratch/mszinte/data)
sys.argv[2]: project name           (e.g. amblyo7T_prf)
sys.argv[3]: subject name           (e.g. sub-17, group-patient, group-control)
sys.argv[4]: server group           (e.g. 327)
-----------------------------------------------------------------------------------------
To run:
cd ~/projects/pRF_analysis/amblyo7T_prf/postproc/prf/postfit/
python make_ecc_comp_tsv.py /scratch/mszinte/data amblyo7T_prf sub-17 327
python make_ecc_comp_tsv.py /scratch/mszinte/data amblyo7T_prf group-patient 327
python make_ecc_comp_tsv.py /scratch/mszinte/data amblyo7T_prf group-control 327
-----------------------------------------------------------------------------------------
Written by Martin Szinte (martin.szinte@gmail.com)
-----------------------------------------------------------------------------------------
"""
import warnings
warnings.filterwarnings("ignore")

import os
import sys
import numpy as np
import pandas as pd
import pingouin as pg
import ipdb
deb = ipdb.set_trace

sys.path.append("{}/../../../../analysis_code/utils".format(os.getcwd()))
from settings_utils import load_settings
from maths_utils import weighted_nan_median, weighted_nan_percentile

# -----------------------------------------------------------------------------------------
# Inputs
# -----------------------------------------------------------------------------------------
main_dir    = sys.argv[1]
project_dir = sys.argv[2]
subject     = sys.argv[3]
group       = sys.argv[4]

# -----------------------------------------------------------------------------------------
# Load settings
# -----------------------------------------------------------------------------------------
base_dir             = os.path.abspath(os.path.join(os.getcwd(), "../../../../"))
settings_path        = os.path.join(base_dir, project_dir, "settings.yml")
prf_settings_path    = os.path.join(base_dir, project_dir, "prf-analysis.yml")
figure_settings_path = os.path.join(base_dir, project_dir, "figure-settings.yml")
settings             = load_settings([settings_path, prf_settings_path, figure_settings_path])
analysis_info        = settings[0]

formats              = analysis_info['formats']
extensions           = analysis_info['extensions']
rois_methods         = analysis_info['rois_methods']
preproc_prep         = analysis_info['preproc_prep']
filtering            = analysis_info['filtering']
normalization        = analysis_info['normalization']
avg_methods          = analysis_info['avg_methods']
prf_tasks_eyes_names = analysis_info['prf_tasks_eyes_names']

ecc_threshold        = analysis_info['ecc_th']
size_threshold       = analysis_info['size_th']
rsqr_threshold       = analysis_info['rsqr_th']
amplitude_threshold  = analysis_info['prf_amp_th']
stats_threshold      = analysis_info['stats_th']
n_threshold          = analysis_info['n_th']
foveal_bound         = analysis_info['foveal_ecc_bound']
params               = analysis_info['ecc_comp_params']

# Load participant info
participants_path = os.path.join(main_dir, project_dir, 'participants.tsv')
participants_df   = pd.read_table(participants_path)
amblyopic_eyes    = dict(zip(participants_df['participant_id'], participants_df['amblyopic_eye']))
subject_groups    = dict(zip(participants_df['participant_id'], participants_df['group']))

patients  = analysis_info['group_patient']
controls  = analysis_info['group_control']

# Determine which subjects to process for individual TSVs
# and whether to run group aggregation + stats
if 'group' not in subject:
    indiv_subjects = [subject]
    run_group      = None
elif subject == 'group-patient':
    indiv_subjects = None
    run_group      = 'patient'
elif subject == 'group-control':
    indiv_subjects = None
    run_group      = 'control'


ecc_categories = ['all', 'foveal', 'peripheral']
alternatives   = ['two-sided']


# -----------------------------------------------------------------------------------------
# Helper: asterisk string from p-value
# -----------------------------------------------------------------------------------------
def pval_to_stars(p):
    if p < 0.001: return '***'
    elif p < 0.01: return '**'
    elif p < 0.05: return '*'
    else: return 'ns'


# -----------------------------------------------------------------------------------------
# Main loop over avg_method / format / rois_method / task
# -----------------------------------------------------------------------------------------
for avg_method in avg_methods:
    rsq2use = 'prf_loo_rsq' if 'loo' in avg_method else 'prf_rsq'

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

                # ------------------------------------------------------------------
                # STEP 1: Per-subject TSV
                # ------------------------------------------------------------------
                # For group-combined we still need all subjects for stats,
                # but skip saving individual TSVs (they were already saved)
                all_subjects_for_stats = patients + controls
                df_all_subjects = []

                for subj in (indiv_subjects if indiv_subjects is not None else []):
                    subject_group  = subject_groups[subj]
                    amblyopic_eye  = amblyopic_eyes[subj]

                    tsv_dir = '{}/{}/derivatives/pp_data/{}/{}/prf/tsv'.format(
                        main_dir, project_dir, subj, format_)
                    if not os.path.isdir(tsv_dir):
                        print(f"[SKIP] {subj}: tsv_dir not found: {tsv_dir}")
                        continue

                    df_eyes = []

                    for prf_task_eye_name in prf_task_eyes_names:
                        fn_spec = "task-{}_{}_{}_{}_{}_{}" .format(
                            prf_task_eye_name, preproc_prep, filtering,
                            normalization, avg_method, rois_method_format)
                        tsv_fn = '{}/{}_{}_prf-css-deriv.tsv'.format(tsv_dir, subj, fn_spec)
                        if not os.path.isfile(tsv_fn):
                            print(f"[SKIP] {subj}: file not found: {tsv_fn}")
                            continue

                        data = pd.read_table(tsv_fn, sep="\t")

                        # Apply same thresholds as make_ecc_size_pcm_tsv.py
                        if stats_threshold == 0.05: stats_col = 'corr_pvalue_5pt'
                        elif stats_threshold == 0.01: stats_col = 'corr_pvalue_1pt'
                        data.loc[
                            (data.amplitude < amplitude_threshold[0]) |
                            (data.prf_ecc   < ecc_threshold[0])  | (data.prf_ecc  > ecc_threshold[1]) |
                            (data.prf_size  < size_threshold[0]) | (data.prf_size > size_threshold[1]) |
                            (data.prf_n     < n_threshold[0])    | (data.prf_n    > n_threshold[1]) |
                            (data[rsq2use]  < rsqr_threshold) |
                            (data[stats_col] > stats_threshold)] = np.nan
                        data = data.dropna().reset_index(drop=True)

                        # Eye labeling
                        if 'RightEye' in prf_task_eye_name:
                            eye_val = 'right eye'
                        else:
                            eye_val = 'left eye'

                        if subject_group == 'patient':
                            amblyopic_side = amblyopic_eye[0]
                            if   eye_val == 'right eye' and amblyopic_side == 'R': eye_type = 'amblyopic eye'
                            elif eye_val == 'right eye' and amblyopic_side == 'L': eye_type = 'fellow eye'
                            elif eye_val == 'left eye'  and amblyopic_side == 'L': eye_type = 'amblyopic eye'
                            elif eye_val == 'left eye'  and amblyopic_side == 'R': eye_type = 'fellow eye'
                        else:
                            eye_type = eye_val

                        eye_condition = 'AE-RE' if eye_type in ['amblyopic eye', 'right eye'] else 'FE-LE'

                        # Eccentricity category
                        data.loc[data.prf_ecc <= foveal_bound, 'ecc_category'] = 'foveal'
                        data.loc[data.prf_ecc >  foveal_bound, 'ecc_category'] = 'peripheral'

                        # Compute R²-weighted median + 2.5/97.5 CI per roi × ecc_category
                        rows_list = []
                        for roi in rois:
                            df_roi = data.loc[data.roi == roi]
                            for ecc_cat in ecc_categories:
                                if ecc_cat == 'all':
                                    df_cat = df_roi
                                elif ecc_cat == 'foveal':
                                    df_cat = df_roi.loc[df_roi.prf_ecc <= foveal_bound]
                                else:
                                    df_cat = df_roi.loc[df_roi.prf_ecc > foveal_bound]

                                w   = df_cat[rsq2use].values if len(df_cat) > 0 else np.array([])
                                has = len(df_cat) > 0

                                row = {
                                    'subject':          subj,
                                    'subject_group':    subject_group,
                                    'roi':              roi,
                                    'eye_condition':    eye_condition,
                                    'ecc_category':     ecc_cat,
                                    'n_vert':           len(df_cat),
                                    'prf_rsq':          weighted_nan_median(df_cat[rsq2use].values, w)          if has else np.nan,
                                    'prf_rsq_ci_lo':    weighted_nan_percentile(df_cat[rsq2use].values, w, 2.5) if has else np.nan,
                                    'prf_rsq_ci_hi':    weighted_nan_percentile(df_cat[rsq2use].values, w, 97.5) if has else np.nan,
                                    'prf_size':         weighted_nan_median(df_cat['prf_size'].values, w)          if has else np.nan,
                                    'prf_size_ci_lo':   weighted_nan_percentile(df_cat['prf_size'].values, w, 2.5) if has else np.nan,
                                    'prf_size_ci_hi':   weighted_nan_percentile(df_cat['prf_size'].values, w, 97.5) if has else np.nan,
                                    'prf_ecc':          weighted_nan_median(df_cat['prf_ecc'].values, w)          if has else np.nan,
                                    'prf_ecc_ci_lo':    weighted_nan_percentile(df_cat['prf_ecc'].values, w, 2.5) if has else np.nan,
                                    'prf_ecc_ci_hi':    weighted_nan_percentile(df_cat['prf_ecc'].values, w, 97.5) if has else np.nan,
                                    'pcm_median':       weighted_nan_median(df_cat['pcm_median'].values, w)          if has else np.nan,
                                    'pcm_median_ci_lo': weighted_nan_percentile(df_cat['pcm_median'].values, w, 2.5) if has else np.nan,
                                    'pcm_median_ci_hi': weighted_nan_percentile(df_cat['pcm_median'].values, w, 97.5) if has else np.nan,
                                }
                                rows_list.append(row)

                        df_eyes.append(pd.DataFrame(rows_list))

                    if not df_eyes:
                        continue

                    df_subj = pd.concat(df_eyes, ignore_index=True)

                    # Add CTRL rows for controls (median of AE-RE and FE-LE per roi/ecc_cat)
                    if subject_group == 'control':
                        ctrl_rows = []
                        for roi in rois:
                            for ecc_cat in ecc_categories:
                                df_rc = df_subj.loc[
                                    (df_subj.roi == roi) & (df_subj.ecc_category == ecc_cat)]
                                ctrl_row = {
                                    'subject':       subj,
                                    'subject_group': subject_group,
                                    'roi':           roi,
                                    'eye_condition': 'CTRL',
                                    'ecc_category':  ecc_cat,
                                    'n_vert':        df_rc['n_vert'].median(),
                                }
                                for param in [p for p in params if p != 'n_vert']:
                                    ctrl_row[param]              = df_rc[param].median()
                                    ctrl_row[f'{param}_ci_lo']   = df_rc[f'{param}_ci_lo'].median()
                                    ctrl_row[f'{param}_ci_hi']   = df_rc[f'{param}_ci_hi'].median()
                                ctrl_rows.append(ctrl_row)
                        df_subj = pd.concat([df_subj, pd.DataFrame(ctrl_rows)], ignore_index=True)

                    # Save per-subject TSV
                    tsv_fn_out = '{}/{}_{}_ecc-comp.tsv'.format(tsv_dir, subj, fn_spec_combined)
                    print(f'Saving per-subject TSV ({subj}): {tsv_fn_out}')
                    df_subj.to_csv(tsv_fn_out, sep='\t', na_rep='NaN', index=False)
                    df_all_subjects.append(df_subj)

                # ------------------------------------------------------------------
                # STEP 2: Group aggregation + stats (only when run_group is set)
                # ------------------------------------------------------------------
                if run_group is None:
                    continue

                # For group-patient, group-control, group-combined:
                # reload individual TSVs from disk
                if not df_all_subjects:
                    subjects_to_load = {
                        'patient': all_subjects_for_stats,  # needs controls for CTRL
                        'control': controls,
                    }[run_group]
                    for subj in subjects_to_load:
                        td = '{}/{}/derivatives/pp_data/{}/{}/prf/tsv'.format(
                            main_dir, project_dir, subj, format_)
                        tf = '{}/{}_{}_ecc-comp.tsv'.format(td, subj, fn_spec_combined)
                        if os.path.isfile(tf):
                            df_all_subjects.append(pd.read_table(tf, sep='\t'))
                        else:
                            print(f"[SKIP] individual TSV not found: {tf}")

                if not df_all_subjects:
                    print(f"[SKIP] No subject data for group stats: {fn_spec_combined}")
                    continue

                df_all = pd.concat(df_all_subjects, ignore_index=True)

                out_dir = '{}/{}/derivatives/pp_data/group-{}/{}/prf/tsv'.format(
                    main_dir, project_dir, run_group, format_)
                os.makedirs(out_dir, exist_ok=True)

                def aggregate_group(df_in, group_cols):
                    """Compute median and 2.5/97.5 CI across subjects using uniform weights."""
                    rows = []
                    for keys, df_grp in df_in.groupby(group_cols):
                        if not isinstance(keys, tuple):
                            keys = (keys,)
                        row = dict(zip(group_cols, keys))
                        for p in params:
                            vals = df_grp[p].values.astype(float)
                            w    = np.ones(len(vals))
                            row[f'{p}_median'] = weighted_nan_median(vals, w)
                            row[f'{p}_ci_lo']  = weighted_nan_percentile(vals, w, 2.5)
                            row[f'{p}_ci_hi']  = weighted_nan_percentile(vals, w, 97.5)
                        rows.append(row)
                    return pd.DataFrame(rows)

                # -- Fig 1: patients AE-RE vs FE-LE + CTRL from controls --
                if run_group == 'patient':
                    df_pat       = df_all.loc[df_all.subject_group == 'patient']
                    df_ctrl      = df_all.loc[df_all.subject_group == 'control']
                    df_ctrl_ctrl = df_ctrl.loc[df_ctrl.eye_condition == 'CTRL'].copy()
                    df_pat_full  = pd.concat([
                        df_pat.loc[df_pat.eye_condition.isin(['AE-RE', 'FE-LE'])],
                        df_ctrl_ctrl], ignore_index=True)
                    df_pat_grp   = aggregate_group(df_pat_full, ['roi', 'ecc_category', 'eye_condition'])
                    fn_pat = '{}/group-patient_{}_ecc-comp.tsv'.format(out_dir, fn_spec_combined)
                    print(f'Saving patient group TSV: {fn_pat}')
                    df_pat_grp.to_csv(fn_pat, sep='\t', na_rep='NaN', index=False)

                    df_pat_indiv = df_pat_full.copy()
                    fn_pat_indiv = '{}/group-patient_{}_ecc-comp_indiv.tsv'.format(out_dir, fn_spec_combined)
                    df_pat_indiv.to_csv(fn_pat_indiv, sep='\t', na_rep='NaN', index=False)

                # -- Fig 2: controls, AE-RE (RE) vs FE-LE (LE) --
                elif run_group == 'control':
                    df_ctrl      = df_all.loc[df_all.subject_group == 'control']
                    df_ctrl_grp  = aggregate_group(
                        df_ctrl.loc[df_ctrl.eye_condition.isin(['AE-RE', 'FE-LE'])],
                        ['roi', 'ecc_category', 'eye_condition'])
                    fn_ctrl = '{}/group-control_{}_ecc-comp.tsv'.format(out_dir, fn_spec_combined)
                    print(f'Saving control group TSV: {fn_ctrl}')
                    df_ctrl_grp.to_csv(fn_ctrl, sep='\t', na_rep='NaN', index=False)

                    df_ctrl_indiv = df_ctrl.loc[df_ctrl.eye_condition.isin(['AE-RE', 'FE-LE'])]
                    fn_ctrl_indiv = '{}/group-control_{}_ecc-comp_indiv.tsv'.format(out_dir, fn_spec_combined)
                    df_ctrl_indiv.to_csv(fn_ctrl_indiv, sep='\t', na_rep='NaN', index=False)



                # ------------------------------------------------------------------
                # STEP 3: Statistics
                # ------------------------------------------------------------------

                def run_stats(df_indiv, comparisons, paired_map, group_col, out_fn):
                    """
                    Run pairwise Wilcoxon/Mann-Whitney tests per roi/ecc_category/param.

                    Parameters
                    ----------
                    df_indiv   : individual-level dataframe
                    comparisons: list of (A, B) tuples of eye_condition values to compare
                    paired_map : dict mapping (A,B) -> bool (True=paired Wilcoxon, False=Mann-Whitney)
                    group_col  : column used as grouping variable (e.g. 'eye_condition')
                    out_fn     : output filename
                    """
                    rows_stats = []
                    for param in params:
                        # Collect all raw p-values for this parameter (family)
                        param_rows = []
                        for roi in rois:
                            for ecc_cat in ecc_categories:
                                df_rc = df_indiv.loc[
                                    (df_indiv.roi == roi) & (df_indiv.ecc_category == ecc_cat)]
                                for (cond_a, cond_b) in comparisons:
                                    paired = paired_map[(cond_a, cond_b)]
                                    df_ab  = df_rc.loc[df_rc[group_col].isin([cond_a, cond_b])].copy()

                                    if paired:
                                        # Matched subjects — pivot and drop NaN pairs
                                        df_piv = df_ab.pivot_table(
                                            index='subject', columns=group_col, values=param)
                                        df_piv = df_piv.dropna()
                                        if len(df_piv) < 2:
                                            p_unc = np.nan
                                        else:
                                            df_long = df_piv.reset_index().melt(
                                                id_vars='subject', var_name=group_col, value_name=param)
                                            try:
                                                res   = pg.pairwise_tests(
                                                    data=df_long, dv=param,
                                                    within=group_col, subject='subject',
                                                    parametric=False,
                                                    alternative='two-sided')
                                                p_unc = res['p_unc'].values[0]
                                            except Exception as e:
                                                print(f"  Stat error (paired) {roi}/{ecc_cat}/{param}/{cond_a}v{cond_b}: {e}")
                                                p_unc = np.nan
                                    else:
                                        # Unpaired Mann-Whitney
                                        data_a = df_ab.loc[df_ab[group_col] == cond_a, param].dropna()
                                        data_b = df_ab.loc[df_ab[group_col] == cond_b, param].dropna()
                                        if len(data_a) < 2 or len(data_b) < 2:
                                            p_unc = np.nan
                                        else:
                                            try:
                                                res   = pg.pairwise_tests(
                                                    data=df_ab, dv=param,
                                                    between=group_col,
                                                    parametric=False,
                                                    alternative='two-sided')
                                                p_unc = res['p_unc'].values[0]
                                            except Exception as e:
                                                print(f"  Stat error (unpaired) {roi}/{ecc_cat}/{param}/{cond_a}v{cond_b}: {e}")
                                                p_unc = np.nan

                                    param_rows.append({
                                        'roi':          roi,
                                        'ecc_category': ecc_cat,
                                        'param':        param,
                                        'cond_A':       cond_a,
                                        'cond_B':       cond_b,
                                        'paired':       paired,
                                        'p_unc':        p_unc,
                                    })

                        # FDR-BH correction per parameter family
                        df_param = pd.DataFrame(param_rows)
                        all_pvals  = df_param['p_unc'].values
                        valid_mask = ~np.isnan(all_pvals)
                        p_corr     = np.full(len(all_pvals), np.nan)
                        if valid_mask.sum() > 0:
                            _, p_corr_valid = pg.multicomp(all_pvals[valid_mask], method='fdr_bh')
                            p_corr[valid_mask] = p_corr_valid
                        df_param['p_corr'] = p_corr
                        df_param['stars']  = df_param['p_corr'].apply(
                            lambda p: pval_to_stars(p) if not np.isnan(p) else 'ns')
                        rows_stats.append(df_param)

                    df_stats = pd.concat(rows_stats, ignore_index=True)
                    df_stats['roi'] = pd.Categorical(df_stats['roi'], categories=rois, ordered=True)
                    df_stats = df_stats.sort_values(['roi', 'ecc_category', 'param'])
                    print(f'Saving stats TSV: {out_fn}')
                    df_stats.to_csv(out_fn, sep='\t', na_rep='NaN', index=False)
                    return df_stats

                # Fig 1 stats: AE vs FE (paired), AE vs CTRL and FE vs CTRL (unpaired)
                if run_group == 'patient':
                    run_stats(
                        df_indiv    = df_pat_indiv,
                        comparisons = [('AE-RE', 'FE-LE'), ('AE-RE', 'CTRL'), ('FE-LE', 'CTRL')],
                        paired_map  = {
                            ('AE-RE', 'FE-LE'): True,
                            ('AE-RE', 'CTRL'):  False,
                            ('FE-LE', 'CTRL'):  False,
                        },
                        group_col   = 'eye_condition',
                        out_fn      = '{}/group-patient_{}_ecc-comp_stats.tsv'.format(out_dir, fn_spec_combined))

                # Fig 2 stats: RE vs LE (paired, within controls)
                elif run_group == 'control':
                    run_stats(
                        df_indiv    = df_ctrl_indiv,
                        comparisons = [('AE-RE', 'FE-LE')],
                        paired_map  = {('AE-RE', 'FE-LE'): True},
                        group_col   = 'eye_condition',
                        out_fn      = '{}/group-control_{}_ecc-comp_stats.tsv'.format(out_dir, fn_spec_combined))



print('Done.')