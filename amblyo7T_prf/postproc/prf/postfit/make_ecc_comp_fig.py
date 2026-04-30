"""
-----------------------------------------------------------------------------------------
make_ecc_comp_fig.py
-----------------------------------------------------------------------------------------
Goal of the script:
Produce figures comparing pRF parameters across eye conditions per subject or group:

  Individual patient  : AE vs FE         (2 conditions, no stats)
  Individual control  : RE vs LE vs CTRL  (3 conditions, no stats)
  group-patient       : AE vs FE vs CTRL  (3 conditions + stats)
  group-control       : RE vs LE vs CTRL  (3 conditions + stats)

Each figure has:
  - 5 rows    : pRF R², pRF size, pRF eccentricity, pRF CM, n_vert
  - 3 columns : All data | Foveal (ecc <= 2.5 dva) | Peripheral (ecc > 2.5 dva)
  - x axis    : ROIs on a fake continuous axis with markers offset per eye condition
  - Legend    : top-right of first subplot, eye condition styles in neutral dark gray
  - Significance lines with asterisks/ns (group only)
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory (e.g. /scratch/mszinte/data)
sys.argv[2]: project name           (e.g. amblyo7T_prf)
sys.argv[3]: subject name           (e.g. sub-17, group-patient, group-control)
sys.argv[4]: server group           (e.g. 327)
-----------------------------------------------------------------------------------------
To run:
cd ~/projects/pRF_analysis/amblyo7T_prf/postproc/prf/postfit/
python make_ecc_comp_fig.py /scratch/mszinte/data amblyo7T_prf sub-17 327
python make_ecc_comp_fig.py /scratch/mszinte/data amblyo7T_prf group-patient 327
python make_ecc_comp_fig.py /scratch/mszinte/data amblyo7T_prf group-control 327
-----------------------------------------------------------------------------------------
Written by Martin Szinte (martin.szinte@gmail.com)
-----------------------------------------------------------------------------------------
"""
import warnings
warnings.filterwarnings("ignore")

import os
import sys
import pandas as pd
import numpy as np

sys.path.insert(0, "{}/../../../utils".format(os.getcwd()))
sys.path.append("{}/../../../../analysis_code/utils".format(os.getcwd()))
from plot_utils import ecc_comp_plot
from settings_utils import load_settings

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
params               = analysis_info['ecc_comp_params']

# Load participant info
participants_path  = os.path.join(main_dir, project_dir, 'participants.tsv')
participants_df    = pd.read_table(participants_path)
subject_groups_map = dict(zip(participants_df['participant_id'], participants_df['group']))

# -----------------------------------------------------------------------------------------
# Main loop
# -----------------------------------------------------------------------------------------
for avg_method in avg_methods:

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

                figure_info = analysis_info.copy()
                figure_info['rois'] = rois

                # ------------------------------------------------------------------
                # Determine run mode
                # ------------------------------------------------------------------
                if 'group' not in subject:
                    indiv_subjects = [subject]
                    run_group      = None
                elif subject == 'group-patient':
                    indiv_subjects = None
                    run_group      = 'patient'
                elif subject == 'group-control':
                    indiv_subjects = None
                    run_group      = 'control'

                # ------------------------------------------------------------------
                # INDIVIDUAL figures
                # ------------------------------------------------------------------
                for subj in (indiv_subjects if indiv_subjects is not None else []):
                    subj_group = subject_groups_map[subj]

                    tsv_dir = '{}/{}/derivatives/pp_data/{}/{}/prf/tsv'.format(
                        main_dir, project_dir, subj, format_)
                    fig_dir = '{}/{}/derivatives/pp_data/{}/{}/prf/figures'.format(
                        main_dir, project_dir, subj, format_)

                    if not os.path.isdir(tsv_dir):
                        print(f"[SKIP] {subj}: tsv_dir not found: {tsv_dir}")
                        continue
                    os.makedirs(fig_dir, exist_ok=True)

                    tsv_fn = '{}/{}_{}_ecc-comp.tsv'.format(tsv_dir, subj, fn_spec_combined)
                    if not os.path.isfile(tsv_fn):
                        print(f"[SKIP] {subj}: TSV not found: {tsv_fn}")
                        continue

                    df_subj = pd.read_table(tsv_fn, sep='\t')
                    # n_vert has no vertex-level CI — add NaN columns for uniform handling
                    df_subj['n_vert_median'] = df_subj['n_vert']
                    df_subj['n_vert_ci_lo']  = np.nan
                    df_subj['n_vert_ci_hi']  = np.nan

                    if subj_group == 'patient':
                        eye_conds        = ['AE-RE', 'FE-LE']
                        subj_group_label = 'patient'
                    else:
                        eye_conds        = ['AE-RE', 'FE-LE', 'CTRL']
                        subj_group_label = 'control'

                    figure_info['subject_group'] = subj_group_label
                    print(f'Plotting individual {subj} ({subj_group_label}) - {fn_spec_combined}')

                    fig = ecc_comp_plot(df=df_subj,
                                        df_stats=None,
                                        figure_info=figure_info,
                                        eye_conditions=eye_conds,
                                        show_stats=False)

                    fig_fn = '{}/{}_{}_ecc-comp.pdf'.format(fig_dir, subj, fn_spec_combined)
                    print(f'Saving: {fig_fn}')
                    fig.write_image(fig_fn)

                # ------------------------------------------------------------------
                # GROUP figure
                # ------------------------------------------------------------------
                if run_group is None:
                    continue

                tsv_dir_grp = '{}/{}/derivatives/pp_data/group-{}/{}/prf/tsv'.format(
                    main_dir, project_dir, run_group, format_)
                fig_dir_grp = '{}/{}/derivatives/pp_data/group-{}/{}/prf/figures'.format(
                    main_dir, project_dir, run_group, format_)
                os.makedirs(fig_dir_grp, exist_ok=True)

                tsv_grp  = '{}/group-{}_{}_ecc-comp.tsv'.format(
                    tsv_dir_grp, run_group, fn_spec_combined)
                tsv_stat = '{}/group-{}_{}_ecc-comp_stats.tsv'.format(
                    tsv_dir_grp, run_group, fn_spec_combined)

                if not os.path.isfile(tsv_grp):
                    print(f"[SKIP] group TSV not found: {tsv_grp}")
                    continue

                df_grp  = pd.read_table(tsv_grp,  sep='\t')
                df_stat = pd.read_table(tsv_stat, sep='\t')

                # Both patient and control group figures include CTRL
                eye_conds = ['AE-RE', 'FE-LE', 'CTRL']
                figure_info['subject_group'] = run_group
                print(f'Plotting group-{run_group} - {fn_spec_combined}')

                fig = ecc_comp_plot(df=df_grp,
                                    df_stats=df_stat,
                                    figure_info=figure_info,
                                    eye_conditions=eye_conds,
                                    show_stats=True)

                fig_fn = '{}/group-{}_{}_ecc-comp.pdf'.format(
                    fig_dir_grp, run_group, fn_spec_combined)
                print(f'Saving: {fig_fn}')
                fig.write_image(fig_fn)

print('Done.')