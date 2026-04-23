"""
-----------------------------------------------------------------------------------------
make_corr_fig.py
-----------------------------------------------------------------------------------------
Goal of the script:
Make per-subject or per-group figures of pRF parameter correlations between eyes
(AE/RE vs FE/LE). One figure per task with rows = ROIs, columns = parameters.
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory (e.g. /home/mszinte/disks/meso_S/data)
sys.argv[2]: project name (e.g. amblyo7T_prf)
sys.argv[3]: subject name (e.g. sub-17, group-patient, group-control)
sys.argv[4]: server group (e.g. 327)
-----------------------------------------------------------------------------------------
Output(s):
- Per-subject/group, per-task correlation figure (PDF):
  {subject}_{fn_spec_combined}_prf-css_corr.pdf
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/pRF_analysis/amblyo7T_prf/postproc/prf/postfit/
2. run python command
python make_corr_fig.py [main directory] [project name] [subject] [group]
-----------------------------------------------------------------------------------------
Example:
cd ~/projects/pRF_analysis/amblyo7T_prf/postproc/prf/postfit/
python make_corr_fig.py /scratch/mszinte/data amblyo7T_prf sub-17 327
python make_corr_fig.py /scratch/mszinte/data amblyo7T_prf group-patient 327
python make_corr_fig.py /scratch/mszinte/data amblyo7T_prf group-control 327
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

# Personal import
sys.path.insert(0, "{}/../../../utils".format(os.getcwd()))
sys.path.append("{}/../../../../analysis_code/utils".format(os.getcwd()))
from plot_utils import corr_plot
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
corr_params = analysis_info['corr_params']
corr_bin_eye = analysis_info['corr_bin_eye']

# Load participant info
import pandas as pd
participants_path = os.path.join(main_dir, project_dir, 'participants.tsv')
participants_df = pd.read_table(participants_path)
subject_groups = dict(zip(participants_df['participant_id'], participants_df['group']))

# Determine subject group label for axis labeling
if 'group' not in subject:
    subject_group = subject_groups[subject]
else:
    subject_group = 'patient' if 'patient' in subject else 'control'

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

                print(f'{fn_spec_combined} - {subject}')

                # Define dirs
                tsv_dir = '{}/{}/derivatives/pp_data/{}/{}/prf/tsv'.format(
                    main_dir, project_dir, subject, format_)
                fig_dir = '{}/{}/derivatives/pp_data/{}/{}/prf/figures'.format(
                    main_dir, project_dir, subject, format_)

                if not os.path.isdir(tsv_dir):
                    print(f"[SKIP] tsv_dir not found: {tsv_dir}")
                    continue

                os.makedirs(fig_dir, exist_ok=True)

                # Build figure_info
                figure_info = analysis_info.copy()
                figure_info['rois'] = rois
                figure_info['corr_params'] = corr_params
                figure_info['subject_group'] = subject_group
                figure_info['corr_bin_eye'] = corr_bin_eye

                # Make and save figure
                fig = corr_plot(tsv_dir=tsv_dir,
                                subject=subject,
                                fn_spec_combined=fn_spec_combined,
                                figure_info=figure_info,
                                rsq2use=rsq2use)

                fig_fn = "{}/{}_{}_prf-css_corr.pdf".format(fig_dir, subject, fn_spec_combined)
                print('Saving pdf: {}'.format(fig_fn))
                fig.write_image(fig_fn)

# Define permission cmd
# print('Changing files permissions in {}/{}'.format(main_dir, project_dir))
# os.system("chmod -Rf 771 {}/{}".format(main_dir, project_dir))
# os.system("chgrp -Rf {} {}/{}".format(group, main_dir, project_dir))