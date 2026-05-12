"""
-----------------------------------------------------------------------------------------
tsv_roi_group.py
-----------------------------------------------------------------------------------------
Goal of the script:
Make PMF figure specific TSV - group level, across 3 fit conditions
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject name (e.g. group)
sys.argv[4]: server group (e.g. 327)
-----------------------------------------------------------------------------------------
Output(s):
PMF analysis tsv (params-median) per condition, grouped across subjects
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/pRF_analysis/RetinoMaps/postproc/pmf/postfit/
2. run python command
python tsv_roi_group.py [main directory] [project name] group [group]
-----------------------------------------------------------------------------------------
Exemple:
python tsv_roi_group.py /scratch/mszinte/data RetinoMaps group 327
-----------------------------------------------------------------------------------------
"""
# Stop warnings
import warnings
warnings.filterwarnings("ignore")

# General imports
import os
import sys
import numpy as np
import pandas as pd

# Personal import
script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.abspath(os.path.join(script_dir, "../../../../analysis_code/utils")))
from settings_utils import load_settings
from maths_utils import weighted_nan_median, weighted_nan_percentile

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]   # should be 'group'
group = sys.argv[4]

# Load settings
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../../"))
settings_path = os.path.join(base_dir, project_dir, "settings.yml")
prf_settings_path = os.path.join(base_dir, project_dir, "pmf-settings.yml")
figure_settings_path = os.path.join(base_dir, project_dir, "figure-settings.yml")
settings = load_settings([settings_path, prf_settings_path, figure_settings_path])
analysis_info = settings[0]

formats = analysis_info['formats']
extensions = analysis_info['extensions']
subjects_to_group = analysis_info['subjects']
preproc_prep = analysis_info['preproc_prep']
filtering = analysis_info['filtering']
normalization = analysis_info['normalization']
avg_methods = ['concat']

# The 3 PMF fit conditions
pmf_task_names = ['pmf-bold-gauss', 'pmf-residuals-gauss', 'prf-control-gauss']
rois_methods = analysis_info['rois_methods']

# Main loop
for avg_method in avg_methods:
    if 'loo' in avg_method:
        rsq2use = 'prf_loo_rsq'
    else:
        rsq2use = 'prf_rsq'

    for format_, extension in zip(formats, extensions):

        # define list of rois methods for each format
        rois_methods_format = rois_methods[format_]
        for rois_method_format in rois_methods_format:

            for pmf_task_name in pmf_task_names:

                print(f'{pmf_task_name} - {avg_method} - {format_} - {rois_method_format}')

                fn_spec = "task-{}_{}_{}_{}_{}_{}".format(
                    'SacLoc', preproc_prep, filtering, normalization, avg_method, rois_method_format)

                for i, subject_to_group in enumerate(subjects_to_group):
                    tsv_dir = '{}/{}/derivatives/pp_data/{}/{}/pmf/tsv'.format(
                        main_dir, project_dir, subject_to_group, format_)

                    # Parameters median
                    # ------------------
                    tsv_params_median_fn = "{}/{}_{}_{}_params-median.tsv".format(
                        tsv_dir, subject_to_group, fn_spec, pmf_task_name)

                    if not os.path.isfile(tsv_params_median_fn):
                        print(f"[SKIP] File not found: {tsv_params_median_fn}")
                        continue

                    df_params_median_indiv = pd.read_table(tsv_params_median_fn, sep="\t")

                    if i == 0:
                        df_params_median = df_params_median_indiv.copy()
                    else:
                        df_params_median = pd.concat([df_params_median, df_params_median_indiv])

                # Group output directory
                tsv_dir = '{}/{}/derivatives/pp_data/group/{}/pmf/tsv'.format(
                    main_dir, project_dir, format_)
                os.makedirs(tsv_dir, exist_ok=True)

                # Parameters median: median across subjects per roi
                # --------------------------------------------------
                numeric_cols = df_params_median.select_dtypes(include=[np.number]).columns.tolist()

                df_params_median_group = (df_params_median
                                          .groupby(['roi'], sort=False)[numeric_cols]
                                          .median()
                                          .reset_index())

                tsv_params_median_fn = "{}/{}_{}_{}_params-median.tsv".format(
                    tsv_dir, subject, fn_spec, pmf_task_name)
                print('Saving tsv: {}'.format(tsv_params_median_fn))
                df_params_median_group.to_csv(tsv_params_median_fn, sep="\t", na_rep='NaN', index=False)

# Define permission cmd
print('Changing files permissions in {}/{}'.format(main_dir, project_dir))
os.system("chmod -Rf 771 {}/{}".format(main_dir, project_dir))
os.system("chgrp -Rf {} {}/{}".format(group, main_dir, project_dir))