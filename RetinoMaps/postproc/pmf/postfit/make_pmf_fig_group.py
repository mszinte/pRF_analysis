"""
-----------------------------------------------------------------------------------------
make_pmf_fig_group.py
-----------------------------------------------------------------------------------------
Goal of the script:
Make PMF figure for params-median - group level, across 3 fit conditions
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject name (group)
sys.argv[4]: server group (e.g. 327)
-----------------------------------------------------------------------------------------
Output(s):
PMF params-median figures (pdf) per format and condition
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/pRF_analysis/RetinoMaps/postproc/pmf/postfit/
2. run python command
python make_pmf_fig_group.py [main directory] [project name] group [group]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/projects/pRF_analysis/RetinoMaps/postproc/pmf/postfit/
python make_pmf_fig_group.py /scratch/mszinte/data RetinoMaps group 327
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
import pandas as pd

# Personal import
script_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.abspath(os.path.join(script_dir, "../../../../analysis_code/utils")))
from plot_utils import *
from settings_utils import load_settings

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]   # should be 'group'
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
preproc_prep = analysis_info['preproc_prep']
filtering = analysis_info['filtering']
normalization = analysis_info['normalization']
avg_methods = analysis_info['avg_methods']
rois_methods = analysis_info['rois_methods']

# The 3 PMF fit conditions
pmf_task_names = ['pmf-bold-gauss', 'pmf-residuals-gauss', 'prf-control-gauss']

# Main loop
for avg_method in avg_methods:
    if 'loo' in avg_method: rsq2use = 'prf_loo_rsq'
    else: rsq2use = 'prf_rsq'

    for format_, extension in zip(formats, extensions):
        # define list of rois for each format
        rois_methods_format = rois_methods[format_]
        for rois_method_format in rois_methods_format:
            if rois_method_format == 'rois-drawn':
                analysis_info['rois'] = analysis_info[rois_method_format]
            elif rois_method_format == 'rois-group-mmp':
                analysis_info['rois'] = list(analysis_info[rois_method_format].keys())

            for pmf_task_name in pmf_task_names:

                print(f'{pmf_task_name} - {avg_method} - {format_} - {rois_method_format}')

                # Define folders
                tsv_dir = '{}/{}/derivatives/pp_data/group/{}/pmf/tsv'.format(
                    main_dir, project_dir, format_)
                fig_dir = '{}/{}/derivatives/pp_data/group/{}/pmf/figures'.format(
                    main_dir, project_dir, format_)
                os.makedirs(fig_dir, exist_ok=True)

                fn_spec = "task-{}_{}_{}_{}_{}_{}".format(
                    'SacLoc', preproc_prep, filtering,
                    normalization, avg_method, rois_method_format)

                # Parameters median plot
                # ----------------------
                tsv_params_median_fn = f"{tsv_dir}/{subject}_{fn_spec}_{pmf_task_name}_params-median.tsv"

                if not os.path.isfile(tsv_params_median_fn):
                    print(f"[SKIP] File not found: {tsv_params_median_fn}")
                    continue

                df_params_median = pd.read_table(tsv_params_median_fn, sep="\t")
                fig = prf_params_median_plot(df=df_params_median, figure_info=analysis_info, rsq2use=rsq2use)
                fig_fn = f"{fig_dir}/{subject}_{fn_spec}_{pmf_task_name}_params-median.pdf"
                print('Saving pdf: {}'.format(fig_fn))
                fig.write_image(fig_fn)
                remove_second_page(fig_fn)

# Define permission cmd
print('Changing files permissions in {}/{}'.format(main_dir, project_dir))
os.system("chmod -Rf 771 {}/{}".format(main_dir, project_dir))
os.system("chgrp -Rf {} {}/{}".format(group, main_dir, project_dir))