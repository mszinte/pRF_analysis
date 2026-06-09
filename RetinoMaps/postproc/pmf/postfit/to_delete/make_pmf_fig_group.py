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
avg_methods = ['concat']
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

                # ROI active vertex
                # -----------------
                tsv_roi_active_vert_fn = f"{tsv_dir}/{subject}_{fn_spec}_{pmf_task_name}_active-vert.tsv"
                if not os.path.isfile(tsv_roi_active_vert_fn):
                    print(f"[SKIP] File not found: {tsv_roi_active_vert_fn}")
                    continue
                df_roi_active_vert = pd.read_table(tsv_roi_active_vert_fn, sep="\t")
                fig = prf_roi_active_vert_plot(df=df_roi_active_vert, figure_info=analysis_info, format=format_)
                fig_fn = f"{fig_dir}/{subject}_{fn_spec}_{pmf_task_name}_active-vert.pdf"
                print('Saving pdf: {}'.format(fig_fn))
                fig.write_image(fig_fn)
                remove_second_page(fig_fn)

                # Violins plot
                # ------------
                tsv_violins_fn = f"{tsv_dir}/{subject}_{fn_spec}_{pmf_task_name}_violins.tsv"
                if not os.path.isfile(tsv_violins_fn):
                    print(f"[SKIP] File not found: {tsv_violins_fn}")
                    continue
                df_violins = pd.read_table(tsv_violins_fn, sep="\t")
                fig = prf_violins_plot(df=df_violins, figure_info=analysis_info, rsq2use=rsq2use)
                fig_fn = f"{fig_dir}/{subject}_{fn_spec}_{pmf_task_name}_violins.pdf"
                print('Saving pdf: {}'.format(fig_fn))
                fig.write_image(fig_fn)
                remove_second_page(fig_fn)

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

                # Ecc.size plot
                # -------------
                tsv_ecc_size_fn = f"{tsv_dir}/{subject}_{fn_spec}_{pmf_task_name}_ecc-size.tsv"
                if not os.path.isfile(tsv_ecc_size_fn):
                    print(f"[SKIP] File not found: {tsv_ecc_size_fn}")
                    continue
                df_ecc_size = pd.read_table(tsv_ecc_size_fn, sep="\t")
                fig = prf_ecc_size_plot(df=df_ecc_size, figure_info=analysis_info, rsq2use=rsq2use)
                fig_fn = f"{fig_dir}/{subject}_{fn_spec}_{pmf_task_name}_ecc-size.pdf"
                print('Saving pdf: {}'.format(fig_fn))
                fig.write_image(fig_fn)
                remove_second_page(fig_fn)


                # Polar angle distributions
                # -------------------------
                tsv_polar_angle_fn = f"{tsv_dir}/{subject}_{fn_spec}_{pmf_task_name}_polar-angle.tsv"
                if not os.path.isfile(tsv_polar_angle_fn):
                    print(f"[SKIP] File not found: {tsv_polar_angle_fn}")
                    continue
                df_polar_angle = pd.read_table(tsv_polar_angle_fn, sep="\t")
                figs, hemis = prf_polar_angle_plot(df=df_polar_angle, figure_info=analysis_info)
                for (fig, hemi) in zip(figs, hemis):
                    if hemi == 'hemi-LR':
                        fig_fn = f"{fig_dir}/{subject}_{fn_spec}_{pmf_task_name}_polar-angle.pdf"
                        print('Saving pdf: {}'.format(fig_fn))
                        fig.write_image(fig_fn)
                        remove_second_page(fig_fn)

                # Contralaterality plot
                # ---------------------
                tsv_contralaterality_fn = f"{tsv_dir}/{subject}_{fn_spec}_{pmf_task_name}_contralaterality.tsv"
                if not os.path.isfile(tsv_contralaterality_fn):
                    print(f"[SKIP] File not found: {tsv_contralaterality_fn}")
                    continue
                df_contralaterality = pd.read_table(tsv_contralaterality_fn, sep="\t")
                fig = prf_contralaterality_plot(df=df_contralaterality, figure_info=analysis_info)
                fig_fn = f"{fig_dir}/{subject}_{fn_spec}_{pmf_task_name}_contralaterality.pdf"
                print('Saving pdf: {}'.format(fig_fn))
                fig.write_image(fig_fn)
                remove_second_page(fig_fn)

                # Spatial distribution plot
                # -------------------------
                tsv_distribution_fn = f"{tsv_dir}/{subject}_{fn_spec}_{pmf_task_name}_distribution.tsv"
                if not os.path.isfile(tsv_distribution_fn):
                    print(f"[SKIP] File not found: {tsv_distribution_fn}")
                    continue
                df_distribution = pd.read_table(tsv_distribution_fn, sep="\t")
                figs, hemis = prf_distribution_plot(df=df_distribution, figure_info=analysis_info)
                for (fig, hemi) in zip(figs, hemis):
                    if hemi == 'hemi-LR':
                        fig_fn = f"{fig_dir}/{subject}_{fn_spec}_{pmf_task_name}_distribution.pdf"
                        print('Saving pdf: {}'.format(fig_fn))
                        fig.write_image(fig_fn)
                        remove_second_page(fig_fn)

                # Spatial distribution barycentre plot
                # ------------------------------------
                tsv_barycentre_fn = f"{tsv_dir}/{subject}_{fn_spec}_{pmf_task_name}_barycentre.tsv"
                if not os.path.isfile(tsv_barycentre_fn):
                    print(f"[SKIP] File not found: {tsv_barycentre_fn}")
                    continue
                df_barycentre = pd.read_table(tsv_barycentre_fn, sep="\t")
                fig = prf_barycentre_plot(df=df_barycentre, figure_info=analysis_info)
                fig_fn = f"{fig_dir}/{subject}_{fn_spec}_{pmf_task_name}_barycentre.pdf"
                print('Saving pdf: {}'.format(fig_fn))
                fig.write_image(fig_fn)
                remove_second_page(fig_fn)

# Define permission cmd
print('Changing files permissions in {}/{}'.format(main_dir, project_dir))
os.system("chmod -Rf 771 {}/{}".format(main_dir, project_dir))
os.system("chgrp -Rf {} {}/{}".format(group, main_dir, project_dir))