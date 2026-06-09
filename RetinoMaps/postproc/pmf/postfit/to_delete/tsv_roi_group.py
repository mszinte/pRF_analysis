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
from maths_utils import weighted_nan_median, weighted_nan_percentile, make_prf_barycentre_df

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
distribution_max_ecc = analysis_info['distribution_max_ecc']
distribution_mesh_grain = analysis_info['distribution_mesh_grain']
hot_zone_percent = analysis_info['hot_zone_percent']

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

        # define list of rois for each format
        rois_methods_format = rois_methods[format_]
        for rois_method_format in rois_methods_format:

            if rois_method_format == 'rois-drawn':
                rois = analysis_info[rois_method_format]
            elif rois_method_format == 'rois-group-mmp':
                rois = list(analysis_info[rois_method_format].keys())

            for pmf_task_name in pmf_task_names:

                print(f'{pmf_task_name} - {avg_method} - {format_} - {rois_method_format}')

                fn_spec = "task-{}_{}_{}_{}_{}_{}".format(
                    'SacLoc', preproc_prep, filtering, normalization, avg_method, rois_method_format)

                for i, subject_to_group in enumerate(subjects_to_group):
                    tsv_dir = '{}/{}/derivatives/pp_data/{}/{}/pmf/tsv'.format(
                        main_dir, project_dir, subject_to_group, format_)

                    # ROI active vertices
                    # -------------------
                    tsv_roi_active_vert_fn = "{}/{}_{}_{}_active-vert.tsv".format(
                        tsv_dir, subject_to_group, fn_spec, pmf_task_name)
                    if not os.path.isfile(tsv_roi_active_vert_fn):
                        print(f"[SKIP] File not found: {tsv_roi_active_vert_fn}")
                        continue
                    df_roi_active_vert_indiv = pd.read_table(tsv_roi_active_vert_fn, sep="\t")
                    if i == 0: df_roi_active_vert = df_roi_active_vert_indiv.copy()
                    else: df_roi_active_vert = pd.concat([df_roi_active_vert, df_roi_active_vert_indiv])

                    # Violins
                    # -------
                    tsv_violins_fn = "{}/{}_{}_{}_violins.tsv".format(
                        tsv_dir, subject_to_group, fn_spec, pmf_task_name)
                    if not os.path.isfile(tsv_violins_fn):
                        print(f"[SKIP] File not found: {tsv_violins_fn}")
                        continue
                    df_violins_indiv = pd.read_table(tsv_violins_fn, sep="\t")
                    if i == 0: df_violins = df_violins_indiv.copy()
                    else: df_violins = pd.concat([df_violins, df_violins_indiv])

                    # Ecc.size
                    # --------
                    tsv_ecc_size_fn = "{}/{}_{}_{}_ecc-size.tsv".format(
                        tsv_dir, subject_to_group, fn_spec, pmf_task_name)
                    if not os.path.isfile(tsv_ecc_size_fn):
                        print(f"[SKIP] File not found: {tsv_ecc_size_fn}")
                        continue
                    df_ecc_size_indiv = pd.read_table(tsv_ecc_size_fn, sep="\t")
                    if i == 0: df_ecc_size = df_ecc_size_indiv.copy()
                    else: df_ecc_size = pd.concat([df_ecc_size, df_ecc_size_indiv])


                    # Polar angle
                    # -----------
                    tsv_polar_angle_fn = "{}/{}_{}_{}_polar-angle.tsv".format(
                        tsv_dir, subject_to_group, fn_spec, pmf_task_name)
                    if not os.path.isfile(tsv_polar_angle_fn):
                        print(f"[SKIP] File not found: {tsv_polar_angle_fn}")
                        continue
                    df_polar_angle_indiv = pd.read_table(tsv_polar_angle_fn, sep="\t")
                    if i == 0: df_polar_angle = df_polar_angle_indiv.copy()
                    else: df_polar_angle = pd.concat([df_polar_angle, df_polar_angle_indiv])

                    # Contralaterality
                    # ----------------
                    tsv_contralaterality_fn = "{}/{}_{}_{}_contralaterality.tsv".format(
                        tsv_dir, subject_to_group, fn_spec, pmf_task_name)
                    if not os.path.isfile(tsv_contralaterality_fn):
                        print(f"[SKIP] File not found: {tsv_contralaterality_fn}")
                        continue
                    df_contralaterality_indiv = pd.read_table(tsv_contralaterality_fn, sep="\t")
                    if i == 0: df_contralaterality = df_contralaterality_indiv.copy()
                    else: df_contralaterality = pd.concat([df_contralaterality, df_contralaterality_indiv])

                    # Parameters median
                    # ------------------
                    tsv_params_median_fn = "{}/{}_{}_{}_params-median.tsv".format(
                        tsv_dir, subject_to_group, fn_spec, pmf_task_name)
                    if not os.path.isfile(tsv_params_median_fn):
                        print(f"[SKIP] File not found: {tsv_params_median_fn}")
                        continue
                    df_params_median_indiv = pd.read_table(tsv_params_median_fn, sep="\t")
                    if i == 0: df_params_median = df_params_median_indiv.copy()
                    else: df_params_median = pd.concat([df_params_median, df_params_median_indiv])

                    # Spatial distribution
                    # --------------------
                    tsv_distribution_fn = "{}/{}_{}_{}_distribution.tsv".format(
                        tsv_dir, subject_to_group, fn_spec, pmf_task_name)
                    if not os.path.isfile(tsv_distribution_fn):
                        print(f"[SKIP] File not found: {tsv_distribution_fn}")
                        continue
                    df_distribution_indiv = pd.read_table(tsv_distribution_fn, sep="\t")
                    mesh_indiv = df_distribution_indiv.drop(columns=['roi', 'x', 'y', 'hemi']).values
                    others_columns = df_distribution_indiv[['roi', 'x', 'y', 'hemi']]
                    if i == 0: mesh_group = np.expand_dims(mesh_indiv, axis=0)
                    else: mesh_group = np.vstack((mesh_group, np.expand_dims(mesh_indiv, axis=0)))

                # Group output directory
                tsv_dir = '{}/{}/derivatives/pp_data/group/{}/pmf/tsv'.format(
                    main_dir, project_dir, format_)
                os.makedirs(tsv_dir, exist_ok=True)

                # ROI active vertices
                # -------------------
                df_roi_active_vert = df_roi_active_vert.groupby(['roi'], sort=False).median().reset_index()
                tsv_roi_active_vert_fn = "{}/{}_{}_{}_active-vert.tsv".format(
                    tsv_dir, subject, fn_spec, pmf_task_name)
                print('Saving tsv: {}'.format(tsv_roi_active_vert_fn))
                df_roi_active_vert.to_csv(tsv_roi_active_vert_fn, sep="\t", na_rep='NaN', index=False)

                # Violins
                # -------
                tsv_violins_fn = "{}/{}_{}_{}_violins.tsv".format(
                    tsv_dir, subject, fn_spec, pmf_task_name)
                print('Saving tsv: {}'.format(tsv_violins_fn))
                df_violins.to_csv(tsv_violins_fn, sep="\t", na_rep='NaN', index=False)

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

                # Ecc.size
                # --------
                df_ecc_size = df_ecc_size.groupby(['roi', 'num_bins'], sort=False).median().reset_index()
                tsv_ecc_size_fn = "{}/{}_{}_{}_ecc-size.tsv".format(
                    tsv_dir, subject, fn_spec, pmf_task_name)
                print('Saving tsv: {}'.format(tsv_ecc_size_fn))
                df_ecc_size.to_csv(tsv_ecc_size_fn, sep="\t", na_rep='NaN', index=False)

                # Polar angle
                # -----------
                df_polar_angle = df_polar_angle.groupby(['roi', 'hemi', 'num_bins'], sort=False).median().reset_index()
                tsv_polar_angle_fn = "{}/{}_{}_{}_polar-angle.tsv".format(
                    tsv_dir, subject, fn_spec, pmf_task_name)
                print('Saving tsv: {}'.format(tsv_polar_angle_fn))
                df_polar_angle.to_csv(tsv_polar_angle_fn, sep="\t", na_rep='NaN', index=False)

                # Contralaterality
                # ----------------
                df_contralaterality = df_contralaterality.groupby(['roi'], sort=False).median().reset_index()
                tsv_contralaterality_fn = "{}/{}_{}_{}_contralaterality.tsv".format(
                    tsv_dir, subject, fn_spec, pmf_task_name)
                print('Saving tsv: {}'.format(tsv_contralaterality_fn))
                df_contralaterality.to_csv(tsv_contralaterality_fn, sep="\t", na_rep='NaN', index=False)

                # Spatial distribution
                # --------------------
                median_mesh = np.median(mesh_group, axis=0)
                df_distribution = pd.DataFrame(median_mesh)
                df_distribution = pd.concat([others_columns.reset_index(drop=True), df_distribution], axis=1)
                tsv_distribution_fn = "{}/{}_{}_{}_distribution.tsv".format(
                    tsv_dir, subject, fn_spec, pmf_task_name)
                print('Saving tsv: {}'.format(tsv_distribution_fn))
                df_distribution.to_csv(tsv_distribution_fn, sep="\t", na_rep='NaN', index=False)

                # Spatial distribution hot zone barycentre
                # ----------------------------------------
                hemis = ['hemi-L', 'hemi-R', 'hemi-LR']
                for j, hemi in enumerate(hemis):
                    hemi_values = ['hemi-L', 'hemi-R'] if hemi == 'hemi-LR' else [hemi]
                    df_distribution_hemi = df_distribution.loc[df_distribution.hemi.isin(hemi_values)]
                    df_barycentre_hemi = make_prf_barycentre_df(
                        df_distribution_hemi, rois, distribution_max_ecc,
                        distribution_mesh_grain, hot_zone_percent=hot_zone_percent)

                    df_barycentre_hemi['hemi'] = [hemi] * len(df_barycentre_hemi)
                    if j == 0: df_barycentre = df_barycentre_hemi
                    else: df_barycentre = pd.concat([df_barycentre, df_barycentre_hemi])

                tsv_barycentre_fn = "{}/{}_{}_{}_barycentre.tsv".format(
                    tsv_dir, subject, fn_spec, pmf_task_name)
                print('Saving tsv: {}'.format(tsv_barycentre_fn))
                df_barycentre.to_csv(tsv_barycentre_fn, sep="\t", na_rep='NaN', index=False)

# Define permission cmd
print('Changing files permissions in {}/{}'.format(main_dir, project_dir))
os.system("chmod -Rf 771 {}/{}".format(main_dir, project_dir))
os.system("chgrp -Rf {} {}/{}".format(group, main_dir, project_dir))