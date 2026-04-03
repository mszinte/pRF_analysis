"""
-----------------------------------------------------------------------------------------
make_intertask_rois_fig.py
-----------------------------------------------------------------------------------------
Goal of the script:
Make ROIs-based CSS figures
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject name (e.g. sub-01)
sys.argv[4]: server group (e.g. 327)
-----------------------------------------------------------------------------------------
Output(s):
CSS analysis figures
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/pRF_analysis/RetinoMaps/postproc/intertask/
2. run python command
python make_rois_fig.py [main directory] [project name] [subject] [group]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/projects/pRF_analysis/RetinoMaps/postproc/intertask/
python make_intertask_rois_fig.py /scratch/mszinte/data RetinoMaps sub-01 327
python make_intertask_rois_fig.py /scratch/mszinte/data RetinoMaps sub-170k 327
python make_intertask_rois_fig.py /scratch/mszinte/data RetinoMaps group 327
-----------------------------------------------------------------------------------------
Written by Uriel Lascombes (uriel.lascombes@laposte.net)
Edited by Martin Szinte (mail@martinszinte.net)
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
sys.path.append("{}/../../../analysis_code/utils".format(os.getcwd()))
from plot_utils import *
from settings_utils import load_settings

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
group = sys.argv[4]

# Load settings
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../"))
settings_path = os.path.join(base_dir, project_dir, "settings.yml")
prf_settings_path = os.path.join(base_dir, project_dir, "prf-analysis.yml")
glm_settings_path = os.path.join(base_dir, project_dir, "glm-analysis.yml")
figure_settings_path = os.path.join(base_dir, project_dir, "figure-settings.yml")
settings = load_settings([settings_path, prf_settings_path, glm_settings_path, figure_settings_path])
analysis_info = settings[0]

formats = analysis_info['formats']
extensions = analysis_info['extensions']
subjects = analysis_info['subjects']
group_tasks = analysis_info['task_intertask']
preproc_prep = analysis_info['preproc_prep']
filtering = analysis_info['filtering']
normalization = analysis_info['normalization']
avg_methods = analysis_info['avg_methods']
prf_task_names = analysis_info['prf_task_names']
rois_methods = analysis_info['rois_methods']
prf_task_name = analysis_info['prf_task_names'][0]
categories_to_plot = analysis_info['categories_to_plot']
avg_method = 'loo-avg'
rsq2use = 'prf_loo_rsq'

for tasks in group_tasks : 
    print(tasks)
    if 'SacVELoc' in tasks: intertask_group = 'SacVE-PurVE-pRF'
    else : intertask_group = 'Sac-Pur-pRF'
    # Format loop
    for format_, extension in zip(formats, extensions):
        print(format_)
        # define list of rois for each format
        rois_methods_format = rois_methods[format_]
        for rois_method_format in rois_methods_format:
            print(rois_method_format)
            if rois_method_format == 'rois-drawn':
                analysis_info['rois'] = analysis_info[rois_method_format]
            elif rois_method_format == 'rois-group-mmp':
                analysis_info['rois'] = list(analysis_info[rois_method_format].keys())
        
        
            # Define/create folders
            tsv_dir = '{}/{}/derivatives/pp_data/{}/{}/intertask/tsv'.format(
                main_dir, project_dir, subject, format_)
            
            # Exception if no data for one format (e.g template subject)
            if not os.path.isdir(tsv_dir):
                print(f"[SKIP] tsv_dir not found for format={format_}: {tsv_dir}")
                continue
            
            fig_dir = '{}/{}/derivatives/pp_data/{}/{}/intertask/figures'.format(
                main_dir, project_dir, subject, format_)
            os.makedirs(fig_dir, exist_ok=True)
    
            fn_spec = "task-{}_{}_{}_{}_{}_{}".format(
                prf_task_name, preproc_prep, filtering, 
                normalization, avg_method, rois_method_format)
                        
            # loop over categories
            for categorie_to_plot in categories_to_plot : 
                print(categorie_to_plot)
                fig_dir_categorie = '{}/{}/derivatives/pp_data/{}/{}/intertask/figures/{}'.format(
                    main_dir, project_dir, subject, format_, categorie_to_plot)
                os.makedirs(fig_dir_categorie, exist_ok=True)
                tsv_dir_categorie = '{}/{}/derivatives/pp_data/{}/{}/intertask/tsv/tsv_{}'.format(
                    main_dir, project_dir, subject, format_, categorie_to_plot)
                fn_spec_cate = "task-{}_{}_{}_{}_{}_{}_{}".format(
                    intertask_group, preproc_prep, filtering,
                    normalization, avg_method, rois_method_format, categorie_to_plot)
                                
                # Violins plot
                if 'vision' in categorie_to_plot:
                    tsv_violins_fn = "{}/{}_{}_intertask_violins.tsv".format(tsv_dir_categorie, subject, fn_spec_cate)
                    df_violins = pd.read_table(tsv_violins_fn, sep="\t")
                    fig = prf_violins_plot(df=df_violins, figure_info=analysis_info, rsq2use=rsq2use)
                    fig_fn = "{}/{}_{}_intertask_violins.pdf".format(fig_dir_categorie, subject, fn_spec_cate)
                    print('Saving pdf: {}'.format(fig_fn))
                    fig.write_image(fig_fn)
                    remove_second_page(fig_fn)
                
                # Parameters median plot
                tsv_params_median_fn = "{}/{}_{}_intertask_params-median.tsv".format(tsv_dir_categorie, subject, fn_spec_cate)
                df_params_median = pd.read_table(tsv_params_median_fn, sep="\t")
                fig = prf_params_median_plot(df=df_params_median, figure_info=analysis_info, rsq2use=rsq2use)
                fig_fn = "{}/{}_{}_intertask_params-median.pdf".format(fig_dir_categorie, subject, fn_spec_cate)
                print('Saving pdf: {}'.format(fig_fn))
                fig.write_image(fig_fn)
                remove_second_page(fig_fn)
    
                # Ecc.size plots
                tsv_ecc_size_fn = "{}/{}_{}_intertask_ecc-size.tsv".format(tsv_dir_categorie, subject, fn_spec_cate)
                df_ecc_size = pd.read_table(tsv_ecc_size_fn, sep="\t")
                fig = prf_ecc_size_plot(df=df_ecc_size, figure_info=analysis_info, rsq2use=rsq2use)
                fig_fn = "{}/{}_{}_intertask_ecc-size.pdf".format(fig_dir_categorie, subject, fn_spec_cate)
                print('Saving pdf: {}'.format(fig_fn))
                fig.write_image(fig_fn)
                remove_second_page(fig_fn)
                
                # Ecc.pCM plot
                tsv_ecc_pcm_fn = "{}/{}_{}_intertask_ecc-pcm.tsv".format(tsv_dir_categorie, subject, fn_spec_cate)
                df_ecc_pcm = pd.read_table(tsv_ecc_pcm_fn, sep="\t")
                fig_fn = "{}/{}_{}_intertask_ecc-pcm.pdf".format(fig_dir_categorie, subject, fn_spec_cate)
                fig = prf_ecc_pcm_plot(df=df_ecc_pcm, figure_info=analysis_info, rsq2use=rsq2use)
                print('Saving pdf: {}'.format(fig_fn))
                fig.write_image(fig_fn)
                remove_second_page(fig_fn)
    
                # Polar angle distributions
                tsv_polar_angle_fn = "{}/{}_{}_intertask_polar-angle.tsv".format(tsv_dir_categorie, subject, fn_spec_cate)
                df_polar_angle = pd.read_table(tsv_polar_angle_fn, sep="\t")
                figs, hemis = prf_polar_angle_plot(df=df_polar_angle, figure_info=analysis_info)
                for (fig, hemi) in zip(figs, hemis):
                    if hemi == 'hemi-LR':
                        fig_fn = "{}/{}_{}_intertask_polar-angle.pdf".format(fig_dir_categorie, subject, fn_spec_cate)
                        print('Saving pdf: {}'.format(fig_fn))
                        fig.write_image(fig_fn)
                        remove_second_page(fig_fn)
    
                # Contralaterality plots
                tsv_contralaterality_fn = "{}/{}_{}_intertask_contralaterality.tsv".format(tsv_dir_categorie, subject, fn_spec_cate)
                df_contralaterality = pd.read_table(tsv_contralaterality_fn, sep="\t")
                fig_fn = "{}/{}_{}_intertask_contralaterality.pdf".format(fig_dir_categorie, subject, fn_spec_cate)
                fig = prf_contralaterality_plot(df=df_contralaterality, figure_info=analysis_info)
                print('Saving pdf: {}'.format(fig_fn))
                fig.write_image(fig_fn)
                remove_second_page(fig_fn)
                
                # Spatial distribution plot
                tsv_distribution_fn = "{}/{}_{}_intertask_distribution.tsv".format(tsv_dir_categorie, subject, fn_spec_cate)
                df_distribution = pd.read_table(tsv_distribution_fn, sep="\t")
                figs, hemis = prf_distribution_plot(df=df_distribution, figure_info=analysis_info)
                for (fig, hemi) in zip(figs, hemis):
                    if hemi == 'hemi-LR':
                        fig_fn = "{}/{}_{}_intertask_distribution.pdf".format(fig_dir_categorie, subject, fn_spec_cate)
                        print('Saving pdf: {}'.format(fig_fn))
                        fig.write_image(fig_fn)
                        remove_second_page(fig_fn)
                
                # Spatial distibution barycentre plot
                tsv_barycentre_fn = "{}/{}_{}_intertask_barycentre.tsv".format(tsv_dir_categorie, subject, fn_spec_cate)
                try:
                    df_barycentre = pd.read_table(tsv_barycentre_fn, sep="\t")
                except FileNotFoundError:
                    print(f"File not found: {tsv_barycentre_fn}, skipping...")
                    continue
                fig_fn = "{}/{}_{}_intertask_barycentre.pdf".format(fig_dir_categorie, subject, fn_spec_cate)
                fig = prf_barycentre_plot(df=df_barycentre, figure_info=analysis_info)
                print('Saving pdf: {}'.format(fig_fn))
                fig.write_image(fig_fn)
                remove_second_page(fig_fn)
            
# # Define permission cmd
# print('Changing files permissions in {}/{}'.format(main_dir, project_dir))
# os.system("chmod -Rf 771 {}/{}".format(main_dir, project_dir))
# os.system("chgrp -Rf {} {}/{}".format(group, main_dir, project_dir))