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
>> cd ~/projects/pRF_analysis/RetinoMaps/intertask/
2. run python command
python make_rois_fig.py [main directory] [project name] [subject] [group]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/projects/pRF_analysis/RetinoMaps/intertask/
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
import json
import pandas as pd

# Personal import
sys.path.append("{}/../../analysis_code/utils".format(os.getcwd()))
from plot_utils import *

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
group = sys.argv[4]

# Load settings
# Define analysis parameters
with open('../settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)
if subject == 'sub-170k': 
    formats = ['170k']
    extensions = ['dtseries.nii']
else: 
    formats = analysis_info['formats']
    extensions = analysis_info['extensions']
rois = analysis_info['rois']
subjects = analysis_info['subjects']
group_tasks = analysis_info['task_intertask']
categories = ['saccade', 'pursuit', 'vision', 'all']

# Figure settings
roi_colors = {'V1': 'rgb(243, 231, 155)', 
              'V2': 'rgb(250, 196, 132)', 
              'V3': 'rgb(248, 160, 126)', 
              'V3AB': 'rgb(235, 127, 134)', 
              'LO': 'rgb(150, 0, 90)',  
              'VO': 'rgb(0, 0, 200)', 
              'hMT+': 'rgb(0, 25, 255)', 
              'iIPS': 'rgb(0, 152, 255)', 
              'sIPS': 'rgb(44, 255, 150)', 
              'iPCS': 'rgb(151, 255, 0)', 
              'sPCS': 'rgb(255, 234, 0)', 
              'mPCS': 'rgb(255, 111, 0)'}

categorie_color_map = {'pursuit': '#E377C2', 
                       'saccade': '#8C564B', 
                       'pursuit_and_saccade': '#9467BD', 
                       'vision': '#D62728', 
                       'vision_and_pursuit': '#2CA02C', 
                       'vision_and_saccade': '#FF7F0E', 
                       'vision_and_pursuit_and_saccade': '#1F77B4'}

with open('../figure_settings.json') as f:
    json_s = f.read()
    figure_info = json.loads(json_s)
num_ecc_size_bins = figure_info['num_ecc_size_bins']
num_ecc_pcm_bins = figure_info['num_ecc_pcm_bins']
num_polar_angle_bins = figure_info['num_polar_angle_bins']
max_ecc = figure_info['max_ecc']
screen_side = figure_info['screen_side']
gaussian_mesh_grain = figure_info['gaussian_mesh_grain']
hot_zone_percent = figure_info['hot_zone_percent']
plot_groups = figure_info['plot_groups']
fig_width = figure_info['fig_width']
categories_to_plot = figure_info['categories_to_plot']
categories_active_vertex_plot = figure_info['categories_active_vertex_plot']

for tasks in group_tasks : 
    if 'SacVELoc' in tasks: suffix = 'SacVE_PurVE'
    else : suffix = 'Sac_Pur'
    # Format loop
    for format_, extension in zip(formats, extensions):
        figures_dict = {}
        figures_titles = []
        # Create folders and fns
        fig_dir = '{}/{}/derivatives/pp_data/{}/{}/intertask/figures'.format(
            main_dir, project_dir, subject, format_)
        os.makedirs(fig_dir, exist_ok=True)
        tsv_dir = '{}/{}/derivatives/pp_data/{}/{}/intertask/tsv'.format(
            main_dir, project_dir, subject, format_)
        
        # Active vertex roi plot
        active_vertex_roi_tsv_fn = '{}/{}_active_vertex_roi_{}.tsv'.format(tsv_dir, subject, suffix)
        df_active_vertex_roi = pd.read_table(active_vertex_roi_tsv_fn)
        fig = active_vertex_roi_plot(
            df_active_vertex_roi=df_active_vertex_roi, fig_height=400, 
            fig_width=fig_width, roi_colors=roi_colors, plot_groups=plot_groups)       
        fig_fn = "{}/{}_active_vertex_roi_{}.pdf".format(fig_dir, subject, suffix)
        print('Saving pdf: {}'.format(fig_fn))
        fig.write_image(fig_fn)
        figures_dict['all_categories'] = []
        figures_dict['all_categories'].append(fig)
        figures_titles.append('Active vertex roi plot')
        
        # Active vertex roi mmp plot
        if format_ == '170k':
            for categorie_active_vertex in categories_active_vertex_plot: 
                figures_dict[categorie_active_vertex] = []
                active_vertex_roi_mmp_tsv_fn = '{}/{}_active_vertex_roi_mmp_{}.tsv'.format(tsv_dir, subject, suffix)
                df_active_vertex_roi_mmp = pd.read_table(active_vertex_roi_mmp_tsv_fn)
                fig = active_vertex_roi_mmp_plot(
                    df_active_vertex_roi_mmp=df_active_vertex_roi_mmp, fig_height=1080, 
                    fig_width=1080, roi_colors=roi_colors, plot_groups=plot_groups, categorie=categorie_active_vertex)
                fig_fn = "{}/{}_active_vertex_{}_roi_mmp_{}.pdf".format(fig_dir, subject, categorie_active_vertex, suffix)
                print('Saving pdf: {}'.format(fig_fn))
                fig.write_image(fig_fn)
                figures_dict[categorie_active_vertex].append(fig)
                figures_titles.append('Active vertex roi mmp plot')
        
        # loop over categories
        for categorie_to_plot in categories_to_plot : 
            fig_dir_categorie = '{}/{}/derivatives/pp_data/{}/{}/intertask/figures/{}'.format(
                main_dir, project_dir, subject, format_, categorie_to_plot)
            os.makedirs(fig_dir_categorie, exist_ok=True)
            tsv_dir_categorie = '{}/{}/derivatives/pp_data/{}/{}/intertask/tsv/tsv_{}'.format(
                main_dir, project_dir, subject, format_, categorie_to_plot)
            
            figures_dict[categorie_to_plot] = []
            
            # Violins plot
            tsv_violins_fn = "{}/{}_{}_prf_violins_{}.tsv".format(tsv_dir_categorie, subject, categorie_to_plot, suffix)
            df_violins = pd.read_table(tsv_violins_fn, sep="\t")
            fig = prf_violins_plot(df_violins=df_violins, fig_width=fig_width, fig_height=600, 
                                    rois=rois, roi_colors=roi_colors)
            fig_fn = "{}/{}_{}_prf_violins_{}.pdf".format(fig_dir_categorie, subject, categorie_to_plot, suffix)
            print('Saving pdf: {}'.format(fig_fn))
            fig.write_image(fig_fn)
            figures_dict[categorie_to_plot].append(fig)
            figures_titles.append('Vilolins plot')
        
            # Parameters median plot
            tsv_params_median_fn = "{}/{}_{}_prf_params_median_{}.tsv".format(tsv_dir_categorie, subject, categorie_to_plot, suffix)
            df_params_median = pd.read_table(tsv_params_median_fn, sep="\t")
            fig = prf_params_median_plot(df_params_avg=df_params_median, fig_width=fig_width, fig_height=600, 
                                      rois=rois, roi_colors=roi_colors)
            fig_fn = "{}/{}_{}_prf_params_median_{}.pdf".format(fig_dir_categorie, subject, categorie_to_plot, suffix)
            print('Saving pdf: {}'.format(fig_fn))
            fig.write_image(fig_fn)
            figures_dict[categorie_to_plot].append(fig)
            figures_titles.append('Parameters median')
            
            # Ecc.size plots
            tsv_ecc_size_fn = "{}/{}_{}_prf_ecc_size_{}.tsv".format(tsv_dir_categorie, subject, categorie_to_plot, suffix)
            df_ecc_size = pd.read_table(tsv_ecc_size_fn, sep="\t")
            fig = prf_ecc_size_plot(df_ecc_size=df_ecc_size, fig_width=fig_width, 
                                    fig_height=400, rois=rois, roi_colors=roi_colors,
                                    plot_groups=plot_groups, max_ecc=max_ecc)
            fig_fn = "{}/{}_{}_prf_ecc_size_{}.pdf".format(fig_dir_categorie, subject, categorie_to_plot, suffix)
            print('Saving pdf: {}'.format(fig_fn))
            fig.write_image(fig_fn)
            figures_dict[categorie_to_plot].append(fig)
            figures_titles.append('Eccentricity size relation')
        
            # Ecc.pCM plot
            tsv_ecc_pcm_fn = "{}/{}_{}_prf_ecc_pcm_{}.tsv".format(tsv_dir_categorie, subject, categorie_to_plot, suffix)
            df_ecc_pcm = pd.read_table(tsv_ecc_pcm_fn, sep="\t")
            fig_fn = "{}/{}_{}_prf_ecc_pcm_{}.pdf".format(fig_dir_categorie, subject, categorie_to_plot, suffix)
            fig = prf_ecc_pcm_plot(df_ecc_pcm=df_ecc_pcm, fig_width=fig_width, fig_height=400, 
                                    rois=rois, roi_colors=roi_colors,
                                    plot_groups=plot_groups, max_ecc=max_ecc)
            print('Saving pdf: {}'.format(fig_fn))
            fig.write_image(fig_fn)
            figures_dict[categorie_to_plot].append(fig)
            figures_titles.append('Eccentricity pCM relation')
            
            # Polar angle distributions
            tsv_contralaterality_fn = "{}/{}_{}_prf_contralaterality_{}.tsv".format(tsv_dir_categorie, subject, categorie_to_plot, suffix)
            df_contralaterality = pd.read_table(tsv_contralaterality_fn, sep="\t")
            tsv_polar_angle_fn = "{}/{}_{}_prf_polar_angle_{}.tsv".format(tsv_dir_categorie, subject, categorie_to_plot, suffix)
            df_polar_angle = pd.read_table(tsv_polar_angle_fn, sep="\t")
            figs, hemis = prf_polar_angle_plot(df_polar_angle=df_polar_angle, fig_width=fig_width, 
                                                fig_height=300, rois=rois, roi_colors=roi_colors,
                                                num_polar_angle_bins=num_polar_angle_bins)
            for (fig, hemi) in zip(figs, hemis):
                if hemi == 'hemi-LR':
                    fig_fn = "{}/{}_{}_prf_polar_angle_{}.pdf".format(fig_dir_categorie, subject, categorie_to_plot, suffix)
                    print('Saving pdf: {}'.format(fig_fn))
                    fig.write_image(fig_fn)
                    figures_dict[categorie_to_plot].append(fig)
                    figures_titles.append('Polar angle distribution')
                    
            # Contralaterality plots
            fig_fn = "{}/{}_{}_contralaterality_{}.pdf".format(fig_dir_categorie, subject, categorie_to_plot, suffix)
            fig = prf_contralaterality_plot(df_contralaterality=df_contralaterality, 
                                            fig_width=fig_width, fig_height=300, 
                                            rois=rois, roi_colors=roi_colors)
            print('Saving pdf: {}'.format(fig_fn))
            fig.write_image(fig_fn)
            figures_dict[categorie_to_plot].append(fig)
            figures_titles.append('Contralaterality')
        
            # # Spatial distibution plot
            # tsv_distribution_fn = "{}/{}_{}_prf_distribution_{}.tsv".format(tsv_dir_categorie, subject, categorie_to_plot, suffix)
            # df_distribution = pd.read_table(tsv_distribution_fn, sep="\t")
            # figs, hemis = prf_distribution_plot(df_distribution=df_distribution, 
            #                                     fig_width=fig_width, fig_height=300, 
            #                                     rois=rois, roi_colors=roi_colors, screen_side=screen_side)
        
            # for (fig, hemi) in zip(figs, hemis):
            #     if hemi == 'hemi-LR':
            #         fig_fn = "{}/{}_{}_distribution_{}.pdf".format(fig_dir_categorie, subject, categorie_to_plot, suffix)
            #         print('Saving pdf: {}'.format(fig_fn))
            #         fig.write_image(fig_fn)
            #         figures_dict[categorie_to_plot].append(fig)
            #         figures_titles.append('pRF distribution')
        
        # Export html with all figures
        subject_html = make_figures_html(subject=subject, figures=figures_dict, figs_title=figures_titles)
        html_fn = '{}/{}_figures_html_{}.html'.format(fig_dir, subject, suffix)
        with open(html_fn, "w") as f:
            f.write(subject_html)
                
# # Define permission cmd
# print('Changing files permissions in {}/{}'.format(main_dir, project_dir))
# os.system("chmod -Rf 771 {}/{}".format(main_dir, project_dir))
# os.system("chgrp -Rf {} {}/{}".format(group, main_dir, project_dir))