"""
-----------------------------------------------------------------------------------------
make_rois_fig.py
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
>> cd ~/projects/[PROJECT]/analysis_code/postproc/prf/postfit/
2. run python command
python make_rois_fig.py [main directory] [project name] [subject] [group]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/projects/pRF_analysis/analysis_code/postproc/prf/postfit/
python make_rois_fig.py /scratch/mszinte/data MotConf sub-01 327
python make_rois_fig.py /scratch/mszinte/data MotConf sub-170k 327
python make_rois_fig.py /scratch/mszinte/data MotConf group 327
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
sys.path.append("{}/../../../utils".format(os.getcwd()))
from plot_utils import *

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
group = sys.argv[4]

# Load settings
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../../"))
settings_path = os.path.join(base_dir, project_dir, "settings.json")

with open(settings_path) as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)
if subject == 'sub-170k': 
    formats = ['170k']
    extensions = ['dtseries.nii']
else: 
    formats = analysis_info['formats']
    extensions = analysis_info['extensions']
rois = analysis_info['rois']

# Threshold settings
ecc_th = analysis_info['ecc_th']
size_th = analysis_info['size_th']
rsqr_th = analysis_info['rsqr_th']
amplitude_th = analysis_info['amplitude_th']
stats_th = analysis_info['stats_th']
n_th = analysis_info['n_th'] 
subjects = analysis_info['subjects']

# Figure settings
colormap_dict = {'V1': (243, 231, 155),
                 'V2': (250, 196, 132),
                 'V3': (248, 160, 126),
                 'V3AB': (235, 127, 134),
                 'LO': (150, 0, 90), 
                 'VO': (0, 0, 200),
                 'hMT+': (0, 25, 255),
                 'iIPS': (0, 152, 255),
                 'sIPS': (44, 255, 150),
                 'iPCS': (151, 255, 0),
                 'sPCS': (255, 234, 0),
                 'mPCS': (255, 111, 0)
                }
roi_colors = ['rgb({},{},{})'.format(*rgb) for rgb in colormap_dict.values()]

fig_settings_path = os.path.join(base_dir, project_dir, "figure_settings.json")

with open(settings_path) as f:
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

# Format loop
for format_, extension in zip(formats, extensions):
    
    # Create folders and fns
    fig_dir = '{}/{}/derivatives/pp_data/{}/{}/prf/figures'.format(
        main_dir, project_dir, subject, format_)
    tsv_dir = '{}/{}/derivatives/pp_data/{}/{}/prf/tsv'.format(
        main_dir, project_dir, subject, format_)
    os.makedirs(fig_dir, exist_ok=True)

    # Roi area and stats plot
    tsv_roi_area_fn = "{}/{}_prf_roi_area.tsv".format(tsv_dir, subject)
    df_roi_area = pd.read_table(tsv_roi_area_fn, sep="\t")
    fig = prf_roi_area(df_roi_area=df_roi_area, fig_width=fig_width, fig_height=300, roi_colors=roi_colors)
    fig_fn = "{}/{}_prf_roi_area.pdf".format(fig_dir, subject)
    print('Saving pdf: {}'.format(fig_fn))
    fig.write_image(fig_fn)
    
    # Violins plot
    tsv_violins_fn = "{}/{}_prf_violins.tsv".format(tsv_dir, subject)
    df_violins = pd.read_table(tsv_violins_fn, sep="\t")
    fig = prf_violins_plot(df_violins=df_violins, fig_width=fig_width, fig_height=600, 
                            rois=rois, roi_colors=roi_colors)
    fig_fn = "{}/{}_prf_violins.pdf".format(fig_dir, subject)
    print('Saving pdf: {}'.format(fig_fn))
    fig.write_image(fig_fn)

    # Parameters median plot
    tsv_params_median_fn = "{}/{}_prf_params_median.tsv".format(tsv_dir, subject)
    df_params_median = pd.read_table(tsv_params_median_fn, sep="\t")
    fig = prf_params_median_plot(df_params_avg=df_params_median, fig_width=fig_width, fig_height=600, 
                              rois=rois, roi_colors=roi_colors)
    fig_fn = "{}/{}_prf_params_median.pdf".format(fig_dir, subject)
    print('Saving pdf: {}'.format(fig_fn))
    fig.write_image(fig_fn)
    
    # Ecc.size plots
    tsv_ecc_size_fn = "{}/{}_prf_ecc_size.tsv".format(tsv_dir, subject)
    df_ecc_size = pd.read_table(tsv_ecc_size_fn, sep="\t")
    fig = prf_ecc_size_plot(df_ecc_size=df_ecc_size, fig_width=fig_width, 
                            fig_height=400, rois=rois, roi_colors=roi_colors,
                            plot_groups=plot_groups, max_ecc=max_ecc)
    fig_fn = "{}/{}_prf_ecc_size.pdf".format(fig_dir, subject)
    print('Saving pdf: {}'.format(fig_fn))
    fig.write_image(fig_fn)

    # Ecc.pCM plot
    tsv_ecc_pcm_fn = "{}/{}_prf_ecc_pcm.tsv".format(tsv_dir, subject)
    df_ecc_pcm = pd.read_table(tsv_ecc_pcm_fn, sep="\t")
    fig_fn = "{}/{}_prf_ecc_pcm.pdf".format(fig_dir, subject)
    fig = prf_ecc_pcm_plot(df_ecc_pcm=df_ecc_pcm, fig_width=fig_width, fig_height=400, 
                            rois=rois, roi_colors=roi_colors,
                            plot_groups=plot_groups, max_ecc=max_ecc)
    print('Saving pdf: {}'.format(fig_fn))
    fig.write_image(fig_fn)
    
    # Polar angle distributions
    tsv_contralaterality_fn = "{}/{}_prf_contralaterality.tsv".format(tsv_dir, subject)
    df_contralaterality = pd.read_table(tsv_contralaterality_fn, sep="\t")
    tsv_polar_angle_fn = "{}/{}_prf_polar_angle.tsv".format(tsv_dir, subject)
    df_polar_angle = pd.read_table(tsv_polar_angle_fn, sep="\t")
    figs, hemis = prf_polar_angle_plot(df_polar_angle=df_polar_angle, fig_width=fig_width, 
                                        fig_height=300, rois=rois, roi_colors=roi_colors,
                                        num_polar_angle_bins=num_polar_angle_bins)
    for (fig, hemi) in zip(figs, hemis):
        if hemi == 'hemi-LR':
            fig_fn = "{}/{}_prf_polar_angle.pdf".format(fig_dir, subject)
            print('Saving pdf: {}'.format(fig_fn))
            fig.write_image(fig_fn)

    # Contralaterality plots
    tsv_distribution_fn = "{}/{}_prf_distribution.tsv".format(tsv_dir, subject)
    df_distribution = pd.read_table(tsv_distribution_fn, sep="\t")
    fig_fn = "{}/{}_contralaterality.pdf".format(fig_dir, subject)
    fig = prf_contralaterality_plot(df_contralaterality=df_contralaterality, 
                                    fig_width=fig_width, fig_height=300, 
                                    rois=rois, roi_colors=roi_colors)
    print('Saving pdf: {}'.format(fig_fn))
    fig.write_image(fig_fn)

    # Spatial distibution plot
    figs, hemis = prf_distribution_plot(df_distribution=df_distribution, 
                                        fig_width=fig_width, fig_height=300, 
                                        rois=rois, roi_colors=roi_colors, screen_side=screen_side)

    for (fig, hemi) in zip(figs, hemis):
        if hemi == 'hemi-LR':
            fig_fn = "{}/{}_distribution.pdf".format(fig_dir, subject)
            print('Saving pdf: {}'.format(fig_fn))
            fig.write_image(fig_fn)

    # # Spatial distibution barycentre plot
    # tsv_barycentre_fn = "{}/{}_prf_barycentre.tsv".format(tsv_dir, subject)
    # df_barycentre = pd.read_table(tsv_barycentre_fn, sep="\t")
    # fig_fn = "{}/{}_barycentre.pdf".format(fig_dir, subject)
    # fig = prf_barycentre_plot(df_barycentre=df_barycentre, 
    #                                 fig_width=fig_width, fig_height=400, 
    #                                 rois=rois, roi_colors=roi_colors, screen_side=screen_side)
    # print('Saving pdf: {}'.format(fig_fn))
    # fig.write_image(fig_fn)
    
# Define permission cmd
print('Changing files permissions in {}/{}'.format(main_dir, project_dir))
os.system("chmod -Rf 771 {}/{}".format(main_dir, project_dir))
os.system("chgrp -Rf {} {}/{}".format(group, main_dir, project_dir))