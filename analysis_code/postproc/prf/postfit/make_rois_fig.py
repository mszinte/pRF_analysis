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
cd ~/projects/RetinoMaps/analysis_code/postproc/prf/postfit/
python make_rois_fig.py /scratch/mszinte/data RetinoMaps sub-01 327
python make_rois_fig.py /scratch/mszinte/data RetinoMaps sub-170k 327
python make_rois_fig.py /scratch/mszinte/data RetinoMaps group 327
-----------------------------------------------------------------------------------------
Written by Uriel Lascombes (uriel.lascombes@laposte.net)
Edited by Martin Szinte (mail@martinszinte.net)
-----------------------------------------------------------------------------------------
"""

# Stop warnings
import warnings
warnings.filterwarnings("ignore")

# General imports
import os
import sys
import json
import pandas as pd
import ipdb
deb = ipdb.set_trace

# Personal import
sys.path.append("{}/../../../utils".format(os.getcwd()))
from plot_utils import *

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
group = sys.argv[4]

# Load settings
with open('../../../settings.json') as f:
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
pcm_th = analysis_info['pcm_th']
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
plot_groups = [['V1', 'V2', 'V3'], ['V3AB', 'LO', 'VO'], ['hMT+', 'iIPS', 'sIPS'], ['iPCS', 'sPCS', 'mPCS']]
num_ecc_size_bins = 6
num_ecc_pcm_bins = 6
num_polar_angle_bins = 9
max_ecc = 15
fig_width = 1440

# Format loop
for format_, extension in zip(formats, extensions):
    
    # Create folders and fns
    fig_dir = '{}/{}/derivatives/pp_data/{}/{}/prf/figures'.format(
        main_dir, project_dir, subject, format_)
    os.makedirs(fig_dir, exist_ok=True)

    # Load data
    df_roi_area, df_violins, df_ecc_size, df_ecc_pcm, df_polar_angle, df_contralaterality, df_params_avg = \
        compute_plot_data(subject=subject,
                          main_dir=main_dir,
                          project_dir=project_dir,
                          format_=format_,
                          rois=rois,
                          amplitude_threshold=amplitude_th, 
                          ecc_threshold=ecc_th, 
                          size_threshold=size_th, 
                          rsqr_threshold=rsqr_th,
                          n_threshold=n_th,
                          stats_threshold=stats_th,
                          pcm_threshold=pcm_th,
                          num_ecc_size_bins=num_ecc_size_bins, 
                          num_ecc_pcm_bins=num_ecc_pcm_bins,
                          num_polar_angle_bins=num_polar_angle_bins,
                          max_ecc=max_ecc,
                          subjects_to_group=subjects)
    
    # Roi area and stats plot
    fig = prf_roi_area(df_roi_area=df_roi_area, fig_width=fig_width, fig_height=300, roi_colors=roi_colors)
    fig_fn = "{}/{}_prf_roi_area.pdf".format(fig_dir, subject)
    print('Saving pdf: {}'.format(fig_fn))
    fig.write_image(fig_fn)
    
    # Violins plot
    fig = prf_violins_plot(df_violins=df_violins, fig_width=fig_width, fig_height=600, 
                           rois=rois, roi_colors=roi_colors)
    fig_fn = "{}/{}_prf_violins.pdf".format(fig_dir, subject)
    print('Saving pdf: {}'.format(fig_fn))
    fig.write_image(fig_fn)

    # Parameters average plot
    fig = prf_params_avg_plot(df_params_avg=df_params_avg, fig_width=fig_width, fig_height=600, 
                              rois=rois, roi_colors=roi_colors)
    fig_fn = "{}/{}_prf_params_avg.pdf".format(fig_dir, subject)
    print('Saving pdf: {}'.format(fig_fn))
    fig.write_image(fig_fn)
    
    # Ecc.size plots
    fig = prf_ecc_size_plot(df_ecc_size=df_ecc_size, fig_width=fig_width, 
                            fig_height=400, rois=rois, roi_colors=roi_colors,
                            plot_groups=plot_groups, max_ecc=max_ecc)
    fig_fn = "{}/{}_prf_ecc_size.pdf".format(fig_dir, subject)
    print('Saving pdf: {}'.format(fig_fn))
    fig.write_image(fig_fn)

    # Ecc.pCM plot
    fig_fn = "{}/{}_prf_ecc_pcm.pdf".format(fig_dir, subject)
    fig = prf_ecc_pcm_plot(df_ecc_pcm=df_ecc_pcm, fig_width=fig_width, fig_height=400, 
                           rois=rois, roi_colors=roi_colors,
                           plot_groups=plot_groups, max_ecc=max_ecc)
    print('Saving pdf: {}'.format(fig_fn))
    fig.write_image(fig_fn)
    
    # Polar angle distributions
    figs, hemis = prf_polar_angle_plot(df_polar_angle=df_polar_angle, fig_width=fig_width, 
                                       fig_height=300, rois=rois, roi_colors=roi_colors,
                                       num_polar_angle_bins=num_polar_angle_bins)
    for (fig, hemi) in zip(figs, hemis):
        fig_fn = "{}/{}_prf_polar_angle_{}.pdf".format(fig_dir, subject, hemi)
        print('Saving pdf: {}'.format(fig_fn))
        fig.write_image(fig_fn)

    # Contralaterality plots
    fig_fn = "{}/{}_contralaterality.pdf".format(fig_dir, subject)
    fig = prf_contralaterality_plot(df_contralaterality=df_contralaterality, 
                                    fig_width=fig_width, fig_height=300, 
                                    rois=rois, roi_colors=roi_colors)
    print('Saving pdf: {}'.format(fig_fn))
    fig.write_image(fig_fn)

    # Spatial distibution plot
    
# Define permission cmd
print('Changing files permissions in {}/{}'.format(main_dir, project_dir))
os.system("chmod -Rf 771 {}/{}".format(main_dir, project_dir))
os.system("chgrp -Rf {} {}/{}".format(group, main_dir, project_dir))