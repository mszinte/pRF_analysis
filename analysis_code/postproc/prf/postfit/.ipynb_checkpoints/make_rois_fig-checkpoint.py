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
>> cd ~/projects/pRF_analysis/analysis_code/postproc/prf/postfit/
2. run python command
python make_rois_fig.py [main directory] [project name] [subject] [group] 
-----------------------------------------------------------------------------------------
Exemple:
cd ~/projects/pRF_analysis/analysis_code/postproc/prf/postfit/
python make_rois_fig.py /scratch/mszinte/data RetinoMaps sub-01 327
python make_rois_fig.py /scratch/mszinte/data RetinoMaps sub-170k 327
python make_rois_fig.py /scratch/mszinte/data RetinoMaps group 327
-----------------------------------------------------------------------------------------
Written by Uriel Lascombes (uriel.lascombes@laposte.net)
Edited by Martin Szinte (martin.szinte@gmail.com)
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
import yaml
import pandas as pd

# Personal import
sys.path.append("{}/../../../utils".format(os.getcwd()))
from plot_utils import *
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

if subject == 'sub-170k': 
    formats = ['170k']
    extensions = ['dtseries.nii']
else: 
    formats = analysis_info['formats']
    extensions = analysis_info['extensions']
rois = analysis_info['rois']
subjects = analysis_info['subjects']
preproc_prep = analysis_info['preproc_prep']
filtering = analysis_info['filtering']
normalization = analysis_info['normalization']
avg_methods = analysis_info['avg_methods']
prf_task_names = analysis_info['prf_task_names']

# Main loop
for avg_method in avg_methods:
    if 'loo' in avg_method: rsq2use = 'prf_loo_rsq'
    else: rsq2use = 'prf_rsq'
        
    for format_, extension in zip(formats, extensions):
        
        for prf_task_name in prf_task_names:
            
            print(f'{prf_task_name} - {avg_method} - {format_}')

            # Define/create folders
            tsv_dir = '{}/{}/derivatives/pp_data/{}/{}/prf/tsv'.format(
                main_dir, project_dir, subject, format_)
            fig_dir = '{}/{}/derivatives/pp_data/{}/{}/prf/figures'.format(
                main_dir, project_dir, subject, format_)
            os.makedirs(fig_dir, exist_ok=True)
            fn_spec = "task-{}_{}_{}_{}_{}".format(
                prf_task_name, preproc_prep, filtering, normalization, avg_method)
            
            # Roi active vertex
            tsv_roi_active_vert_fn = "{}/{}_{}_prf-css_active-vert.tsv".format(tsv_dir, subject, fn_spec)
            df_roi_active_vert = pd.read_table(tsv_roi_active_vert_fn, sep="\t")
            fig = prf_roi_active_vert_plot(df=df_roi_active_vert, figure_info=analysis_info, format=format_, )
            fig_fn = "{}/{}_{}_prf-css_active-vert.pdf".format(fig_dir, subject, fn_spec)
            print('Saving pdf: {}'.format(fig_fn))
            fig.write_image(fig_fn)
            remove_second_page(fig_fn)
            
            # Violins plot
            tsv_violins_fn = "{}/{}_{}_prf-css_violins.tsv".format(tsv_dir, subject, fn_spec)
            df_violins = pd.read_table(tsv_violins_fn, sep="\t")
            fig = prf_violins_plot(df=df_violins, figure_info=analysis_info, rsq2use=rsq2use)
            fig_fn = "{}/{}_{}_prf-css_violins.pdf".format(fig_dir, subject, fn_spec)
            print('Saving pdf: {}'.format(fig_fn))
            fig.write_image(fig_fn)
            remove_second_page(fig_fn)
            
            # Parameters median plot
            tsv_params_median_fn = "{}/{}_{}_prf-css_params-median.tsv".format(tsv_dir, subject, fn_spec)
            df_params_median = pd.read_table(tsv_params_median_fn, sep="\t")
            fig = prf_params_median_plot(df=df_params_median, figure_info=analysis_info, rsq2use=rsq2use)
            fig_fn = "{}/{}_{}_prf-css_params-median.pdf".format(fig_dir, subject, fn_spec)
            print('Saving pdf: {}'.format(fig_fn))
            fig.write_image(fig_fn)
            remove_second_page(fig_fn)

            # Ecc.size plots
            tsv_ecc_size_fn = "{}/{}_{}_prf-css_ecc-size.tsv".format(tsv_dir, subject, fn_spec)
            df_ecc_size = pd.read_table(tsv_ecc_size_fn, sep="\t")
            fig = prf_ecc_size_plot(df=df_ecc_size, figure_info=analysis_info, rsq2use=rsq2use)
            fig_fn = "{}/{}_{}_prf-css_ecc-size.pdf".format(fig_dir, subject, fn_spec)
            print('Saving pdf: {}'.format(fig_fn))
            fig.write_image(fig_fn)
            remove_second_page(fig_fn)
            
            # Ecc.pCM plot
            tsv_ecc_pcm_fn = "{}/{}_{}_prf-css_ecc-pcm.tsv".format(tsv_dir, subject, fn_spec)
            df_ecc_pcm = pd.read_table(tsv_ecc_pcm_fn, sep="\t")
            fig_fn = "{}/{}_{}_prf-css_ecc-pcm.pdf".format(fig_dir, subject, fn_spec)
            fig = prf_ecc_pcm_plot(df=df_ecc_pcm, figure_info=analysis_info, rsq2use=rsq2use)
            print('Saving pdf: {}'.format(fig_fn))
            fig.write_image(fig_fn)
            remove_second_page(fig_fn)

            # Polar angle distributions
            tsv_polar_angle_fn = "{}/{}_{}_prf-css_polar-angle.tsv".format(tsv_dir, subject, fn_spec)
            df_polar_angle = pd.read_table(tsv_polar_angle_fn, sep="\t")
            figs, hemis = prf_polar_angle_plot(df=df_polar_angle, figure_info=analysis_info)
            for (fig, hemi) in zip(figs, hemis):
                if hemi == 'hemi-LR':
                    fig_fn = "{}/{}_{}_prf-css_polar-angle.pdf".format(fig_dir, subject, fn_spec)
                    print('Saving pdf: {}'.format(fig_fn))
                    fig.write_image(fig_fn)
                    remove_second_page(fig_fn)

            # Contralaterality plots
            tsv_contralaterality_fn = "{}/{}_{}_prf-css_contralaterality.tsv".format(tsv_dir, subject, fn_spec)
            df_contralaterality = pd.read_table(tsv_contralaterality_fn, sep="\t")
            fig_fn = "{}/{}_{}_prf-css_contralaterality.pdf".format(fig_dir, subject, fn_spec)
            fig = prf_contralaterality_plot(df=df_contralaterality, figure_info=analysis_info)
            print('Saving pdf: {}'.format(fig_fn))
            fig.write_image(fig_fn)
            remove_second_page(fig_fn)
            
            # Spatial distribution plot
            tsv_distribution_fn = "{}/{}_{}_prf-css_distribution.tsv".format(tsv_dir, subject,fn_spec )
            df_distribution = pd.read_table(tsv_distribution_fn, sep="\t")
            figs, hemis = prf_distribution_plot(df=df_distribution, figure_info=analysis_info)
            for (fig, hemi) in zip(figs, hemis):
                if hemi == 'hemi-LR':
                    fig_fn = "{}/{}_{}_prf-css_distribution.pdf".format(fig_dir, subject, fn_spec)
                    print('Saving pdf: {}'.format(fig_fn))
                    fig.write_image(fig_fn)
                    remove_second_page(fig_fn)
            
            # Spatial distibution barycentre plot
            tsv_barycentre_fn = "{}/{}_{}_prf-css_barycentre.tsv".format(tsv_dir, subject, fn_spec)
            df_barycentre = pd.read_table(tsv_barycentre_fn, sep="\t")
            fig_fn = "{}/{}_{}_prf-css_barycentre.pdf".format(fig_dir, subject, fn_spec)
            fig = prf_barycentre_plot(df=df_barycentre, figure_info=analysis_info)
            print('Saving pdf: {}'.format(fig_fn))
            fig.write_image(fig_fn)
            remove_second_page(fig_fn)
    
# Define permission cmd
print('Changing files permissions in {}/{}'.format(main_dir, project_dir))
os.system("chmod -Rf 771 {}/{}".format(main_dir, project_dir))
os.system("chgrp -Rf {} {}/{}".format(group, main_dir, project_dir))