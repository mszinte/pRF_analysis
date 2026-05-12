"""
-----------------------------------------------------------------------------------------
make_rois_fig_gauss.py
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
>> cd ~/projects/pRF_analysis/RetinoMaps/postproc/pmf/postfit/
2. run python command
python make_rois_fig.py [main directory] [project name] [subject] [group] 
-----------------------------------------------------------------------------------------
Exemple:
cd ~/projects/pRF_analysis/analysis_code/postproc/pmf/postfit/
python make_rois_fig_gauss.py /scratch/mszinte/data RetinoMaps sub-01 327 bold 
python make_rois_fig_gauss.py /scratch/mszinte/data RetinoMaps sub-01 327 residuals
python make_rois_fig_gauss.py /scratch/mszinte/data RetinoMaps sub-hcp1.6mm 327
python make_rois_fig_gauss.py /scratch/mszinte/data RetinoMaps group 327
python make_rois_fig_gauss.py /scratch/mszinte/data amblyo7T_prf sub-02 327
-----------------------------------------------------------------------------------------
Written by Uriel Lascombes (uriel.lascombes@laposte.net)
Edited by Martin Szinte (martin.szinte@gmail.com) 
Adapted by Sina Kling (sina.kling@outlook.de)
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
subject = sys.argv[3]
group = sys.argv[4]
deriv_condition = sys.argv[5]


# Load settings
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../../"))
settings_path = os.path.join(base_dir, project_dir, "settings.yml")
prf_settings_path = os.path.join(base_dir, project_dir, "pmf-settings.yml")
figure_settings_path = os.path.join(base_dir, project_dir, "figure-settings.yml")
settings = load_settings([settings_path, prf_settings_path, figure_settings_path])
analysis_info = settings[0]

formats = analysis_info['formats']
extensions = analysis_info['extensions']
subjects = analysis_info['subjects']
preproc_prep = analysis_info['preproc_prep']
filtering = analysis_info['filtering']
normalization = analysis_info['normalization']
avg_methods = ['concat']
prf_task_names = ['SacLoc']
rois_methods = analysis_info['rois_methods']

# Define fit_type-specific filename patterns
fit_type_config = {
    'residuals': {
        'deriv_tag':  'pmf-residuals-gauss_deriv',
        'output_tag': 'pmf-residuals',
    },
    'bold': {
        'deriv_tag':  'pmf2-gauss_deriv',
        'output_tag': 'pmf-bold',
    },
    'prf': {
        'deriv_tag':  'prf-gauss_deriv',
        'output_tag': 'prf-control',
    },
}

cfg        = fit_type_config[deriv_condition]
deriv_tag  = cfg['deriv_tag']
output_tag = cfg['output_tag']


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

            for prf_task_name in prf_task_names:
                
                print(f'{prf_task_name} - {avg_method} - {format_}')
    
                # Define/create folders
                tsv_dir = '{}/{}/derivatives/pp_data/{}/{}/pmf/tsv'.format(
                    main_dir, project_dir, subject, format_)
                
                # Exception if no data for one format (e.g template subject)
                if not os.path.isdir(tsv_dir):
                    print(f"[SKIP] tsv_dir not found for format={format_}: {tsv_dir}")
                    continue
                
                fig_dir = '{}/{}/derivatives/pp_data/{}/{}/pmf/figures'.format(
                    main_dir, project_dir, subject, format_)
                os.makedirs(fig_dir, exist_ok=True)

                fn_spec = "task-{}_{}_{}_{}_{}_{}".format(
                    prf_task_name, preproc_prep, filtering, 
                    normalization, avg_method, rois_method_format)

                # Roi active vertex
                tsv_roi_active_vert_fn = "{}/{}_{}_{}-gauss_active-vert.tsv".format(tsv_dir, subject, fn_spec, output_tag)
                df_roi_active_vert = pd.read_table(tsv_roi_active_vert_fn, sep="\t")
                fig = prf_roi_active_vert_plot(df=df_roi_active_vert, figure_info=analysis_info, format=format_, )
                fig_fn = "{}/{}_{}_{}-gauss_active-vert.pdf".format(fig_dir, subject, fn_spec, output_tag)
                print('Saving pdf: {}'.format(fig_fn))
                fig.write_image(fig_fn)
                remove_second_page(fig_fn)
                
                # Violins plot
                tsv_violins_fn = "{}/{}_{}_{}-gauss_violins.tsv".format(tsv_dir, subject, fn_spec,output_tag)
                df_violins = pd.read_table(tsv_violins_fn, sep="\t")
                fig = prf_violins_plot(df=df_violins, figure_info=analysis_info, rsq2use=rsq2use)
                fig_fn = "{}/{}_{}_{}-gauss_violins.pdf".format(fig_dir, subject, fn_spec, output_tag)
                print('Saving pdf: {}'.format(fig_fn))
                fig.write_image(fig_fn)
                remove_second_page(fig_fn)
                
                # Parameters median plot
                tsv_params_median_fn = "{}/{}_{}_{}-gauss_params-median.tsv".format(tsv_dir, subject, fn_spec, output_tag)
                df_params_median = pd.read_table(tsv_params_median_fn, sep="\t")
                fig = prf_params_median_plot(df=df_params_median, figure_info=analysis_info, rsq2use=rsq2use)
                fig_fn = "{}/{}_{}_{}-gauss_params-median.pdf".format(fig_dir, subject, fn_spec, output_tag)
                print('Saving pdf: {}'.format(fig_fn))
                fig.write_image(fig_fn)
                remove_second_page(fig_fn)
    
                # Ecc.size plots
                tsv_ecc_size_fn = "{}/{}_{}_{}-gauss_ecc-size.tsv".format(tsv_dir, subject, fn_spec, output_tag)
                df_ecc_size = pd.read_table(tsv_ecc_size_fn, sep="\t")
                fig = prf_ecc_size_plot(df=df_ecc_size, figure_info=analysis_info, rsq2use=rsq2use)
                fig_fn = "{}/{}_{}_{}-gauss_ecc-size.pdf".format(fig_dir, subject, fn_spec, output_tag)
                print('Saving pdf: {}'.format(fig_fn))
                fig.write_image(fig_fn)
                remove_second_page(fig_fn)
                
    
                # Polar angle distributions
                tsv_polar_angle_fn = "{}/{}_{}_{}-gauss_polar-angle.tsv".format(tsv_dir, subject, fn_spec, output_tag)
                df_polar_angle = pd.read_table(tsv_polar_angle_fn, sep="\t")
                figs, hemis = prf_polar_angle_plot(df=df_polar_angle, figure_info=analysis_info)
                for (fig, hemi) in zip(figs, hemis):
                    if hemi == 'hemi-LR':
                        fig_fn = "{}/{}_{}_{}-gauss_polar-angle.pdf".format(fig_dir, subject, fn_spec, output_tag)
                        print('Saving pdf: {}'.format(fig_fn))
                        fig.write_image(fig_fn)
                        remove_second_page(fig_fn)
    
                # Contralaterality plots
                tsv_contralaterality_fn = "{}/{}_{}_{}-gauss_contralaterality.tsv".format(tsv_dir, subject, fn_spec, output_tag)
                df_contralaterality = pd.read_table(tsv_contralaterality_fn, sep="\t")
                fig_fn = "{}/{}_{}_{}-gauss_contralaterality.pdf".format(fig_dir, subject, fn_spec, output_tag)
                fig = prf_contralaterality_plot(df=df_contralaterality, figure_info=analysis_info)
                print('Saving pdf: {}'.format(fig_fn))
                fig.write_image(fig_fn)
                remove_second_page(fig_fn)
                
                # Spatial distribution plot
                tsv_distribution_fn = "{}/{}_{}_{}-gauss_distribution.tsv".format(tsv_dir, subject,fn_spec, output_tag)
                df_distribution = pd.read_table(tsv_distribution_fn, sep="\t")
                figs, hemis = prf_distribution_plot(df=df_distribution, figure_info=analysis_info)
                for (fig, hemi) in zip(figs, hemis):
                    if hemi == 'hemi-LR':
                        fig_fn = "{}/{}_{}_{}-gauss_distribution.pdf".format(fig_dir, subject, fn_spec, output_tag)
                        print('Saving pdf: {}'.format(fig_fn))
                        fig.write_image(fig_fn)
                        remove_second_page(fig_fn)
                
                # Spatial distibution barycentre plot
                tsv_barycentre_fn = "{}/{}_{}_{}-gauss_barycentre.tsv".format(tsv_dir, subject, fn_spec, output_tag)
                try:
                    df_barycentre = pd.read_table(tsv_barycentre_fn, sep="\t")
                except FileNotFoundError:
                    print(f"File not found: {tsv_barycentre_fn}, skipping...")
                    continue
                fig_fn = "{}/{}_{}_{}-gauss_barycentre.pdf".format(fig_dir, subject, fn_spec, output_tag)
                fig = prf_barycentre_plot(df=df_barycentre, figure_info=analysis_info)
                print('Saving pdf: {}'.format(fig_fn))
                fig.write_image(fig_fn)
                remove_second_page(fig_fn)
    
# Define permission cmd
print('Changing files permissions in {}/{}'.format(main_dir, project_dir))
os.system("chmod -Rf 771 {}/{}".format(main_dir, project_dir))
os.system("chgrp -Rf {} {}/{}".format(group, main_dir, project_dir))