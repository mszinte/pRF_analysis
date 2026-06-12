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
sys.argv[4]: analysis name 9e.g prf
sys.argv[5]: server group (e.g. 327)
-----------------------------------------------------------------------------------------
Output(s):
CSS analysis figures
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/pRF_analysis/analysis_code/postproc/prf/postfit/
2. run python command
python make_rois_fig.py [main directory] [project name] [subject] [analysis name] [group] 
-----------------------------------------------------------------------------------------
Exemple:
cd ~/projects/pRF_analysis/analysis_code/postproc/prf/postfit/
python make_rois_fig.py /scratch/mszinte/data RetinoMaps sub-01 prf 327
python make_rois_fig.py /scratch/mszinte/data RetinoMaps sub-hcp1.6mm prf 327
python make_rois_fig.py /scratch/mszinte/data RetinoMaps group prf 327
python make_rois_fig.py /scratch/mszinte/data amblyo7T_prf sub-02 prf 327
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
import pandas as pd

# Personal import
sys.path.append("{}/../../../utils".format(os.getcwd()))
from plot_utils import *
from settings_utils import load_settings

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
analysis_name = sys.argv[4]
group = sys.argv[5]

# Load settings
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../../"))
general_settings_path = os.path.join(base_dir, project_dir, "settings.yml")
analysis_settings_path = os.path.join(base_dir, project_dir, f"{analysis_name}-analysis.yml")
figure_settings_path = os.path.join(base_dir, project_dir, "figure-settings.yml")
settings = load_settings([general_settings_path, analysis_settings_path, figure_settings_path])
analysis_info = settings[0]

formats = analysis_info['formats']
extensions = analysis_info['extensions']
subjects = analysis_info['subjects']
preproc_prep = analysis_info['preproc_prep']
filtering = analysis_info['filtering']
normalization = analysis_info['normalization']
avg_methods = analysis_info['avg_methods']
task_names = analysis_info['analysis_task_names']
rois_methods = analysis_info['rois_methods']

output_folder = analysis_info["output_folder"]
dm_name = analysis_info["dm_name"]

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
                # rois_dict = analysis_info['rois-group-mmp']
                # analysis_info['rois'] = [item for sublist in rois_dict.values() for item in sublist]

            for task_name in task_names:
                
                print(f'{task_name} - {avg_method} - {format_}')
    
                # Define/create folders
                tsv_dir = f'{main_dir}/{project_dir}/derivatives/pp_data/{subject}/{format_}/{output_folder}/tsv'
                
                # Exception if no data for one format (e.g template subject)
                if not os.path.isdir(tsv_dir):
                    print(f"[SKIP] tsv_dir not found for format={format_}: {tsv_dir}")
                    continue
                
                fig_dir = f'{main_dir}/{project_dir}/derivatives/pp_data/{subject}/{format_}/{output_folder}/figures'
                os.makedirs(fig_dir, exist_ok=True)

                fn_spec = f"task-{task_name}_{preproc_prep}_{filtering}_{normalization}_{avg_method}_{rois_method_format}"

                # Roi active vertex
                tsv_roi_active_vert_fn = f"{tsv_dir}/{subject}_{fn_spec}_{analysis_name}-css{dm_name}_active-vert.tsv"
                df_roi_active_vert = pd.read_table(tsv_roi_active_vert_fn, sep="\t")
                fig = prf_roi_active_vert_plot(df=df_roi_active_vert, figure_info=analysis_info, format=format_, )
                fig_fn = f"{fig_dir}/{subject}_{fn_spec}_{analysis_name}-css{dm_name}_active-vert.pdf"
                print('Saving pdf: {}'.format(fig_fn))
                fig.write_image(fig_fn)
                remove_second_page(fig_fn)
                
                # Violins plot
                tsv_violins_fn = f"{tsv_dir}/{subject}_{fn_spec}_{analysis_name}-css{dm_name}_violins.tsv"
                df_violins = pd.read_table(tsv_violins_fn, sep="\t")
                fig = prf_violins_plot(df=df_violins, figure_info=analysis_info, rsq2use=rsq2use)
                fig_fn = f"{fig_dir}/{subject}_{fn_spec}_{analysis_name}-css{dm_name}_violins.pdf"
                print('Saving pdf: {}'.format(fig_fn))
                fig.write_image(fig_fn)
                remove_second_page(fig_fn)
                
                # Parameters median plot
                tsv_params_median_fn = f"{tsv_dir}/{subject}_{fn_spec}_{analysis_name}-css{dm_name}_params-median.tsv"
                df_params_median = pd.read_table(tsv_params_median_fn, sep="\t")
                fig = prf_params_median_plot(df=df_params_median, figure_info=analysis_info, rsq2use=rsq2use)
                fig_fn = f"{fig_dir}/{subject}_{fn_spec}_{analysis_name}-css{dm_name}_params-median.pdf"
                print('Saving pdf: {}'.format(fig_fn))
                fig.write_image(fig_fn)
                remove_second_page(fig_fn)
    
                # Ecc.size plots
                tsv_ecc_size_fn = f"{tsv_dir}/{subject}_{fn_spec}_{analysis_name}-css{dm_name}_ecc-size.tsv"
                df_ecc_size = pd.read_table(tsv_ecc_size_fn, sep="\t")
                fig = prf_ecc_size_plot(df=df_ecc_size, figure_info=analysis_info, rsq2use=rsq2use)
                fig_fn = f"{fig_dir}/{subject}_{fn_spec}_{analysis_name}-css{dm_name}_ecc-size.pdf"
                print('Saving pdf: {}'.format(fig_fn))
                fig.write_image(fig_fn)
                remove_second_page(fig_fn)
                
                # Ecc.pCM plot
                tsv_ecc_pcm_fn = f"{tsv_dir}/{subject}_{fn_spec}_{analysis_name}-css{dm_name}_ecc-pcm.tsv"
                df_ecc_pcm = pd.read_table(tsv_ecc_pcm_fn, sep="\t")
                fig_fn = f"{fig_dir}/{subject}_{fn_spec}_{analysis_name}-css{dm_name}_ecc-pcm.pdf"
                fig = prf_ecc_pcm_plot(df=df_ecc_pcm, figure_info=analysis_info, rsq2use=rsq2use)
                print('Saving pdf: {}'.format(fig_fn))
                fig.write_image(fig_fn)
                remove_second_page(fig_fn)
    
                # Polar angle distributions
                tsv_polar_angle_fn = f"{tsv_dir}/{subject}_{fn_spec}_{analysis_name}-css{dm_name}_polar-angle.tsv"
                df_polar_angle = pd.read_table(tsv_polar_angle_fn, sep="\t")
                figs, hemis = prf_polar_angle_plot(df=df_polar_angle, figure_info=analysis_info)
                for (fig, hemi) in zip(figs, hemis):
                    if hemi == 'hemi-LR':
                        fig_fn = f"{fig_dir}/{subject}_{fn_spec}_{analysis_name}-css{dm_name}_polar-angle.pdf"
                        print('Saving pdf: {}'.format(fig_fn))
                        fig.write_image(fig_fn)
                        remove_second_page(fig_fn)
    
                # Contralaterality plots
                tsv_contralaterality_fn = f"{tsv_dir}/{subject}_{fn_spec}_{analysis_name}-css{dm_name}_contralaterality.tsv"
                df_contralaterality = pd.read_table(tsv_contralaterality_fn, sep="\t")
                fig_fn = f"{fig_dir}/{subject}_{fn_spec}_{analysis_name}-css{dm_name}_contralaterality.pdf"
                fig = prf_contralaterality_plot(df=df_contralaterality, figure_info=analysis_info)
                print('Saving pdf: {}'.format(fig_fn))
                fig.write_image(fig_fn)
                remove_second_page(fig_fn)
                
                # Spatial distribution plot
                tsv_distribution_fn = f"{tsv_dir}/{subject}_{fn_spec}_{analysis_name}-css{dm_name}_distribution.tsv"
                df_distribution = pd.read_table(tsv_distribution_fn, sep="\t")
                figs, hemis = prf_distribution_plot(df=df_distribution, figure_info=analysis_info)
                for (fig, hemi) in zip(figs, hemis):
                    if hemi == 'hemi-LR':
                        fig_fn = f"{fig_dir}/{subject}_{fn_spec}_{analysis_name}-css{dm_name}_distribution.pdf"
                        print('Saving pdf: {}'.format(fig_fn))
                        fig.write_image(fig_fn)
                        remove_second_page(fig_fn)
                
                # Spatial distibution barycentre plot
                tsv_barycentre_fn = f"{tsv_dir}/{subject}_{fn_spec}_{analysis_name}-css{dm_name}_barycentre.tsv"
                try:
                    df_barycentre = pd.read_table(tsv_barycentre_fn, sep="\t")
                except FileNotFoundError:
                    print(f"File not found: {tsv_barycentre_fn}, skipping...")
                    continue
                fig_fn = f"{fig_dir}/{subject}_{fn_spec}_{analysis_name}-css{dm_name}_barycentre.pdf"
                fig = prf_barycentre_plot(df=df_barycentre, figure_info=analysis_info)
                print('Saving pdf: {}'.format(fig_fn))
                fig.write_image(fig_fn)
                remove_second_page(fig_fn)
    
# Define permission cmd
print('Changing files permissions in {}/{}'.format(main_dir, project_dir))
os.system("chmod -Rf 771 {}/{}".format(main_dir, project_dir))
os.system("chgrp -Rf {} {}/{}".format(group, main_dir, project_dir))