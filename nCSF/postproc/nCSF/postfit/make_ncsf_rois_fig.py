"""
-----------------------------------------------------------------------------------------
make_ncsf_rois_fig.py
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
nCSF analysis figures
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/pRF_analysis/nCSF/postproc/nCSF/postfit/
2. run python command
python make_ncsf_rois_fig.py [main directory] [project name] [subject] [group] 
-----------------------------------------------------------------------------------------
Exemple:
d ~/projects/pRF_analysis/nCSF/postproc/nCSF/postfit/
python make_ncsf_rois_fig.py /scratch/mszinte/data nCSF sub-01 327
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
sys.path.append("{}/../../../../analysis_code/utils".format(os.getcwd()))
from settings_utils import load_settings

sys.path.append("{}/../../../utils".format(os.getcwd()))
from nCSF_plot_utils import *

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
group = sys.argv[4]

# Load settings
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../../"))
settings_path = os.path.join(base_dir, project_dir, "settings.yml")
prf_settings_path = os.path.join(base_dir, project_dir, "prf-analysis.yml")
ncsf_settings_path = os.path.join(base_dir, project_dir, "nCSF-analysis.yml")
figure_settings_path = os.path.join(base_dir, project_dir, "figure-settings.yml")
settings = load_settings([settings_path, prf_settings_path, figure_settings_path, ncsf_settings_path])
analysis_info = settings[0]

# formats = analysis_info['formats']
# extensions = analysis_info['extensions']


formats = [analysis_info['formats'][0]] #######################################
extensions = [analysis_info['extensions'][0]]

subjects = analysis_info['subjects']
preproc_prep = analysis_info['preproc_prep']
filtering = analysis_info['filtering']
normalization = analysis_info['normalization']
avg_methods = analysis_info['avg_methods']
ncsf_task_name = analysis_info['nCSF_task_name']
rois_methods = analysis_info['rois_methods']

# Main loop
for avg_method in avg_methods:
    if 'loo' in avg_method: rsq2use = 'ncsf_loo_rsq'
    else: rsq2use = 'ncsf_rsq'
    rsq2use_median_task = 'ncsf_ncsf_median_{}'.format(rsq2use)

    for format_, extension in zip(formats, extensions):
        # define list of rois for each format
        rois_methods_format = rois_methods[format_]
        for rois_method_format in rois_methods_format:
            if rois_method_format == 'rois-drawn':
                analysis_info['rois'] = analysis_info[rois_method_format]
            elif rois_method_format == 'rois-group-mmp':
                analysis_info['rois'] = list(analysis_info[rois_method_format].keys())

          
                
            print(f'{ncsf_task_name} - {avg_method} - {format_}')

            # Define/create folders
            tsv_dir = '{}/{}/derivatives/pp_data/{}/{}/ncsf/tsv'.format(
                main_dir, project_dir, subject, format_)
            
            # Exception if no data for one format (e.g template subject)
            if not os.path.isdir(tsv_dir):
                print(f"[SKIP] tsv_dir not found for format={format_}: {tsv_dir}")
                continue
            
            fig_dir = '{}/{}/derivatives/pp_data/{}/{}/ncsf/figures'.format(
                main_dir, project_dir, subject, format_)
            os.makedirs(fig_dir, exist_ok=True)

            fn_spec = "task-{}_{}_{}_{}_{}_{}".format(
                ncsf_task_name, preproc_prep, filtering, 
                normalization, avg_method, rois_method_format)

            # Roi active vertex
            tsv_roi_active_vert_fn = "{}/{}_{}_ncsf_active-vert.tsv".format(tsv_dir, subject, fn_spec)
            df_roi_active_vert = pd.read_table(tsv_roi_active_vert_fn, sep="\t")
            fig = nCSF_roi_active_vert_plot(df=df_roi_active_vert, figure_info=analysis_info, format=format_, )
            fig_fn = "{}/{}_{}_ncsf_active-vert.pdf".format(fig_dir, subject, fn_spec)
            print('Saving pdf: {}'.format(fig_fn))
            fig.write_image(fig_fn)
            remove_second_page(fig_fn)
            
            # Violins plot
            tsv_violins_fn = "{}/{}_{}_ncsf_violins.tsv".format(tsv_dir, subject, fn_spec)
            df_violins = pd.read_table(tsv_violins_fn, sep="\t")
            fig = ncsf_violins_plot(df=df_violins, figure_info=analysis_info, rsq2use=rsq2use)
            fig_fn = "{}/{}_{}_ncsf_violins.pdf".format(fig_dir, subject, fn_spec)
            print('Saving pdf: {}'.format(fig_fn))
            fig.write_image(fig_fn)
            remove_second_page(fig_fn)
            
            # Parameters median plot
            tsv_params_median_fn = "{}/{}_{}_ncsf_params-median.tsv".format(tsv_dir, subject, fn_spec)
            df_params_median = pd.read_table(tsv_params_median_fn, sep="\t")
            fig = ncsf_params_median_plot(df=df_params_median, figure_info=analysis_info, rsq2use=rsq2use)
            fig_fn = "{}/{}_{}_ncsf_params-median.pdf".format(fig_dir, subject, fn_spec)
            print('Saving pdf: {}'.format(fig_fn))
            fig.write_image(fig_fn)
            remove_second_page(fig_fn)

            # Ecc.SFp
            tsv_ecc_SFp_fn = "{}/{}_{}_ncsf_ecc-SFp.tsv".format(tsv_dir, subject, fn_spec)
            df_ecc_SFp = pd.read_table(tsv_ecc_SFp_fn, sep="\t")
            fig = ecc_SFp_plot(df=df_ecc_SFp, figure_info=analysis_info, rsq2use=rsq2use_median_task)
            fig_fn = "{}/{}_{}_ecc-SFp.pdf".format(fig_dir, subject, fn_spec)
            print('Saving pdf: {}'.format(fig_fn))
            fig.write_image(fig_fn)
            remove_second_page(fig_fn)
            
            # Ecc.auc
            tsv_ecc_auc_fn = "{}/{}_{}_ncsf_ecc-auc.tsv".format(tsv_dir, subject, fn_spec)
            df_ecc_auc = pd.read_table(tsv_ecc_auc_fn, sep="\t")
            fig = ecc_SFp_plot(df=df_ecc_SFp, figure_info=analysis_info, rsq2use=rsq2use_median_task)
            fig_fn = "{}/{}_{}_ecc-auc.pdf".format(fig_dir, subject, fn_spec)
            print('Saving pdf: {}'.format(fig_fn))
            fig.write_image(fig_fn)
            remove_second_page(fig_fn)
            
            # Size.SFp
            tsv_size_sfp_fn = "{}/{}_{}_ncsf_size-SFp.tsv".format(tsv_dir, subject, fn_spec)
            df_size_sfp = pd.read_table(tsv_size_sfp_fn, sep="\t")
            fig = size_SFp_plot(df=df_size_sfp, figure_info=analysis_info, rsq2use=rsq2use_median_task)
            fig_fn = "{}/{}_{}_size-SFp.pdf".format(fig_dir, subject, fn_spec)
            print('Saving pdf: {}'.format(fig_fn))
            fig.write_image(fig_fn)
            remove_second_page(fig_fn)
            
            # Size.auc
            tsv_size_auc_fn = "{}/{}_{}_ncsf_size-auc.tsv".format(tsv_dir, subject, fn_spec)
            df_size_auc = pd.read_table(tsv_size_auc_fn, sep="\t")
            fig = size_auc_plot(df=df_size_auc, figure_info=analysis_info, rsq2use=rsq2use_median_task)
            fig_fn = "{}/{}_{}_size-auc.pdf".format(fig_dir, subject, fn_spec)
            print('Saving pdf: {}'.format(fig_fn))
            fig.write_image(fig_fn)
            remove_second_page(fig_fn)
            


# Define permission cmd
print('Changing files permissions in {}/{}'.format(main_dir, project_dir))
os.system("chmod -Rf 771 {}/{}".format(main_dir, project_dir))
os.system("chgrp -Rf {} {}/{}".format(group, main_dir, project_dir))