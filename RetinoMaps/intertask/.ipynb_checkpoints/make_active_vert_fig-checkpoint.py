"""
-----------------------------------------------------------------------------------------
make_active_vert_fig.py
-----------------------------------------------------------------------------------------
Goal of the script:
Make tsv for active vertex figures
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject name (e.g. sub-01)
sys.argv[4]: server group (e.g. 327)
-----------------------------------------------------------------------------------------
Output(s):
TSV for active vertex figures
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/pRF_analysis/RetinoMaps/intertask/
2. run python command
python make_rois_fig.py [main directory] [project name] [subject] [group]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/projects/pRF_analysis/RetinoMaps/intertask/
python make_active_vert_fig.py /scratch/mszinte/data RetinoMaps sub-01 327
python make_active_vert_fig.py /scratch/mszinte/data RetinoMaps sub-170k 327
python make_active_vert_fig.py /scratch/mszinte/data RetinoMaps group 327
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

# Imports
import os
import sys
import json
import pandas as pd

# Plotly imports
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# Personal import
sys.path.append("{}/../../analysis_code/utils".format(os.getcwd()))
from plot_utils import plotly_template

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
subjects = analysis_info['subjects']
group_tasks = analysis_info['task_intertask']
categories = ['saccade', 'pursuit', 'vision', 'all']


# Threshold settings
ecc_th = analysis_info['ecc_th']
size_th = analysis_info['size_th']
rsqr_th = analysis_info['rsqr_th']
amplitude_th = analysis_info['amplitude_th']
stats_th = analysis_info['stats_th']
n_th = analysis_info['n_th'] 

# Figure settings
with open('../figure_settings.json') as f:
    json_s = f.read()
    figure_info = json.loads(json_s)
    
plot_groups = figure_info['plot_groups']
fig_width = figure_info['fig_width']

roi_colors_dict = {'V1': 'rgb(243, 231, 155)', 
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

# Template settings
template_specs = dict(axes_color="rgba(0, 0, 0, 1)",
                      axes_width=2,
                      axes_font_size=15,
                      bg_col="rgba(255, 255, 255, 1)",
                      font='Arial',
                      title_font_size=15,
                      plot_width=1.5)
fig_template = plotly_template(template_specs)
standoff = 8

for tasks in group_tasks : 
    if 'SacVELoc' in tasks: suffix = 'SacVE_PurVE'
    else : suffix = 'Sac_Pur'
    # Format loop
    for format_, extension in zip(formats, extensions):
        fig_dir = '{}/{}/derivatives/pp_data/{}/{}/intertask/figures'.format(
            main_dir, project_dir, subject, format_)
        os.makedirs(fig_dir, exist_ok=True)
        tsv_dir = '{}/{}/derivatives/pp_data/{}/{}/intertask/tsv'.format(
            main_dir, project_dir, subject, format_)
        
        # Load TSV 
        active_vertex_roi_tsv_fn = '{}/{}_active_vertex_roi_{}.tsv'.format(tsv_dir, subject, suffix)
        df_active_vertex_roi = pd.read_table(active_vertex_roi_tsv_fn)
        
        # Individual subject analysis
        if 'group' not in subject:
            # Active vertex roi figure
            # ------------------------
            rows, cols = 1, len(plot_groups)
            fig = make_subplots(rows=rows, cols=cols, print_grid=False)
            for l, line_label in enumerate(plot_groups):
                for j, roi in enumerate(line_label):
                    df_roi = df_active_vertex_roi.loc[df_active_vertex_roi['roi']==roi]
                    fig.add_trace(go.Bar(x=df_roi['categorie'], 
                                         y=df_roi['percentage_active'], 
                                         name=roi,  
                                         marker=dict(color=roi_colors_dict[roi]), 
                                         showlegend=True), 
                                  row=1, col=l + 1)
            
            # Set axes
            fig.update_xaxes(showline=True)
            fig.update_yaxes(title='Active vertex (%)', 
                             range=[0,100], 
                             showline=True)
            
            # Update layout of the figure
            fig.update_layout(template=fig_template, 
                              barmode='group', 
                              height=400, 
                              width=fig_width
                             )
            # Export figure
            fig_fn = "{}/{}_active_vertex_roi_{}.pdf".format(fig_dir, subject, suffix)
            print('Saving pdf: {}'.format(fig_fn))
            fig.write_image(fig_fn)
            
            # Active vertex roi mmp figure
            # ----------------------------
            if format_ == '170k' :
                # Load DF 
                active_vertex_roi_mmp_tsv_fn = '{}/{}_active_vertex_roi_mmp_{}.tsv'.format(tsv_dir, subject, suffix)
                df_active_vertex_roi_mmp = pd.read_table(active_vertex_roi_mmp_tsv_fn)
                
                rows, cols = len(plot_groups), 3 
                for n_cat, categorie in enumerate(categories):
                    fig = make_subplots(rows=rows, 
                                        cols=cols,  
                                        horizontal_spacing=0.15, 
                                        vertical_spacing=0.15,
                                        print_grid=False
                                       )
                    
                    for row_idx, line_label in enumerate(plot_groups):  
                        for col_idx, roi in enumerate(line_label):      
                            df_roi = df_active_vertex_roi_mmp.loc[(df_active_vertex_roi_mmp['roi'] == roi) & 
                                                                  (df_active_vertex_roi_mmp['categorie'] == categorie)].sort_values(by='percentage_active', 
                                                                                                                                    ascending=True)
                
                            rois_mmp = df_roi['roi_mmp'].unique()
                            for n_roi_mmp, roi_mmp in enumerate(rois_mmp):
                                showlegend = (n_roi_mmp == 0)
                                df_roi_mmp = df_roi[df_roi['roi_mmp'] == roi_mmp]
                
                                fig.add_trace(
                                    go.Bar(x=df_roi_mmp['percentage_active'], 
                                           y=df_roi_mmp['roi_mmp'], 
                                           orientation='h', 
                                           name=roi, 
                                           marker=dict(color=roi_colors_dict[roi]),  
                                           width=0.9, 
                                           showlegend=showlegend), 
                                    row=row_idx + 1, 
                                    col=col_idx + 1
                                )
                            
                            # Set axes
                            fig.update_xaxes(title=dict(text='Active vertex (%)', standoff=standoff), 
                                             range=[0, 100], 
                                             tickvals=[25, 50, 75, 100],  
                                             ticktext=[str(val) for val in [25, 50, 75, 100]], 
                                             showline=True, 
                                             row=row_idx + 1, 
                                             col=col_idx + 1
                                            )
                
                            y_title = 'Glasser parcellation' if col_idx == 0 else ''
                            fig.update_yaxes(title=dict(text=y_title, standoff=standoff), 
                                             showline=True, 
                                             row=row_idx + 1, 
                                             col=col_idx + 1)
                    
                    # Update layout of the figure
                    fig.update_layout(title='{}'.format(categorie), 
                                      template=fig_template, 
                                      height=1080, 
                                      width=1080, 
                                      margin_l=100, 
                                      margin_r=100, 
                                      margin_t=100, 
                                      margin_b=100)
                    # Export figure
                    fig_fn = "{}/{}_active_vertex_{}_roi_mmp_{}.pdf".format(fig_dir, subject,categorie, suffix)
                    print('Saving pdf: {}'.format(fig_fn))
                    fig.write_image(fig_fn)
                                        
        # Group Analysis    
        else :
            # Active vertex roi figure
            # ------------------------
            rows, cols = 1, len(plot_groups)
            fig = make_subplots(rows=rows, cols=cols, print_grid=False)
            for l, line_label in enumerate(plot_groups):
                for j, roi in enumerate(line_label):
                    df_roi = df_active_vertex_roi.loc[df_active_vertex_roi['roi']==roi]
                    fig.add_trace(go.Bar(x=df_roi['categorie'], 
                                         y=df_roi['median'], 
                                         name=roi,  
                                         marker=dict(color=roi_colors_dict[roi]), 
                                         error_y=dict(type='data', 
                                                      array=df_roi['ci_high'] - df_roi['median'], 
                                                      arrayminus=df_roi['median'] - df_roi['ci_low'], 
                                                      visible=True, 
                                                      width=0, 
                                                      color='black'), 
                                         showlegend=True), 
                                  row=1, 
                                  col=l + 1)
            
            # Set axes 
            fig.update_xaxes(showline=True)
            fig.update_yaxes(title='Active vertex (%)', 
                             range=[0,100], 
                             showline=True)
            
            # Update layout of the figure
            fig.update_layout(template=fig_template, 
                              barmode='group', 
                              height=400, 
                              width=fig_width
                             )
            # Export figure
            fig_fn = "{}/{}_active_vertex_roi_{}.pdf".format(fig_dir, subject, suffix)
            print('Saving pdf: {}'.format(fig_fn))
            fig.write_image(fig_fn)
            
            # Active vertex roi mmp figure
            # ----------------------------
            if format_ == '170k' :
                # Load DF 
                active_vertex_roi_mmp_tsv_fn = '{}/{}_active_vertex_roi_mmp_{}.tsv'.format(tsv_dir, subject, suffix)
                df_active_vertex_roi_mmp = pd.read_table(active_vertex_roi_mmp_tsv_fn)
                
                rows, cols = len(plot_groups), 3 
                for n_cat, categorie in enumerate(categories):
                    fig = make_subplots(rows=rows, 
                                        cols=cols,  
                                        horizontal_spacing=0.15, 
                                        vertical_spacing=0.15,
                                        print_grid=False
                                       )
                    
                    for row_idx, line_label in enumerate(plot_groups):  
                        for col_idx, roi in enumerate(line_label):      
                            df_roi = df_active_vertex_roi_mmp.loc[(df_active_vertex_roi_mmp['roi'] == roi) & 
                                                                  (df_active_vertex_roi_mmp['categorie'] == categorie)].sort_values(by='median', 
                                                                                                                                    ascending=True)
                
                            rois_mmp = df_roi['roi_mmp'].unique()
                            for n_roi_mmp, roi_mmp in enumerate(rois_mmp):
                                showlegend = (n_roi_mmp == 0)
                                df_roi_mmp = df_roi[df_roi['roi_mmp'] == roi_mmp]
                
                                fig.add_trace(go.Bar(x=df_roi_mmp['median'], 
                                                     y=df_roi_mmp['roi_mmp'], 
                                                     orientation='h', 
                                                     name=roi, 
                                                     marker=dict(color=roi_colors_dict[roi]), 
                                                     error_x=dict(type='data', 
                                                                  array=(df_roi_mmp['ci_high'] - df_roi_mmp['median']).values, 
                                                                  arrayminus=(df_roi_mmp['median'] - df_roi_mmp['ci_low']).values, 
                                                                  visible=True,  
                                                                  width=0, 
                                                                  color='black'), 
                                                     width=0.9, 
                                                     showlegend=showlegend), 
                                              row=row_idx + 1, 
                                              col=col_idx + 1
                                             )
                            
                            # Set axes 
                            fig.update_xaxes(title=dict(text='Active vertex (%)', standoff=standoff), 
                                             range=[0, 100], 
                                             tickvals=[25, 50, 75, 100],  
                                             ticktext=[str(val) for val in [25, 50, 75, 100]], 
                                             showline=True, 
                                             row=row_idx + 1, 
                                             col=col_idx + 1
                                            )
                
                
                            y_title = 'Glasser parcellation' if col_idx == 0 else ''
                            fig.update_yaxes(title=dict(text=y_title, standoff=standoff), 
                                             showline=True, 
                                             row=row_idx + 1, 
                                             col=col_idx + 1)
                    
                    # Update layout of the figure
                    fig.update_layout(title='{}'.format(categorie), 
                                      template=fig_template, 
                                      height=1080, 
                                      width=1080, 
                                      margin_l=100, 
                                      margin_r=100, 
                                      margin_t=100, 
                                      margin_b=100)
                    
                    # Export figure
                    fig_fn = "{}/{}_active_vertex_{}_roi_mmp_{}.pdf".format(fig_dir, subject,categorie, suffix)
                    print('Saving pdf: {}'.format(fig_fn))
                    fig.write_image(fig_fn)
        
# # Define permission cmd
# print('Changing files permissions in {}/{}'.format(main_dir, project_dir))
# os.system("chmod -Rf 771 {}/{}".format(main_dir, project_dir))
# os.system("chgrp -Rf {} {}/{}".format(group, main_dir, project_dir))