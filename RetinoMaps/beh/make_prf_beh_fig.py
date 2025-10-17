"""
-----------------------------------------------------------------------------------------
make_prf_beh_fig.py
-----------------------------------------------------------------------------------------
Goal of the script:
make pRF behaviour figures 
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject name (e.g. sub-01)
sys.argv[4]: server group (e.g. 327)
-----------------------------------------------------------------------------------------
Output(s):
figure for pRF behaviour analysis 
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/pRF_analysis/RetinoMaps/pRF_beh/
2. run python command
python make_rois_fig.py [main directory] [project name] [subject] [group]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/projects/pRF_analysis/RetinoMaps/beh/
python make_prf_beh_fig.py /scratch/mszinte/data RetinoMaps sub-01 327
python make_prf_beh_fig.py /scratch/mszinte/data RetinoMaps group 327
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
import numpy as np
import pandas as pd
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
with open('../settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)
if 'group' not in subject: subjects=[subject]
else: subjects = analysis_info['subjects']
subject_excluded = analysis_info['outlier_prf_beh']

# Load figure settings
with open('../figure_settings.json') as f:
    json_s = f.read()
    figure_info = json.loads(json_s)
fig_width = figure_info['fig_width']

# Template settings
template_specs = dict(axes_color="rgba(0, 0, 0, 1)",
                      axes_width=2,
                      axes_font_size=15,
                      bg_col="rgba(255, 255, 255, 1)",
                      font='Arial',
                      title_font_size=15,
                      plot_width=1.5)
fig_template = plotly_template(template_specs)

# color_map
colormap_subject_dict = {
    'sub-01': '#AA0DFE',
    'sub-02': '#3283FE',
    'sub-03': '#85660D',
    'sub-04': '#782AB6',
    'sub-05': '#565656',
    'sub-06': '#1C8356',
    'sub-07': '#16FF32',
    'sub-08': '#F7E1A0',
    'sub-09': '#E2E2E2',
    'sub-11': '#1CBE4F',
    'sub-12': '#C4451C',
    'sub-13': '#DEA0FD',
    'sub-14': '#FBE426',
    'sub-15': '#325A9B',
    'sub-16': '#FEAF16',
    'sub-17': '#F8A19F',
    'sub-18': '#90AD1C',
    'sub-20': '#F6222E',
    'sub-21': '#1CFFCE',
    'sub-22': '#2ED9FF',
    'sub-23': '#B10DA1',
    'sub-24': '#C075A6',
    'sub-25': '#FC1CBF'
}

# settinhs 
step = 128
x_max = 640

# import data 
pRF_beh_tsv_dir = '{}/{}/derivatives/pp_data/{}/pRF_beh/tsv'.format(
    main_dir, project_dir, subject)
pRF_beh_trials_fn = '{}/{}_pRF_beh_trials.tsv'.format(
    pRF_beh_tsv_dir, subject)
pRF_beh_trials = pd.read_table(pRF_beh_trials_fn)

# Defind sublot 
# fig = make_subplots(rows=2, cols=1, vertical_spacing=0.20, print_grid=False)
fig = make_subplots(rows=2, cols=2, column_widths=[0.75, 0.25], vertical_spacing=0.2, print_grid=False)
print('{} is processing...'.format(subject))
for n_subject, subject_to_plot in enumerate(subjects):
    pRF_beh_trials_subject = pRF_beh_trials.loc[pRF_beh_trials['subject']==subject_to_plot]

    # Perf 
    # plot perf across trials 
    fig.add_trace(go.Scatter(y=pRF_beh_trials_subject.perf, 
                             line=dict(color=colormap_subject_dict[subject_to_plot]), 
                             name=subject_to_plot, 
                             showlegend=False), 
                  row=1, col=1)

    # Staircase 
    # plot perf across trials 
    fig.add_trace(go.Scatter(y=pRF_beh_trials_subject.stair, 
                             line=dict(color=colormap_subject_dict[subject_to_plot]), 
                             name=subject_to_plot), 
                  row=2, col=1)

for i in range(0, x_max, step):
    fig.add_shape(type="rect", 
                  x0=i, 
                  x1=i + step, 
                  y0=0, 
                  y1=15, 
                  fillcolor="gray" if (i // step) % 2 == 0 else "lightgray", 
                  opacity=0.2,
                  layer="below", 
                  line_width=0,
                  row="all",
                  col=1) 
# Update axes
# perf 
fig.update_xaxes(showline=True, 
                 range=[0,x_max], 
                 tickvals=[0, 128, 256, 384, 512, 640], 
                 ticktext=["0", "128", "256", "384", "512", "640"], 
                 row=1, col=1)
fig.update_yaxes(showline=True, 
                 range=[0,1], 
                 nticks=6, 
                 title='Perf (%)', 
                 row=1, col=1)
# stair 
fig.update_xaxes(showline=True, range=[0,x_max], 
                 tickvals=[0, 128, 256, 384, 512, 640], 
                 ticktext=["0", "128", "256", "384", "512", "640"], 
                 title='Trials', 
                 row=2, col=1)
fig.update_yaxes(showline=True, range=[0,10], title='Kappa (a.u.)', row=2, col=1)
# Update layout
fig.update_layout(template=fig_template, 
                  width=fig_width, 
                  height=800, 
                  showlegend=True, 
                  legend=dict(orientation="v", 
                              yanchor="top", 
                              y=1, 
                              xanchor="left", 
                              x=1.05))
# Export figure
prf_beh_fig_dir = '{}/{}/derivatives/pp_data/{}/pRF_beh/figures'.format(
    main_dir, project_dir, subject)
os.makedirs(prf_beh_fig_dir, exist_ok=True)

prf_beh_fig_fn = '{}/{}_prf_beh_trials.pdf'.format(prf_beh_fig_dir, subject)
print('Export {}'.format(prf_beh_fig_fn))
fig.write_image(prf_beh_fig_fn)

if 'group' in subject:
    # import data 
    subject_excluded = analysis_info['outlier_prf_beh']
    subject = 'group'
    pRF_beh_tsv_dir = '{}/{}/derivatives/pp_data/{}/pRF_beh/tsv'.format(main_dir, project_dir, subject)
    pRF_beh_trials_group_median_fn = '{}/{}_pRF_beh_trials_group_median.tsv'.format(pRF_beh_tsv_dir, subject)
    pRF_beh_median_group_median_fn = '{}/{}_pRF_beh_median_group_median.tsv'.format(pRF_beh_tsv_dir, subject)
    
    pRF_beh_trials_group_median = pd.read_table(pRF_beh_trials_group_median_fn)
    pRF_beh_median_group_median = pd.read_table(pRF_beh_median_group_median_fn)
    # Filter data to get group data 
    pRF_beh_trials_group = pRF_beh_trials_group_median.loc[pRF_beh_trials_group_median['subject']=='group']
    pRF_beh_median_group = pRF_beh_median_group_median.loc[pRF_beh_median_group_median['subject']=='group']
    
    # Defind sublot 
    fig = make_subplots(rows=2, cols=2, column_widths=[0.75, 0.25], vertical_spacing=0.2, print_grid=False)
    
    
    # plot perf across trials 
    fig.add_trace(go.Scatter(y=pRF_beh_trials_group['median_perf'], 
                             line=dict(color='black'), 
                             name='group', 
                             showlegend=True), 
                  row=1, col=1)
    
    # Error area
    fig.add_trace(go.Scatter(y=np.concatenate([pRF_beh_trials_group['ci_up_perf'], pRF_beh_trials_group['ci_low_perf'][::-1]]), 
                             x=np.concatenate([np.arange(len(pRF_beh_trials_group['ci_up_perf'])), np.arange(len(pRF_beh_trials_group['ci_low_perf']))[::-1]]), 
                             mode='lines', 
                             fill='toself', 
                             fillcolor='rgba(169, 169, 169, 0.4)', 
                             line=dict(width=0), 
                             showlegend=False), 
                  row=1, col=1)
    
    # make the bar plot 
    fig.add_trace(go.Bar(y=pRF_beh_median_group['perf_median'], 
                         marker=dict(color='black'),
                         error_y=dict(type='data', 
                                      symmetric=False, 
                                      array=pRF_beh_median_group['perf_ci_up'] - pRF_beh_median_group['perf_median'], 
                                      arrayminus=pRF_beh_median_group['perf_median'] - pRF_beh_median_group['perf_ci_low'], 
                                      thickness=5, 
                                      width=0,
                                      color='gray'), 
                         width=0.5, 
                         showlegend=False),
    
                  row=1, col=2)
    
    # plot stair across trials 
    fig.add_trace(go.Scatter(y=pRF_beh_trials_group['median_stair'], 
                             line=dict(color='black'), 
                             name='group', 
                             showlegend=True), 
                  row=2, col=1)
    
    # Error area
    fig.add_trace(go.Scatter(y=np.concatenate([pRF_beh_trials_group['ci_up_stair'], pRF_beh_trials_group['ci_low_stair'][::-1]]), 
                             x=np.concatenate([np.arange(len(pRF_beh_trials_group['ci_up_stair'])), np.arange(len(pRF_beh_trials_group['ci_low_stair']))[::-1]]), 
                             mode='lines', 
                             fill='toself', 
                             fillcolor='rgba(169, 169, 169, 0.4)', 
                             line=dict(width=0), 
                             showlegend=False), 
                  row=2, col=1)
    
    # make the bar plot 
    fig.add_trace(go.Bar(y=pRF_beh_median_group['stair_median'], 
                         marker=dict(color='black'),
                         error_y=dict(type='data', 
                                      symmetric=False, 
                                      array=pRF_beh_median_group['stair_ci_up'] - pRF_beh_median_group['stair_median'], 
                                      arrayminus=pRF_beh_median_group['stair_median'] - pRF_beh_median_group['stair_ci_low'], 
                                      thickness=5, 
                                      width=0,
                                      color='gray'), 
                         width=0.5, 
                         showlegend=False),
    
                  row=2, col=2)
    
    # Plot Outliers
    for n_outlier, outlier in enumerate(subject_excluded):
        # showlegend = True if n_outlier == 0 else False
        pRF_beh_trials_outliers_indiv = pRF_beh_trials_group_median.loc[pRF_beh_trials_group_median['subject']==outlier]
        pRF_beh_median_outliers_indiv = pRF_beh_median_group_median.loc[pRF_beh_median_group_median['subject']==outlier]
        
        # plot perf for outliers 
        fig.add_trace(go.Scatter(y=pRF_beh_trials_outliers_indiv['median_perf'], 
                                 line=dict(color=colormap_subject_dict[outlier], 
                                           dash='dot'), 
                                 opacity=0.5, 
                                 name=outlier, 
                                 showlegend=True),  
                      row=1, col=1)
    
        # plot stair for outliers 
        fig.add_trace(go.Scatter(y=pRF_beh_trials_outliers_indiv['median_stair'], 
                                 line=dict(color=colormap_subject_dict[outlier], 
                                           dash='dot'), 
                                 opacity=0.5, 
                                 name=outlier, 
                                 showlegend=True),  
                      row=2, col=1)
        
        # add point in boxplot 
        # Perf
        fig.add_trace(go.Scatter(x=[0, 0],
                                 y=pRF_beh_median_outliers_indiv['perf_median'], 
                                 mode='markers', 
                                 marker=dict(color='rgba(255, 255, 255, 0)',  
                                             line=dict(color=colormap_subject_dict[outlier], width=2), 
                                             size=8),  
                                 showlegend=False),
                      row=1, col=2)
        
        # Stair
        fig.add_trace(go.Scatter(x=[0, 0],
                                 y=pRF_beh_median_outliers_indiv['stair_median'], 
                                 mode='markers', 
                                 marker=dict(color='rgba(255, 255, 255, 0)',  
                                             line=dict(color=colormap_subject_dict[outlier], width=2),  
                                             size=8),  
                                 showlegend=False),
                      row=2, col=2)
    
    # add runs square
    for i in range(0, x_max, step):
        fig.add_shape(type="rect", 
                      x0=i, 
                      x1=i + step, 
                      y0=0, 
                      y1=15, 
                      fillcolor="gray" if (i // step) % 2 == 0 else "lightgray", 
                      opacity=0.2,
                      layer="below", 
                      line_width=0,
                      row="all",
                      col=1) 
    # Update axes
    # perf 
    # left plot
    fig.update_xaxes(showline=True, 
                     range=[0,x_max], 
                     tickvals=[0, 128, 256, 384, 512, 640], 
                     ticktext=["0", "128", "256", "384", "512", "640"], 
                     row=1, col=1)
    fig.update_yaxes(showline=True, 
                     range=[0,1], 
                     nticks=6, 
                     title='Perf (%)', 
                     row=1, col=1)
    # right plot
    fig.update_xaxes(showline=False, 
                     ticks='', 
                     showticklabels=False, 
                     title='Median performances', 
                     row=1, col=2)
    fig.update_yaxes(showline=True,range=[0,1], 
                     nticks=6, 
                     row=1, col=2)
    
    # stair 
    # left plot
    fig.update_xaxes(showline=True, 
                     range=[0,x_max], 
                     tickvals=[0, 128, 256, 384, 512, 640], 
                     ticktext=["0", "128", "256", "384", "512", "640"], 
                     title='Trials', 
                     row=2, col=1)
    fig.update_yaxes(showline=True, 
                     range=[0,15], 
                     title='Kappa (a.u.)', 
                     row=2, col=1)
    # right plot
    fig.update_xaxes(showline=False, 
                     ticks='', 
                     showticklabels=False, 
                     title='Median kappa', 
                     row=2, col=2)
    fig.update_yaxes(showline=True, 
                     range=[0,15], 
                     row=2, col=2)
    
    # Update layout
    fig.update_layout(template=fig_template, 
                      width=fig_width, 
                      height=800, 
                      showlegend=True)

    # Export figure
    prf_beh_fig_fn = '{}/{}_prf_beh.pdf'.format(prf_beh_fig_dir, subject)
    print('Export {}'.format(prf_beh_fig_fn))
    fig.write_image(prf_beh_fig_fn)

# Define permission cmd
print('Changing files permissions in {}/{}'.format(main_dir, project_dir))
os.system("chmod -Rf 771 {}/{}".format(main_dir, project_dir))
os.system("chgrp -Rf {} {}/{}".format(group, main_dir, project_dir))    
    
