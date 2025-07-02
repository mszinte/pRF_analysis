"""
-----------------------------------------------------------------------------------------
stats_figures.py
-----------------------------------------------------------------------------------------
Goal of the script:
Create figures for all subjects together showing the percentage of amount of data of the 
euclidean error under each threshold (precision) as well as under one specific 
threshold 
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject name
sys.argv[4]: group of shared data (e.g. 327)
-----------------------------------------------------------------------------------------
Output(s):
{subject}_{task}_threshold_precision.pdf
group_{task}_threshold_precision.pdf
{subject}_{task}_threshold_ranking.pdf
group_{task}_threshold_ranking.pdf
{subject}_{task}_stats_figure.pdf
group_{task}_stats_figure.pdf
-----------------------------------------------------------------------------------------
To run:
cd ~/projects/pRF_analysis/RetinoMaps/eyetracking/
python stats_figures.py /scratch/mszinte/data RetinoMaps group 327
-----------------------------------------------------------------------------------------
Written by Sina Kling
Edited by Uriel Lascombes (uriel.lascombes@laposte.net)
-----------------------------------------------------------------------------------------
"""
# Stop warnings
import warnings
warnings.filterwarnings("ignore")

# Debug
import ipdb
deb = ipdb.set_trace

# Imports
import re
import os
import sys
import glob 
import json
import numpy as np
import pandas as pd
import plotly.graph_objects as go 

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
group = sys.argv[4]

# Personal imports
sys.path.append("{}/../../analysis_code/utils".format(os.getcwd()))
from plot_utils import plotly_template

# Load general settings
with open('../settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)
tasks = analysis_info['task_intertask'][0]
prf_task_name = analysis_info['prf_task_name']
subjects = analysis_info['subjects']


# General figure settings
template_specs = dict(axes_color="rgba(0, 0, 0, 1)",
                      axes_width=2,
                      axes_font_size=15,
                      bg_col="rgba(255, 255, 255, 1)",
                      font='Arial',
                      title_font_size=15,
                      plot_width=1.5)

fig_template = plotly_template(template_specs)

colormap_subject_dict = {'sub-01': '#AA0DFE', 
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
                         'sub-25': '#FC1CBF'}


for task in tasks :
    base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../"))
    task_settings_path = os.path.join(base_dir, project_dir, '{}_settings.json'.format(task))
    with open(task_settings_path) as f:
        settings = json.load(f)
    
    eye = settings['eye']
    num_run = settings['num_run']
    num_seq = settings['num_seq']
    seq_trs = settings['seq_trs']
    eye_mov_seq = settings['eye_mov_seq']
    trials_seq = settings['trials_seq']
    rads = settings['rads']
    pursuits_tr = np.arange(0,seq_trs,2)
    saccades_tr = np.arange(1,seq_trs,2)
    eyetracking_sampling = settings['eyetrack_sampling']
    screen_size = settings['screen_size']
    ppd = settings['ppd']
    
    threshold = settings['threshold']
    

       
    precision_data = {}
    for subject_to_group in subjects:
        # Defind directories
        eye_tracking_dir = '{}/{}/derivatives/pp_data/{}/eyetracking'.format(
            main_dir, project_dir, subject_to_group)
        
        # Lad data
        precision_summary = pd.read_csv("{}/stats/{}_task-{}_precision_summary.tsv".format(
            eye_tracking_dir, subject_to_group, task), delimiter='\t')
        precision_one_summary = pd.read_csv("{}/stats/{}_task-{}_precision_one_threshold_summary.tsv".format(
            eye_tracking_dir, subject_to_group, task), delimiter='\t')
        
        #Concat subjects
        precision_data[subject_to_group] = {
            "precision_mean": precision_summary["precision_mean"],
            "precision_one_thrs_mean": precision_one_summary["precision_one_thrs_mean"].item()
        }
        
    fig = go.Figure()
    sorted_precision_data = dict(sorted(precision_data.items(),
                                        key=lambda item: item[1]["precision_one_thrs_mean"],
                                        reverse=True))
    
    thresholds=np.linspace(0, 9, 100)
    for subject_to_group, data in sorted_precision_data.items():
        fig.add_trace(go.Scatter(x=thresholds, 
                                  y=data["precision_mean"], 
                                  mode='lines', 
                                  name='{}'.format(subject_to_group), 
                                  line=dict(color=colormap_subject_dict[subject_to_group])))
        
                
    fig.update_xaxes(showline=True, range=[0,6])
    fig.update_yaxes(showline=True, range=[0,1])
    
    fig.update_layout(template=fig_template, 
                      height=800, 
                      width=400)
                
                
    
    
    # Save figure 
    fig_dir = "{}/{}/derivatives/pp_data/{}/eyetracking/figures".format(main_dir, project_dir, subject)
    os.makedirs(fig_dir, exist_ok=True)
    fig_fn = "{}/{}_{}_threshold_precision.pdf".format(fig_dir, subject, task)
    
    print(f"Saving {fig_fn}")
    fig.write_image(fig_fn)
    
# # Define permission cmd
# print('Changing files permissions in {}/{}'.format(main_dir, project_dir))
# os.system("chmod -Rf 771 {}/{}".format(main_dir, project_dir))
# os.system("chgrp -Rf {} {}/{}".format(group, main_dir, project_dir))    
