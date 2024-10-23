
#%%
import numpy as np
import h5py
import json
import os
import plotly 
import sys
import cortex
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import os
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import pandas as pd
import matplotlib.pyplot as plt
import cortex
       

from sac_utils import plotly_layout_template

def ensure_save_dir(base_dir, subject):
    save_dir = f"{base_dir}/{subject}/eyetracking"
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    return save_dir


subjects = ['sub-01']
task = "SacLoc"


with open('/Users/sinakling/projects/pRF_analysis/RetinoMaps/eyetracking/dev/PurLoc_SacLoc/behavior_settings.json') as f:
    json_s = f.read()
    settings = json.loads(json_s)

runs = np.arange(0,settings['num_run'],1)
sequences = np.arange(0,settings['num_seq'],1)
trials_seq = settings['trials_seq']
rads = settings['rads']
eye_mov_seq = settings['eye_mov_seq']
seq_trs = settings['seq_trs']

for subject in subjects: 

    # Load data
    # ---------

    main_dir = settings['main_dir_mac']
    file_dir_save = ensure_save_dir(f'{main_dir}/derivatives/pp_data', subject)
    fig_dir_save = f'{file_dir_save}/figures'
    if not os.path.exists(fig_dir_save):
        os.makedirs(fig_dir_save)
    h5_filename = '{file_dir}/stats/{sub}_task-{task}_eyedata_sac_stats.h5'.format(file_dir = file_dir_save, sub = subject, task = task)
    h5_file = h5py.File(h5_filename,'r')
    time_start_seq = np.array(h5_file['time_start_seq'])
    time_end_seq = np.array(h5_file['time_end_seq'])
    time_start_trial = np.array(h5_file['time_start_trial'])
    time_end_trial = np.array(h5_file['time_end_trial'])
    amp_sequence = np.array(h5_file['amp_sequence'])
    saccades_output = np.array(h5_file['saccades_output'])

    

    eye_data_run_1 = pd.read_csv(f"{file_dir_save}/timeseries/{subject}_task-{task}_run_01_eyedata.tsv.gz", compression='gzip', delimiter='\t')
    eye_data_run_1 = eye_data_run_1[['timestamp', 'x', 'y', 'pupil_size']].to_numpy()
    eye_data_run_2 = pd.read_csv(f"{file_dir_save}/timeseries/{subject}_task-{task}_run_01_eyedata.tsv.gz", compression='gzip', delimiter='\t')
    eye_data_run_2 = eye_data_run_2[['timestamp', 'x', 'y', 'pupil_size']].to_numpy()
    eye_data_all_runs = [eye_data_run_1,eye_data_run_2]

    if task == 'SacLoc': 
        # plot eyetraces per run
        for run, eye_data_run in enumerate(eye_data_all_runs):
            
            dur_run = (eye_data_run[-1,0]-eye_data_run[0,0])
            time_prct = (eye_data_run[:,0]- eye_data_run[0,0])/dur_run

            fig = plotly_layout_template(task,run)

            fig.add_trace(go.Scatter(x = time_prct, y= eye_data_run[:,1],showlegend=False,line=dict(color='black', width=2)), row = 1, col = 1)
            fig.add_trace(go.Scatter(x = time_prct, y= eye_data_run[:,2],showlegend=False,line=dict(color='black', width=2)), row = 2, col = 1)
            fig.add_trace(go.Scatter(x = eye_data_run[:,1], y= eye_data_run[:,2],showlegend=False,line=dict(color='black', width=2)), row = 1, col = 2)

            fig_fn = f"{fig_dir_save}/{subject}_task-{task}_run-0{run+1}_eyetraces.pdf"
            print('Saving {}'.format(fig_fn))
            fig.write_image(fig_fn)

        # Plot per trial with saccades highlighted 
        for run, eye_data_run in enumerate(eye_data_all_runs):
            run_saccade_logic = saccades_output[:,0] == run
            if run >= 9:run_txt = '{}'.format(run+1)
            else:run_txt = '0{}'.format(run+1)

            for sequence in sequences:
                trials = np.arange(0,trials_seq[sequence],1)
                seq_eye_data_logic = np.logical_and(eye_data_run[:,0] >= time_start_seq[sequence,run],\
                                                    eye_data_run[:,0] <= time_end_seq[sequence,run])
                sequence_saccade_logic = saccades_output[:,1] == sequence

                if sequence >= 9:sequence_txt = '{}'.format(sequence+1)
                else:sequence_txt = '0{}'.format(sequence+1)

                # trial start and end
                # define trial start and trial end
                for trial in trials:
                    if trial >= 9:trial_txt = '{}'.format(trial+1)
                    else:trial_txt = '0{}'.format(trial+1)

                    t_trial_start = time_start_trial[trial,sequence,run]
                    t_trial_end = time_end_trial[trial,sequence,run]
                    trial_eye_data_logic = np.logical_and(  eye_data_run[:,0] >= t_trial_start,\
                                                            eye_data_run[:,0] <= t_trial_end)
                    trial_saccade_logic = saccades_output[:,2] == trial
                    
                    eye_data_logic = np.logical_and.reduce(np.array((seq_eye_data_logic,trial_eye_data_logic)))
                    saccade_logic = np.logical_and.reduce(np.array((run_saccade_logic,sequence_saccade_logic,trial_saccade_logic)))

                    
                    time_prct = ((eye_data_run[eye_data_logic][:,0]- t_trial_start)/(t_trial_end - t_trial_start))
                    t, p, x, y = eye_data_run[eye_data_logic,0],time_prct,eye_data_run[eye_data_logic,1],eye_data_run[eye_data_logic,2]
                    #xnb, ynb = eye_data_int_blink[eye_data_logic,1],eye_data_int_blink[eye_data_logic,2]

                    # make basic figure layout
                    fig = plotly_layout_template(task,run)

                    # add data
                    fig.add_trace(go.Scatter(x = time_prct, y= x,showlegend=False,line=dict(color='black', width=2)), row = 1, col = 1)
                    fig.add_trace(go.Scatter(x = time_prct, y= y,showlegend=False,line=dict(color='black', width=2)), row = 2, col = 1)
                    fig.add_trace(go.Scatter(x = x, y= y,showlegend=False,line=dict(color='black', width=2)), row = 1, col = 2)

                    # add text annotation
                    fig.add_annotation(x=8, y=-10,
                    text=f"Seq {sequence + 1}, Trial {trial + 1}",
                    showarrow=False, row = 1, col = 2)

                    # highlight saccades 
        
                    saccades_trial = saccades_output[saccade_logic]
                    if np.isnan(saccades_trial[:,3][0]) == False:
                        for saccade_num in saccades_trial[:,3]:

                            saccade_num_logic = saccades_trial[:,3] == saccade_num
                            saccade_output = saccades_trial[saccade_num_logic,:]

                            sac_t_onset = saccade_output[0,8]
                            sac_t_offset = saccade_output[0,9]
                            sac_p_onset = saccade_output[0,10]
                            sac_p_offset = saccade_output[0,11]
                            blink_saccade = saccade_output[0,25]

                            if blink_saccade == 0:
                                x_sac = x[np.logical_and(t > sac_t_onset,t < sac_t_offset)]
                                y_sac = y[np.logical_and(t > sac_t_onset,t < sac_t_offset)]
                                t_sac = t[np.logical_and(t > sac_t_onset,t < sac_t_offset)]
                                p_sac = p[np.logical_and(t > sac_t_onset,t < sac_t_offset)]


                                fig.add_vrect(x0=sac_p_onset,x1=sac_p_offset, row="all", col=1,
                                fillcolor="darkgrey", opacity=0.45, line_width=0)


                    #fig.show()


                    #fig_fn = f"{file_dir_save}/figures/{task}/per_trial/{subject}_task-{task}_run-0{run+1}_seq-0{sequence+ 1}_trial-0{trial+1}_eyetrace.pdf"
                    fig_fn = f"{fig_dir_save}/{subject}_task-{task}_run-0{run+1}_seq-0{sequence+ 1}_trial-0{trial+1}_eyetrace.pdf"
                    print('Saving {}'.format(fig_fn))

                    fig.write_image(fig_fn)


    elif task == 'PurLoc': 

        # plot eyetraces per run
        for run, eye_data_run in enumerate(eye_data_all_runs):
            
            dur_run = (eye_data_run[-1,0]-eye_data_run[0,0])
            time_prct = (eye_data_run[:,0]- eye_data_run[0,0])/dur_run

            fig = plotly_layout_template(task,run)

            fig.add_trace(go.Scatter(x = time_prct, y= eye_data_run[:,1],showlegend=False,line=dict(color='black', width=2)), row = 1, col = 1)
            fig.add_trace(go.Scatter(x = time_prct, y= eye_data_run[:,2],showlegend=False,line=dict(color='black', width=2)), row = 2, col = 1)
            fig.add_trace(go.Scatter(x = eye_data_run[:,1], y= eye_data_run[:,2],showlegend=False,line=dict(color='black', width=2)), row = 1, col = 2)

            fig_fn = f"{fig_dir_save}/{subject}_task-{task}_run-0{run+1}_eyetraces.pdf"
            print('Saving {}'.format(fig_fn))
            fig.write_image(fig_fn)


        # Plot Saccades per Sequence with contious color scale 
        pursuits_tr = np.arange(0,seq_trs,2)
        saccades_tr = np.arange(1,seq_trs,2)
        trials_seq = saccades_tr


        # Define colors
        # -------------
        cmap = 'hsv'
        col_offset = 0.5


        # Draw figure
        # -----------
        for run, eye_data_run in enumerate(eye_data_all_runs):

            if run > 10: run_txt = '{}'.format(run+1)
            else: run_txt = '0{}'.format(run+1)


            for sequence in eye_mov_seq:

                trials = np.arange(0,trials_seq[sequence],1)
                seq_eye_data_logic = np.logical_and(eye_data_run[:,0] >= time_start_seq[sequence,run],\
                                                    eye_data_run[:,0] <= time_end_seq[sequence,run])

                dur_seq = time_end_seq[sequence,run]-time_start_seq[sequence,run]

                if sequence > 10: sequence_txt = '{}'.format(sequence+1)
                else: sequence_txt = '0{}'.format(sequence+1)
                
                # make basic figure layout
                fig = plotly_layout_template(task,run)
                
                # create legend 
                t = np.linspace(0,2*np.pi,360,endpoint=True)
                r = np.ones_like(t)
                
                fig.add_trace(go.Barpolar(
                    r= r, 
                    theta=np.degrees(t) + 180,
                    marker_line_width=0.001,
                    width=np.pi / np.pi,
                    marker=dict(
                        color=np.degrees(t),  # Color based on theta grid
                        colorscale=cmap,  
                        showscale=False  
                    

                    ), 
                ),row=1, col=3)

                # Add lines dividing the circle into 8 bins
                theta_bins = np.linspace(0, 2*np.pi, 9)  
                for theta in theta_bins:
                    fig.add_trace(go.Scatterpolar(
                        r=[0, 1],
                        theta=[np.degrees(theta), np.degrees(theta)],
                        mode='lines',
                        line=dict(color='rgb(204,204,204)', width=1)
                    ),row=1, col=3)

                # Add two smaller circles
                r_small = [0.25, 0.5, 0.75]
                for r_val in r_small:
                    fig.add_trace(go.Scatterpolar(
                        r=np.full_like(t, r_val),
                        theta=np.degrees(t),
                        mode='lines',
                        line=dict(color='rgb(204,204,204)', width=1)
                    ),row=1, col=3)
                            
                
                fig.update_polars(angularaxis=dict(visible=False), 
                                radialaxis=dict(visible=False))
                
                fig.update_layout(polar=dict(domain=dict(y=[0.85, 1], x= [0.84,0.9])))

                # Add text annotation for current sequence 
                text_annotation = f"Seq {sequence_txt}"
                fig.add_annotation(x=10, y=-10,
                    text=text_annotation,
                    showarrow=False, row=1, col=2)

                for trial_num, trial_plot in enumerate(trials_seq):

                    trial_eye_data_logic = np.logical_and(eye_data_run[:,0] >= time_start_trial[trial_plot,sequence,run],
                                                            eye_data_run[:,0] <= time_end_trial[trial_plot,sequence,run])

                    data_logic = np.logical_and.reduce(np.array((seq_eye_data_logic,trial_eye_data_logic)))

                    time_prct = ((eye_data_run[data_logic][:,0]- time_start_seq[sequence,run])/dur_seq)
                    
                    plot_color = colors.hsv_to_rgb([(trial_num / len(trials_seq) + col_offset) % 1.0, 1.0, 1.0])
                    plot_color_hex = colors.to_hex(plot_color)
            
                    # add data
                    fig.add_trace(go.Scatter(x = time_prct, y= eye_data_run[data_logic,1],showlegend=False,line=dict(color=plot_color_hex, width=2)), row = 1, col = 1)
                    fig.add_trace(go.Scatter(x = time_prct, y= eye_data_run[data_logic,2],showlegend=False,line=dict(color=plot_color_hex, width=2)), row = 2, col = 1)
                    fig.add_trace(go.Scatter(x = eye_data_run[data_logic,1], y= eye_data_run[data_logic,2],showlegend=False,line=dict(color=plot_color_hex, width=2)), row = 1, col = 2)


                fig_fn = f"{fig_dir_save}/{subject}_task-{task}_run-0{run+1}_seq-0{sequence+ 1}_eyetrace.png"
                print('Saving {}'.format(fig_fn))
                fig.write_image(fig_fn)

                
