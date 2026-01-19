"""
-----------------------------------------------------------------------------------------
extract_saccades.py
-----------------------------------------------------------------------------------------
Goal of the script:
Extract saccade metrics
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject name
sys.argv[4]: group of shared data (e.g. 327)
-----------------------------------------------------------------------------------------
Output(s):
h5 files with vals_all
vals_all[:,00]:	run number
vals_all[:,01]:	sequence number
vals_all[:,02]:	trial number
vals_all[:,03]:	saccade number detected (nan if no),
vals_all[:,04]:	saccade onset x coordinate (dva from screen center)
vals_all[:,05]:	saccade offset x coordinate (dva from screen center))
vals_all[:,06]:	saccade onset y coordinate (dva from screen center)
vals_all[:,07]:	saccade offset y coordinate (dva from screen center))
vals_all[:,08]:	saccade onset time trigger
vals_all[:,09]:	saccade offset time trigger
vals_all[:,10]:	saccade onset time relative to trial onset (proportion of the trial)
vals_all[:,11]:	saccade offset time relative to trial onset (proportion of the trial)
vals_all[:,12]:	saccade duration (ms)
vals_all[:,13]:	saccade velocity peak (dva/sec)
vals_all[:,14]:	saccade distance (dva)
vals_all[:,15]:	saccade amplitude (dva)
vals_all[:,16]:	saccade distance angle (degrees),
vals_all[:,17]:	saccade amplitude angle (degrees),
vals_all[:,18]:	saccade trial with correct fixation (start within boundary)
vals_all[:,19]:	saccade trial with correct saccade (end within boundary)
vals_all[:,20]:	saccade_task trial
vals_all[:,21]:	trial with missing time stamps
vals_all[:,22]:	saccade task with accurate saccade outwards
vals_all[:,23]:	saccade task with accurate saccade inswards
vals_all[:,24]:	trial with no_saccade detected,
vals_all[:,25]:	microsaccade detected (<1 dva)
-----------------------------------------------------------------------------------------
To run:
cd ~/projects/pRF_analysis/RetinoMaps/eyetracking
python extract_saccades.py /scratch/mszinte/data RetinoMaps sub-01 327
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

# General imports
# ---------------
import os
import sys
import json
import h5py
import ipdb
import numpy as np
import pandas as pd


# Personal imports
sys.path.append("{}/../../analysis_code/utils".format(os.getcwd()))
from sac_utils import vecvel, microsacc_merge, saccpar, isincircle

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
group = sys.argv[4]

# Load general settings
with open('../settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)
tasks = analysis_info['task_intertask'][0]
prf_task_name = analysis_info['prf_task_name']

if prf_task_name in tasks:
    tasks.remove(prf_task_name)


for task in tasks :
    print('Processing {} ...'.format(task))
    base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../"))
    settings_path = os.path.join(base_dir, project_dir, f'{task}_settings.json')
    with open(settings_path) as f:
    	settings = json.load(f)
    if subject == 'sub-01':
    	if task == 'pRF': ses = 'ses-02'
    	else: ses = 'ses-01'
    else: ses = settings['session']
    
    # Load main experiment settings 
    runs = np.arange(0, settings['num_run'],1)
    sequences = np.arange(0, settings['num_seq'],1)
    trials_seq = settings['trials_seq']
    rads = settings['rads']
    polar_ang = np.deg2rad(np.arange(0,360,settings['ang_steps']))
    pursuits_tr = np.arange(0,settings['seq_trs'],2)
    saccades_tr = np.arange(1,settings['seq_trs'],2)
    seq_type = settings['seq_type']
    
    #-------------------- Load data --------------------------------------
    file_dir_save = '{}/{}/derivatives/pp_data/{}/eyetracking/'.format(
        main_dir, project_dir, subject)
    # os.makedirs(file_dir_save, exist_ok=True)
    h5_filename = '{file_dir}/stats/{sub}_task-{task}_eyedata_sac_stats.h5'.format(
        file_dir = file_dir_save, sub = subject, task = task)
    h5_file = h5py.File(h5_filename,'r') 
    time_start_seq = np.array(h5_file['time_start_seq'])
    time_end_seq = np.array(h5_file['time_end_seq'])
    time_start_trial = np.array(h5_file['time_start_trial'])
    time_end_trial = np.array(h5_file['time_end_trial'])
    amp_sequence = np.array(h5_file['amp_sequence'])
    eye_data_run_01_nan_blink_interpol = pd.read_csv(
        f"{file_dir_save}/timeseries/{subject}_task-{task}_run_01_eyedata.tsv.gz", compression='gzip', delimiter='\t')
    eye_data_run_01_nan_blink_interpol = eye_data_run_01_nan_blink_interpol[['timestamp', 
                                                                             'x', 
                                                                             'y', 
                                                                             'pupil_size']].to_numpy()
    eye_data_run_02_nan_blink_interpol = pd.read_csv(
        f"{file_dir_save}/timeseries/{subject}_task-{task}_run_02_eyedata.tsv.gz", compression='gzip', delimiter='\t')
    eye_data_run_02_nan_blink_interpol = eye_data_run_02_nan_blink_interpol[['timestamp', 'x', 'y', 'pupil_size']].to_numpy()
    eye_data_all_runs = [eye_data_run_01_nan_blink_interpol,eye_data_run_02_nan_blink_interpol]
    
    
    # Get saccade model
    sampling_rate = settings['sampling_rate']
    velocity_th = settings['velocity_th']
    min_dur = settings['min_dur']
    merge_interval = settings['merge_interval']
    tolerance_ratio = settings['tolerance_ratio']
    
    #----------------------- Main loop ------------------------------------
    mat = 0
    for run in runs:
    	print(f'--extracting saccades from run: {run+1}--')
    	eye_data_run = pd.read_csv(
            f"{file_dir_save}/timeseries/{subject}_task-{task}_run_0{run+1}_eyedata.tsv.gz", compression='gzip', delimiter='\t')
    	eye_data_run = eye_data_run[['timestamp', 'x', 'y', 'pupil_size']].to_numpy()
    
    		
    	for sequence in sequences: 
    		print('sequence: {}'.format(sequence))
    		trials = np.arange(0,trials_seq[sequence],1)
    		seq_data_logic = np.logical_and(eye_data_run[:,0] >= time_start_seq[sequence,run],eye_data_run[:,0] <= time_end_seq[sequence,run])	
    		trial_with_sac = 0
    		for trial in trials:
    			#print('trial: {}'.format(trial))
    			trial_data_logic = np.logical_and(eye_data_run[:,0] >= time_start_trial[trial,sequence,run],\
    											  eye_data_run[:,0] <= time_end_trial[trial,sequence,run])
    
    			data_logic = np.logical_and.reduce(np.array((seq_data_logic,trial_data_logic)))
    			
    			# fixation target position
    			if (amp_sequence[sequence] == 5) :
    				amp_sac = 0
    				fix_pos_x, fix_pos_y = 0,0
    				sac_pos_x,sac_pos_y = 0,0
    				
    			else:
    				amp_sac = rads[int(amp_sequence[sequence])]
    				fix_pos_x,fix_pos_y = 0,0
    				sac_pos_x, sac_pos_y = np.round(np.cos(polar_ang[trial_with_sac])*amp_sac,decimals=3),\
    								   	   np.round(np.sin(polar_ang[trial_with_sac])*amp_sac,decimals=3)
    			
    			# trial start and end
    			# define trial start and trial end
    			t_trial_start = time_start_trial[trial,sequence,run]
    			t_trial_end = time_end_trial[trial,sequence,run]
    			dur_trial = t_trial_end - t_trial_start
    			time_prct = ((eye_data_run[data_logic][:,0]- t_trial_start)/dur_trial)
    			
    
    			# Indicators
    			saccade_task = 0	  #0 saccade task
    			miss_time = 0         #1 missing data
    			sac_out_accuracy = 0  #2 saccade out accuracy 
    			sac_in_accuracy = 0   #3 saccade out accuracy 
    			microsaccade = 0      #4 microsaccade
    			no_saccade = 0        #5 no saccade detected
    			blink_saccade = 0	  #6 blink saccade
    
    			# Saccade analysis parameters
    			num_res = 26
    			sac_fix_rad = tolerance_ratio*amp_sac
    
    			#0 Saccade task
    			if seq_type[sequence] == 1:
    				if np.sum(saccades_tr==trial):
    					saccade_task = 1
    					trial_with_sac += 1
    
    			#1 Missing data point
       			
    			if np.sum(np.diff(eye_data_run[trial_data_logic,0])>1000/sampling_rate) > 0:
    				miss_time = 1
    
    				
    			#2 saccade detection
    			if not miss_time:
    				t, p, x, y = eye_data_run[trial_data_logic,0],time_prct,eye_data_run[trial_data_logic,1],eye_data_run[trial_data_logic,2]
    				vx, vy = vecvel(x,y,sampling_rate)
    				sac = microsacc_merge(x,y,vx,vy,velocity_th,min_dur,merge_interval)
    				ms = saccpar(sac)
    
    				if np.isnan(ms[0,0]):
    					#4 no saccade
    					no_saccade = 1
    					s1 = 0
    
    					fix_cor = 0
    					sac_cor = 0
    
    					if mat == 0:
    						vals_all = np.array([	run,			sequence,		trial,			    np.nan,			  np.nan,
    												np.nan,			np.nan,			np.nan,			    np.nan,			  np.nan,
    												np.nan,			np.nan,			np.nan,			    np.nan,			  np.nan,
    												np.nan,			np.nan,			np.nan,			    fix_cor,		  sac_cor,
    												saccade_task, 	miss_time,		sac_out_accuracy,	sac_out_accuracy, no_saccade,		
    												microsaccade,   blink_saccade])
    
    						mat = 1
    					else:
    						vals_all = np.vstack((vals_all,np.array([	run,			sequence,		trial,			    np.nan,			 np.nan,
    																	np.nan,			np.nan,			np.nan,			    np.nan,			 np.nan,
    																	np.nan,			np.nan,			np.nan,			    np.nan,			 np.nan,
    																	np.nan,			np.nan,			np.nan,			    fix_cor,		 sac_cor,
    																	saccade_task, 	miss_time,		sac_out_accuracy,	sac_in_accuracy, no_saccade,		
    																	microsaccade, blink_saccade])))
    
    				else:
    					n_sac_tot = ms.shape[0]
    					s1 = 0
    					fix_cor = 1
    					sac_cor = 1
    					while s1 < n_sac_tot:
    						sac_x_onset,	sac_x_offset	=	x[int(ms[s1,0])], 	x[int(ms[s1,1])]
    						sac_y_onset,	sac_y_offset	=	y[int(ms[s1,0])], 	y[int(ms[s1,1])]
    						sac_t_onset,	sac_t_offset	=	t[int(ms[s1,0])], 	t[int(ms[s1,1])]
    						sac_p_onset,	sac_p_offset	=	p[int(ms[s1,0])], 	p[int(ms[s1,1])]
    						sac_dur,		sac_vpeak		=	ms[s1,2],			ms[s1,3]
    						sac_dist,		sac_amp			= 	ms[s1,4],			ms[s1,6]
    						sac_dist_ang,	sac_amp_ang		= 	ms[s1,5],			ms[s1,7]
    
    						fix_out = isincircle(sac_x_onset,sac_y_onset,fix_pos_x,fix_pos_y,sac_fix_rad)
    						sac_out = isincircle(sac_x_offset,sac_y_offset,sac_pos_x,sac_pos_y,sac_fix_rad)
    						# only if saccade onset is within fixation radius and saccade offset is within saccade radius (saccade outwards)
    
    						fix_in = isincircle(sac_x_offset,sac_y_offset,fix_pos_x,fix_pos_y,sac_fix_rad)
    						sac_in = isincircle(sac_x_onset,sac_y_onset,sac_pos_x,sac_pos_y,sac_fix_rad)
    						# only if saccade offset is within fixation radius and saccade onset is within saccade radius (saccade inwards)
    
    
    						if np.logical_and(fix_out,sac_out):sac_out_accuracy = 1
    						else:sac_out_accuracy = 0
    
    						if np.logical_and(fix_in,sac_in):sac_in_accuracy = 1
    						else:sac_in_accuracy = 0
    
    						#3 microsaccade
    						if sac_amp <= 1.0:microsaccade = 1
    
    						# extract metrics
    						if mat == 0:
    							vals_all = np.array([	run,			sequence,		trial,			    s1,				 sac_x_onset,
    													sac_x_offset,	sac_y_onset,	sac_y_offset,	    sac_t_onset,	 sac_t_offset,
    													sac_p_onset,	sac_p_offset,	sac_dur,		    sac_vpeak,		 sac_dist,
    													sac_amp,		sac_dist_ang,	sac_amp_ang,	    fix_cor,		 sac_cor,
    													saccade_task, 	miss_time,		sac_out_accuracy,	sac_in_accuracy, no_saccade,		
    													microsaccade,   blink_saccade])
    							mat = 1
    						else:
    						    vals_all = np.vstack((vals_all,np.array([	run,			sequence,		trial,			    s1,				 sac_x_onset,
    																		sac_x_offset,	sac_y_onset,	sac_y_offset,	    sac_t_onset,	 sac_t_offset,
    																		sac_p_onset,	sac_p_offset,	sac_dur,		    sac_vpeak,		 sac_dist,
    																		sac_amp,		sac_dist_ang,	sac_amp_ang,	    fix_cor,		 sac_cor,
    																		saccade_task, 	miss_time,		sac_out_accuracy,	sac_in_accuracy, no_saccade,		
    																		microsaccade,   blink_saccade])))
    							
    						s1 += 1
    			else:
    				if mat == 0:
    					vals_all = np.array([	run,			sequence,		trial,			    np.nan,			 np.nan,
    											np.nan,			np.nan,			np.nan,			    np.nan,			 np.nan,
    											np.nan,			np.nan,			np.nan,			    np.nan,			 np.nan,
    											np.nan,			np.nan,			np.nan,			    fix_cor,		 sac_cor,
    											saccade_task, 	miss_time,		sac_out_accuracy,	sac_in_accuracy, no_saccade,		
    											microsaccade,   blink_saccade])
    					mat = 1
    				else:
    					vals_all = np.vstack((vals_all,np.array([	run,			sequence,		trial,			    np.nan,			 np.nan,
    																np.nan,			np.nan,			np.nan,			    np.nan,			 np.nan,
    																np.nan,			np.nan,			np.nan,			    np.nan,			 np.nan,
    																np.nan,			np.nan,			np.nan,			    fix_cor,		 sac_cor,
    																saccade_task, 	miss_time,		sac_out_accuracy,	sac_in_accuracy, no_saccade,		
    																microsaccade,   blink_saccade])))
    
    
    
    # Save all
    # --------
    h5_file = "{file_dir}/stats/{sub}_task-{task}_eyedata_sac_stats.h5".format(file_dir = file_dir_save, sub = subject, task = task)
    h5file = h5py.File(h5_file, "a")
    h5file.create_dataset('saccades_output',data = vals_all,dtype ='float32')
    h5file.close()
    
# # Define permission cmd
# print('Changing files permissions in {}/{}'.format(main_dir, project_dir))
# os.system("chmod -Rf 771 {}/{}".format(main_dir, project_dir))
# os.system("chgrp -Rf {} {}/{}".format(group, main_dir, project_dir))

    