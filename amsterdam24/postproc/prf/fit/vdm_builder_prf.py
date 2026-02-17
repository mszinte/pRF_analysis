"""
-----------------------------------------------------------------------------------------
vdm_builder_prf.py
-----------------------------------------------------------------------------------------
Goal of the script:
This code takes a video of a pRF visual stimulus with dimensions width × height 
(e.g., 1920×1080 pixels) as input. It removes n pixels (e.g., 420 pixels) 
from both the left and right sides to produce a square height × height 
(e.g., 1080×1080 pixels) VDM for pRF fitting.
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: group (e.g. 327)
-----------------------------------------------------------------------------------------
Output(s):
npy of vdm 
mp4 video of vdm to visually inspect 
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/pRF_analysis/analysis_code/postproc/prf/fit
2. run python command
>> python vdm_builder_prf.py [main directory] [project name] [group] 
-----------------------------------------------------------------------------------------
Exemple:
cd ~/projects/pRF_analysis/amsterdam24/postproc/prf/fit
python vdm_builder_prf.py /scratch/mszinte/data RetinoMaps 327
-----------------------------------------------------------------------------------------
Written by Martin Szinte (martin.szinte@gmail.com)
edited by Uriel Lascombes (uriel.lascombes@laposte.net)
edited by Sina Kling (sina.kling@outlook.de)
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
import cv2
import yaml
import numpy as np

# Personal imports
sys.path.append("{}/../../../../analysis_code/utils".format(os.getcwd()))
from settings_utils import load_settings

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
group = sys.argv[3]

# Load settings
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../../"))
settings_path = os.path.join(base_dir, project_dir, "settings.yml")
prf_settings_path = os.path.join(base_dir, project_dir, "prf-analysis.yml")
settings = load_settings([settings_path, prf_settings_path])
analysis_info = settings[0]

prf_task_names = analysis_info['prf_task_names']
TR = analysis_info['TR']
vdm_width = analysis_info['vdm_size_pix'][0]
vdm_height = analysis_info['vdm_size_pix'][1]

deb()

for prf_task_name in prf_task_names:
    
    # Video dir
    stim_video_fn = '{}/{}/derivatives/vdm/task-{}_vid.mp4'.format(
        main_dir, project_dir, prf_task_name)
    
    # open video
    cap = cv2.VideoCapture(stim_video_fn)
    
    # get video settings
    vid_width, vid_height = cap.get(3), cap.get(4)
    vid_fps = cap.get(cv2.CAP_PROP_FPS)
    vid_frames = cap.get(cv2.CAP_PROP_FRAME_COUNT)
    
    # set VDM settings
    vdm_frames_sampled = vid_fps * TR
    vdm_frames = vid_frames / vdm_frames_sampled
    vdm_mat = np.zeros((vdm_width, vdm_height, int(vdm_frames)))
    
    # Defind the portion of width to remove
    width_to_remove = (vid_width - vid_height) / 2 
    
    # create VDM matrix
    vid_idx, vdm_idx = 0, 0
    while(cap.isOpened()):
        
        # read video frames
        ret, vid_frame = cap.read()
        if not ret: break
    
        # process sampled frames
        if np.mod(vid_idx, int(vdm_frames_sampled)) == 0:
    
            # convert frame to grayscale
            gray_mat = cv2.cvtColor(vid_frame, cv2.COLOR_BGR2GRAY)
    
            # convert to binary
            binary_mat = (gray_mat > 5).astype(np.uint8)
        
            # Remove pixels in left and right to reach height x height
            binary_reshape_mat = binary_mat[:, int(width_to_remove):-int(width_to_remove)]
            
            # resize to create VDM
            binary_reshape_resize_mat = cv2.resize(binary_reshape_mat, dsize=(vdm_width, vdm_height),
                                                   interpolation=cv2.INTER_NEAREST)
            
            # fill VDM matrix
            vdm_mat[...,vdm_idx] = binary_reshape_resize_mat
            vdm_idx += 1
            
        vid_idx += 1
    
    cap.release()
    
    # save VDM as numpy matrix
    vdm_numpy_fn = '{}/{}/derivatives/vdm/task-{}_vdm.npy'.format(
        main_dir, project_dir, prf_task_name, vdm_width, vdm_height)
    print('Saving {}'.format(vdm_numpy_fn))
    np.save(vdm_numpy_fn, vdm_mat)
    
    # save VDM as video
    vdm_video_fn = '{}/{}/derivatives/vdm/task-{}_vdm.mp4'.format(
        main_dir, project_dir, prf_task_name, vdm_width, vdm_height)
    print('Saving {}'.format(vdm_video_fn))
    fourcc = cv2.VideoWriter_fourcc(*'mp4v')
    out = cv2.VideoWriter(vdm_video_fn, fourcc, 1/TR, (vdm_width, vdm_height), False)
    [out.write(np.uint8(frame*255)) for frame in np.split(vdm_mat, vdm_frames, axis=2)]
    out.release()
    
# Define permission cmd
print('Changing files permissions in {}/{}'.format(main_dir, project_dir))
os.system("chmod -Rf 771 {}/{}".format(main_dir, project_dir))
os.system("chgrp -Rf {} {}/{}".format(group, main_dir, project_dir))