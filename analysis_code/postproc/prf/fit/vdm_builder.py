"""
-----------------------------------------------------------------------------------------
vdm_builder.py
-----------------------------------------------------------------------------------------
Goal of the script:
Create visual design matrix from stimulus video
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: task 
-----------------------------------------------------------------------------------------
Output(s):
npy of vdm 
mp4 video of vdm to visually inspect 
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/pRF_analysis/analysis_code/postproc/prf/fit
2. run python command
python vdm_builder.py [main directory] [project name] [task] 
>> python vdm_builder.py /scratch/mszinte/data/ centbids pRF
-----------------------------------------------------------------------------------------
Written by Martin Szinte (martin.szinte@gmail.com)
Edited by Uriel Lascombes (uriel.lascombes@laposte.net) & Sina Kling (sina.kling@outlook.de)
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
import cv2
import json
import numpy as np
import sys

# Define directories 
main_dir = sys.argv[1]
project_name = sys.argv[2]
project_dir = os.path.join(main_dir,  project_name)
print(project_dir)

# Analysis parameters from settings
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../../"))
settings_path = os.path.join(base_dir, project_name, "settings.json")

with open(settings_path) as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)


task = sys.argv[3]
stim_video_fn = f'{project_dir}/derivatives/vdm/{task}_vid.mp4'

# open video
cap = cv2.VideoCapture(stim_video_fn)
if not cap.isOpened():
    raise FileNotFoundError(f"Could not open video file {stim_video_fn}")

# get video settings
vid_width, vid_height = cap.get(3), cap.get(4)
vid_fps = cap.get(cv2.CAP_PROP_FPS)
vid_frames = cap.get(cv2.CAP_PROP_FRAME_COUNT)

print(f"Video properties= Width: {vid_width}, Height: {vid_height}, FPS: {vid_fps}, Frames: {vid_frames}")

# set VDM settings
TR = analysis_info['TR']

vdm_width = analysis_info['vdm_size_pix'][0] 
vdm_height = analysis_info['vdm_size_pix'][1]

vdm_video_fn = f'{project_dir}/derivatives/vdm/vdm_{task}_{vdm_width}_{vdm_height}.mp4'
vdm_numpy_fn = f'{project_dir}/derivatives/vdm/vdm_{task}_{vdm_width}_{vdm_height}.npy'

vdm_frames_sampled = vid_fps * TR
vdm_frames = vid_frames / vdm_frames_sampled
vdm_mat = np.zeros((vdm_width, vdm_height, int(vdm_frames)))

# define top and bottom part as 16/9 screen experiment and VDM should be a square: example visualization for 7T CRMBM setup

#=============== Screen ===============================
#               1920 px  (~73.3 cm)                                                      ============ Converting rectangular screen to square ============
# +--------------------------------------------------+												
# |                                                  |
# |                                                  |													1920 px  (~73.3 cm)
# |                                                  |									  +--------------------------------------------------+
# |                                                  |	1080 px  (~43.5 cm)  ->			  |                  MARGIN                          |
# +--------------------------------------------------+					        		  |              (420 px each)                       |
#																						  +--------------------------------------------------+
#																						  |                                                  |
#																						  |                SQUARE AREA                       |
#																						  |               1920 x 1920 px                     |
#																						  |                                                  |
#																						  +--------------------------------------------------+
#																						  |                  MARGIN                          |
#																						  |              (420 px each)                       |
#																						  +--------------------------------------------------+
# Added margins = (1920 - 1080) / 2 = 420 px (top & bottom).  Final dimensions, to be taken into account when giving dimensions to prfpy!
# Square size   = 1920 px x 1920 px

height_to_add = (vid_width - vid_height) / 2 
add_mat =  np.zeros((int(height_to_add),int(vid_width)))

# create VDM matrix
vid_idx, vdm_idx = 0, 0
print("\nStarting frame processing...")
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
    
        # add top and bottom blank    
        binary_reshape_mat = np.concatenate((add_mat,binary_mat,add_mat), axis=0)
        
        # resize to create VDM
        binary_reshape_resize_mat = cv2.resize(binary_reshape_mat, dsize=(vdm_width, vdm_height),
                                               interpolation=cv2.INTER_NEAREST)
        
        # fill VDM matrix
        vdm_mat[...,vdm_idx] = binary_reshape_resize_mat
        vdm_idx += 1


    vid_idx += 1
    
cap.release()
print("Finished processing all frames")


import numpy as np

# find coordinates of white pixels
white_coords = np.argwhere(vdm_mat == 1)
white_heights = white_coords[:, 0]
white_widths = white_coords[:, 1]
h = max(white_heights) - min(white_heights) + 1
w = max(white_widths) - min(white_widths) + 1
print(f"Vdm Stimulus grid in pixels: {h} x {w}")

# ===== Stimulus grid ==========
#
#   0        18       32       50   (px)
#   |--------|--------|--------|
#   |        |        |        |
#   |        +--------+        |
#   |        | 15x15  |        |
#   |        |  box   |        |
#   |        +--------+        |
#   |                          |
#   +--------------------------+
#  0                          50
#        (50 x 50 px matrix)

# save VDM as numpy matrix
np.save(vdm_numpy_fn, vdm_mat)

# save VDM as video
print(f"Saving VDM video: {vdm_video_fn}")
fourcc = cv2.VideoWriter_fourcc(*'mp4v')
out = cv2.VideoWriter(vdm_video_fn, fourcc, 1/TR, (vdm_width, vdm_height), False)
[out.write(np.uint8(frame*255)) for frame in np.split(vdm_mat, vdm_frames, axis=2)]
out.release()