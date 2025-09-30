"""
-----------------------------------------------------------------------------------------
vdm_modifier.py
-----------------------------------------------------------------------------------------
Goal of the script :
Create a slack of images generating the visual design of a pRF experiment and shift it 
using a set of eyetracking data position
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: server group (e.g. 327)
sys.argv[4]: eye movement dataset file root (e.g. 'prf_em_ctrl')
-----------------------------------------------------------------------------------------
Output(s) :
Video of the vdm
np array of the vdm
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/amblyo_prf/analysis_code/postproc/prf/fit
2. run python command
python vdm_builder.py  [main directory] [project name] [group] [eye movement file]
-----------------------------------------------------------------------------------------
Written by Martin Szinte (martin.szinte@gmail.com)
Edited by Adrien Chopin (adrien.chopin@gmail.com)
-----------------------------------------------------------------------------------------
"""

#General imports
import cv2 #import json
import math
import numpy as np
import os
import sys
import csv
from scipy.ndimage import shift
sys.path.append("{}/../../../utils".format(os.getcwd()))
from conversion_utils import conversion

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
group = sys.argv[3]
em_ctrl_fn = sys.argv[4]

# Settings
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../../"))
settings_path = os.path.join(base_dir, project_dir, "settings.json")

with open(settings_path) as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)

# Parameters
screen_converter = conversion(screen_size_pix=analysis_info['screen_size_pix'], 
                              screen_size_cm=analysis_info['screen_size_cm'],
                              screen_distance_cm=analysis_info['screen_distance_cm'])

TR = analysis_info['TR']
n = analysis_info['screen_size_pix'][1]
vdm_size_pix = analysis_info['vdm_size_pix']
TRs = analysis_info['TRs']
apperture_rad_dva = analysis_info['apperture_rad_dva']
apperture_rad_pix = round(np.array([screen_converter.dva2pix(apperture_rad_dva)[0],screen_converter.dva2pix(apperture_rad_dva)[1]]).mean()) 
bar_length_pix = apperture_rad_pix
bar_width_dva = analysis_info['bar_width_dva']
bar_width_pix = round(np.array(screen_converter.dva2pix(bar_width_dva)[0],screen_converter.dva2pix(bar_width_dva)[1]).mean())

pass_duration = 13
blank_duration = 6.5
delays = [5, 6.5]
bar_list = np.array([1, 2, 0, 3, 4, 0, 5, 6, 0, 7, 8])
list_angles = np.array([np.nan, 45, 0, 90, 135, 225, 180, 270, 315])

# Eye movement file
em_filepath = '{}/{}/derivatives/vdm/{}.csv'.format(main_dir, project_dir, em_ctrl_fn)
with open(em_filepath, 'r') as file: 
    csv_reader = csv.reader(file)
    data = list(csv_reader)
    eye_movements_data = np.array(data)
    eye_movements_data = eye_movements_data[1:, :]
    eye_movements = np.array(eye_movements_data, dtype=float).astype(int)

# Define directories
fileName_em_ctrl = 'vdm_{}_{}_{}'.format(vdm_size_pix[0], vdm_size_pix[1], em_ctrl_fn)
filepath_em_ctrl = os.path.join(rootpath,fileName_em_ctrl + '.npy')
videopath_shift = os.path.join(rootpath,fileName_em_ctrl + '.mp4')

# Create a meshgrid of image coordinates x, y, a list of angles by TR
x, y = np.meshgrid(range(0,n), range(0,n))
angle_list = list_angles[bar_list];
angle_halfTR = np.empty((1,2*TRs)); angle_halfTR.fill(np.nan)
head = 0; newhead = 2*delays[0];
angle_halfTR[0,head:newhead]=np.nan
for i in angle_list:
    head = newhead
    if np.isnan(i):
        newhead=int(head+2*blank_duration)
    else:
        newhead=int(head+2*pass_duration) 
    angle_halfTR[0,head:newhead]=i
angle_halfTR[0,newhead:]=np.nan

# Define the draw_frame function
def draw_frame(x,y,position,n,bar_length_pix,bar_width_pix,angle,apperture_rad_pix):
    frame = np.zeros((n,n))
    center_x = round(n/2)
    center_y = round(n/2)
    if ((position[0]-center_x)!=0) & ((position[1]-center_y)!=0):
        position_to_center_line_slope = (position[1]-center_y)/(position[0]-center_x)
        a = -1/position_to_center_line_slope
        b_low = (position[1]-np.sin(np.radians(angle))*bar_width_pix/2)-a*(position[0]-np.cos(np.radians(angle))*bar_width_pix/2)
        b_up = (position[1]+np.sin(np.radians(angle))*bar_width_pix/2)-a*(position[0]+np.cos(np.radians(angle))*bar_width_pix/2)
        frame[(y>(a*x+min(b_low,b_up)))&(y<(a*x+max(b_low,b_up)))]=1
    else:
        if (position[0]-center_x)==0:
            frame[(y>(position[1]-bar_width_pix/2))&(y<(position[1]+bar_width_pix/2))]=1
        elif (position[1]-center_y)==0:
            frame[(x>(position[0]-bar_width_pix/2))&(x<(position[0]+bar_width_pix/2))]=1
        else:
            print('oops!')
    # apply aperture
    frame[((x-center_x)**2+(y-center_y)**2)>apperture_rad_pix**2]=0
    return frame

# Create the frame through the list of angles
current_angle = np.nan
center_x = round(n/2)
center_y = round(n/2)
list_im=np.array([])
frames = np.zeros((n,n,2*TRs));

for i in range(0,np.size(angle_halfTR)):
    angle = angle_halfTR[0,i]
    if ~np.isnan(angle):
        if angle!=current_angle:
            current_angle = angle
            start_position = np.array([center_x, center_y])+(apperture_rad_pix+bar_width_pix/2)*np.array([math.cos(math.radians(current_angle+180)),math.sin(math.radians(current_angle+180))]) 
            end_position = np.array([center_x, center_y])+(apperture_rad_pix+bar_width_pix/2)*np.array([math.cos(math.radians(current_angle)),math.sin(math.radians(current_angle))])
            distance = end_position - start_position
            step = distance/(2*pass_duration-1)
            position = start_position
            list_im=np.append(list_im,i)
            print('Done direction angle '+str(angle))
        else:
            position = position + step
        position_rnd = position.round() 
        frames[:,:,i]=draw_frame(x,y,position_rnd,n,bar_length_pix,bar_width_pix,angle,apperture_rad_pix)

frames = frames[:,:,0::2] # only save the full TR, not the half-TR

# shit the data in frames by the amounts documented in eye_movements (in a new frame)
def shift_frames_by_eye_movements(frames, eye_movements):
    shifted_frames = np.zeros_like(frames)
    for i in range(frames.shape[2]):  # Iterate through each frame
        x_shift = -eye_movements[i, 0]  # sign does not matter because we are only doing a simulation - if matters, should be at least negative (opposite to eye movements)
        y_shift = -eye_movements[i, 1]  # sign does not matter because we are only doing a simulation - if matters, should be at least negative (opposite to eye movements)
        
        # Shift the frame
        
        shifted_frame = shift(frames[:, :, i], (y_shift, x_shift), mode='constant', cval=0)
        
        # Add the shifted frame to our new array
        shifted_frames[:, :, i] = shifted_frame
    
    return shifted_frames

# apply the shift
frames_shift = shift_frames_by_eye_movements(frames, eye_movements)

# Downsampling : resize the frames and inverse y axis
frames_reshape = np.zeros((vdm_size_pix[0],vdm_size_pix[1],TRs))
for k in range(frames_reshape.shape[-1]):
    frames_reshape[:,:,k] = cv2.resize(frames[:,:,k], dsize=(vdm_size_pix[0], vdm_size_pix[1]), interpolation=cv2.INTER_NEAREST)
frames = frames_reshape
frames_reshape_shift = np.zeros((vdm_size_pix[0],vdm_size_pix[1],TRs))
for k in range(frames_reshape_shift.shape[-1]):
    frames_reshape_shift[:,:,k] = cv2.resize(frames_shift[:,:,k], dsize=(vdm_size_pix[0], vdm_size_pix[1]), interpolation=cv2.INTER_NEAREST)
frames_shift = frames_reshape_shift

frames_rotate = np.zeros((vdm_size_pix[0],vdm_size_pix[1],TRs))
for num, frame in enumerate(np.split(frames, TRs, axis=2)):
    frames_rotate[:,:,num] = frame[-1::-1,:,0]
frames = frames_rotate
frames_rotate_shift = np.zeros((vdm_size_pix[0],vdm_size_pix[1],TRs))
for num, frame in enumerate(np.split(frames_shift, TRs, axis=2)):
    frames_rotate_shift[:,:,num] = frame[-1::-1,:,0]
frames_shift = frames_rotate_shift

# Export videos
fourcc_shift = cv2.VideoWriter_fourcc(*'mp4v')
out_shift = cv2.VideoWriter(videopath_shift, fourcc_shift, 1/TR, (vdm_size_pix[0], vdm_size_pix[1]), False)
[out_shift.write(np.uint8(frame*255)) for frame in np.split(frames_shift, TRs, axis=2)]
out_shift.release()
print('Video conversion done (eye-movement-shifted version), saved to:'+videopath_shift)

# Save numpy array
np.save(filepath_shift, frames_shift)
print('Data (eye-movement-shifted version) saved to :'+filepath_shift)

# Define permission cmd
print('Changing files permissions in {}/{}'.format(main_dir, project_dir))
os.system("chmod -Rf 771 {}/{}".format(main_dir, project_dir))
os.system("chgrp -Rf {} {}/{}".format(group, main_dir, project_dir))