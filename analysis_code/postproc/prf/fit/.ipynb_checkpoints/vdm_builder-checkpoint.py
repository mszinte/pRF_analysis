"""
-----------------------------------------------------------------------------------------
vdm_builder.py
-----------------------------------------------------------------------------------------
Goal of the script :
Create a slack of images generating the visual design of a pRF experiment
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: group of shared data (e.g. 327)
-----------------------------------------------------------------------------------------
Output(s) :
Video of the vdm
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/stereo_prf/analysis_code/postproc/prf
2. run python command
python vdm_builder.py [main directory] [project name] [group]
-----------------------------------------------------------------------------------------
Exemple:
python vdm_builder.py /scratch/mszinte/data amblyo_prf 327
-----------------------------------------------------------------------------------------
Written by Martin Szinte (mail@martinszinte.net)
-----------------------------------------------------------------------------------------
"""

#General imports
import cv2
import json
import math
import numpy as np
import os
import sys
sys.path.append("{}/../../utils".format(os.getcwd()))
from conversion_utils import conversion

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
group = sys.argv[3]

# Load settings
with open('../../settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)

# Define directories
rootpath = "{}/{}/derivatives/vdm/".format(main_dir, project_dir)
fileName = 'vdm'
filepath = os.path.join(rootpath,fileName+'.npy')
videopath = os.path.join(rootpath,fileName+'.mp4')
os.makedirs(rootpath, exist_ok=True)

# Parameters
screen_converter = conversion(screen_size_pix = analysis_info['screen_size_pix'], 
                      screen_size_cm = analysis_info['screen_size_cm'],
                      screen_distance_cm = analysis_info['screen_distance_cm'])

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

# Downsampling : resize the frames and inverse y axis
frames_reshape = np.zeros((vdm_size_pix[0],vdm_size_pix[1],TRs))
for k in range(frames_reshape.shape[-1]):
    frames_reshape[:,:,k] = cv2.resize(frames[:,:,k], dsize=(vdm_size_pix[0], vdm_size_pix[1]), interpolation=cv2.INTER_NEAREST)
frames = frames_reshape

frames_rotate = np.zeros((vdm_size_pix[0],vdm_size_pix[1],TRs))
for num, frame in enumerate(np.split(frames, TRs, axis=2)):
    frames_rotate[:,:,num] = frame[-1::-1,:,0]
frames = frames_rotate

# Export video
fourcc = cv2.VideoWriter_fourcc(*'mp4v')
out = cv2.VideoWriter(videopath, fourcc, 1/TR, (vdm_size_pix[0], vdm_size_pix[1]), False)
[out.write(np.uint8(frame*255)) for frame in np.split(frames, TRs, axis=2)]
out.release()
print('Video conversion done, save to:'+videopath)

# Save numpy array
np.save(filepath,frames)
print('Data saved to :'+filepath)

# Define permission cmd
os.system("chmod -Rf 771 {main_dir}/{project_dir}".format(main_dir=main_dir, project_dir=project_dir))
os.system("chgrp -Rf {group} {main_dir}/{project_dir}".format(main_dir=main_dir, project_dir=project_dir, group=group))