"""
-----------------------------------------------------------------------------------------
sagital_view.py
-----------------------------------------------------------------------------------------
Goal of the script:
Make freeview sagital video of segmentation (to run before and after manual edit)
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject name (e.g. sub-01)
sys.argv[4]: video name ('before_edit', 'after_edit')
sys.argv[5]: session name (optional, e.g. ses-01)
-----------------------------------------------------------------------------------------
Output(s):
Sagital view video and images of the brain segmentation
-----------------------------------------------------------------------------------------
To run:
0. TO RUN LOCALLY WITH FREEWIEW INSTALLED (not on server)
1. cd to function
>> cd ~/disks/meso_H/projects/pRF_analysis/analysis_code/preproc/anatomical/
2. run python command
python sagital_view.py [main directory] [project name] [subject num] [video name] [session (optional)]
-----------------------------------------------------------------------------------------
Exemple:
python sagital_view.py ~/disks/meso_S/data/ RetinoMaps sub-01 before_edit
python sagital_view.py ~/disks/meso_S/data/ RetinoMaps sub-01 before_edit ses-01
-----------------------------------------------------------------------------------------
# Written by Martin Szinte (martin.szinte@gmail.com)
Edited by Uriel Lascombes (uriel.lascombes@laposte.net)
-----------------------------------------------------------------------------------------
"""
# imports modules
import os
import sys
import ipdb
import numpy as np
import subprocess as sb
deb = ipdb.set_trace

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
vid_name = sys.argv[4]
session = sys.argv[5] if len(sys.argv) > 5 else None

# define directory
if session:
    subject_name = f"{subject}_{session}"
else:
    subject_name = subject

fs_dir = "{}/{}/derivatives/fmriprep/freesurfer/{}".format(main_dir, project_dir, subject_name)
vid_dir = "{}/vid/{}".format(fs_dir, vid_name)
os.makedirs(vid_dir, exist_ok=True)
image_dir = "{}/img".format(vid_dir)
os.makedirs(image_dir, exist_ok=True)

# list commands
anat_cmd = '-v {}:grayscale=10,100'.format('{}/mri/T1.mgz'.format(fs_dir))
volumes_cmd = '-f {fs_dir}/surf/lh.white:color=red:edgecolor=red \
-f {fs_dir}/surf/rh.white:color=red:edgecolor=red \
-f {fs_dir}/surf/lh.pial:color=white:edgecolor=white \
-f {fs_dir}/surf/rh.pial:color=white:edgecolor=white \
-layout 1'.format(fs_dir = fs_dir)

slice_cmd = ''
x_start, x_end = 50, 210
for x in np.arange(x_start,x_end):
    if x < 10:x_name = '00{}'.format(x)
    elif x >= 10 and x < 100:x_name = '0{}'.format(x)
    else: x_name = '{}'.format(x)

    slice_cmd += ' -slice {} 127 127 \n -ss {}/{}.png \n'.format(x,image_dir,x_name)

# main command
freeview_cmd = '{} {} -viewport sagittal {}-quit '.format(anat_cmd, volumes_cmd, slice_cmd)

# save command in shell file
sh_fn = "{}/{}_{}.sh".format(vid_dir, subject, vid_name)
of = open(sh_fn, 'w')
of.write(freeview_cmd)
of.close()

# run freeview cmd
sb.call('freeview -cmd {}'.format(sh_fn), shell=True)

# convert images in video
mk_vid_cmd = 'ffmpeg -framerate 5 -pattern_type glob -i "{}/*.png" -b:v 2M -c:v mpeg4 {}/{}_{}.mp4'.format(
    image_dir, vid_dir, subject, vid_name)
os.system(mk_vid_cmd)
