
"""
-----------------------------------------------------------------------------------------
save_design_coords.py
-----------------------------------------------------------------------------------------
Goal of the script:
save general screen design coordinates of experiment as usable numpy arrays 
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject name
sys.argv[4]: group of shared data (e.g. 327)
-----------------------------------------------------------------------------------------
Output(s):

-----------------------------------------------------------------------------------------
To run : 
cd ~/projects/pRF_analysis/RetinoMaps/eyetracking/preproc/
python save_design_coords.py /scratch/mszinte/data RetinoMaps 327
-----------------------------------------------------------------------------------------
Written by Sina Kling
-----------------------------------------------------------------------------------------
"""
# Stop warnings
import warnings
warnings.filterwarnings("ignore")

# Debug
import ipdb
deb = ipdb.set_trace

# General imports
import numpy as np 
import scipy.io

# path of utils folder  
sys.path.append("{}/../../../analysis_code/utils".format(os.getcwd()))
from settings_utils import load_settings
from eyetrack_utils import load_event_files

# inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
group = sys.argv[3]


mat = scipy.io.loadmat(f'{main_dir}/{project_dir}/sourcedata/sub-01/ses-01/add/sub-01_ses-01_task-SacLoc_run-01_matFile.mat') 
# all subjects have same coordinates saved: just amplitude reading is depending on design! 
config = mat['config'] 
print(config.shape) 
print(config.dtype) 
cfg0 = config[0][0] 
const = cfg0['const'][0][0]


design_coordinates_X = const['saccade_matX'][0]
design_coordinates_Y = const['saccade_matY'][0]

np.save(f"{main_dir}/{project_dir}/exp_design/design_coordinates_x", design_coordinates_X)
np.save(f"{main_dir}/{project_dir}/exp_design/design_coordinates_y", design_coordinates_Y)


# Define permission cmd
print('Changing files permissions in {}/{}'.format(main_dir, project_dir))
os.system("chmod -Rf 771 {}/{}".format(main_dir, project_dir))
os.system("chgrp -Rf {} {}/{}".format(group, main_dir, project_dir))

