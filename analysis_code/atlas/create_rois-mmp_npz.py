"""
-----------------------------------------------------------------------------------------
create_rois-mmp_npz.py
-----------------------------------------------------------------------------------------
Goal of the script:
Create 170k npz mmp atlas files and a tsv of rois numbers and names
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: data project directory
sys.argv[1]: codee project directory
sys.argv[2]: project name (correspond to directory)
-----------------------------------------------------------------------------------------
Output(s):
.npz masks for roi-mmp and roi-group-mmp for both 91k and 170k template 
-----------------------------------------------------------------------------------------
To run:
1. cd to function
>> cd ~/projects/[PROJECT]/analysis_code/atlas/
2. run python command
python create_rois-mmp_npz.py [main directory] [project name] [subject name] [group]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/projects/pRF_analysis/analysis_code/atlas/
python create_rois-mmp_npz.py /scratch/mszinte/data RetinoMaps 327
-----------------------------------------------------------------------------------------
Written by Uriel Lascombes (uriel.lascombes@laposte.net)
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
import yaml
import cortex
import numpy as np
import pandas as pd

# Personal imports
sys.path.append("{}/../../analysis_code/utils".format(os.getcwd()))
from settings_utils import load_settings
from cifti_utils import from_91k_to_32k, from_170k_to_59k
from surface_utils import load_surface, make_surface_image
from pycortex_utils import load_surface_pycortex, set_pycortex_config_file, draw_cortex


