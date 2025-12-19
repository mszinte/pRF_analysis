"""
-----------------------------------------------------------------------------------------
pycortex_maps_rois.py
-----------------------------------------------------------------------------------------
Goal of the script:
Create flatmap plots and dataset of ROIs
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject name (e.g. sub-01)
sys.argv[4]: save in svg (e.g. no)
-----------------------------------------------------------------------------------------
Output(s):
Pycortex flatmaps figures and dataset
-----------------------------------------------------------------------------------------
To run:
0. TO RUN ON INVIBE SERVER (with Inkscape)
1. cd to function
>> cd ~/disks/meso_H/projects/pRF_analysis/analysis_code/postproc/prf/postfit/
2. run python command
>> python pycortex_maps_rois.py [main directory] [project name] [subject num] 
                                [save_in_svg]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/disks/meso_H/projects/pRF_analysis/analysis_code/postproc/prf/postfit/
python pycortex_maps_rois.py ~/disks/meso_S/data MotConf sub-01 n
-----------------------------------------------------------------------------------------
Written by Uriel Lascombes (uriel.lascombes@laposte.net)
Edited by Martin Szinte (martin.szinte@gmail.com)
-----------------------------------------------------------------------------------------
"""
# Stop warnings
import warnings
warnings.filterwarnings("ignore")

# Debug import 
import ipdb
deb = ipdb.set_trace

# General imports
import os
import sys
import json
import cortex
import matplotlib.pyplot as plt

# Personal imports
sys.path.append("{}/../../../utils".format(os.getcwd()))
from pycortex_utils import draw_cortex, set_pycortex_config_file, load_surface_pycortex, create_colormap

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]
subject = sys.argv[3]
save_svg_in = sys.argv[4]

try:
    if save_svg_in == 'yes' or save_svg_in == 'y':
        save_svg = True
    elif save_svg_in == 'no' or save_svg_in == 'n':
        save_svg = False
    else:
        raise ValueError
except ValueError:
    sys.exit('Error: incorrect input (Yes, yes, y or No, no, n)')
if subject == 'sub-170k': save_svg = False
else: save_svg = save_svg

# Define analysis parameters
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../../"))
settings_path = os.path.join(base_dir, project_dir, "settings.json")

with open(settings_path) as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)
if subject == 'sub-170k': formats = ['170k']
else: formats = analysis_info['formats']

# Set pycortex db and colormaps
cortex_dir = "{}/{}/derivatives/pp_data/cortex".format(main_dir, project_dir)
set_pycortex_config_file(cortex_dir)

# Define/create colormap
colormap_name = 'rois_colors'
colormap_dict = {'n/a': (255, 255, 255),
                 'V1': (243, 231, 155),
                 'V2': (250, 196, 132),
                 'V3': (248, 160, 126),
                 'V3AB': (235, 127, 134),
                 'LO': (150, 0, 90), 
                 'VO': (0, 0, 200),
                 'hMT+': (0, 25, 255),
                 'iIPS': (0, 152, 255),
                 'sIPS': (44, 255, 150),
                 'iPCS': (151, 255, 0),
                 'sPCS': (255, 234, 0),
                 'mPCS': (255, 111, 0)}

create_colormap(cortex_dir=cortex_dir, 
                colormap_name=colormap_name, 
                colormap_dict=colormap_dict,
                recreate=True
               )

# Create flatmaps
for format_, pycortex_subject in zip(formats, [subject, 'sub-170k']):
    # Define directories and fn
    rois_dir = "{}/{}/derivatives/pp_data/{}/{}/rois".format(main_dir, project_dir, subject, format_)
    flatmaps_dir = '{}/pycortex/flatmaps_rois'.format(rois_dir)
    datasets_dir = '{}/pycortex/datasets_rois'.format(rois_dir)
    
    os.makedirs(flatmaps_dir, exist_ok=True)
    os.makedirs(datasets_dir, exist_ok=True)
    
    if format_ == 'fsnative':
        roi_fn_L = '{}/{}_hemi-L_rois.func.gii'.format(rois_dir, subject)
        roi_fn_R = '{}/{}_hemi-R_rois.func.gii'.format(rois_dir, subject)
        results = load_surface_pycortex(L_fn=roi_fn_L, 
                                        R_fn=roi_fn_R)
        roi_mat = results['data_concat']
        
    elif format_ == '170k':
        roi_fn = '{}/{}_rois.dtseries.nii'.format(rois_dir, subject)
        results = load_surface_pycortex(brain_fn=roi_fn)
        roi_mat = results['data_concat']
        if subject == 'sub-170k': save_svg = save_svg
        else: save_svg = False
    rois_opacity = 0.5
    alpha_mat = roi_mat*0+rois_opacity
    alpha_mat[roi_mat==0]=0
    print('Creating flatmaps...')

    # Rois
    roi_name = 'rois'
    param_rois = {'subject': pycortex_subject,
                  'data': roi_mat, 
                  'cmap': colormap_name,
                  'alpha': alpha_mat,
                  'cbar': 'discrete_personalized', 
                  'vmin': 0,
                  'vmax': 12,
                  'cmap_steps': len(colormap_dict),
                  'cmap_dict': colormap_dict,
                  'cortex_type': 'VertexRGB',
                  'description': 'ROIs',
                  'curv_brightness': 1, 
                  'curv_contrast': 0.25,
                  'add_roi': save_svg,
                  'with_labels': True,
                  'roi_name': roi_name}
                  
    # Draw flatmaps
    volume_roi = draw_cortex(**param_rois)
    plt.savefig('{}/{}_rois.pdf'.format(flatmaps_dir, subject))
    plt.close()

    # Save flatmap as dataset
    volumes = {}
    vol_description = param_rois["description"]
    volumes.update({vol_description:volume_roi})

    # Save dataset
    dataset_file = "{}/{}_rois.hdf".format(datasets_dir, subject)
    if os.path.exists(dataset_file): os.system("rm -fv {}".format(dataset_file))
    dataset = cortex.Dataset(data=volumes)
    dataset.save(dataset_file)