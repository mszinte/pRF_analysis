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
python pycortex_maps_rois.py ~/disks/meso_S/data RetinoMaps sub-01 n
python pycortex_maps_rois.py ~/disks/meso_S/data RetinoMaps sub-170k n
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
import yaml
import cortex
import matplotlib.pyplot as plt

# Personal imports
sys.path.append("{}/../../../utils".format(os.getcwd()))
from pycortex_utils import draw_cortex, set_pycortex_config_file, load_surface_pycortex, create_colormap
from settings_utils import load_settings

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

# Load settings
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../../"))
settings_path = os.path.join(base_dir, project_dir, "settings.yml")
prf_settings_path = os.path.join(base_dir, project_dir, "prf-analysis.yml")
figure_settings_path = os.path.join(base_dir, project_dir, "figure-settings.yml")
settings = load_settings([settings_path, prf_settings_path, figure_settings_path])
analysis_info = settings[0]
if subject == 'sub-170k': formats = ['170k']
else: formats = analysis_info['formats']
preproc_prep = analysis_info['preproc_prep']
filtering = analysis_info['filtering']
normalization = analysis_info['normalization']
rois_methods = analysis_info['rois_methods']
pycortex_subject_template = analysis_info['pycortex_subject_template']

# Set pycortex db and colormaps
cortex_dir = "{}/{}/derivatives/pp_data/cortex".format(main_dir, project_dir)
set_pycortex_config_file(cortex_dir)

# Define/create colormap
rois_colors = analysis_info['rois_colors']
colormap_name = 'rois_colors'
colormap_dict = {key: tuple(value) for key, value in rois_colors.items()}
create_colormap(cortex_dir=cortex_dir, 
                colormap_name=colormap_name, 
                colormap_dict=colormap_dict,
                recreate=True
               )

# Create flatmaps
for format_, pycortex_subject in zip(formats, [subject, pycortex_subject_template]):
    
    # Define directories and fn
    rois_dir = "{}/{}/derivatives/pp_data/{}/{}/rois".format(main_dir, project_dir, subject, format_)
    flatmaps_dir = '{}/pycortex/flatmaps_rois'.format(rois_dir)
    datasets_dir = '{}/pycortex/datasets_rois'.format(rois_dir)
    
    os.makedirs(flatmaps_dir, exist_ok=True)
    os.makedirs(datasets_dir, exist_ok=True)

    rois_methods_format = rois_methods[format_]
    
    for rois_method_format in rois_methods_format:
    
        if format_ == 'fsnative':
            roi_fn_L = '{}/{}_hemi-L_{}_{}_{}_{}.func.gii'.format(rois_dir, subject, 
                                                                  preproc_prep, filtering, 
                                                                  normalization, rois_method_format)
            roi_fn_R = '{}/{}_hemi-R_{}_{}_{}_{}.func.gii'.format(rois_dir, subject,
                                                                  preproc_prep, filtering, 
                                                                  normalization, rois_method_format)
            results = load_surface_pycortex(L_fn=roi_fn_L, R_fn=roi_fn_R)
            roi_mat = results['data_concat']
            
        elif format_ == '170k':
            roi_fn = '{}/{}_{}_{}_{}_{}.dtseries.nii'.format(rois_dir, subject,
                                                      preproc_prep, filtering, 
                                                      normalization, rois_method_format)
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
                      'vmax': len(colormap_dict),
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
        plt.savefig('{}/{}_{}_{}_{}_{}_rois.pdf'.format(flatmaps_dir, subject,
                                                        preproc_prep, filtering, 
                                                        normalization, rois_method_format))
        plt.close()
    
        # Save flatmap as dataset
        volumes = {}
        vol_description = param_rois["description"]
        volumes.update({vol_description:volume_roi})
    
        # Save dataset
        dataset_file = "{}/{}_{}_{}_{}_{}.hdf".format(datasets_dir, subject,
                                                      preproc_prep, filtering, 
                                                      normalization, rois_method_format)
        if os.path.exists(dataset_file): os.system("rm -fv {}".format(dataset_file))
        dataset = cortex.Dataset(data=volumes)
        dataset.save(dataset_file)