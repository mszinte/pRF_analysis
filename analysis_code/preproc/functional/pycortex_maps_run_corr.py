"""
-----------------------------------------------------------------------------------------
pycortex_maps_run_corr.py
-----------------------------------------------------------------------------------------
Goal of the script:
Create flatmap plots of inter-run correlations
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
sys.argv[3]: subject name (e.g. sub-01)
sys.argv[4]: save map in svg (y/n)
-----------------------------------------------------------------------------------------
Output(s):
Pycortex flatmaps figures
-----------------------------------------------------------------------------------------
To run:
0. TO RUN ON INVIBE SERVER (with Inkscape)
1. cd to function
>> cd ~/disks/meso_H/projects/pRF_analysis/analysis_code/preproc/functional/
2. run python command
>> python pycortex_maps_run_cor.py [main directory] [project name] [subject num] [save_svg]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/disks/meso_H/projects/pRF_analysis/analysis_code/preproc/functional/
python pycortex_maps_run_corr.py ~/disks/meso_S/data RetinoMaps sub-01 n
python pycortex_maps_run_corr.py ~/disks/meso_S/data RetinoMaps sub-170k n
-----------------------------------------------------------------------------------------
Written by Martin Szinte (martin.szinte@gmail.com)
and Uriel Lascombes (uriel.lascombes@laposte.net)
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
import copy
import cortex
import matplotlib.pyplot as plt
import numpy as np

# Personal imports
sys.path.append("{}/../../utils".format(os.getcwd()))
from pycortex_utils import draw_cortex, set_pycortex_config_file, load_surface_pycortex
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

# Define analysis parameters
base_dir = os.path.abspath(os.path.join(os.getcwd(), "../../../"))
settings_path = os.path.join(base_dir, project_dir, "settings.yml")
figure_settings_path = os.path.join(base_dir, project_dir, "figure-settings.yml")
settings = load_settings([settings_path, figure_settings_path])
analysis_info = settings[0]

if subject == 'sub-170k': formats = ['170k']
else: formats = analysis_info['formats']
extensions = analysis_info['extensions']
tasks = analysis_info['task_names']
maps_names_corr = analysis_info['maps_names_corr']
preproc_prep = analysis_info['preproc_prep']
filtering = analysis_info['filtering']
normalization = analysis_info['normalization']
corr_scale = analysis_info['flatmap_corr_scale']
pycortex_subject_template = analysis_info['pycortex_subject_template']

# Set pycortex db and colormaps
cortex_dir = "{}/{}/derivatives/pp_data/cortex".format(main_dir, project_dir)
set_pycortex_config_file(cortex_dir)

# Maps settings
cmap_corr = 'BuBkRd'

# Index
for idx, col_name in enumerate(maps_names_corr):
    exec("{}_idx = idx".format(col_name))

# Plot scales

for format_, pycortex_subject in zip(formats, [subject, pycortex_subject_template]):
    corr_dir = "{}/{}/derivatives/pp_data/{}/{}/corr/{}_{}_{}_corr".format(
        main_dir, project_dir, subject, format_,
        preproc_prep, filtering, normalization)
    flatmaps_dir = '{}/{}/derivatives/pp_data/{}/{}/corr/pycortex/flatmaps_corr'.format(
        main_dir, project_dir, subject, format_)
    datasets_dir = '{}/{}/derivatives/pp_data/{}/{}/corr/pycortex/datasets_corr'.format(
        main_dir, project_dir, subject, format_)

    os.makedirs(flatmaps_dir, exist_ok=True)
    os.makedirs(datasets_dir, exist_ok=True)

    for task in tasks : 
        
        if format_ == 'fsnative':
            corr_fn_L = "{}/{}_task-{}_hemi-L_{}_{}_{}_corr_bold.func.gii".format(
                corr_dir, subject, task, preproc_prep, filtering, normalization)
            corr_fn_R = "{}/{}_task-{}_hemi-R_{}_{}_{}_corr_bold.func.gii".format(
                corr_dir, subject, task, preproc_prep, filtering, normalization)
            results = load_surface_pycortex(L_fn=corr_fn_L, R_fn=corr_fn_R)
        elif format_ == '170k':
            cor_fn = '{}/{}_task-{}_{}_{}_{}_corr_bold.dtseries.nii'.format(
                corr_dir, subject, task, preproc_prep, filtering, normalization)
            results = load_surface_pycortex(brain_fn=cor_fn)
            if subject == 'sub-170k': save_svg = save_svg
            else: save_svg = False
        corr_mat = results['data_concat']
        maps_names = []
        
        # Correlation uncorrected
        corr_mat_uncorrected = corr_mat[rvalue_idx, :]
        
        # correlation uncorrected
        param_corr = {'data': corr_mat_uncorrected, 
                      'alpha': corr_mat_uncorrected*0+1, # no alpha correction
                      'cmap': cmap_corr,
                      'vmin': corr_scale[0], 
                      'vmax': corr_scale[1], 
                      'cbar': 'discrete', 
                      'cortex_type': 'VertexRGB', 
                      'description': 'Inter-run correlation (uncorrected): task-{}'.format(task), 
                      'curv_brightness': 0.1, 
                      'curv_contrast': 0.25, 
                      'add_roi': save_svg, 
                      'cbar_label': 'Pearson coefficient',
                      'with_labels': True,
                      'overlay_fn': 'overlay.svg'
                     }
        maps_names.append('corr')

        # Correlation corrected mat
        corr_mat_corrected = copy.copy(corr_mat)
        corr_mat_corrected_th = corr_mat_corrected
        if analysis_info['corr_stats_th'] == 0.05: stats_th_down = corr_mat_corrected_th[corr_pvalue_5pt_idx,...] <= 0.05
        elif analysis_info['corr_stats_th'] == 0.01: stats_th_down = corr_mat_corrected_th[corr_pvalue_1pt_idx,...] <= 0.01
        corr_mat_corrected[rvalue_idx, stats_th_down==False]=0 # put this to zero to not plot it
        corr_mat_corrected = corr_mat_corrected[rvalue_idx, :]
        
        # correlation corrected
        param_corr_stats = {'data': corr_mat_corrected, 
                            'alpha': corr_mat_corrected*0+1,
                            'cmap': cmap_corr ,
                            'vmin': corr_scale[0],
                            'vmax': corr_scale[1], 
                            'cbar': 'discrete', 
                            'cortex_type': 'VertexRGB', 
                            'description': 'Inter-run correlation (corrected): task-{}'.format(task),
                            'curv_brightness': 0.1, 
                            'curv_contrast': 0.25, 
                            'add_roi': save_svg, 
                            'cbar_label': 'Pearson coefficient',
                            'with_labels': True,
                            'overlay_fn': 'overlay.svg'
                           }
        maps_names.append('corr_stats')

        # draw flatmaps
        volumes = {}
        for maps_name in maps_names:
            
            # create flatmap
            roi_name = '{}_{}'.format(task, maps_name)
            roi_param = {'subject': pycortex_subject, 
                         'roi_name': roi_name}
            print(roi_name)
            exec('param_{}.update(roi_param)'.format(maps_name))
            exec('volume_{maps_name} = draw_cortex(**param_{maps_name})'.format(maps_name=maps_name))
            exec("plt.savefig('{}/{}_task-{}_{}_{}_{}_{}.pdf')".format(
                flatmaps_dir, subject, task, 
                preproc_prep, filtering, normalization, maps_name))
            plt.close()
        
            # save flatmap as dataset
            exec('vol_description = param_{}["description"]'.format(maps_name))
            exec('volume = volume_{}'.format(maps_name))
            volumes.update({vol_description:volume})

        # save dataset
        dataset_file = "{}/{}_task-{}_{}_{}_{}_corr.hdf".format(datasets_dir, subject, task, 
                                                                preproc_prep, filtering, normalization)
        if os.path.exists(dataset_file): os.system("rm -fv {}".format(dataset_file))
        dataset = cortex.Dataset(data=volumes)
        dataset.save(dataset_file)