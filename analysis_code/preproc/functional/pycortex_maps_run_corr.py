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
>> cd ~/disks/meso_H/projects/[PROJECT]/analysis_code/preproc/functional/
2. run python command
>> python pycortex_maps_run_cor.py [main directory] [project name] [subject num] [save_svg]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/disks/meso_H/projects/pRF_analysis/analysis_code/preproc/functional/
python pycortex_maps_run_corr.py ~/disks/meso_S/data MotConf sub-01 n
python pycortex_maps_run_corr.py ~/disks/meso_S/data MotConf sub-170k n
-----------------------------------------------------------------------------------------
Written by Martin Szinte (mail@martinszinte.net)
Edited by Uriel Lascombes (uriel.lascombes@laposte.net)
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
import numpy as np
import copy
import cortex
import matplotlib.pyplot as plt

# Personal imports
sys.path.append("{}/../../utils".format(os.getcwd()))
from pycortex_utils import draw_cortex, set_pycortex_config_file, load_surface_pycortex

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

# Define analysis parameters
with open('../../settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)
if subject == 'sub-170k': formats = ['170k']
else: formats = analysis_info['formats']
extensions = analysis_info['extensions']
tasks = analysis_info['task_names']
alpha_range = analysis_info["alpha_range"]
maps_names_corr = analysis_info['maps_names_corr']

# Set pycortex db and colormaps
cortex_dir = "{}/{}/derivatives/pp_data/cortex".format(main_dir, project_dir)
set_pycortex_config_file(cortex_dir)

# Maps settings
cmap_corr = 'BuBkRd'

# Index
for idx, col_name in enumerate(maps_names_corr):
    exec("{}_idx = idx".format(col_name))

# Plot scales
corr_scale = [-1, 1]
 
for format_, pycortex_subject in zip(formats, [subject, 'sub-170k']):

    corr_dir = "{}/{}/derivatives/pp_data/{}/{}/corr/fmriprep_dct_corr".format(main_dir, project_dir, subject, format_)
    flatmaps_dir = '{}/{}/derivatives/pp_data/{}/{}/corr/pycortex/flatmaps_inter-run-corr'.format(main_dir, project_dir, subject, format_)
    datasets_dir = '{}/{}/derivatives/pp_data/{}/{}/corr/pycortex/datasets_inter-run-corr'.format(main_dir, project_dir, subject, format_)

    os.makedirs(flatmaps_dir, exist_ok=True)
    os.makedirs(datasets_dir, exist_ok=True)

    for task in tasks : 
        print(task)
        
        if format_ == 'fsnative':
            corr_fn_L = "{}/{}_task-{}_hemi-L_fmriprep_dct_corr_bold.func.gii".format(corr_dir, subject, task)
            corr_fn_R = "{}/{}_task-{}_hemi-R_fmriprep_dct_corr_bold.func.gii".format(corr_dir, subject, task)
            results = load_surface_pycortex(L_fn=corr_fn_L, R_fn=corr_fn_R)
        elif format_ == '170k':
            cor_fn = '{}/{}_task-{}_fmriprep_dct_corr_bold.dtseries.nii'.format(corr_dir, subject, task)
            results = load_surface_pycortex(brain_fn=cor_fn)
            if subject == 'sub-170k':
                save_svg = save_svg
            else: 
                save_svg = False
        corr_mat = results['data_concat']
        maps_names = []        
        
        # Correlation uncorrected
        corr_mat_uncorrected = corr_mat[rvalue_idx, :]
        
        # Compute alpha
        alpha_uncorrected = np.abs(corr_mat_uncorrected)
        alpha_uncorrected = (alpha_uncorrected - alpha_range[0]) / (alpha_range[1] - alpha_range[0])
        alpha_uncorrected[alpha_uncorrected>1] = 1
        
        # correlation uncorrected
        param_run_corr = {'data': corr_mat_uncorrected, 
                          'alpha': alpha_uncorrected,
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
                          'with_labels': True}
        maps_names.append('run_corr')

        # Correlation corrected mat
        corr_mat_corrected = copy.copy(corr_mat)
        corr_mat_corrected_th = corr_mat_corrected
        if analysis_info['stats_th'] == 0.05: stats_th_down = corr_mat_corrected_th[corr_pvalue_5pt_idx,...] <= 0.05
        elif analysis_info['stats_th'] == 0.01: stats_th_down = corr_mat_corrected_th[corr_pvalue_1pt_idx,...] <= 0.01
        corr_mat_corrected[rvalue_idx, stats_th_down==False]=0 # put this to zero to not plot it
        corr_mat_corrected = corr_mat_corrected[rvalue_idx, :]

        # Compute alpha
        alpha_corrected = np.abs(corr_mat_corrected)
        alpha_corrected = (alpha_corrected - alpha_range[0]) / (alpha_range[1] - alpha_range[0])
        alpha_corrected[alpha_corrected>1]=1
        
        # correlation corrected
        param_run_corr_stats = {'data': corr_mat_corrected, 
                                'alpha': alpha_corrected,
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
                                'with_labels': True}
        maps_names.append('run_corr_stats')

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
            exec("plt.savefig('{}/{}_task-{}_{}.pdf')".format(flatmaps_dir, subject, task, maps_name))
            plt.close()
        
            # save flatmap as dataset
            exec('vol_description = param_{}["description"]'.format(maps_name))
            exec('volume = volume_{}'.format(maps_name))
            volumes.update({vol_description:volume})
        
        # save dataset
        dataset_file = "{}/{}_task-{}_inter-run-corr.hdf".format(datasets_dir, subject, task)
        if os.path.exists(dataset_file): os.system("rm -fv {}".format(dataset_file))
        dataset = cortex.Dataset(data=volumes)
        dataset.save(dataset_file)