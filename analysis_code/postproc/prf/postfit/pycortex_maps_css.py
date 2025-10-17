"""
-----------------------------------------------------------------------------------------
pycortex_maps_css.py
-----------------------------------------------------------------------------------------
Goal of the script:
Create flatmap plots and dataset of css results
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
>> cd ~/disks/meso_H/projects/[PROJECT]/analysis_code/postproc/prf/postfit/
2. run python command
>> python pycortex_maps_css.py [main directory] [project] [subject] [save_svg_in]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/disks/meso_H/projects/pRF_analysis/analysis_code/postproc/prf/postfit/

python pycortex_maps_css.py ~/disks/meso_S/data MotConf sub-01 n
python pycortex_maps_css.py ~/disks/meso_S/data MotConf sub-170k n

python pycortex_maps_css.py ~/disks/meso_S/data RetinoMaps sub-01 n
python pycortex_maps_css.py ~/disks/meso_S/data RetinoMaps sub-170k n

python pycortex_maps_css.py ~/disks/meso_S/data amblyo_prf sub-01 n
python pycortex_maps_css.py ~/disks/meso_S/data amblyo_prf sub-170k n

python pycortex_maps_css.py ~/disks/meso_S/data centbids sub-2100247523 n
-----------------------------------------------------------------------------------------
Written by Martin Szinte (mail@martinszinte.net)
Edited by Uriel Lascombes (uriel.lascombes@laposte.net)
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
import json
import cortex
import numpy as np
import matplotlib.pyplot as plt
import glob

# Personal import
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
extensions = analysis_info['extensions']
prf_task_name = analysis_info['prf_task_name']
maps_names_css = analysis_info['maps_names_css']
maps_names_pcm = analysis_info['maps_names_pcm']
maps_names_css_stats = analysis_info['maps_names_css_stats']

# Set pycortex db and colormaps
cortex_dir = "{}/{}/derivatives/pp_data/cortex".format(main_dir, project_dir)
set_pycortex_config_file(cortex_dir)

# Maps settings 
for idx, col_name in enumerate(maps_names_css + maps_names_css_stats + maps_names_pcm):
    exec("{}_idx = idx".format(col_name))

cmap_polar, cmap_uni, cmap_ecc_size = 'hsv', 'Reds', 'Spectral'
col_offset = 1.0/14.0
cmap_steps = 255

# plot scales
rsq_scale = [0, 1]
ecc_scale = [0, 8]
size_scale = [0, 5]
n_scale = [0, 2]
pcm_scale = [0, 10]

for format_ in formats:
    
    # Define directories and fn
    prf_dir = "{}/{}/derivatives/pp_data/{}/{}/prf".format(main_dir, project_dir, subject, format_)
    prf_deriv_dir = "{}/prf_derivatives".format(prf_dir)
    flatmaps_dir = '{}/pycortex/flatmaps_css_loo-median'.format(prf_dir)
    datasets_dir = '{}/pycortex/datasets_css_loo-median'.format(prf_dir)
    
    os.makedirs(flatmaps_dir, exist_ok=True)
    os.makedirs(datasets_dir, exist_ok=True)

    # Load all data
    if format_ == 'fsnative':
        pycortex_subject = subject
        # Derivatives
        deriv_median_fn_L = '{}/{}_task-{}_hemi-L_fmriprep_dct_avg_prf-deriv_css_loo-median.func.gii'.format(
            prf_deriv_dir, subject, prf_task_name)
        deriv_median_fn_R = '{}/{}_task-{}_hemi-R_fmriprep_dct_avg_prf-deriv_css_loo-median.func.gii'.format(
            prf_deriv_dir, subject, prf_task_name)
        deriv_results = load_surface_pycortex(L_fn=deriv_median_fn_L, R_fn=deriv_median_fn_R)
        deriv_mat = deriv_results['data_concat']
        
        # Stats
        stats_median_fn_L = '{}/{}_task-{}_hemi-L_fmriprep_dct_avg_prf-stats_loo-median.func.gii'.format(
            prf_deriv_dir, subject, prf_task_name)
        stats_median_fn_R = '{}/{}_task-{}_hemi-R_fmriprep_dct_avg_prf-stats_loo-median.func.gii'.format(
            prf_deriv_dir, subject, prf_task_name)
        stats_results = load_surface_pycortex(L_fn=stats_median_fn_L, R_fn=stats_median_fn_R)
        stats_mat = stats_results['data_concat']

        # pRF CM
        pcm_median_fn_L = '{}/{}_task-{}_hemi-L_fmriprep_dct_avg_prf-pcm_css_loo-median.func.gii'.format(
            prf_deriv_dir, subject, prf_task_name)
        pcm_median_fn_R = '{}/{}_task-{}_hemi-R_fmriprep_dct_avg_prf-pcm_css_loo-median.func.gii'.format(
            prf_deriv_dir, subject, prf_task_name)
        pcm_results = load_surface_pycortex(L_fn=pcm_median_fn_L, R_fn=pcm_median_fn_R)
        pcm_mat = pcm_results['data_concat']
        
    elif format_ == '170k':
        pycortex_subject = 'sub-170k'
        # Derivatives
        deriv_median_fn = '{}/{}_task-{}_fmriprep_dct_avg_prf-deriv_css_loo-median.dtseries.nii'.format(
            prf_deriv_dir, subject, prf_task_name)
        deriv_results = load_surface_pycortex(brain_fn=deriv_median_fn)
        deriv_mat = deriv_results['data_concat']

        # Stats
        stats_median_fn = '{}/{}_task-{}_fmriprep_dct_avg_prf-stats_loo-median.dtseries.nii'.format(
            prf_deriv_dir, subject, prf_task_name)
        stats_results = load_surface_pycortex(brain_fn=stats_median_fn)
        stats_mat = stats_results['data_concat']

        # pRF CM
        pcm_median_fn = '{}/{}_task-{}_fmriprep_dct_avg_prf-pcm_css_loo-median.dtseries.nii'.format(
            prf_deriv_dir, subject, prf_task_name)
        pcm_results = load_surface_pycortex(brain_fn=pcm_median_fn)
        pcm_mat = pcm_results['data_concat']
        
        if subject == 'sub-170k': save_svg = save_svg
        else: save_svg = False

    # Combine mat
    all_deriv_mat = np.concatenate((deriv_mat, stats_mat, pcm_mat))
    # Threshold mat
    all_deriv_mat_th = all_deriv_mat
    amp_down = all_deriv_mat_th[amplitude_idx,...] > 0
    rsq_down = all_deriv_mat_th[prf_loo_r2_idx,...] >= analysis_info['rsqr_th']
    size_th_down = all_deriv_mat_th[prf_size_idx,...] >= analysis_info['size_th'][0]
    size_th_up = all_deriv_mat_th[prf_size_idx,...] <= analysis_info['size_th'][1]
    ecc_th_down = all_deriv_mat_th[prf_ecc_idx,...] >= analysis_info['ecc_th'][0]
    ecc_th_up = all_deriv_mat_th[prf_ecc_idx,...] <= analysis_info['ecc_th'][1]
    n_th_down = all_deriv_mat_th[prf_n_idx,...] >= analysis_info['n_th'][0]
    n_th_up = all_deriv_mat_th[prf_n_idx,...] <= analysis_info['n_th'][1]
    if analysis_info['stats_th'] == 0.05: stats_th_down = all_deriv_mat_th[corr_pvalue_5pt_idx,...] <= 0.05
    elif analysis_info['stats_th'] == 0.01: stats_th_down = all_deriv_mat_th[corr_pvalue_1pt_idx,...] <= 0.01
    all_th = np.array((amp_down,
                       rsq_down,
                       size_th_down,size_th_up, 
                       ecc_th_down, ecc_th_up,
                       n_th_down, n_th_up#,stats_th_down
                      )) 
    all_deriv_mat[prf_loo_r2_idx, np.logical_and.reduce(all_th)==False]=0 # put this to zero to not plot it
    
    # Create flatmaps
    maps_names = []
    
    # Loo r-square and alpha
    loo_rsq_data = all_deriv_mat[prf_loo_r2_idx,...]
    alpha = loo_rsq_data
    alpha_range = analysis_info["alpha_range"]
    alpha = (alpha - alpha_range[0]) / (alpha_range[1] - alpha_range[0])
    alpha[alpha>1]=1

    param_loo_rsq = {'data': loo_rsq_data, 
                     'cmap': cmap_uni, 
                     'alpha': alpha, 
                     'vmin': rsq_scale[0], 
                     'vmax': rsq_scale[1], 
                     'cbar': 'discrete', 
                     'cortex_type': 'VertexRGB',
                     'description': 'CSS pRF LOO R2',
                     'curv_brightness':1,
                     'curv_contrast': 0.1,
                     'add_roi': save_svg,
                     'cbar_label': 'pRF LOO R2', 
                     'with_labels': True}
    maps_names.append('loo_rsq')
    
    # Polar angle
    pol_comp_num = all_deriv_mat[polar_real_idx,...] + 1j * all_deriv_mat[polar_imag_idx,...]
    polar_ang = np.angle(pol_comp_num)
    ang_norm = (polar_ang + np.pi) / (np.pi * 2.0)
    ang_norm = np.fmod(ang_norm + col_offset,1)
    param_polar = {'data': ang_norm, 
                   'cmap': cmap_polar, 
                   'alpha': alpha, 
                   'vmin': 0, 
                   'vmax': 1, 
                   'cmap_steps': cmap_steps, 
                   'cortex_type': 'VertexRGB',
                   'cbar': 'polar', 
                   'col_offset': col_offset, 
                   'description': 'CSS pRF polar angle', 
                   'curv_brightness': 0.1, 
                   'curv_contrast': 0.25, 
                   'add_roi': save_svg, 
                   'with_labels': True}
    exec('param_polar_{cmap_steps} = param_polar'.format(cmap_steps = int(cmap_steps)))
    exec('maps_names.append("polar_{cmap_steps}")'.format(cmap_steps = int(cmap_steps)))
    
    # Eccentricity
    ecc_data = all_deriv_mat[prf_ecc_idx,...]
    param_ecc = {'data': ecc_data, 
                 'cmap': cmap_ecc_size, 
                 'alpha': alpha,
                 'vmin': ecc_scale[0], 
                 'vmax': ecc_scale[1], 
                 'cbar': 'ecc', 
                 'cortex_type': 'VertexRGB',
                 'description': 'CSS pRF eccentricity', 
                 'curv_brightness': 1,
                 'curv_contrast': 0.1, 
                 'add_roi': save_svg, 
                 'with_labels': True}
    maps_names.append('ecc')
    
    # Size
    size_data = all_deriv_mat[prf_size_idx,...]
    param_size = {'data': size_data, 
                  'cmap': cmap_ecc_size, 
                  'alpha': alpha, 
                  'vmin': size_scale[0], 
                  'vmax': size_scale[1], 
                  'cbar': 'discrete', 
                  'cortex_type': 'VertexRGB', 
                  'description': 'CSS pRF size', 
                  'curv_brightness': 1, 
                  'curv_contrast': 0.1, 
                  'add_roi': False, 
                  'cbar_label': 'pRF size (dva)',
                  'with_labels': True}
    maps_names.append('size')
    
    # n
    n_data = all_deriv_mat[prf_n_idx,...]
    param_n = {'data': n_data, 
               'cmap': cmap_ecc_size, 
               'alpha': alpha, 
               'vmin': n_scale[0], 
               'vmax': n_scale[1], 
               'cbar': 'discrete', 
               'cortex_type': 'VertexRGB', 
               'description': 'CSS pRF n',
               'curv_brightness': 1, 
               'curv_contrast': 0.1, 
               'add_roi': False, 
               'cbar_label': 'pRF n',
               'with_labels': True}
    maps_names.append('n')
    
    # pcm
    pcm_data = all_deriv_mat[pcm_median_idx,...]
    param_pcm = {'data': pcm_data, 
                 'cmap': cmap_ecc_size, 
                 'alpha': alpha, 
                 'vmin': pcm_scale[0], 
                 'vmax': pcm_scale[1], 
                 'cbar': 'discrete', 
                 'cortex_type': 'VertexRGB', 
                 'description': 'CSS pRF CM',
                 'curv_brightness': 1, 
                 'curv_contrast': 0.1, 
                 'add_roi': False, 
                 'cbar_label': 'pRF CM (mm/dva)',
                 'with_labels': True}
    maps_names.append('pcm')
        
    # draw flatmaps
    volumes = {}
    for maps_name in maps_names:
    
        # create flatmap
        roi_name = 'pRF_{}'.format(maps_name)
        roi_param = {'subject': pycortex_subject, 'xfmname': None, 'roi_name': roi_name}
        print(roi_name)
        exec('param_{}.update(roi_param)'.format(maps_name))
        exec('volume_{maps_name} = draw_cortex(**param_{maps_name})'.format(maps_name=maps_name))
        exec("plt.savefig('{}/{}_task-{}_{}_css_transparency_0.pdf')".format(flatmaps_dir, subject, prf_task_name, maps_name))
        plt.close()
    
        # save flatmap as dataset
        exec('vol_description = param_{}["description"]'.format(maps_name))
        exec('volume = volume_{}'.format(maps_name))
        volumes.update({vol_description:volume})
    
    # save dataset
    dataset_file = "{}/{}_task-{}_css_loo-median.hdf".format(datasets_dir, subject, prf_task_name)
    if os.path.exists(dataset_file): os.system("rm -fv {}".format(dataset_file))
    dataset = cortex.Dataset(data=volumes)
    dataset.save(dataset_file)